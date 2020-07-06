# Aerospike Nozzle Design Solver
# Matthew Vernacchia - MIT Rocket Team
# 2014 Jan 15

from math import atan, sin, tan, pi, asin
import numpy as np
from matplotlib import pyplot as plt
from collections import namedtuple

plt.ion()

### Warning logging ###
import logging
formatter = logging.Formatter('%(levelname)s {%(pathname)s:%(lineno)d}: %(message)s')
ch = logging.StreamHandler()
ch.setFormatter(formatter)
logger = logging.getLogger('simple_example')
logger.addHandler(ch)

### Physical Constants ###
# Molar Gas Constant [J mol^-1 K^-1]
R_mol = 8.314
# Acceleration due to gravity at surface [m s^-2]
g = 9.81


### Set up the engine parameters tuple ###
EngineParameters = namedtuple('EngineParameters', 'Tc Pc molar_m gamma Re er')
default_params = EngineParameters(Tc=1900, Pc=5e6, molar_m=20.18, gamma=1.27, Re=0.015, er=12)
zero_params = EngineParameters(Tc=0, Pc=0, molar_m=0, gamma=0, Re=0, er=0)

class NozzleSolver:
    def __init__(self):
        # Number of points along spike to solve for
        self.N = 10000
        # the mach numbers to examine
        self.M = np.zeros((self.N,))
        # the radius of the plug at some point x, normalized by the outer radius
        self.Rx_over_Re = np.zeros((self.N,))
        # The axial distance from the lip to point x, normalized by the outer
        # radius
        self.X_over_Re = np.zeros((self.N,))
        # The flow velocity to mach wave angle [rad]
        self.mu = np.zeros((self.N,))
        # the turning angle [rad]
        self.nu = np.zeros((self.N,))
        # the pressure at point x [Pa]
        self.P = np.zeros((self.N,))
        # the temperature at point x [K]
        self.T = np.zeros((self.N,))
        # the cumulative Isp up to point x [s]
        self.Isp = np.zeros((self.N,))
        # Thrust force [N]
        self.F = 0
        # Throat mass flow [kg s^-1]
        self.m_dot = 0
        # Throat Area [m^2]
        self.At = 0
        # Throat gap width [m]
        self. ht = 0
        # Angle between shroud surface and a plane normal to the nozzle axis [rad]
        self.delta = 0
        # Ambient pressure [Pa]
        self.Pa = 0
        # Nozzle Thrust Coefficent [-]
        self.Cf = 0
        self.last_params = zero_params
        self.results_string = ''

        # Figure and axes for plotting
        self.set_up_plot()

    def solve(self, params, Pa):
        ''' Solve for the spike contour and nozzle conditions given the engine parameters.
            params     Engine Parameters namedtuple
            Pa         Ambient pressure [Pa]
        '''
        self.last_params = params
        self.Pa = Pa
        y = params.gamma
        # Find the Specific Gas Constant
        R = R_mol / params.molar_m * 1000

        #### Compute the Exit Mach Number ####
        Me = get_Mach( params )
        # use this to find the exit (tip-plane) pressure  (assumung isetropic expansion)
        Pe = params.Pc * (1 + (y-1)/2*Me**2)**(-y/(y-1))
        
        #### Find the total flow turning angle [rad] ####
        nu_e = ((y+1)/(y-1))**0.5 * atan( ((y-1)/(y+1)*(Me**2-1))**0.5 ) \
            - atan( (Me**2-1)**0.5 )
        # use this to find the angle of the shroud edge with the normal to
        # the axial direction [rad]
        self.delta = pi/2 - nu_e
        
        #### Find the throat gap width / outer radius ratio [] ####
        htRe = (params.er - (params.er*(params.er-sin(self.delta)))**0.5) \
            / (params.er*sin(self.delta))
        
        #### Find the velocity and temperature at the throat (assuming M=1) ####
        # throat thermo temperature [K]
        Tt = params.Tc / (1 + (y-1)/2)
        # throat pressure [Pa]
        Pt = params.Pc * (1 + (y-1)/2) ** (-y/(y-1))
        # throat fluid velocity [m s^-1]
        vt = np.sqrt( y*R*Tt )

        #### Examine a range of Mach numbers to determine the shape of the nozzle
        # the Mach numbers to examine
        self.M = np.linspace(1, Me, self.N)
        # the Isp due to momentum flux and pressure at the throat [sec]
        self.Isp[0] = vt*sin(self.delta)*( 1 + 1/y*(1-((y+1)/2)**(y/(y-1))*Pa/params.Pc) ) / g
        for x in range(self.N):
            self.mu[x] = asin( 1/self.M[x] ) # See Lozano's Fall2012 Lec17 Notes
            # use the Prandtl-Meyer equation to find nu[x]
            self.nu[x] = ((y+1)/(y-1))**0.5 * atan( ((y-1)/(y+1)*(self.M[x]**2-1))**0.5 ) \
            - atan( (self.M[x]**2-1)**0.5 )   
            # use CC Lee Eqn (26) to find Rx/Re for the point
            RxRe2 = ( 1 - ( 2/(y+1)*(1+(y-1)/2*self.M[x]**2))**((y+1)/(2*y-2)) \
                * sin( nu_e - self.nu[x] + self.mu[x] ) / params.er )
            if RxRe2 > 0:
                self.Rx_over_Re[x] = RxRe2**0.5
            else: self.Rx_over_Re[x] = 0
            # Compute the phi angle for the point. See CC Lee Fig 2 for the definition of
            # the phi angle.
            phi_x = nu_e - self.nu[x] + self.mu[x]
            # Find the X (axial) coordinate of the point.
            # This uses a re-arrangement of CC Lee Eqn (19). Lee's Eqn (19) is:
            #   tan(phi_x) = (R_e - R_x) / X_x
            # Solve for X_x:
            #   X_x = (R_e - R_x) / tan()
            # Then divide both sides by R_e:
            #   X_x / R_e = (1 - R_x / R_e) / tan(phi)
            # Note that X_over_Re[x] is (X_x / R_e) and Rx_over_Re[x] is (R_x / R_e).
            self.X_over_Re[x] = (1 - self.Rx_over_Re[x]) / tan( phi_x )
            # find the pressure
            self.P[x] = params.Pc * (1 + (y-1)/2*self.M[x]**2)**(-y/(y-1))
            if x > 0:
                self.Isp[x] = self.Isp[x-1] + ( vt/y*((y+1)/2)**(y/(y-1)) * \
                    params.er/2 * ( (self.P[x-1]-Pa)/params.Pc + (self.P[x]-Pa)/params.Pc ) * \
                    (self.Rx_over_Re[x-1]**2 - self.Rx_over_Re[x]**2) ) / g
            #find the temperature
            self.T[x] = params.Tc / (1+ (y-1)/2 * self.M[x]**2)
        # throat width [meter]
        self.ht = htRe * params.Re
        # throat area [meter^2]
        self.At = pi*self.ht*(2*params.Re - self.ht*sin(self.delta))
        # fluid density at the throat [kg m^-3]
        rho_t = Pt / (R*Tt)
        # throat mass flow [kg sec^-1]
        self.m_dot = rho_t*self.At*vt
        # Find the thrust [N]
        self.F = self.Isp[-1] * self.m_dot * g
        # Find the thrust coefficient [-]
        self.Cf = self.F / (params.Pc*self.At)

        self.results_string = \
        'Engine Geometry:' \
        + '\n' + '\tShroud angle,     delta = %.1f degrees'%(self.delta*180/pi) \
        + '\n' + '\tShroud lip radius,   Re = %.1f mm'%(params.Re*1000) \
        + '\n' + '\tThroat width,        ht = %.2f mm'%(self.ht*1000) \
        + '\n' + '\tThroat area,         At = %.8f m^2'%(self.At) \
        + '\n' + '\tExpansion ratio,     er = %.2f'%(params.er) \
        + '\n' + 'Chamber Conditions:' \
        + '\n' + '\tChamber pressure,    Pc = %.3f MPa'%(params.Pc/1.0e6) \
        + '\n' + '\tChamber temperature, Tc = %.0f K'%(params.Tc) \
        + '\n' + '\tExhaust Avg Molar Mass  = %.1f g mol^-1'%(params.molar_m) \
        + '\n' + 'Engine Performance:' \
        + '\n' + '\tMass flow rate,   m_dot = %.3f kg sec^-1'%(self.m_dot) \
        + '\n' + '\tThrust force,         F = %.1f N'%(self.F) \
        + '\n' + '\tSpecific impulse,   Isp = %.1f sec'%(self.Isp[-1]) \
        + '\n' + '\tExit pressure,       Pe = %.2f Pa'%(Pe) \
        + '\n' + '\tExit Mach number,    Me = %.2f'%(Me) \
        + '\n' + '\tThrust coefficient,  Cf = %.2f'%(self.Cf)
        return

    def set_up_plot(self):
        '''Set up the figure and axes for plotting.'''
        self.fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(10, 6))
        self.axes = {
            'radius': axes[0][0],
            'pressure': axes[0][1],
            'Isp': axes[0][2],
            'Mach': axes[1][0],
            'temperature': axes[1][1],
            'text': axes[1][2],
        }
        self.axes['radius'].set_ylabel('$R_x / R_e$ [-]')
        self.axes['radius'].grid(True)

        self.axes['pressure'].set_ylabel('Pressure [MPa]')
        self.axes['pressure'].grid(True)

        self.axes['Isp'].set_ylabel(r'Cumulative $I_{sp}$ [s]')
        self.axes['Isp'].set_xlabel(r'$X_x \, / \, R_e$')
        self.axes['Isp'].grid(True)

        self.axes['Mach'].set_ylabel('Mach Number [-]')
        self.axes['Mach'].set_xlabel(r'$X_x \, / \, R_e$')
        self.axes['Mach'].grid(True)

        self.axes['temperature'].set_ylabel('Temperature [K]')
        self.axes['temperature'].set_xlabel(r'$X_x \, / \, R_e$')
        self.axes['temperature'].grid(True)

        self.fig.suptitle(
            'Simulated Nozzle Conditions vs Axial Distance from Throat Normalized by Exit Radius')

        self.fig.tight_layout()
        self.fig.subplots_adjust(top=0.93)

    def plot(self):
        '''Plot the current solution results for the nozzle conditions.'''
        for name, ax in self.axes.items():
            for artist in ax.lines + ax.texts:
                artist.remove()

        color = 'tab:blue'

        self.axes['radius'].plot( self.X_over_Re, self.Rx_over_Re, color=color )
  
        self.axes['Mach'].plot( self.X_over_Re, self.M, color=color )
        
 
        self.axes['pressure'].plot( self.X_over_Re, self.P/1.0e6, color=color)
        self.axes['pressure'].plot( self.X_over_Re, np.ones((self.N,))*self.Pa/1.0e6, 'r' )

        self.axes['temperature'].plot( self.X_over_Re, self.T, color=color )

        

        self.axes['Isp'].plot( self.X_over_Re, self.Isp, color=color )

        self.axes['text'].text(0.1,0.8, r'$A_e/A_t$ = %.1f'%(self.last_params.er))
        self.axes['text'].text(0.1,0.7, r'$P_c =$ %.2f MPa'%(self.last_params.Pc/1e6))
        self.axes['text'].text(0.1,0.6, r'$T_c =$ %.0f K'%(self.last_params.Tc))
        self.axes['text'].text(0.1,0.5, r'$M_{molar} =$ %.1f g mol^-1'%(self.last_params.molar_m))
        self.axes['text'].text(0.1,0.4, r'$\gamma =$ %.2f'%(self.last_params.gamma))

        plt.show()
        plt.pause(0.01)

### Utility Thermofluids Functions ###
def get_Mach( params ):
    '''
    Find the exit Mach number given the engine parameters.
    This function uses the params.er (expansion ratio) and params.gamma (gas specific heat ratio)

    Explicit Inversion of Stodola's Area-Mach Equation
    Source: J. Majdalani and B. A. Maickie
    http://maji.utsi.edu/publications/pdf/HT02_11.pdf
    '''    
    n = 5 # order of the aproximation
    X = np.zeros((n,))
    M = np.zeros((n,))

    e = 1/float(params.er) # expansion ratio
    y = params.gamma # ratio of specific heats
    B = (y+1)/(y-1)
    k = np.sqrt( 0.5*(y-1) )
    u = e**(1/B) / np.sqrt( 1+k**2 )
    X[0] = (u*k)**(B/(1-B))
    M[0] = X[0]

    for i in range(1,n):
        lamb = 1/( 2*M[i-1]**(2/B)*(B-2) + M[i-1]**2 *B**2*k**2*u**2 )
        X[i] = lamb*M[i-1]*B*( M[i-1]**(2/B) - M[i-1]**2*B*k**2*u**2 \
            + ( M[i-1]**(2+2/B)*k**2*u**2*(B**2-4*B+4) \
            - M[i-1]**2*B**2*k**2*u**4 + M[i-1]**(4/B)*(2*B-3) \
            + 2*M[i-1]**(2/B)*u**2*(2-B) )**0.5 )
        M[i] = M[i-1] + X[i]
    if abs( np.imag( M[n-1] ) ) > 1e-5:
        logger.warning('Exit Mach Number has nonzero imaginary part!')
    Me = np.real( M[n-1] )
    return Me

def get_Pe( params ):
    ''' Find the nozzle exit static pressure, given the engine parameters.
        returns       nozzle exit pressure [Pa]
    '''
    # Find the exit Mach number
    Me = get_Mach(params)
    y = params.gamma
    # Use this to find the exit (tip-plane) pressure  (assumung isetropic expansion)
    Pe = params.Pc * (1 + (y-1)/2*Me**2)**(-y/(y-1))
    return Pe

def get_Pa_from_alt( alt ):
    ''' Find the ambient atmospheric pressure, given the altitude above mean sea level in Earth's atmosphere.
        This functionnis only valid for input altitudes up to 44 km.
        alt       Altitiude above mean sea level [m]
        returns   Atmospheric pressure [Pa]
    '''
    if alt < 44331:
        # Pressure (http://psas.pdx.edu/RocketScience/PressureAltitude_Derived.pdf, eqn 9)
        Pa = 100 * ((44331.514 - alt)/(11880.516))**(1/0.1902632)
        return Pa
    else:
        logger.warning('Atmospheric pressure model does not extend to the altitude requested. returning 0 Pa pressure')
        return 0

def get_alt_from_Pa( Pa ):
    ''' Given the ambient atmospheric pressure, find the altitude above mean sea level in Earth's atmosphere.
        This functionnis only valid for input altitudes up to 44 km.
        alt       Altitiude above mean sea level [m]
        returns   Atmospheric pressure [Pa]
    '''
    if Pa > 0:
        # Altitude (http://psas.pdx.edu/RocketScience/PressureAltitude_Derived.pdf, eqn 10)
        alt = 44331.5 - 4946.62*Pa**0.190263
        return alt
    else:
        logger.warning('Negative pressure input. returning zero altitude')
        return 0

def get_Re_from_thrust( params, F ):
    ''' Find the shroud radius, given the thrust force and the other parameters.
        F          Desired pressure-matched thrust force [N]
        returns    Needed shroud radius [m]
    '''
    y = params.gamma
    Pe = get_Pe(params)
    # Find the pressure-matched thrust coefficient, using
    # Rocket Propulsion Elements 7th Ed, Equation 3-30
    C_F = ( (2*y**2)/(y-1) * (2/(y+1))**((y+1)/(y-1)) \
        * (1-(Pe/params.Pc)**((y-1)/y)) )**0.5
    # Use the thrust coefficient to find the throat area required
    At = F / (C_F * params.Pc)
    Ae = At * params.er
    # The Exit Area is a circle of radius Re
    Re = (Ae/pi)**0.5
    return Re

def get_thrust( params ):
    ''' Find the pressure-matched thrust of the engine, given the engine parameters.
        returns    Thrust force at matched pressure (Pa=Pe) [N]
    '''
    y = params.gamma
    Pe = get_Pe(params)
    # Find the pressure-matched thrust coefficient, using
    # Rocket Propulsion Elements 7th Ed, Equation 3-30
    C_F = ( (2*y**2)/(y-1) * (2/(y+1))**((y+1)/(y-1)) \
        * (1-(Pe/params.Pc)**((y-1)/y)) )**0.5
    # Find the throat area
    Ae = pi*params.Re**2
    At = Ae/params.er
    F = C_F*At*params.Pc
    return F

def get_Re_from_mass_flow( params, m_dot ):
    ''' Find the shroud radius, given the mass flow and the other parameters.
        m_dot      Desired mass flow [kg / s]
        returns    Needed shroud radius [m]
    '''
    y = params.gamma
    # Find the Specific Gas Constant
    R = R_mol / params.molar_m * 1000
    # Find the Throat Area require for the specified mass flow, using
    # Eocket Propul;sion Equations  7th Ed, Equation 3-24
    At = m_dot / ( params.Pc * y * (2/(y+1))**((y+1)/(2*y-2)) \
        / (y*R*params.Tc)**0.5)
    # Use the Throat Area and the Expansion Ratio to find the Exit Area
    Ae = At * params.er
    # The Exit Area is a circle of radius Re
    Re = (Ae/pi)**0.5
    return Re

def get_er_from_Pe( params, Pe ):
    ''' Find the needed expansion ratio, given a desired exit pressure and the other engine parameters.
        Pe        Desired exit pressure [Pa]
        returns   Needed expansion ratio [-]
    '''
    y  = params.gamma
    # Rocket Propulsion Elements 7th Ed, Equation 3-25
    AtAe = ((y+1)/2)**(1/(y-1)) * (Pe/params.Pc)**(1/y) \
        * ( (y+1)/(y-1)*( 1 - (Pe/params.Pc)**((y-1)/y) ) )**0.5
    er = 1/AtAe
    return er
