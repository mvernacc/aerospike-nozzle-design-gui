# Aerospike Nozzle Design Solver Tests
# Matthew Vernacchia - MIT Rocket Team
# 2014 Jan 17

from nozzle_solver import *
import time


for Pa in [101e3, 50e3, 25e3, 10e3, 1e3]:
    alt = get_alt_from_Pa(Pa)
    print 'Atmo pressure of %.0f KPa at %.0f m altitude'%(Pa/1e3, alt)
    
test_params = EngineParameters(Tc=1900, Pc=5e6, molar_m=20.18, gamma=1.27, Re=0.015, er=12)

Me = get_Mach(test_params)
print 'Test params exit Mach number = %.2f (should be 3.64)'%(Me)

Pe = get_Pe(test_params)
print 'Test params exit pressure = %.2f (should be 40283.20 Pa)'%(Pe)

Re = get_Re_from_mass_flow(test_params, 0.2196)
print 'Test params exit radius = %.4f (should be 0.0150 m)'%(Re)

Re = get_Re_from_thrust(test_params, 470)
print 'Test params exit radius = %.4f (should be about 0.0150 m)'%(Re)

ns = NozzleSolver()
ns.solve(test_params, 57e3)
ns.plot()
print 'plotted once'
ns.plot()
print 'plotted twice'

time.sleep(10)