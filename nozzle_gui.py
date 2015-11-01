# Aerospike Nozzle Design GUI
# Matthew Vernacchia - MIT Rocket Team
# 2014 Jan 15
import Tkinter, tkFileDialog, tkMessageBox
from collections import namedtuple
import pickle
import numpy as np
from math import tan
import nozzle_solver

class NozzleGUI(Tkinter.Tk):
    def __init__(self, parent):
        # Engine Parameters storage
        self.params = nozzle_solver.default_params
        self.Pa = 57e3
        # set up file saving/loading options
        self.pkl_file_opt = {}
        self.pkl_file_opt['filetypes'] = [('all files', '.*'), ('pickle files', '.pkl')]
        self.pkl_file_opt['defaultextension'] = '.pkl'
        self.pkl_file_opt['parent'] = self
        # set up export file options
        self.export_file_opt = {}
        self.export_file_opt['filetypes'] = [('all files', '.*')]
        self.export_file_opt['defaultextension'] = ''
        self.export_file_opt['parent'] = self
        # init TKinter GUI
        Tkinter.Tk.__init__(self, parent)
        self.parent = parent
        self.initialize()
        # Solver
        self.ns = nozzle_solver.NozzleSolver()
        

    def initialize(self):
        # Use grid layout
        self.grid()

        # Chamber conditions
        chamb_con_label = Tkinter.Label(self, text='Chamber Conditions', anchor='w')
        chamb_con_label.grid(column=0, row=0, columnspan=2, sticky='EW')
        # Chamber Temp label
        Tc_label = Tkinter.Label(self, text='Chamber Temperature, Tc [K]', anchor='w')
        Tc_label.grid(column=0, row=1, sticky='EW')
        # Chamber Temp input box
        vcmd = (self.register(self.validate_input_float_nonneg), 
                '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')
        self.Tc_str = Tkinter.StringVar()
        self.Tc_str.set(str(self.params.Tc))
        self.Tc_box = Tkinter.Spinbox(from_=0, to=10000, textvariable=self.Tc_str, validate='all', validatecommand=vcmd, justify='right')
        self.Tc_box.grid(column=1, row=1, sticky='EW')
        self.Tc_box.bind('<Return>', self.On_Tc_change)
        self.Tc_box.bind('<FocusOut>', self.On_Tc_change)
        # Chamber pressure label
        Pc_label = Tkinter.Label(self, text='Chamber Pressure, Pc [MPa]', anchor='w')
        Pc_label.grid(column=0, row=2, sticky='EW')
        # Chamber Pressure input box
        vcmd = (self.register(self.validate_input_float_nonneg), 
                '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')
        self.Pc_str = Tkinter.StringVar()
        self.Pc_str.set(str(self.params.Pc/1e6))
        self.Pc_box = Tkinter.Spinbox(from_=0, to=100, textvariable=self.Pc_str, validate='all', validatecommand=vcmd, justify='right', format='%.2f', increment=0.1)
        self.Pc_box.grid(column=1, row=2, sticky='EW')
        self.Pc_box.bind('<Return>', self.On_Pc_change)
        self.Pc_box.bind('<FocusOut>', self.On_Pc_change)
        
        # Gas properties
        gas_label = Tkinter.Label(self, text='Gas Properties', anchor='w')
        gas_label.grid(column=0, row=3, columnspan=2, sticky='EW')
        # Molar mass label
        molar_m_label = Tkinter.Label(self, text='Molar Mass [g mol^-1]', anchor='w')
        molar_m_label.grid(column=0, row=4, sticky='EW')
        # Molar Mass input box
        vcmd = (self.register(self.validate_input_float_nonneg), 
                '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')
        self.molar_m_str = Tkinter.StringVar()
        self.molar_m_str.set(str(self.params.molar_m))
        self.molar_m_box = Tkinter.Spinbox(self, from_=0, to=1000, textvariable=self.molar_m_str, validate='all', validatecommand=vcmd, format='%.1f', justify='right')
        self.molar_m_box.grid(column=1, row=4, sticky='EW')
        self.molar_m_box.bind('<Return>', self.On_molar_m_change)
        self.molar_m_box.bind('<FocusOut>', self.On_molar_m_change)
        # Gamma label
        gamma_label = Tkinter.Label(self, text='Ratio of Specific Heats, gamma = cp/cv [-]', anchor='w')
        gamma_label.grid(column=0, row=5, sticky='EW')
        # Gamma input box
        vcmd = (self.register(self.validate_input_float_nonneg), 
                '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')
        self.gamma_str = Tkinter.StringVar()
        self.gamma_str.set(str(self.params.gamma))
        self.gamma_box = Tkinter.Spinbox(self, from_=1, to=2, textvariable=self.gamma_str, validate='all', validatecommand=vcmd, justify='right', format='%.2f', increment=0.1)
        self.gamma_box.grid(column=1, row=5, sticky='EW')
        self.gamma_box.bind('<Return>', self.On_gamma_change)
        self.gamma_box.bind('<FocusOut>', self.On_gamma_change)

        # Ambient conditions
        amb_label = Tkinter.Label(self, text='Ambient Pressure', anchor='w')
        amb_label.grid(column=0, row=6, columnspan=2, sticky='EW')
        self.which_amb = Tkinter.IntVar()
        self.Pa_button = Tkinter.Radiobutton(self, text = 'Ambient Pressure, Pa [KPa]', variable=self.which_amb, value='0', command=self.On_amb_change, anchor='w')
        self.Pa_button.grid(column=0, row=7, sticky='EW')
        # ambient Pressure input box
        vcmd = (self.register(self.validate_input_float_nonneg), 
                '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')
        self.Pa_str = Tkinter.StringVar()
        self.Pa_str.set(str(self.Pa/1e3))
        self.Pa_box = Tkinter.Spinbox(from_=0, to=110, textvariable=self.Pa_str, validate='all', validatecommand=vcmd, justify='right', format='%.1f', increment=1)
        self.Pa_box.grid(column=1, row=7, sticky='EW')
        self.Pa_box.bind('<Return>', self.On_Pa_change)
        self.Pa_box.bind('<FocusOut>', self.On_Pa_change)
        # altitude button
        self.alt_button = Tkinter.Radiobutton(self, text = 'Altitude AMSL [m]', variable=self.which_amb, value='1', command=self.On_amb_change, anchor='w')
        self.alt_button.grid(column=0, row=8, sticky='EW')
        # Altitude input box
        vcmd = (self.register(self.validate_input_float_nonneg), 
                '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')
        self.alt_str = Tkinter.StringVar()
        self.alt_str.set(str(nozzle_solver.get_alt_from_Pa(self.Pa)))
        self.alt_box = Tkinter.Spinbox(from_=0, to=44e3, textvariable=self.alt_str, validate='all', validatecommand=vcmd, justify='right', format='%.0f', increment=100)
        self.alt_box.grid(column=1, row=8, sticky='EW')
        self.alt_box.bind('<Return>', self.On_alt_change)
        self.alt_box.bind('<FocusOut>', self.On_alt_change)
        self.alt_box.configure(state='disabled')

        # Nozzle Expansion
        exp_label = Tkinter.Label(self, text='Nozzle Expansion', anchor='w')
        exp_label.grid(column=0, row=9, columnspan=2, sticky='EW')
        self.which_exp = Tkinter.IntVar()
        # Expansion ratio button
        self.er_button = Tkinter.Radiobutton(self, text = 'Expansion Area Ratio, er = Ae/At[-]', variable=self.which_exp, value='0', command=self.On_exp_change, anchor='w')
        self.er_button.grid(column=0, row=10, sticky='EW')
        # Expansion ratio input box
        vcmd = (self.register(self.validate_input_float_nonneg), 
                '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')
        self.er_str = Tkinter.StringVar()
        self.er_str.set(str(self.params.er))
        self.er_box = Tkinter.Spinbox(from_=0, to=100, textvariable=self.er_str, validate='all', validatecommand=vcmd, justify='right', format='%.1f', increment=1)
        self.er_box.grid(column=1, row=10, sticky='EW')
        self.er_box.bind('<Return>', self.On_er_change)
        self.er_box.bind('<FocusOut>', self.On_er_change)
        # Exit pressure button
        self.Pe_button = Tkinter.Radiobutton(self, text = 'Exit Pressure, Pe [KPa]', variable=self.which_exp, value='1', command=self.On_exp_change, anchor='w')
        self.Pe_button.grid(column=0, row=11, sticky='EW')
        # Exit Pressure input box
        vcmd = (self.register(self.validate_input_float_nonneg), 
                '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')
        self.Pe_str = Tkinter.StringVar()
        self.Pe_str.set(str( nozzle_solver.get_Pe(self.params)/1e3 ))
        self.Pe_box = Tkinter.Spinbox(from_=0, to=500, textvariable=self.Pe_str, validate='all', validatecommand=vcmd, justify='right', format='%.1f', increment=1)
        self.Pe_box.grid(column=1, row=11, sticky='EW')
        self.Pe_box.bind('<Return>', self.On_Pe_change)
        self.Pe_box.bind('<FocusOut>', self.On_Pe_change)
        self.Pe_box.configure(state='disabled')

        # Sizing
        size_label = Tkinter.Label(self, text='Sizing', anchor='w')
        size_label.grid(column=0, row=12, columnspan=2, sticky='EW')
        self.which_size = Tkinter.IntVar()
        # Trust button
        self.F_button = Tkinter.Radiobutton(self, text = 'Desired Thrust, F [N]', variable=self.which_size, value='0', command=self.On_size_change, anchor='w')
        self.F_button.grid(column=0, row=13, sticky='EW')
        # Thrust input box
        vcmd = (self.register(self.validate_input_float_nonneg), 
                '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')
        self.F_str = Tkinter.StringVar()
        self.F_str.set(str( nozzle_solver.get_thrust(self.params) ))
        self.F_box = Tkinter.Spinbox(from_=0, to=1e6, textvariable=self.F_str, validate='all', validatecommand=vcmd, justify='right', format='%.0f', increment=100)
        self.F_box.grid(column=1, row=13, sticky='EW')
        self.F_box.bind('<Return>', self.On_F_change)
        self.F_box.bind('<FocusOut>', self.On_F_change)
        # Exit radius button
        self.Re_button = Tkinter.Radiobutton(self, text = 'Exit Radius [mm]', variable=self.which_size, value='1', command=self.On_size_change, anchor='w')
        self.Re_button.grid(column=0, row=14, sticky='EW')
        # Exit Pressure input box
        vcmd = (self.register(self.validate_input_float_nonneg), 
                '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')
        self.Re_str = Tkinter.StringVar()
        self.Re_str.set(str( self.params.Re/1e3 ))
        self.Re_box = Tkinter.Spinbox(from_=0, to=500, textvariable=self.Re_str, validate='all', validatecommand=vcmd, justify='right', format='%.1f', increment=1)
        self.Re_box.grid(column=1, row=14, sticky='EW')
        self.Re_box.bind('<Return>', self.On_Re_change)
        self.Re_box.bind('<FocusOut>', self.On_Re_change)
        self.Re_box.configure(state='disabled')

        # Results Text Box
        self.results_text = Tkinter.Text(self)       
        self.results_text.config(state='disabled')
        self.results_text.grid(column=2,row=0, columnspan=1, rowspan=11)

        # Export Contour Button
        self.export_button = Tkinter.Button(self, text='Export Spike+Shroud Shape to CAD', command=self.On_export_button)
        self.export_button.grid(column=2, row=11)

        # Save Button
        self.save_button = Tkinter.Button(self, text='Save Design', command=self.On_save_button)
        self.save_button.grid(column=2, row=13)

        # Load Button
        self.load_button = Tkinter.Button(self, text='Load Design', command=self.On_load_button)
        self.load_button.grid(column=2, row=14)

    def resolve(self):
        # re-slove
        self.ns.solve(self.params, self.Pa)
        # re-plot
        self.ns.plot()
        # update results text with new results
        self.results_text.config(state='normal')
        self.results_text.delete('1.0','end')
        self.results_text.insert('1.0', self.ns.results_string)
        self.results_text.config(state='disabled')

    def On_Tc_change(self, event):
        try:
            Tc_new = float(self.Tc_box.get())
            if Tc_new >= 0:
                self.params = self.params._replace(Tc=Tc_new)
                print 'Tc = %f'%(self.params.Tc)
                self.resolve()
            return
        except ValueError:
            return

    def On_Pc_change(self, event):
        try:
            Pc_new = 1e6*float(self.Pc_box.get())
            if Pc_new >= 0:
                self.params = self.params._replace(Pc=Pc_new)
                print 'Pc = %f'%(self.params.Pc)
                self.resolve()
            return
        except ValueError:
            return

    def On_molar_m_change(self, event):
        try:
            molar_m_new = float(self.molar_m_box.get())
            if molar_m_new > 0:
                self.params = self.params._replace(molar_m=molar_m_new)
                print 'molar_m = %f'%(self.params.molar_m)
                self.resolve()
            return
        except ValueError:
            return

    def On_gamma_change(self, event):
        try:
            gamma_new = float(self.gamma_box.get())
            if gamma_new >= 0:
                self.params = self.params._replace(gamma=gamma_new)
                print 'gamma = %f'%(self.params.gamma)
                self.resolve()
            return
        except ValueError:
            return

    def On_amb_change(self):
        if self.which_amb.get() == 0:
            self.Pa_box.configure(state='normal')
            self.Pa_box.update()
            self.alt_box.configure(state='disabled')
            self.alt_box.update()
        else:
            self.Pa_box.configure(state='disabled')
            self.Pa_box.update()
            self.alt_box.configure(state='normal')
            self.alt_box.update()

    def On_Pa_change(self, event):
        try:
            Pa_new = 1e3*float(self.Pa_box.get())
            if Pa_new >= 0:
                self.Pa =Pa_new
                self.alt_str.set('%.0f'%( nozzle_solver.get_alt_from_Pa(self.Pa) ))
                print 'Pa = %f'%(self.Pa)
                self.resolve()
            return
        except ValueError:
            return

    def On_alt_change(self, event):
        try:
            alt_new = float(self.alt_box.get())
            self.Pa = nozzle_solver.get_Pa_from_alt(alt_new)
            self.Pa_str.set('%.1f'%( self.Pa/1e3 ))
            print 'Pa = %f'%(self.Pa)
            self.resolve()
            return
        except ValueError:
            return

    def On_exp_change(self):
        if self.which_exp.get() == 0:
            self.er_box.configure(state='normal')
            self.er_box.update()
            self.Pe_box.configure(state='disabled')
            self.Pe_box.update()
        else:
            self.er_box.configure(state='disabled')
            self.er_box.update()
            self.Pe_box.configure(state='normal')
            self.Pe_box.update()

    def On_er_change(self, event):
        try:
            er_new = float(self.er_box.get())
            if er_new > 0:
                self.params = self.params._replace(er=er_new)
                self.Pe_str.set('%.1f'%( nozzle_solver.get_Pe(self.params)/1e3 ))
                print 'er = %f'%(self.params.er)
                self.resolve()
            return
        except ValueError:
            return

    def On_Pe_change(self, event):
        try:
            Pe_new = 1e3*float(self.Pe_box.get())
            if Pe_new >= 0:
                er_new = nozzle_solver.get_er_from_Pe(self.params, Pe_new)
                self.params = self.params._replace(er=er_new)
                self.er_str.set('%.1f'%(self.params.er))
                print 'er = %f'%(self.params.er)
                self.resolve()
            return
        except ValueError:
            return

    def On_size_change(self):
        if self.which_size.get() == 0:
            self.F_box.configure(state='normal')
            self.F_box.update()
            self.Re_box.configure(state='disabled')
            self.Re_box.update()
        else:
            self.F_box.configure(state='disabled')
            self.F_box.update()
            self.Re_box.configure(state='normal')
            self.Re_box.update()

    def On_F_change(self, event):
        try:
            F_new = float(self.F_box.get())
            if F_new > 0:
                Re_new = nozzle_solver.get_Re_from_thrust(self.params, F_new)
                self.params = self.params._replace(Re=Re_new)
                self.Re_str.set( '%.1f'%(Re_new*1e3) )
                print 'Re = %f'%(self.params.Re)
                self.resolve()
            return
        except ValueError:
            return

    def On_Re_change(self, event):
        try:
            Re_new = 1e-3*float(self.Re_box.get())
            if Re_new > 0:
                self.params = self.params._replace(Re=Re_new)
                self.F_str.set('%.0f'%( nozzle_solver.get_thrust(self.params) ))
                print 'Re = %f'%(self.params.Re)
                self.resolve()
            return
        except ValueError:
            return

    def On_save_button(self):
        # Open the file via dialog
        save_file = tkFileDialog.asksaveasfile( mode='wb', title='Save Nozzle Design Parameters', **self.pkl_file_opt)
        # pickle the engine parameters to the file
        pickle.dump(self.params, save_file)
        # close the file
        save_file.close()

    def On_load_button(self):
        # open the file via dialog
        load_file = tkFileDialog.askopenfile(mode='rb', title='Load Nozzle Design Parameters', **self.pkl_file_opt)
        # try to unpickle the data
        try:
            self.params = pickle.load(load_file)
        except (PickleError, AttributeError, EOFError, ImportError):
            tkMessageBox.showerror('Load File Failed', 'Cannot load data from the selected file')
        self.Tc_str.set(str(self.params.Tc))
        self.Pc_str.set(str(self.params.Pc/1e6))
        self.molar_m_str.set(str(self.params.molar_m))
        self.gamma_str.set(str(self.params.gamma))
        self.er_str.set(str(self.params.er))
        self.On_er_change(None)
        self.Re_str.set(str(self.params.Re*1e3))
        self.On_Re_change(None)
        self.resolve()

    def On_export_button(self):
        file_name = tkFileDialog.asksaveasfilename(title='Export Spike and Shroud Curves to CAD Software', **self.export_file_opt)
        # trim any file extension
        a = file_name.partition('.')
        file_name = a[0]
        # Save the spike and shroud curve files
        #try:
        f = open(file_name+'_spike_curve.txt', 'w')
        for x in xrange(self.ns.N):
            f.write('%.4f,%.4f,%.4f\r\n'%(0, self.ns.X[x]*self.params.Re*-1000, self.ns.RxRe[x]*self.params.Re*1000))
        f.close()
        f = open(file_name+'_shroud_curve.txt', 'w')
        for r in np.arange(0, 0.003, 0.0001):
            f.write('%.4f,%.4f,%.4f\r\n'%( 0, r*tan(self.ns.delta)*1000, (self.params.Re+r)*1000))
        f.close()
        #except:
        #    tkMessageBox.showerror('Export Failed', \
        #        'Cannot export spike and shroud curves with the given filename:\n%s'%(file_name))


    def validate_input_float_nonneg(self, d, i, P, s, S, v, V, W):
        # Validate that a decimal number greater than zero was input into the box
        # %d = Type of action (1=insert, 0=delete, -1 for others)
        # %i = index of char string to be inserted/deleted, or -1
        # %P = value of the entry if the edit is allowed
        # %s = value of entry prior to editing
        # %S = the text string being inserted or deleted, if any
        # %v = the type of validation that is currently set
        # %V = the type of validation that triggered the callback
        #      (key, focusin, focusout, forced)
        # %W = the tk name of the widget
        if S in '0123456789.-+':
            try:
                x = float(P)
                #print 'Valid numeric input!'
                if x >= 0:
                    return True
                else:
                    return False
            except ValueError:
                #print 'Invalid non-numeric input'
                return False
        else:
            #print 'Invalid non-numeric input'
            return False

if __name__ == "__main__":
    app = NozzleGUI(None)
    app.title('Aerospike Nozzle Design')
    app.mainloop()