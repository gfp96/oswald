#from asyncore import close_all
import tkinter as tk
from tkinter import ttk, filedialog
from tkinter.messagebox import showinfo
from Functions import browseFiles, browseFolders, Set_graph_params, Read_BE_file, Update_df, Read_terratek_signal
import pandas as pd
import numpy as np
import seaborn as sns
import os
from datetime import datetime
from signal_interp import Find_start, Get_Max_AIC_velocity, Get_STALTA_AIC_velocity, Get_filtered_signal, Find_start_burst
from functools import partial
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backend_bases import MouseButton
matplotlib.use("TkAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg,
    NavigationToolbar2Tk
)
import ctypes #for stopping window resizing
 

#-------------------------------------------------------------------------------------------------------------
#----------------------------------------------------Launch software------------------------------------------
#-------------------------------------------------------------------------------------------------------------

class App(tk.Tk):
    def __init__(self):
        # Define general size
        super().__init__()
        # tk.Frame.__init__(self, parent)
        # self.parent = parent

        #turn off matplotlib interactive mode
        plt.ioff() #to switch OFF the interactive mode

        self.vratio = 9/10
        self.hratio = 1 # 9/10
        self.screen_width = self.winfo_screenwidth()
        self.screen_height = self.winfo_screenheight()
        self.window_width = int(self.screen_width*self.hratio)
        self.window_height = int(self.screen_height*self.vratio)
        # find the center point
        self.center_x = int(self.screen_width/2 - self.window_width / 2)
        self.center_y = int(self.screen_height/2 - self.window_height / 2)
        self.geometry(f'{self.window_width}x{self.window_height}+{self.center_x}+{self.center_y}')

        #Set folder location and data df
        self.folder = filedialog.askdirectory(initialdir = "C:\\Users\\gfloo\\Documents\\Professionnel\\PhD\\8_Experimental\\Exp_results\\Hostun_HN31\\HN31_ID60_P_1", title = "Select folder containing bender element files")
        self.datafile = filedialog.askopenfilename(initialdir = self.folder, title = "Select file containing experimental database",filetypes = (("Excel files", "*.xlsx*"), ("all files","*.*")))
        if self.folder=="" or self.datafile=="": #stop execution if use does not select anything
            quit()
        self.data_original = pd.read_excel(self.datafile)
        self.data_original = self.data_original.replace([np.inf, -np.inf], np.nan) # replace inf values by nan
        self.data_original = self.data_original[(self.data_original.Sweep==False)&(self.data_original.freqlev>=0)] #(self.data_original.valid==True)
        self.data = self.data_original.drop_duplicates(subset={"isvp", "freqlev", "stageno"}, keep="first").copy(deep=True) #"pelev", "pflev", 
        #ensure that all filtering columns are integers
        self.data.stageno = self.data.stageno.astype(int)
        self.data.pelev = self.data.pelev.astype(int)
        self.data.pflev = self.data.pflev.astype(int)
        self.data.freqlev = self.data.freqlev.astype(float)
        #ensure bools are bool
        self.data.isvp = self.data.isvp.astype(bool)
        self.data.Sweep = self.data.Sweep.astype(bool) # tells if signals are sweeps
        self.data.Burst = self.data.Burst.astype(bool) # tells if signals are bursts
        self.data.valid = self.data.valid.astype(bool)
        #self.data.dry = self.data.dry.astype(bool)
        #self.data.saturated = self.data.saturated.astype(bool)
        #self.data.drainage = self.data.drainage.astype(bool)
        #self.data.consolidated = self.data.consolidated.astype(bool)
        self.data.valid_aicmax = self.data.valid_aicmax.astype(bool)
        self.data.valid_SLA = self.data.valid_SLA.astype(bool)
        #make lists
        self.pe_list = np.unique(self.data.pelev)
        self.pf_list = np.unique(self.data.pflev)
        self.stageno_list = np.unique(self.data.stageno)
        #Array for storing manual vels
        #self.manual_vels = -np.ones(self.data.shape[0], dtype = float) #filled with -1 as default value
        self.maxaic_vels = -np.ones(self.data.shape[0], dtype = float) #filled with -1 as default value
        self.stalta_aic_vels = -np.ones(self.data.shape[0], dtype = float) #filled with -1 as default value
        #Default interpretation method : defines shape of zoom plots
        self.method = "Max_AIC" #default

        #check if new columns have been made, otherwise create them
        if "Grade" in self.data:
            pass
        else:
            Update_df(self.data, "Grade", -1)

        # Check if the column has already been created
        if "vel_manual" in self.data:
            pass
        else:
            Update_df(self.data, "vel_manual", -1.0)

        #Set grid
        # configure the grid
        self.cols = 10
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.columnconfigure(2, weight=1)
        self.columnconfigure(3, weight=1)
        self.columnconfigure(4, weight=1)
        self.columnconfigure(5, weight=1)
        self.columnconfigure(6, weight=1)
        self.columnconfigure(7, weight=1)
        self.columnconfigure(8, weight=1)
        self.columnconfigure(9, weight=1)
        self.rows = 20
        for i in range(self.rows): # configure rows
            self.rowconfigure(i, weight=1)
        

        # Create a File Explorer : this is original start point of code
        self.label_folder = ttk.Label(self, text = f"Folder: {os.path.basename(os.path.normpath(self.folder))}",width = 50, foreground = "black")
        self.label_folder.grid(column = 8, row = 0, columnspan=3)

        self.button_explore = ttk.Button(self,text = "Folder",command = browseFolders)
        self.button_explore.grid(column = 9,  row = 8)
        # set the disabled flag
        self.button_explore.state(['disabled'])

        #Place textboxes for details of possible datasets
        #Make combobox of effective stress
        ttk.Label(self, text="p' [kPa]").grid(column=0, row = 0)
        nlist = self.pe_list.tolist()
        nlist.append("Any")
        self.pe_cbox = ttk.Combobox(self, values = nlist, state = "readonly") #np.append(self.pe_list, -1)
        self.pe_cbox.grid(column=1, row=0)
        self.pe_cbox.bind("<<ComboboxSelected>>", self.Update_pe_choice)

        #Make combobox of fluid pressure
        ttk.Label(self, text="pf [kPa]").grid(column=2, row = 0)
        nlist = self.pf_list.tolist()
        nlist.append("Any")
        self.pf_cbox = ttk.Combobox(self, values = nlist, state = "readonly")
        self.pf_cbox.grid(column=3, row=0)
        self.pf_cbox.bind("<<ComboboxSelected>>", self.Update_pf_choice)

        #Make combobox of stageno
        nlist = np.unique(self.data.stageno).tolist()
        nlist.append("Any")
        ttk.Label(self, text="Stage n°").grid(column=4, row = 0)
        self.stage_cbox = ttk.Combobox(self, values = nlist, state = "readonly")
        self.stage_cbox.grid(column=5, row=0)
        self.stage_cbox.bind("<<ComboboxSelected>>", self.Update_stage_choice)

        #Make method for LTT
        self.Ltt = 0.193 #tip to tip length of BE
        ttk.Label(self, text="Ltt [m]").grid(column=6, row = 0)
        self.Ltt_entry = ttk.Entry(text=f"{self.Ltt}")
        self.Ltt_entry.grid(column=7, row=0)
        self.Ltt_entry.insert(0, f"{self.Ltt}")
        self.Ltt_entry.bind("<1>", self.Update_Ltt)


        # Configure event : window size change
        self.first_change_time = datetime.now() #time at which window was changed
        self.configure_time = np.timedelta64(5, 's')
        self.bind("<Configure>", self.On_window_resize)

        #make sure exit "x" actually closes software, not just window
        self.protocol("WM_DELETE_WINDOW", self._quit)
    

        #-------------------------------------------------------------------------------------------------------------
        #----------------------------------------------------Make Plots----------------------------------------------
        #-------------------------------------------------------------------------------------------------------------

        #Make multiplot
        Set_graph_params()
        self.px = 1/plt.rcParams['figure.dpi']  # pixel in inches
        self.Make_figmain()
        self.Make_figfreq()
        self.Make_figzoom()
        
        #Check is anything has been plotted
        self.virginmain = True
        self.virginfreq = True
        self.virginzoom = True

        #Check is manual Vp has been clicked
        self.Pclicked = False
        #Check is manual Vs has been clicked
        self.Sclicked = False
        #Store which is last click, defult : none
        self.Lastclick = "None"
        #Setting txtboxes
        self.props = dict(boxstyle='round', facecolor='white', alpha=0)

        #-------------------------------------------------------------------------------------------------------------
        #----------------------------------------------------Main button----------------------------------------------
        #-------------------------------------------------------------------------------------------------------------
        #Button to plot signals
        self.button_plotsig = ttk.Button(self, text = "Get Data", command = self.Plot_all_signals)
        self.button_plotsig.grid(column = 9,  row = 2)

        #-------------------------------------------------------------------------------------------------------------
        #----------------------------------------------------Choose interpretation------------------------------------
        #-------------------------------------------------------------------------------------------------------------
        #Button to plot signals
        self.methods = ["None", "Filter", "Max_AIC", "STA/LTA_AIC", "Show CC"]
        ttk.Label(self, text="Method").grid(column=9, row = 3)
        self.method_cbox = ttk.Combobox(self, values = self.methods) #, state = "disabled")
        self.method_cbox.grid(column=9, row=4)
        self.method_cbox.bind("<<ComboboxSelected>>", self.Choose_method)

        #Button to choose signal format
        self.formats = ["Navier_BE", "Terratek"]
        self.format = "Navier_BE" #default format
        ttk.Label(self, text="File format").grid(column=9, row = 6)
        self.format_cbox = ttk.Combobox(self, values = self.formats) #, state = "disabled")
        self.format_cbox.grid(column=9, row=7)
        self.format_cbox.bind("<<ComboboxSelected>>", self.Choose_format)

        #If Terratek data format, indicate encap to encap travel time
        self.ee_timeP = 4.6e-6 #encap to encap travel time P-wave
        ttk.Label(self, text="E-E time (P) [µs]").grid(column=9, row = 8)
        self.ee_timeP_entry = ttk.Entry(text=f"{int(np.round(1e6*self.ee_timeP))}")
        self.ee_timeP_entry.grid(column=9, row=9)
        self.ee_timeP_entry.insert(0, f"{np.round(1e6*self.ee_timeP, 1)}")
        self.ee_timeP_entry.bind("<1>", self.Update_ee_timeP)
        self.ee_timeP_entry.state(["disabled"])

        self.ee_timeS = 7.1e-6 #encap to encap travel time S-wave
        ttk.Label(self, text="E-E time (S) [µs]").grid(column=9, row = 10)
        self.ee_timeS_entry = ttk.Entry(text=f"{int(np.round(1e6*self.ee_timeS))}")
        self.ee_timeS_entry.grid(column=9, row=11)
        self.ee_timeS_entry.insert(0, f"{np.round(1e6*self.ee_timeS, 1)}")
        self.ee_timeS_entry.bind("<1>", self.Update_ee_timeS)
        self.ee_timeS_entry.state(["disabled"])

        #frequency ranges
        self.freq_rangeP = [100, 200e3] #acceptable frequency range for P waves
        self.freq_rangeS = [100, 100e3] #acceptable frequency range for S waves
        self.min_frange = 50e3 #default width of filtering bandpass

        

        #-------------------------------------------------------------------------------------------------------------
        #----------------------------------------------------Interpretation----------------------------------------------
        #-------------------------------------------------------------------------------------------------------------
        #Button to interpret signals
        self.button_interpret = ttk.Button(self, text = "Analyse", command = self.Interpret_signals)
        self.button_interpret.grid(column = 9,  row = 5)
        self.button_interpret.state(["disabled"])
        
        
        #partial(Plot_all_signals, figure_canvas, figure, ax, data, pe_cbox, pf_cbox, folder, button_interpret)
        #Button to Refresh plots
        # self.button_refresh = ttk.Button(self, text = "Refresh")#, command = partial(Plot_all_signals, figure_canvas, figure, ax, data, pe_cbox, pf_cbox, folder))
        # self.button_refresh.grid(column = 9,  row = 8)
        # self.button_refresh.state(["disabled"])

        #For future Cbox freq
        ttk.Label(self, text="Frequency [kHz]").grid(column=9, row = 12)
        self.freq_cbox = ttk.Combobox(self, values = self.data.freqlev, state = "disabled")
        self.freq_cbox.grid(column=9, row=13)
        self.freq_cbox.bind("<<ComboboxSelected>>", self.Zoom_on_freq)
        self.freq = -1 #-1 is default value : indicates no data has been selected

        #Make method for Burst type signals
        self.nbperiods_burst = 21 #tip to tip length of BE
        ttk.Label(self, text="Periods (burst)").grid(column=9, row = 14)
        self.nbper_entry = ttk.Entry(text=f"{self.nbperiods_burst}")
        self.nbper_entry.grid(column=9, row=15)
        self.nbper_entry.insert(0, f"{self.nbperiods_burst}")
        self.nbper_entry.bind("<1>", self.Update_nbperiods_burst)

        #bind tab key
        self.bind("<Tab>", self.Press_tab)

        #bind number keys to grading the signal
        self.bind("<Key>", self.Grading)
        # self.bind("<KP_0>", self.Grading)
        # self.bind("<KP_1>", self.Grading)
        # self.bind("<KP_2>", self.Grading)
        # self.bind("<KP_3>", self.Grading)
        # self.bind("<KP_4>", self.Grading)
        # self.bind("<KP_5>", self.Grading)
        
        #self.freq_cbox.bind("<Return>", self.Zoom_on_freq)
        
        #Button to plot manual vel calculation on freqplot
        self.button_plotmanual = ttk.Button(self, text = "Plot Manual", command = self.Plot_manual)
        self.button_plotmanual.grid(column = 9,  row = 16)
        self.button_plotmanual.state(["disabled"])

        #Button to Save manual vel calculation (all freqs)
        self.button_save = ttk.Button(self, text = "Save", command = self.Save)#,
        self.button_save.grid(column = 9,  row = 18)
        self.button_save.state(["disabled"])



        
        # self.button_interpret = ttk.Button(self, text = "Interpretation", command = Interpret_signals)
    #-------------------------------------------------------------------------------------------------------------
    #----------------------------------------------------Functions------------------------------------------------
    #------------------------------------------------------------------------------------------------------------
    def On_window_resize(self, event):
        if datetime.now()-self.first_change_time>self.configure_time:
            self.screen_width = self.winfo_width()
            self.screen_height = self.winfo_height()
            self.window_width = int(self.screen_width)
            self.window_height = int(self.screen_height)
            # find the center point
            self.center_x = int(self.screen_width/2 - self.window_width / 2)
            self.center_y = int(self.screen_height/2 - self.window_height / 2)
            self.geometry(f'{self.window_width}x{self.window_height}+{self.center_x}+{self.center_y}')

            #replot figure
            # self.px = 1/plt.rcParams['figure.dpi']  # pixel in inches
            # self.Make_figmain()
            # self.Make_figfreq()
            # self.Make_figzoom()
            self.first_change_time = datetime.now()

    #Function for buuilding array and plotting 
    def Plot_all_signals(self): #, canvas, ax, data, sel_pe, sel_pf, folder):
        #Get filter p' and pf
        # selpe = self.pe_cbox.get()
        # selpf = self.pf_cbox.get()
        # if selpe=="" and selpf=="":
        #     self.sel_data = self.data.copy(deep=True)
        # elif selpe=="":
        #     self.sel_data = self.data[self.data.pflev==int(selpf)].copy(deep=True)
        # elif selpf=="":
        #     self.sel_data = self.data[self.data.pelev==int(selpe)].copy(deep=True)
        # else:
        #     self.sel_data = self.data[(self.data.pflev==int(selpf))&(self.data.pelev==int(selpe))].copy(deep=True)
        #look at one file for array sizes

        #Filter data according to stage number instead of pe/pf
        stage = self.stage_cbox.get()
        if stage=="" or stage=="Any":
            self.sel_data = self.data.copy(deep=True)
        else:
            self.sel_data = self.data[self.data.stageno==int(stage)].copy(deep=True)

        #Make arrays to atore all signal data in them (empty zeros at end opf each array)
        self.sel_size = self.sel_data.shape[0]
        self.indexes = np.arange(self.sel_size, dtype = int)
        self.rough = [None] * self.sel_size #np.zeros((size, self.selsize), dtype = float)
        self.ins = [None] * self.sel_size #np.zeros_like(self.rough)
        self.times = [None] * self.sel_size # np.zeros_like(self.rough)
        if self.format=="Navier_BE":
            for i, filename in enumerate(self.sel_data.filename):
                self.times[i], self.ins[i], self.rough[i] = Read_BE_file(self.folder, filename) #time, input signal, unfiltered output signal
        elif self.format=="Terratek":
            self.instimeP, self.eePsig = Read_terratek_signal(self.encapencapP) #self.instimesP records times for encap encap file
            self.instimeS, self.eeSsig = Read_terratek_signal(self.encapencapS) #self.instimesS records times for encap encap file
            for i, filename in enumerate(self.sel_data.filename):
                self.times[i], self.rough[i] = Read_terratek_signal(os.path.join(self.folder, filename)) #time, input signal, unfiltered output signal
                # if self.sel_data.isvp.iloc[i]: #P wave
                #     self.ins[i] = eePsig
                # else: #S wave
                #     self.ins[i] = eeSsig
        #Get max values to normalize plot+
        if self.format=="Terratek":
            self.maxamp = max(np.nanmax(self.eeSsig), np.nanmax(self.eePsig))
        else: #if self.format=="Navier_BE":
            self.maxamp = max([np.nanmax(self.ins[i])+np.abs(np.nanmin(self.ins[i])) for i in self.indexes])     
        
        self.amp_rat = self.maxamp/(max([np.nanmax(self.rough[i])+np.abs(np.nanmin(self.rough[i])) for i in self.indexes]))
        self.colours = sns.color_palette(n_colors = 2, as_cmap=False)  # set color scheme
        #Setting txtboxes
        props = dict(boxstyle='round', facecolor='white', alpha=0)
        # place a text box in upper left in axes coords
        self.unique_freqs = np.unique(self.sel_data.freqlev)
        lfreqs = self.unique_freqs.size
        full_size = (lfreqs+5)
        offset_labels = 1/lfreqs
        
        starttime = datetime.now()
        
        #self.Clear_plots(self.figmain, self.virginmain)
        if self.virginmain==False:
            plt.close(self.figmain)
            self.Make_figmain()
        #Loop while making plots
        for i, freq in enumerate(self.unique_freqs):
            #Vp
            name = self.sel_data[(self.sel_data.freqlev==freq)&(self.sel_data.isvp)].filename
            if name.size!=0:
                if name.size==1:
                    ind = self.sel_data.index.get_loc(self.sel_data[(self.sel_data.filename==name.item())].index.item())
                else:
                    name = name.iloc[0]
                    ind = self.sel_data.index.get_loc(self.sel_data[(self.sel_data.filename==name)].index.item())
                if self.format=="Terratek":
                    self.ax.plot(self.eePsig+i*self.maxamp, self.instimeP, color= self.colours[0], alpha=0.5)
                else:#all other cases : BE_navier format 
                    self.ax.plot(self.ins[ind]+i*self.maxamp, self.times[ind], color= self.colours[0], alpha=0.5)
                self.ax.plot(self.rough[ind]*self.amp_rat+i*self.maxamp, self.times[ind], color= self.colours[0], label="Vp")
                self.virginmain=False #indicates that something has been plotted
            #Vs
            name = self.sel_data[(self.sel_data.freqlev==freq)&(self.sel_data.isvp==False)].filename
            if name.size!=0:
                if name.size==1:
                    ind = self.sel_data.index.get_loc(self.sel_data[(self.sel_data.filename==name.item())].index.item())
                else:
                    name = name.iloc[0]
                    ind = self.sel_data.index.get_loc(self.sel_data[(self.sel_data.filename==name)].index.item())
                if self.format=="Terratek":
                    self.ax.plot(self.eeSsig+i*self.maxamp, self.instimeS, color= self.colours[1], alpha=0.5)
                else:
                    self.ax.plot(self.ins[ind]+i*self.maxamp, self.times[ind], color= self.colours[1], alpha=0.5)
                self.ax.plot(self.rough[ind]*self.amp_rat+i*self.maxamp, self.times[ind], color= self.colours[1], label = "Vs")
                self.virginmain=False #indicates that something has been plotted
            #Frequency txt box
            self.ax.text((1+i)/full_size, 0.02, f"{np.round(freq)}", transform=self.ax.transAxes, fontsize=10,va='bottom', ha="center", bbox=props)
        #ax.legend(loc= "center right")
        self.ax.set_xlim(-self.maxamp, (4+lfreqs)*self.maxamp)
        handles, labels = self.ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        self.ax.legend(by_label.values(), by_label.keys(), loc = "center right")
        self.figmain_canvas.draw()
        print(datetime.now()-starttime)
        #Modify buttons
        self.button_interpret.state(["!disabled"])
        #Make combobox of frequency choice
        self.freq_cbox['values']=self.unique_freqs.tolist()
        self.Lastclick=="None" #no click has been made
        



    def Interpret_signals(self):
        if self.virginfreq==False:
            plt.close(self.figfreq)
            self.Make_figfreq()

        self.s_inds = np.zeros(self.sel_size, dtype = int)
        for i in self.indexes:
            if self.format=="Terratek":
                if self.sel_data.isvp.iloc[i]: #P wave
                    self.s_inds[i] = int(np.round(self.ee_timeP/(self.times[i][1]-self.times[i][0]))) # index corresponding to arrival
                else: #S wave
                    self.s_inds[i] = int(np.round(self.ee_timeS/(self.times[i][1]-self.times[i][0]))) # index corresponding to arrival
            elif self.sel_data.Burst.iloc[i]: #signal is burst type
                self.s_inds[i] = Find_start_burst(self.ins[i], self.nbperiods_burst)
            else:
                self.s_inds[i] = Find_start(self.ins[i])
        #Array for storing manual vels at given filtering
        #self.manual_vels_sel = -100*np.ones(self.sel_data.shape[0], dtype = float) #filled with -100 as default value
        #TODO Make the filtering depend on the wavetype : 0:50 kHz for Vs
        if self.method in ["Filter", "Show CC"]: #only filter the signal
            self.filtered = Get_filtered_signal(self.times, self.s_inds, self.rough, self.sel_data.isvp, self.sel_data.infreq, self.freq_rangeP, self.freq_rangeS)
        if self.method in ["Max_AIC", "STA/LTA_AIC"]:    
            #AIC interpretation
            if self.method=="Max_AIC":
                self.interpdata, self.aics, self.filtered = Get_Max_AIC_velocity(self.times, self.s_inds, self.rough, self.Ltt, self.sel_data.infreq, self.freq_rangeP) #
            elif self.method=="STA/LTA_AIC":
                self.interpdata, self.aic1s, self.aic2s, self.filtered = Get_STALTA_AIC_velocity(self.times, self.s_inds, self.rough, self.Ltt, self.sel_data.infreq, self.freq_rangeP)
                vels_calc = -1*np.ones(self.s_inds.size, dtype = float)
            #PLots
            #Plot on main plot
            for i, freq in enumerate(self.unique_freqs):
                #Vp
                name = self.sel_data[(self.sel_data.freqlev==freq)&(self.sel_data.isvp)].filename
                if name.size!=0:
                    ind = self.sel_data.index.get_loc(self.sel_data[(self.sel_data.filename==name.item())].index.item())
                    ind_big = self.data.index.get_loc(self.data[(self.data.filename==name.item())].index.item()) #index in non filtered array
                    if self.method=="Max_AIC":
                        arrival = self.interpdata.iloc[ind].a_time
                        self.maxaic_vels[ind_big] = self.interpdata.iloc[ind].vel
                    elif self.method=="STA/LTA_AIC": #For P waves : first arrival should be the right one
                        arrival = self.interpdata.iloc[ind].a1_time
                        self.stalta_aic_vels[ind_big] = self.interpdata.iloc[ind].vel1
                        vels_calc[ind] = self.interpdata.iloc[ind].vel1
                    self.ax.plot(np.array([(i-0.5)*self.maxamp, (i+0.5)*self.maxamp], dtype = float), np.array([arrival, arrival], dtype = float), color= self.colours[0], linestyle = "dotted")
                #Vs
                name = self.sel_data[(self.sel_data.freqlev==freq)&(self.sel_data.isvp==False)].filename
                if name.size!=0:
                    ind = self.sel_data.index.get_loc(self.sel_data[(self.sel_data.filename==name.item())].index.item())
                    ind_big = self.data.index.get_loc(self.data[(self.data.filename==name.item())].index.item()) #index in non filtered array
                    if self.method=="Max_AIC":
                        arrival = self.interpdata.iloc[ind].a_time
                        self.maxaic_vels[ind_big] = self.interpdata.iloc[ind].vel
                    elif self.method=="STA/LTA_AIC": #For S waves : second arrival should be the right one
                        arrival = self.interpdata.iloc[ind].a2_time
                        self.stalta_aic_vels[ind_big] = self.interpdata.iloc[ind].vel2
                        vels_calc[ind] = self.interpdata.iloc[ind].vel2
                    self.ax.plot(np.array([(i-0.5)*self.maxamp, (i+0.5)*self.maxamp], dtype = float), np.array([arrival, arrival], dtype = float), color= self.colours[1], linestyle = "dotted")
            self.figmain_canvas.draw()
            #Adding columns to dataframe
            # Check if the column interpretation has already been created
            if "vel_aic" in self.sel_data:
                if self.method=="Max_AIC":
                    self.sel_data.vel_aic = np.array(self.interpdata.vel) # if yes : replace the data
                elif self.method=="STA/LTA_AIC":
                    self.sel_data.vel_aic = vels_calc
            else:
                if self.method=="Max_AIC":
                    self.sel_data.insert(loc=len(self.sel_data.columns), column = "vel_aic", value = np.array(self.interpdata.vel)) #if no : create
                elif self.method=="STA/LTA_AIC":
                    self.sel_data.insert(loc=len(self.sel_data.columns), column = "vel_aic", value = vels_calc)
            #Plot 
            sns.scatterplot(ax= self.axf, x=self.sel_data[self.sel_data.isvp].freqlev, y=self.sel_data[self.sel_data.isvp].vel_aic, label = "Vp", color = self.colours[0], marker = "+")
            sns.scatterplot(ax= self.axf, x=self.sel_data[self.sel_data.isvp==False].freqlev, y=self.sel_data[self.sel_data.isvp==False].vel_aic, label = "Vs", color = self.colours[1], marker = "+")
        elif self.method in ["None", "Filter"]:
            if "vel_aicmax" in self.sel_data:
                sns.scatterplot(ax= self.axf, x=self.sel_data[self.sel_data.isvp].freqlev, y=self.sel_data[self.sel_data.isvp].vel_aicmax, label = "Vp-MAIC", color = self.colours[0], marker = "1")
                sns.scatterplot(ax= self.axf, x=self.sel_data[self.sel_data.isvp==False].freqlev, y=self.sel_data[self.sel_data.isvp==False].vel_aicmax, label = "Vs-MAIC", color = self.colours[1], marker = "1")
            if "vel_SLA" in self.sel_data:
                sns.scatterplot(ax= self.axf, x=self.sel_data[self.sel_data.isvp].freqlev, y=self.sel_data[self.sel_data.isvp].vel_SLA, label = "Vp-SLA", color = self.colours[0], marker = "x")
                sns.scatterplot(ax= self.axf, x=self.sel_data[self.sel_data.isvp==False].freqlev, y=self.sel_data[self.sel_data.isvp==False].vel_SLA, label = "Vs-SLA", color = self.colours[1], marker = "x")
        elif self.method in ["Show CC"]: #Mode to show Cross correlation
            if "vel_CC" in self.sel_data:
                sns.scatterplot(ax= self.axf, x=self.sel_data[self.sel_data.isvp].freqlev, y=self.sel_data[self.sel_data.isvp].vel_CC, label = "Vp-CC", color = self.colours[0], marker = "1")
                sns.scatterplot(ax= self.axf, x=self.sel_data[self.sel_data.isvp==False].freqlev, y=self.sel_data[self.sel_data.isvp==False].vel_CC, label = "Vs-CC", color = self.colours[1], marker = "1")
            if "vel_SLA" in self.sel_data:
                sns.scatterplot(ax= self.axf, x=self.sel_data[self.sel_data.isvp].freqlev, y=self.sel_data[self.sel_data.isvp].vel_SLA, label = "Vp-SLA", color = self.colours[0], marker = "x")
                sns.scatterplot(ax= self.axf, x=self.sel_data[self.sel_data.isvp==False].freqlev, y=self.sel_data[self.sel_data.isvp==False].vel_SLA, label = "Vs-SLA", color = self.colours[1], marker = "x")
        
        #For all cases
        if "vel_manual" in self.sel_data:
            sns.scatterplot(ax= self.axf, x=self.sel_data[self.sel_data.isvp].freqlev, y=self.sel_data[self.sel_data.isvp].vel_manual, label = "Vp-Manual", color = self.colours[0], marker = ".")
            sns.scatterplot(ax= self.axf, x=self.sel_data[self.sel_data.isvp==False].freqlev, y=self.sel_data[self.sel_data.isvp==False].vel_manual, label = "Vs-Manual", color = self.colours[1], marker = ".")

        #xlim, ylim and save
        self.axf.set_ylim(bottom=0, top=2000)
        self.axf.set_xlim(0, np.nanmax(self.unique_freqs))
        self.axf.legend(loc = "lower right")
        #Setting txtboxes
        props = dict(boxstyle='round', facecolor='white', alpha=0)
        #Frequency txt box
        self.axf.text(0.05, 0.05, f"p' = {np.round(np.nanmean(self.sel_data.pelev))} kPa\npf = {np.round(np.nanmean(self.sel_data.pflev))} kPa", transform=self.ax.transAxes, fontsize=10,va='bottom', ha="left", bbox=props)
        self.figfreq_canvas.draw()
        self.virginfreq=False
        self.freq_cbox.state(["!disabled"])
        self.button_plotmanual.state(["!disabled"])
        self.button_save.state(["!disabled"])
        self.freq = -2 #-2 indicates no freq  values have yet been chosen

    def Zoom_on_freq(self, event):
        freq = float(self.freq_cbox.get()) #int(self.freq_cbox.get())
        self.freq = freq
        #if self.virginzoom==False:
        plt.close(self.figzoom)
        self.Make_figzoom()
        # click lines are always erased by previous command:
        self.Pclicked=False
        self.Sclicked=False
        self.Lastclick=="None" #no click has been made
        #self.Clear_plots(self.figzoom, self.virginzoom)
        lwt = 1 #special linewidth for aic plots
        #Vp
        name = self.sel_data[(self.sel_data.freqlev==freq)&(self.sel_data.isvp)].filename
        if name.size!=0:
            ind = self.sel_data.index.get_loc(self.sel_data[(self.sel_data.filename==name.item())].index.item())
            time = self.times[ind]*1e3
            if self.method in ["Max_AIC", "STA/LTA_AIC"]: # if aic analysis -> color on input signal
                if self.format=="Terratek":
                    self.axz1.plot(self.instimeP*1e3, self.eePsig, color= self.colours[0], alpha = 0.2)
                else:
                    self.axz1.plot(time, self.ins[ind], color= self.colours[0], alpha = 0.5)
                self.axz1.axvline(x=time[self.s_inds[ind]], color=self.colours[0], linestyle = "dashed", linewidth=lwt)
            elif self.method in ["Filter", "None", "Show CC"]: # otherwise show black input and coloured output
                if self.format=="Terratek":
                    self.axz1.plot(self.instimeP*1e3, self.eePsig, color= "black", alpha = 0.2)
                else:
                    self.axz1.plot(time, self.ins[ind], color= "black", alpha = 0.5)
                self.axz1.axvline(x=time[self.s_inds[ind]], color="black", linestyle = "dashed", linewidth=lwt)
            if self.method in ["None"]: #show rough signal
                self.axz1b.plot(time, self.rough[ind]*self.amp_rat, color= self.colours[0], label = "Vp", alpha = 0.5)
            elif self.method in ["Filter", "Show CC", "Max_AIC", "STA/LTA_AIC"]:
                self.axz1b.plot(time, self.filtered[ind]*self.amp_rat, color= self.colours[0], label = "Vp", alpha = 1)
            #self.axz[0].plot(self.times[:, ind], self.filtered[:, ind]*self.amp_rat, color= self.colours[0], label = "Vp")
            # self.axz1.axvline(x=time[self.s_inds[ind]], color=self.colours[0], linestyle = "dashed", linewidth=lwt)
            if self.method=="Max_AIC":
                #Aic
                self.axz2.plot(time[self.interpdata.iloc[ind].start_aic:self.interpdata.iloc[ind].max_ind+1], self.aics[ind], color= self.colours[0], linewidth=lwt)
                #Vlines
                self.axz1.axvline(x=self.interpdata.iloc[ind].a_time*1e3, color=self.colours[0], linestyle = "dotted")
                self.axz2.axvline(x=self.interpdata.iloc[ind].a_time*1e3, color=self.colours[0], linestyle = "dotted", linewidth=lwt)
                self.axz2.axvline(x=time[self.interpdata.iloc[ind].start_aic], color=self.colours[0], linestyle = "dashdot", linewidth=lwt)
                self.axz2.axvline(x=time[self.interpdata.iloc[ind].max_ind], color=self.colours[0], linestyle = "dashdot", linewidth=lwt)
            elif self.method=="STA/LTA_AIC":
                #Aic
                self.axz2.plot(time[self.interpdata.iloc[ind].saic1_ind:self.interpdata.iloc[ind].eaic1_ind], self.aic1s[ind], color= self.colours[0], linestyle = "dashdot", linewidth=lwt)
                self.axz2b.plot(time[self.interpdata.iloc[ind].saic2_ind:self.interpdata.iloc[ind].eaic2_ind], self.aic2s[ind], color= self.colours[0], linestyle = "dashed", linewidth=lwt)
                #Vlines
                self.axz1.axvline(x=self.interpdata.iloc[ind].a1_time*1e3, color=self.colours[0], linestyle = "dotted")#, linewidth=lwt)
                self.axz2.axvline(x=self.interpdata.iloc[ind].a1_time*1e3, color=self.colours[0], linestyle = "dashdot", linewidth=lwt)
                self.axz2.axvline(x=time[self.interpdata.iloc[ind].saic1_ind], color=self.colours[0], linestyle = "dashdot", linewidth=lwt)
                self.axz2.axvline(x=time[self.interpdata.iloc[ind].eaic1_ind], color=self.colours[0], linestyle = "dashdot", linewidth=lwt)
                self.axz2.axvline(x=self.interpdata.iloc[ind].a2_time*1e3, color=self.colours[0], linestyle = "dashed", linewidth=lwt)
                self.axz2.axvline(x=time[self.interpdata.iloc[ind].saic2_ind], color=self.colours[0], linestyle = "dashed", linewidth=lwt)
                self.axz2.axvline(x=time[self.interpdata.iloc[ind].eaic2_ind], color=self.colours[0], linestyle = "dashed", linewidth=lwt)
            elif self.method in ["None", "Filter"]:
                if "vel_aicmax" in self.sel_data:
                    if self.sel_data.vel_aicmax.iloc[ind]!=-1:
                        arrival = time[self.s_inds[ind]]+self.Ltt/self.sel_data.vel_aicmax.iloc[ind]*1e3
                        self.axz1.axvline(x=arrival, color=self.colours[0], linestyle = "dotted", label = "Vp-MAIC", linewidth=lwt)
                if "vel_SLA" in self.sel_data:
                    if self.sel_data.vel_SLA.iloc[ind]!=-1:
                        arrival = time[self.s_inds[ind]]+self.Ltt/self.sel_data.vel_SLA.iloc[ind]*1e3
                        self.axz1.axvline(x=arrival, color=self.colours[0], linestyle = "dashdot", label = "Vp-SLA", linewidth=lwt)
            elif self.method in ["Show CC"]:
                if "vel_CC" in self.sel_data:
                    if self.sel_data.vel_CC.iloc[ind]!=-1:
                        arrival = time[self.s_inds[ind]]+self.Ltt/self.sel_data.vel_CC.iloc[ind]*1e3
                        self.axz1.axvline(x=arrival, color=self.colours[0], linestyle = "dotted", label = "Vp-CC", linewidth=lwt)
                if "vel_SLA" in self.sel_data:
                    if self.sel_data.vel_SLA.iloc[ind]!=-1:
                        arrival = time[self.s_inds[ind]]+self.Ltt/self.sel_data.vel_SLA.iloc[ind]*1e3
                        self.axz1.axvline(x=arrival, color=self.colours[0], linestyle = "dashdot", label = "Vp-SLA", linewidth=lwt)
            self.virginzoom=False

        #Vs
        name = self.sel_data[(self.sel_data.freqlev==freq)&(self.sel_data.isvp==False)].filename
        if name.size!=0:
            ind = self.sel_data.index.get_loc(self.sel_data[(self.sel_data.filename==name.item())].index.item())
            time = self.times[ind]*1e3 #pass to ms
            if self.method in ["None"]: #show rough signal
                self.axz2b.plot(time, self.rough[ind]*self.amp_rat, color= self.colours[1], label = "Vs", alpha = 0.5)
            elif self.method in ["Filter", "Show CC"]: # show filtered signal
                self.axz2b.plot(time, self.filtered[ind]*self.amp_rat, color= self.colours[1], label = "Vs", alpha = 1)
            elif self.method in ["Max_AIC", "STA/LTA_AIC"]: # show aic calculatiion
                self.axz1c.plot(time, self.filtered[ind]*self.amp_rat, color= self.colours[1], label = "Vs", alpha = 1)
            #self.axz[0].plot(self.times[ind], self.filtered[ind]*self.amp_rat, color= self.colours[1], label = "Vs")
            
            if self.method in ["None", "Filter", "Show CC"]: #Condition necessary as None and filter separate by wavetype
                if self.format=="Terratek":
                    self.axz2.plot(self.instimeS*1e3, self.eeSsig, color= "black", alpha = 0.2)
                else:
                    self.axz2.plot(time, self.ins[ind], color= "black", alpha = 0.5)
                self.axz2.axvline(x=time[self.s_inds[ind]], color="black", linestyle = "dashed", linewidth=lwt)
            elif self.method in ["Max_AIC", "STA/LTA_AIC"]:
                if self.format=="Terratek":
                    self.axz1.plot(self.instimeS*1e3, self.eeSsig, color= self.colours[1], alpha = 0.2)
                else:
                    self.axz1.plot(time, self.ins[ind], color= self.colours[1], alpha = 0.5)
                self.axz1.axvline(x=time[self.s_inds[ind]], color=self.colours[1], linestyle = "dashed")

            if self.method=="Max_AIC":
                #Aic
                self.axz2b.plot(time[self.interpdata.iloc[ind].start_aic:self.interpdata.iloc[ind].max_ind+1], self.aics[ind], color= self.colours[1], linewidth=lwt)
                #Vlines
                self.axz1.axvline(x=self.interpdata.iloc[ind].a_time*1e3, color=self.colours[1], linestyle = "dotted")
                self.axz2.axvline(x=self.interpdata.iloc[ind].a_time*1e3, color=self.colours[1], linestyle = "dotted", linewidth=lwt)
                self.axz2.axvline(x=time[self.interpdata.iloc[ind].start_aic], color=self.colours[1], linestyle = "dashdot", linewidth=lwt)
                self.axz2.axvline(x=time[self.interpdata.iloc[ind].max_ind], color=self.colours[1], linestyle = "dashdot", linewidth=lwt)
            elif self.method=="STA/LTA_AIC":
                #Aic
                self.axz2c.plot(time[self.interpdata.iloc[ind].saic1_ind:self.interpdata.iloc[ind].eaic1_ind], self.aic1s[ind], color= self.colours[1], linestyle = "dashdot", linewidth=lwt)
                self.axz2d.plot(time[self.interpdata.iloc[ind].saic2_ind:self.interpdata.iloc[ind].eaic2_ind], self.aic2s[ind], color= self.colours[1], linestyle = "dashed", linewidth=lwt)
                #Vlines
                self.axz1.axvline(x=self.interpdata.iloc[ind].a2_time*1e3, color=self.colours[1], linestyle = "dotted")
                self.axz2.axvline(x=self.interpdata.iloc[ind].a1_time*1e3, color=self.colours[1], linestyle = "dashdot", linewidth=lwt)
                self.axz2.axvline(x=time[self.interpdata.iloc[ind].saic1_ind], color=self.colours[1], linestyle = "dashdot", linewidth=lwt)
                self.axz2.axvline(x=time[self.interpdata.iloc[ind].eaic1_ind], color=self.colours[1], linestyle = "dashdot", linewidth=lwt)
                self.axz2.axvline(x=self.interpdata.iloc[ind].a2_time*1e3, color=self.colours[1], linestyle = "dashed", linewidth=lwt)
                self.axz2.axvline(x=time[self.interpdata.iloc[ind].saic2_ind], color=self.colours[1], linestyle = "dashed", linewidth=lwt)
                self.axz2.axvline(x=time[self.interpdata.iloc[ind].eaic2_ind], color=self.colours[1], linestyle = "dashed", linewidth=lwt)
            elif self.method in ["None", "Filter"]:
                if "vel_aicmax" in self.sel_data:
                    if self.sel_data.vel_aicmax.iloc[ind]!=-1:
                        arrival = time[self.s_inds[ind]]+self.Ltt/self.sel_data.vel_aicmax.iloc[ind]*1e3
                        self.axz2.axvline(x=arrival, color=self.colours[1], linestyle = "dotted", label = "Vs-MAIC", linewidth=lwt)
                if "vel_SLA" in self.sel_data:
                    if self.sel_data.vel_SLA.iloc[ind]!=-1:
                        arrival = time[self.s_inds[ind]]+self.Ltt/self.sel_data.vel_SLA.iloc[ind]*1e3
                        self.axz2.axvline(x=arrival, color=self.colours[1], linestyle = "dashdot", label = "Vs-SLA", linewidth=lwt)
            elif self.method in ["Show CC"]:
                if "vel_CC" in self.sel_data:
                    if self.sel_data.vel_CC.iloc[ind]!=-1:
                        arrival = time[self.s_inds[ind]]+self.Ltt/self.sel_data.vel_CC.iloc[ind]*1e3
                        self.axz2.axvline(x=arrival, color=self.colours[1], linestyle = "dotted", label = "Vs-CC", linewidth=lwt)
                if "vel_SLA" in self.sel_data:
                    if self.sel_data.vel_SLA.iloc[ind]!=-1:
                        arrival = time[self.s_inds[ind]]+self.Ltt/self.sel_data.vel_SLA.iloc[ind]*1e3
                        self.axz2.axvline(x=arrival, color=self.colours[1], linestyle = "dashdot", label = "Vs-SLA", linewidth=lwt)
                    
            self.virginzoom=False
        self.axz1.legend(loc="lower right")
        if self.method in ["None", "Filter", "Show CC"]: # show legend for Vs plot
            self.axz2.legend(loc="lower right")
        self.figzoom_canvas.draw()

    def Onclick(self, event):
        freq = float(self.freq_cbox.get())
        x= event.xdata
        if event.button is MouseButton.LEFT:
            #Pwave for left click
            name = self.sel_data[(self.sel_data.freqlev==freq)&(self.sel_data.isvp)].filename
            if name.size!=0:
                if self.Pclicked:
                    self.Pclickline.remove()
                ind = self.sel_data.index.get_loc(self.sel_data[(self.sel_data.filename==name.item())].index.item())
                start = self.times[ind][self.s_inds[ind]] #self.interpdata.iloc[ind].s_time
                vel = self.Ltt/(x/1e3-start) 
                ind_big = self.data.index.get_loc(self.data[(self.data.filename==name.item())].index.item()) #index in non filtered array
                #self.manual_vels[ind_big] = vel
                #self.manual_vels_sel[ind] = vel
                self.data.loc[self.data.filename==name.item(), "vel_manual"] = vel #send straight to dataframe
                self.sel_data.loc[self.sel_data.filename==name.item(), "vel_manual"] = vel #send straight to selection dataframe
                self.Pclickline = self.axz1.axvline(x=x, color=self.colours[0], linestyle = "solid")
                self.figzoom_canvas.draw()
                self.Pclicked=True #manual Vp has been clicked 
                self.Lastclick = "Pclick" #last click is on Pwave signal
        elif event.button is MouseButton.RIGHT:
            #Swave for right click
            name = self.sel_data[(self.sel_data.freqlev==freq)&(self.sel_data.isvp==False)].filename
            if name.size!=0:
                if self.Sclicked:
                    self.Sclickline.remove()
                    # if self.method in ["None", "Filter", "Show CC"]: #Condition necessary as None and filter separate by wavetype
                    #     self.axz2.lines.remove(self.Sclickline)
                    # elif self.method in ["Max_AIC", "STA/LTA_AIC"]:
                    #     self.axz1.lines.remove(self.Sclickline)
                ind = self.sel_data.index.get_loc(self.sel_data[(self.sel_data.filename==name.item())].index.item())
                start = self.times[ind][self.s_inds[ind]] #self.interpdata.iloc[ind].s_time
                vel = self.Ltt/(x/1e3-start) 
                ind_big = self.data.index.get_loc(self.data[(self.data.filename==name.item())].index.item()) #index in non filtered array
                #self.manual_vels[ind_big] = vel
                #self.manual_vels_sel[ind] = vel
                self.data.loc[self.data.filename==name.item(), "vel_manual"] = vel #send straight to dataframe
                self.sel_data.loc[self.sel_data.filename==name.item(), "vel_manual"] = vel #send straight to selection dataframe
                if self.method in ["None", "Filter", "Show CC"]: #Condition necessary as None and filter separate by wavetype
                    self.Sclickline = self.axz2.axvline(x=x, color=self.colours[1], linestyle = "solid")
                elif self.method in ["Max_AIC", "STA/LTA_AIC"]:
                    self.Sclickline = self.axz1.axvline(x=x, color=self.colours[1], linestyle = "solid")
                self.figzoom_canvas.draw()
                self.Sclicked=True #manual Vs has been clicked
                self.Lastclick = "Sclick" #last click is on Swave signal

    def Update_Ltt(self, event):
        self.Ltt = float(self.Ltt_entry.get())

    def Update_nbperiods_burst(self, event):
        self.nbperiods_burst = float(self.nbper_entry.get())
    
    def Update_ee_timeP(self, event):
        """Reads value and updates ee_timeP"""
        self.ee_timeP = float(self.ee_timeP_entry.get())/1e6 #pass from micro s to s

    def Update_ee_timeS(self, event):
        """Reads value and updates ee_timeS"""
        self.ee_timeS = float(self.ee_timeS_entry.get())/1e6 #pass from micro s to s

    def Clear_plots(self, fig, virgin):
        if virgin==False: #if virgin, no lines to delete
            for ax in fig.axes:
                for l in ax.lines:
                    ax.lines.remove(l)

    def Choose_method(self, event):
        self.method = str(self.method_cbox.get())
        self.Interpret_signals() #Launch interpretation

    def Choose_format(self, event):
        """Changes format of input wave signals"""
        self.format = str(self.format_cbox.get())
        if self.format=="Terratek":
            self.encapencapP = filedialog.askopenfilename(initialdir = self.folder, title = "Select file P-wave encap-encap recording",filetypes = (("TRC files", "*.TRC*"), ("all files","*.*")))
            self.ee_timeP_entry.state(["!disabled"])
            self.encapencapS = filedialog.askopenfilename(initialdir = self.folder, title = "Select file S-wave encap-encap recording",filetypes = (("TRC files", "*.TRC*"), ("all files","*.*")))
            self.ee_timeS_entry.state(["!disabled"])
            #frequency ranges
            self.freq_rangeP = [100, 5e6] #acceptable frequency range for P waves
            self.freq_rangeS = [100, 5e6] #acceptable frequency range for S waves
            self.min_frange = 1e6 #default width of filtering bandpass
        elif self.format=="Navier_BE":
            #frequency ranges
            self.freq_rangeP = [100, 100e3] #acceptable frequency range for P waves
            self.freq_rangeS = [100, 100e3] #acceptable frequency range for S waves
            self.min_frange = 50e3 #default width of filtering bandpass

    def Update_pe_choice(self, event):
        selpe = self.pe_cbox.get()
        #Make combobox of frequency choice
        if selpe=="" or selpe=="Any": #nothing selected
            pflist = np.unique(self.data.pflev).tolist()
            stagelist=np.unique(self.data.stageno).tolist()
        else:
            pflist = np.unique(self.data[self.data.pelev==int(selpe)].pflev).tolist()
            stagelist=np.unique(self.data[self.data.pelev==int(selpe)].stageno).tolist()
        pflist.append("Any")
        stagelist.append("Any")
        self.pf_cbox['values'] = pflist
        self.stage_cbox['values'] = stagelist


    def Update_pf_choice(self, event):
        selpf = self.pf_cbox.get()
        #Make combobox of frequency choice
        if selpf=="" or selpf=="Any": #nothing selected
            pelist = np.unique(self.data.pelev).tolist()
            stagelist=np.unique(self.data.stageno).tolist()
        else:
            pelist = np.unique(self.data[self.data.pflev==int(selpf)].pelev).tolist()
            self.stage_cbox['values']=np.unique(self.data[self.data.pflev==int(selpf)].stageno).tolist()
        pelist.append("Any")
        stagelist.append("Any")
        self.pe_cbox['values']=pelist
        self.stage_cbox['values'] = stagelist


    def Update_stage_choice(self, event):
        stage = self.stage_cbox.get()
        #Make combobox of frequency choice
        if stage=="" or "Any": #nothing selected
            pelist = np.unique(self.data.pelev).tolist()
            pflist = np.unique(self.data.pflev).tolist()
        else:
            pelist = np.unique(self.data[self.data.stageno==int(stage)].pelev).tolist()
            pflist = np.unique(self.data[self.data.stageno==int(stage)].pflev).tolist()
        pelist.append("Any")
        pflist.append("Any")
        self.pe_cbox['values']=pelist
        self.pf_cbox['values']=pflist
    
    def Press_tab(self, event):
        if self.freq==-1: #no data has been selected
            pass
        else: #okay to proceed
            if self.freq==-2: #indicates data selected but no frequencies : choose first one
                self.freq = self.unique_freqs[0]
                self.freq_cbox.set(self.unique_freqs[0])
            elif self.freq==self.unique_freqs[-1]: #last value of array choose first value
                self.freq = self.unique_freqs[0]
                self.freq_cbox.set(self.unique_freqs[0])
            else : # any other value
                ind = np.argmin(np.abs(self.freq-self.unique_freqs))
                self.freq = self.unique_freqs[ind+1]
                self.freq_cbox.set(self.unique_freqs[ind+1])
            self.Zoom_on_freq(event)

    def Grading(self, event):
        try:
            grade = int(event.char)
            if grade in [0, 1, 2, 3, 4, 5]:
                freq = float(self.freq_cbox.get())
                if self.Lastclick=="Pclick": # grading Pwave signal
                    name = self.sel_data[(self.sel_data.freqlev==freq)&(self.sel_data.isvp)].filename.item()
                    self.data.loc[self.data.filename==name, "Grade"] = grade
                    try:
                        self.txtPgrade.remove()
                    except:
                        pass
                    self.txtPgrade = self.axz1.text(0.05, 0.95, f"Pwave grade = {grade}", transform=self.axz1.transAxes, fontsize=10,va='top', ha="left", bbox=self.props)
                    self.figzoom_canvas.draw()
                    if grade==0: #means signal is unreadable
                        self.data.loc[self.data.filename==name, "valid"] = False
                elif self.Lastclick=="Sclick": # grading Swave signal
                    name = self.sel_data[(self.sel_data.freqlev==freq)&(self.sel_data.isvp==False)].filename.item()
                    self.data.loc[self.data.filename==name, "Grade"] = grade
                    try:
                        self.txtSgrade.remove()
                    except:
                        pass
                    self.txtSgrade = self.axz2.text(0.05, 0.95, f"Swave grade = {grade}", transform=self.axz2.transAxes, fontsize=10,va='top', ha="left", bbox=self.props)
                    self.figzoom_canvas.draw()
                    if grade==0: #means signal is unreadable
                        self.data.loc[self.data.filename==name, "valid"] = False
                #if other than Pclick or Sclick : do nothing
            return
        except:
            pass




            
    #-------------------------------------------------------------------------------------------------------------
    #----------------------------------------------------Make Plots----------------------------------------------
    #-------------------------------------------------------------------------------------------------------------

    #Make multiplot
    def Make_figmain(self):
        sns.set_style("white")
        self.figmain, self.ax = plt.subplots(1,1, figsize=(6/10*self.window_width*self.px, 5/10*self.window_height*self.px), tight_layout=True)
        # create FigureCanvasTkAgg object
        self.figmain_canvas = FigureCanvasTkAgg(self.figmain, self)
        # create the toolbar
        #NavigationToolbar2Tk(self.figmain_canvas, self, pack_toolbar=False)
        #Plotting
        self.ax.axes.yaxis.set_ticklabels([])
        self.ax.axes.xaxis.set_ticklabels([])
        self.ax.invert_yaxis()
        self.ax.set_ylabel("Time")
        self.ax.set_xlabel("Frequency [kHz]")
        self.figmain_canvas.draw()
        self.figmain_canvas.get_tk_widget().grid(column = 0, row = 1, columnspan=6, rowspan=10, sticky="W")

    def Make_figfreq(self):
        #Make Freqplot
        sns.set_style("white")
        self.figfreq, self.axf = plt.subplots(1, 1, sharex=True, figsize=(3/10*self.window_width*self.px, 5/10*self.window_height*self.px), tight_layout=True)
        # create FigureCanvasTkAgg object
        self.figfreq_canvas = FigureCanvasTkAgg(self.figfreq, self)
        #Plotting
        # axf.axes.yaxis.set_ticklabels([])
        # axf.axes.xaxis.set_ticklabels([])
        self.axf.set_ylabel("Velocity [m/s]")
        # self.axf[0].set_ylabel("Vp [m/s]")
        self.axf.set_xlabel("Frequency [kHz]")
        self.figfreq_canvas.draw()
        self.figfreq_canvas.get_tk_widget().grid(column = 6, row = 1, columnspan=3, rowspan=10, sticky = "W")
    
    def Make_figzoom(self):
        #Make zoomed plot
        sns.set_style("white") #("ticks")#, {"ytick.left":False, "ytick.right":False})
        #plt.tick_params(left = False, right=False)
        if self.method in ["", "None", "Filter", "Show CC"]:
            self.figzoom, self.axz = plt.subplots(2,1, sharex = True, figsize=(9/10*self.window_width*self.px, 4/10*self.window_height*self.px), tight_layout=True)
            # create FigureCanvasTkAgg object
            self.figzoom_canvas = FigureCanvasTkAgg(self.figzoom, self)
            #Plotting
            self.axz1 = self.axz[0]
            self.axz2 = self.axz[1]
            self.axz1.axes.yaxis.set_ticklabels([])
            self.axz2.axes.yaxis.set_ticklabels([])
            #self.axz1.axes.yaxis.set_tick_params(left=False, right=False)
            self.axz2.set_xlabel("Time [ms]")
            self.axz1.set_ylabel("P waves")
            self.axz2.set_ylabel("S waves")
            self.axz1b =self.axz1.twinx()
            self.axz2b =self.axz2.twinx() 
            self.axz1b.axes.yaxis.set_ticklabels([])
            self.axz2b.axes.yaxis.set_ticklabels([])
            
        elif self.method in ["Max_AIC", "STA/LTA_AIC"]: #=="Max_AIC":
            self.figzoom, self.axz = plt.subplots(2,1, sharex = True, figsize=(9/10*self.window_width*self.px, 4/10*self.window_height*self.px), tight_layout=True)
            # create FigureCanvasTkAgg object
            self.figzoom_canvas = FigureCanvasTkAgg(self.figzoom, self)
            # create the toolbar
            #NavigationToolbar2Tk(self.figzoom_canvas, self, pack_toolbar=False)
            #Plotting
            self.axz1 = self.axz[0]
            self.axz2 = self.axz[1]
            self.axz1.axes.yaxis.set_ticklabels([])
            self.axz2.axes.yaxis.set_ticklabels([])
            self.axz2.set_xlabel("Time [ms]")
            self.axz1.set_ylabel("Amplitude")
            self.axz2.set_ylabel("AIC")
            self.axz1b =self.axz1.twinx()
            self.axz1c =self.axz1.twinx()
            self.axz2b =self.axz2.twinx() 
            self.axz1b.axes.yaxis.set_ticklabels([])
            self.axz2b.axes.yaxis.set_ticklabels([])
            self.axz1c.axes.yaxis.set_ticklabels([])
        if self.method == "STA/LTA_AIC": #add twinxs for second aic plots
            self.axz2c =self.axz2.twinx() 
            self.axz2c.axes.yaxis.set_ticklabels([])
            self.axz2d =self.axz2.twinx() 
            self.axz2d.axes.yaxis.set_ticklabels([])
        self.figzoom_canvas.draw()
        self.figzoom_canvas.get_tk_widget().grid(column = 0, row = 11, columnspan=9, rowspan=9, sticky = "W")
        plt.connect('button_press_event', self.Onclick)

    def Plot_manual(self):
        # Check if the column has already been created
        # if "vel_manual" in self.sel_data:
        #     return_vels = np.nanmax([self.manual_vels_sel, np.array(self.sel_data.vel_manual)], axis = 0)
        #     self.sel_data.vel_manual = return_vels # if yes : replace the data
        # else:
        #     self.sel_data.insert(loc=len(self.sel_data.columns), column = "vel_manual", value = self.manual_vels_sel) #if no : create
        #Remove existing plot
        if self.virginfreq==False:
            plt.close(self.figfreq)
            self.Make_figfreq()
        #AIC vels
        if self.method in ["Max_AIC", "STA/LTA_AIC"]:
            sns.scatterplot(ax= self.axf, x=self.sel_data[self.sel_data.isvp].freqlev, y=self.sel_data[self.sel_data.isvp].vel_aic, label = "Vp", color = self.colours[0], marker = "+")
            sns.scatterplot(ax= self.axf, x=self.sel_data[self.sel_data.isvp==False].freqlev, y=self.sel_data[self.sel_data.isvp==False].vel_aic, label = "Vs", color = self.colours[1], marker = "+")
        elif self.method in ["None", "Filter"]:
            if "vel_aicmax" in self.sel_data:
                sns.scatterplot(ax= self.axf, x=self.sel_data[self.sel_data.isvp].freqlev, y=self.sel_data[self.sel_data.isvp].vel_aicmax, label = "Vp-MAIC", color = self.colours[0], marker = "1")
                sns.scatterplot(ax= self.axf, x=self.sel_data[self.sel_data.isvp==False].freqlev, y=self.sel_data[self.sel_data.isvp==False].vel_aicmax, label = "Vs-MAIC", color = self.colours[1], marker = "1")
            if "vel_SLA" in self.sel_data:
                sns.scatterplot(ax= self.axf, x=self.sel_data[self.sel_data.isvp].freqlev, y=self.sel_data[self.sel_data.isvp].vel_SLA, label = "Vp-SLA", color = self.colours[0], marker = "x")
                sns.scatterplot(ax= self.axf, x=self.sel_data[self.sel_data.isvp==False].freqlev, y=self.sel_data[self.sel_data.isvp==False].vel_SLA, label = "Vs-SLA", color = self.colours[1], marker = "x")

        #Manual vels
        sns.scatterplot(ax= self.axf, x=self.sel_data[self.sel_data.isvp].freqlev, y=self.sel_data[self.sel_data.isvp].vel_manual, label = "Vp-Manual", color = self.colours[0], marker = ".")
        sns.scatterplot(ax= self.axf, x=self.sel_data[self.sel_data.isvp==False].freqlev, y=self.sel_data[self.sel_data.isvp==False].vel_manual, label = "Vs-Manual", color = self.colours[1], marker = ".")
        #self.axf.set_ylim(0, max(np.nanmax(self.sel_data.vel_manual), np.nanmax(self.sel_data.vel_aic))+10)
        self.axf.set_ylim(bottom=0, top=None)
        self.figfreq_canvas.draw()
    
    def Save(self):
        #Manual
        # if "vel_manual" in self.data:
        #     return_vels = np.nanmax([self.manual_vels, np.array(self.data.vel_manual)], axis = 0) #to avoid replacing good value with -100
        #     self.data.vel_manual = return_vels # if yes : replace the data
        # else:
        #     self.data.insert(loc=len(self.data.columns), column = "vel_manual", value = self.manual_vels)
        #Max aic
        if "vel_maxaic" in self.data:
            return_vels = np.nanmax([self.maxaic_vels, np.array(self.data.vel_maxaic)], axis = 0) #to avoid replacing good value with -100
            self.data.vel_maxaic = return_vels # if yes : replace the data
        else:
            self.data.insert(loc=len(self.data.columns), column = "vel_maxaic", value = self.maxaic_vels)
        #Max aic
        if "vel_stalta_aic" in self.data:
            return_vels = np.nanmax([self.stalta_aic_vels, np.array(self.data.vel_stalta_aic)], axis = 0) #to avoid replacing good value with -100
            self.data.vel_stalta_aic = return_vels # if yes : replace the data
        else:
            self.data.insert(loc=len(self.data.columns), column = "vel_stalta_aic", value = self.stalta_aic_vels)
        self.data.to_excel(self.datafile, index=False) #f"{self.datafile[0:-5]}_manual_analysis.xlsx")


    def _quit(self):
        quit()
        #destroy()

if __name__ == "__main__":
    ctypes.windll.shcore.SetProcessDpiAwareness(1)
    app = App()
    app.mainloop()
root = tk.Tk()
#root.protocol("WM_DELETE_WINDOW", root.destroy)
my_gui = App()
root.mainloop()

