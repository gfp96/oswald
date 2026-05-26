# import filedialog module
from tkinter import filedialog
import tkinter as tk
from tkinter import ttk
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
from datetime import datetime
import pandas.errors  # for errors
from signal_interp import Find_start

def browseFiles():
    filename = filedialog.askopenfilename(initialdir = "/",
                                          title = "Select a File",
                                          filetypes = (("Text files",
                                                        "*.txt*"),
                                                       ("all files",
                                                        "*.*")))
      
    # Change label contents
    ttk.label_file_explorer.configure(text="File Opened: "+filename)



def browseFolders():
    folder = filedialog.askdirectory(initialdir = "C:\\Users\\g.flood-page\\Documents\\Professionnel\\PhD\\8_Experimental\\Exp_results", title = "Select a Folder")
    datafile = filedialog.askopenfilename(initialdir = folder, title = "Select file containing experiment data",filetypes = (("Excel files", "*.xlsx*"), ("all files","*.*")))
    data = pd.read_excel(datafile)
    data = data[data.valid==True]
    pe_list = np.unique(data.pelev)
    pf_list = np.unique(data.pflev)

def Set_graph_params():
    graph_font = {'family': 'DejaVu Sans',
                  'weight': 'normal',
                  'size': 18}
    plt.rcParams['figure.figsize'] = [16, 10]
    plt.rc('font', **graph_font)  # set font parameters
    plt.rcParams['lines.linewidth'] = 3  # set width of lines
    plt.rcParams['figure.dpi'] = 100  # 200 e.g. is really fine, but slower
    plt.rcParams['figure.titlesize'] = 22
    plt.rcParams['axes.titlesize'] = 20
    plt.rcParams['axes.labelsize'] = 18
    plt.rcParams['legend.fontsize'] = 18
    plt.rcParams['xtick.labelsize'] = 18
    plt.rcParams['ytick.labelsize'] = 18
    #plt.rcParams["figure.autolayout"] = True
    #plt.rcParams['figure.constrained_layout.use'] = True
    plt.rcParams['lines.markersize'] = 10

    custom_params = {
        'figure.figsize': [16, 10],
        #'font': graph_font,  # set font parameters
        'lines.linewidth':0.5,  # set width of lines
        'figure.dpi':100,  # 200 e.g. is really fine, but slower
        'figure.titlesize':14,
        'axes.titlesize':12,
        'axes.labelsize':12,
        'legend.fontsize':12,
        'xtick.labelsize':12,
        'ytick.labelsize':12,
        'lines.markersize':14,
        "axes.xmargin" : 0,
        "axes.ymargin" : 0
    }
    sns.set_theme(rc=custom_params)
    sns.set_style("dark")

repE = lambda x: (x.replace('E+0','')) #replaces bugged E+0
def Read_BE_file(folder_path, filename):
    """This function reads a Bender Element file and returns measurement time, dt, input signal array, output signal array"""
    stime = datetime.now()
    try:
        read = pd.read_csv(os.path.join(folder_path, filename), sep="\t", header=3, index_col=False,
                                    skiprows=0, usecols = [1, 3], converters = {"Y[0]": repE,"Y[1]": repE},  #dtype={"Y[0]": np.float,"Y[1]": np.float}
                                    encoding_errors = 'ignore', decimal = ',', engine = "python", on_bad_lines="skip") #decimal = ',',
    except pandas.errors.EmptyDataError:
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan  # time, dt, input time, input signal, output time, output signal

    if read.size==0: #in case of bugged file
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan  # time, dt, input time, input signal, output time, output signal

    input_signal = pd.to_numeric(read.iloc[:, 0].str.replace(',', '.'), errors = "coerce")
    output_signal = pd.to_numeric(read.iloc[:, 1].str.replace(',', '.'), errors = "coerce")
    top = pd.read_csv(os.path.join(folder_path, filename), sep="\t", header=None, index_col=False,
                                    dtype=str, usecols=[1], skiprows=1, nrows=2)  # time increment in between 2 data points
    dt = float(top.iloc[1].item().replace(',', '.'))
    time = np.arange(input_signal.size)*dt
    print(datetime.now()-stime)
    return time, input_signal, output_signal


def Read_terratek_signal(file, dayFirst_be_files=False):
    """This function reads a wave signal file and returns measurement measurement datetime, time resolution, time, signal arrays"""

    with open(file) as f:
        try:
            read = pd.read_csv(file, sep="\t", skiprows=2)#, header=2)#, index_col=False, dayfirst=dayFirst_be_files,
                                        #dtype={"time[0]": str, "Y[0]": str, "time[1]": str, "T[1]": str}, skiprows=2))
        except pandas.errors.EmptyDataError:
            return np.nan, np.nan  # time, signal

        if read.size==0: #in case of bugged file
            return np.nan, np.nan  # measurement time, dt, time, signal

        times = (np.array(read.iloc[:,0], dtype = float)-read.iloc[0,0])*1e-6 #pass from micros to s
        #dt = times[1]-times[0]
        signal = np.array(read.iloc[:,1], dtype = float)
        # date = pd.read_csv(os.path.join(folder_path, filename), sep="\t", header=None, index_col=False, dtype=str, usecols=[1], skiprows=1, nrows=1).iloc[0].item()  # time increment in between 2 data points
        # time = pd.read_csv(os.path.join(folder_path, filename), sep="\t", header=None, index_col=False, dtype=str, usecols=[0], skiprows=1, nrows=1).iloc[0].item()  # time increment in between 2 data points
        # ini_time = date+" "+time
        # if dayFirst_be_files: #code dates normally
        #     ini_time = datetime.strptime(ini_time, "%d/%m/%y %H:%M:%S")
        # else: #weird american date encoding
        #     ini_time = datetime.strptime(ini_time, "%m/%d/%y %H:%M:%S")
        return times, signal




def Update_df(df, colname, values, position = -1):
    """Checks if a column already exists in df before inserting new column - column is either the last one (-1) or given by position"""
    if colname in df:
        df.loc[:, colname]=values
    else:
        if position==-1:
            df.insert(loc=len(df.columns), column = colname, value = values)
        else:
            df.insert(loc=position, column = colname, value = values)

