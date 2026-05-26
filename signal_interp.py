from tracemalloc import start
from scipy.signal import butter, sosfiltfilt, find_peaks
import numpy as np
import pandas as pd


#---------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------AIC functions-----------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------

def AIC(k, N, input_signal):  # compares variance right and left of index k :  signal "input_signal" and length of signal= N
    return (k+1)*np.log(np.nanvar(input_signal[0:k]))+(N-k-2)*np.log(np.nanvar(input_signal[k+1:N-1]))

# def AIC_matrix(k, N, input_signal):  # compares variance right and left of index k :  signal "input_signal" and length of signal= N
#     return (k+1)*np.log(np.nanvar(input_signal[:, 0:k], axis=1))+(N-k-2)*np.log(np.nanvar(input_signal[:, k+1:N-1], axis = 1))

#---------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------SNR functions-----------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------

# def Get_snr_full_signal(signal, index_start):
#     """Returns noise, signal and snr considering noise before emission and signal after emission - for single signal"""
#     #Compares RMS amplitudes of start of signal with useful section
#     noise2 = np.mean(signal[0:index_start]**2)
#     sig2 = np.mean(signal[index_start:-1]**2)
#     return noise2, sig2, 20*np.log10(sig2/noise2)

# def Get_snr_full_signal_multiple(signals, starts):
#     """Returns noise, signal and snr considering noise before emission and signal after emission - for array containing multuple signals"""
#     #Compares RMS amplitudes of start of signal with useful section
#     noise = np.zeros(starts.size)
#     sig = np.zeros_like(noise)
#     snrs = np.zeros_like(noise)
#     for i, start in enumerate(starts):
#         noise[i] = np.nanmean(signals[0:start, i]**2)
#         sig[i] = np.nanmean(signals[start:-1, i]**2)
#     return noise, sig, 20*np.log10(sig/noise)

#---------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------Find start functions----------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------

def Find_start(signal):  # finds positive and negative peaks of input and deduces the start of the wave
    """Find start of input BE signal"""
    ind_min = np.argmin(signal, axis=0) 
    ind_max = np.argmax(signal, axis=0) 
    inis = min(ind_min, ind_max)
    ends = max(ind_min, ind_max)
    index = inis-(ends-inis)//2  # averaged to higher int
    index = max(index, 0)
    return index


def Find_start_burst(input_signal,nbperiods):  # finds positive and negative peaks of input and deduces the start of the wave
    """
    Find start of input BE signal for burst-type signals\n
    - nbperiods: integer number of periods in the used burst
    """
    if nbperiods//2==0: #even number : even code is as of yet untested
        ind_min = np.argmin(input_signal)
        ind_max = np.argmax(input_signal)
        ind_center = int(np.mean([ind_min, ind_max]))
        ind_start = max(0, ind_center-int(nbperiods/2)*2*np.abs(ind_max-ind_min)) #2*abs(max-min) is a period
    else: #uneven number
        ind_center = np.argmax((np.abs(input_signal))) #max of the abs has to be center of signal, no matter if positive or negative peaks
        ind_min = np.argmin(input_signal)
        ind_max = np.argmax(input_signal)
        ind_start=max(0, ind_center-int(nbperiods/2)*2*np.abs(ind_max-ind_min)) #2*abs(max-min) is a period
    return ind_start

#---------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------GFP Max-AIC method for velocity calculation-------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------

def Get_Max_AIC_velocity(times, s_inds, rough, propagation_distance, input_freqs,
    freq_range = np.array([100, 200000]),
    input_amplitude = [10, 30], #min and max acceptable amplitude of input signal [V]
    butter_order = 1,
    Bad_grounding = False #If bad grounding, the input signal is visible in output signal : need to delay start of AIC to end of Input signal
    ):

    sel_size = s_inds.size
    indexes = range(sel_size)
    filtered = [None] * sel_size #for filtered data
    aics = [None] * sel_size
    start_aics = np.zeros(sel_size, dtype = int)
    max_inds = np.zeros_like(start_aics)
    valids = np.ones(sel_size, dtype = bool)
    a_inds = np.zeros(sel_size, dtype = int)
    arrivals = np.zeros(sel_size, dtype = float)
    starts = np.zeros(sel_size, dtype = float)
    vels = np.zeros(sel_size, dtype = float)


    for i in indexes:
        dt = times[i][1]-times[i][0] #timestep supposed identical
        # Set Butterworth filter
        Nyquist = 0.5*1/dt  # Nyquist freq = half the sampling frequency
        butter_order = int(butter_order) #in case user does not provide int
        sos = butter(butter_order, [freq_range[0]/Nyquist, freq_range[1]/Nyquist], btype="bandpass",  output = "sos")
        filtered[i] = sosfiltfilt(sos, rough[i])
        #adjust signal for median to be null
        median_sig = np.nanmedian(filtered[i])
        signal = filtered[i]-median_sig
        filtered[i] = signal

        #Choose start point for aic according to quality of grounding and signal : if electric current 
        #of input is significant on receiving BE signal, then choose Bad_grounding
        if Bad_grounding:
            start_aics[i] = s_inds[i]+np.array((1/input_freqs.iloc[i]/dt), dtype = int) #start after end of input
        else:
            start_aics[i] = s_inds[i]#+int(1/input_freq/dt) #start after start of input

        #Get index for end of aic
        max_inds[i] = np.argmax(np.abs(signal))
        N = max_inds[i]-start_aics[i]
        aici = np.zeros(N+1, dtype = float)
        #Calculate AIC and velocity
        if N>0 : #valid signal
            for k in np.arange(N+1, dtype = int): #this loop is compulsory as k cannot be an array
                aici[k] = AIC(k, N, signal[start_aics[i]:max_inds[i]])
            #Arrival 
            a_inds[i] = start_aics[i]+np.argmin(np.ma.masked_invalid(aici)) 
            aics[i] = aici #store in bigger array
        else:
            a_inds[i] = start_aics[i]+N #bugged signal
            valids[i]=False

        #Arrival time and noise
        arrivals[i] = times[i][a_inds[i]] #arrival time
        starts[i] = times[i][s_inds[i]] #start time

        #noises, sig_amps, SNRs = Get_snr_full_signal_multiple(signals, s_inds) 
        vels[i] =  propagation_distance/(arrivals[i]-starts[i])
        #min_snr = np.mean(SNRs)-np.std(SNRs)
        # for i, snr in enumerate(SNRs): 
        #     if snr<min_snr:
        #         valids[i] = False
    
    # signals = pd.DataFrame({
    #     "aic":AICs,
    #     "signal":coefficients
    # })
    data = pd.DataFrame({
        "s_ind":s_inds,
        "a_ind":a_inds,
        "valid":valids,
        "s_time":starts,
        "a_time":arrivals,
        "start_aic":start_aics,
        "vel":vels,
        "max_ind":max_inds
    })
    return data, aics, filtered

#---------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------GFP Max-AIC method for velocity calculation-------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------

def Get_STALTA_AIC_velocity(times, s_inds, rough, propagation_distance, input_freqs,
    freq_range = np.array([100, 50e3]),
    input_amplitude = [10, 30], #min and max acceptable amplitude of input signal [V]
    butter_order = 1,
    max_vel = 2000, #max possible velocity : anything above is considered to be a bug
    min_frange = 30e3, #minimum range for filtering
    Bad_grounding = False #If bad grounding, the input signal is visible in output signal : need to delay start of AIC to end of Input signal
    ):

    sel_size = s_inds.size
    indexes = range(sel_size)
    filtered = [None] * sel_size #for filtered data
    aic1s = [None] * sel_size #first round of aic
    aic2s = [None] * sel_size #second round of aic
    saic1_inds = np.zeros(sel_size, dtype = int)
    saic2_inds = np.zeros(sel_size, dtype = int)
    eaic1_inds = np.zeros_like(saic1_inds)
    eaic2_inds= np.zeros_like(saic1_inds)
    valids = np.ones(sel_size, dtype = bool)
    a1_inds = np.zeros(sel_size, dtype = int)
    a2_inds = np.zeros(sel_size, dtype = int)
    arrivals1 = np.zeros(sel_size, dtype = float)
    arrivals2 = np.zeros(sel_size, dtype = float)
    starts = np.zeros(sel_size, dtype = float)
    fpeaks = np.zeros(sel_size, dtype = float)
    vels1 = np.zeros(sel_size, dtype = float)
    vels2 = np.zeros(sel_size, dtype = float)


    for i in indexes:
        dt = times[i][1]-times[i][0] #timestep supposed identical

        #FFT to find peak freaquency
        factor = 1 #for spacing of frequencies in result : *2 means twice as many frequencies calculated
        fourier = np.fft.rfft(rough[i], factor*rough[i].size)
        freqs = np.fft.rfftfreq(rough[i].size, dt*factor)
        fft_inds = np.where(np.all([freqs<freq_range[1], freqs>freq_range[0]], axis = 0))
        freqs = freqs[fft_inds]
        fourier = np.abs(fourier[fft_inds])
        freq_ind = np.argmax(fourier)
        fpeaks[i] = freqs[freq_ind]

        # Set Butterworth filter
        Nyquist = 0.5*1/dt  # Nyquist freq = half the sampling frequency
        butter_order = int(butter_order) #in case user does not provide int
        sos = butter(butter_order, [max(freq_range[0], input_freqs.iloc[i]-min_frange)/Nyquist, min(freq_range[1], input_freqs.iloc[i]+min_frange)/Nyquist], btype="bandpass",  output = "sos")
        filtered[i] = sosfiltfilt(sos, rough[i])
        #adjust signal for median to be null
        signal = filtered[i]-np.nanmedian(filtered[i])
        filtered[i] = signal

        #STALTA parameters
        Nl = times[i].size
        Ns = int(2/fpeaks[i]/dt) #1
        Ns2 = Ns//2 #Ns2 = Ns//2 #
        
        #Characteristic function
        char_func = signal**2
        stalta = np.zeros_like(times[i])
        lta = np.mean(char_func)

        for k, cf in enumerate(char_func[Ns2:Nl-Ns2]):
            stalta[k+Ns2] = np.mean(char_func[k:k+2*Ns2+1])#np.mean(char_func[k-Ns2:k+Ns2+1])
        stalta/=lta

        ind_max_vel = int(s_inds[i] + propagation_distance/max_vel/dt)

        #Determine where the signal is
        mean_stalta = np.mean(stalta)
        is_signal = np.where(stalta[ind_max_vel:]>1)[0]+ind_max_vel
        # is_signal = np.where(stalta[ind_max_vel:]>2*mean_stalta)[0]+ind_max_vel # to avoid taking signal before emission or during cross-talk
        # if is_signal.size==0:
        #     is_signal = np.where(stalta[ind_max_vel:]>mean_stalta)[0]+ind_max_vel # to avoid taking signal before emission or during cross-talk
        #     if is_signal.size==0:
        #         is_signal = np.where(stalta[ind_max_vel:]>0.5*mean_stalta)[0]+ind_max_vel
        #         if is_signal.size==0:
        #             continue #problem with signal : amplitude is never large
        if is_signal.size<1/input_freqs.iloc[i]/dt: #Bug in interpretation : if is_signal smaller than input period, must be cross talk
            #valid = False
            continue #bugged
        
        eaic1_inds[i] = is_signal[0]
        eaic2_inds[i] = np.nanmin(np.where((stalta[eaic1_inds[i]:]<=1)))+eaic1_inds[i]
        #end_aic = times[eaic_ind]

        #Choose start point for aic according to quality of grounding and signal : if electric current 
        #of input is significant on receiving BE signal, then choose Bad_grounding
        if Bad_grounding:
            saic1_inds[i] = s_inds[i]+np.array((1/input_freqs.iloc[i]/dt), dtype = int) #start after end of input
        else:
            saic1_inds[i] = s_inds[i]#+int(1/input_freq/dt) #start after start of input

        #Get index for end of aic
        N = eaic1_inds[i]-saic1_inds[i]
        aic1 = np.zeros(N, dtype = float)
        #Calculate AIC and velocity
        if N>0 : #valid signal
            for k in np.arange(N, dtype = int): #this loop is compulsory as k cannot be an array
                aic1[k] = AIC(k, N, signal[saic1_inds[i]:eaic1_inds[i]])
            #Arrival 
            a1_inds[i] = saic1_inds[i]+np.argmin(np.ma.masked_invalid(aic1)) 
        else:
            a1_inds[i] = saic1_inds[i]+N #bugged signal
            valid=False
        aic1s[i] = aic1
        #Get index for end of aic : second aic loop
        saic2_inds[i] = a1_inds[i]-int(1/fpeaks[i]/dt)
        N = eaic2_inds[i]-saic2_inds[i] #from 1 period before arrival to end of stalta
        aic2 = np.zeros(N, dtype = float)
        #Calculate AIC and velocity
        if N>0 : #valid signal
            for k in np.arange(N, dtype = int): #this loop is compulsory as k cannot be an array
                aic2[k] = AIC(k, N, signal[saic2_inds[i]:eaic2_inds[i]])
            #Arrival 
            a2_inds[i] = saic2_inds[i]+np.argmin(np.ma.masked_invalid(aic2)) 
        else:
            a2_inds[i] = saic2_inds[i]+N #bugged signal
            valid=False
        aic2s[i] = aic2 #store aic2

        arrivals1[i] = times[i][a1_inds[i]]
        arrivals2[i] = times[i][a2_inds[i]]
        starts[i] = times[i][s_inds[i]] #start time
        vels1[i] =  (propagation_distance/(arrivals1[i]-starts[i]))
        vels2[i] =  (propagation_distance/(arrivals2[i]-starts[i]))

    data = pd.DataFrame({
        "s_ind":s_inds,
        "a1_ind":a1_inds,
        "a2_ind":a2_inds,
        "valid":valids,
        "s_time":starts,
        "a1_time":arrivals1,
        "a2_time":arrivals2,
        "eaic1_ind":eaic1_inds,
        "eaic2_ind":eaic2_inds,
        "saic1_ind":saic1_inds,
        "saic2_ind":saic2_inds,
        "fpeak":fpeaks,
        "vel1":vels1,
        "vel2":vels2
    })
    return data, aic1s, aic2s, filtered

#---------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------Just filtering signal------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------

def Get_filtered_signal(times, s_inds, rough, isvp, input_freqs,
    freq_rangeP = np.array([100, 200000]),
    freq_rangeS = np.array([100, 50000]),
    min_frange = 50e3, 
    butter_order = 1,
    ):

    sel_size = s_inds.size
    indexes = range(sel_size)
    filtered = [None] * sel_size #for filtered data

    for i in indexes:
        dt = times[i][1]-times[i][0] #timestep supposed identical
        # Set Butterworth filter
        if isvp.iloc[i]:
            freq_range=freq_rangeP
        else:
            freq_range=freq_rangeS
        Nyquist = 0.5*1/dt  # Nyquist freq = half the sampling frequency
        butter_order = int(butter_order) #in case user does not provide int
        sos = butter(butter_order, [max(freq_range[0], input_freqs.iloc[i]-min_frange)/Nyquist, min(freq_range[1], input_freqs.iloc[i]+min_frange)/Nyquist], btype="bandpass",  output = "sos")
        filtered[i] = sosfiltfilt(sos, rough[i])
        #adjust signal for median to be null
        filtered[i] = filtered[i]-np.nanmedian(filtered[i])

    return filtered