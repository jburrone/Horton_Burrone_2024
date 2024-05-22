# -*- coding: utf-8 -*-
"""
Created on Tue May 21 11:49:46 2024

@author: Vincenzo Mastrolia (vincenzo.mastrolia@kcl.ac.uk)
Research associate at the Burrone lab, Centre for Developmental Neurobiology, King's College London
"""

from scipy import signal
from scipy import ndimage
import numpy as np
import pandas as pd
import statistics
import math
from scipy.optimize import curve_fit
from sklearn.metrics import auc
import pyabf


class signal_processing:
    """
    This class contains functions useful to filter and process the signal from
    electrophysiological recordings.
    """
    
    def power_spectrum_analysis(time_seriesY,time_seriesX):
        """Performs a power spectrum analysis of the time series, using Fourier transform"""
        
        x = time_seriesY                                    # Relabel the data variable
        dt = time_seriesX[1] - time_seriesX[0]              # Define the sampling interval
        N = len(x)                                          # Define the total number of data points
        T = N * dt                                          # Define the total duration of the data
        
        xf = np.fft.fft(x - x.mean())                              # Compute Fourier transform of x
        Sxx = 2 * dt ** 2 / T * (xf * xf.conj())            # Compute spectrum
        Sxx = Sxx[:int(len(x) / 2)]
                
        df = 1 / T.max()                                    # Determine frequency resolution
        fNQ = 1 / dt / 2                                    # Determine Nyquist frequency
        faxis = np.arange(0,fNQ,df)                            # Construct frequency axis

        return(Sxx,faxis)
    
    def find_noise(Sxx,faxis):
        """Finds the strongest frequencies from the power spectrum analysis
        from 10 Hz to 500 Hz, at 10 Hz""" 
        
        max_freq = 500  # expressed in Hz
        # freq is a list of frequencies to bin through
        freqs = np.linspace(10,max_freq,num=int(max_freq/10)+1,dtype=int)
        frequencies = []
        for f in range(len(freqs)-1): # Finds the max frequencies to eliminate

            for fi in range(len(Sxx)):
                if Sxx[fi] < 0.2:
                    Sxx[fi] = 0
            
            freq = faxis[np.argmax(Sxx[int(len(faxis)/faxis.max()*freqs[f]):
                                        int(len(faxis)/faxis.max()*freqs[f+1])])+
                                         int(len(faxis)/faxis.max()*freqs[f])]
            frequencies.append(freq)
            
        return(frequencies)
    
    def bessel_filter_combo(time_series,sampling_rate,cutoff_frequency=500):
        """Applies a low-pass filter to remove any noise above the cutoff frequency"""
        
        cutoff_frequency = cutoff_frequency
        b, a = signal.bessel(1,  # Order of the filter
                              cutoff_frequency,  # Cutoff frequency
                              'low', # Type of filter
                              analog=False,  # Analog or digital filter
                              norm='phase', 
                              fs=sampling_rate)  # fs: sampling frequency

        filtered_signal = signal.lfilter(b, a, time_series)
        
        return(filtered_signal)
    
    
    def denoise(time_seriesY,time_seriesX,sampling_rate, Q = 500,graph=False,ext_frequencies = False,frequencies=[]):
        """ Uses the frequencies found by find_noise and eliminates them using a 
            notch filter"""
            
        if ext_frequencies == False: # ext_frequencies allows to use a list of 
            Sxx,faxis = signal_processing.power_spectrum_analysis(time_seriesY,time_seriesX,graph=graph)
            frequencies = signal_processing.find_noise(Sxx, faxis)
        
        else:
            frequencies = frequencies
            
        filtered_signal = time_seriesY
        for frequency in frequencies:
            c, d = signal.iirnotch(frequency, Q, sampling_rate)
            filtered_signal = signal.filtfilt(c, d,filtered_signal)
            
        return(filtered_signal,frequencies)


    def load_sweep(abf,sampling_rate,sweep,cursors,baseline_subtraction = True):
        """
        This function allows to load a specific sweep from the .abf file and 
        subtract the baseline, between cursors[0] and cursors[1].
        """

        abf.setSweep(sweep)
        
        # sweepY is the signal containing the recorded values
        sweepY = abf.sweepY[int(cursors[0]*abf.dataRate):int(cursors[1]*abf.dataRate)]
        # sweepX contains the time positions (in seconds) relative to sweepY
        sweepX = abf.sweepX[int(cursors[0]*abf.dataRate):int(cursors[1]*abf.dataRate)]
        
        if baseline_subtraction == True:
            x = ndimage.uniform_filter(sweepY, sampling_rate)
            # baseline subtraction condition applies
            sweep_sub=sweepY-x
            
            return(sweep_sub,sweepX)
        else:
            return(sweepY,sweepX)


    def peak_finder(time_seriesY,parameters,sampling_rate,std_threshold = 3):
        """Uses a rolling(browsing) window approach to screen for events with a
            higher mean than the previous and follwing browsed intervals. This 
            is a preprocessing step before the real selection, which happens in
            peak typer"""
            
        sweep_mean = ndimage.uniform_filter(time_seriesY,int(sampling_rate/100)) # Highly smoothed trace for peak finding
        sweep_filt = ndimage.uniform_filter(time_seriesY, int(sampling_rate/2000)) # the filtered trace to use to extract the parameters
        
        browsing_window = int(parameters[0]/1000*sampling_rate) # Rolling(browsing) window width
        bw_25 = int(browsing_window*0.25)   # 25% of the browsing window

        pre_peaks = []
        peaks = []

        mean_0 = np.mean(sweep_mean[0:browsing_window])
        
        for m in range(browsing_window,len(sweep_mean),browsing_window): #browsing in the sweep
            mean_next = np.mean(sweep_mean[m:m+browsing_window])
            mean_next_next = np.mean(sweep_mean[m+browsing_window:m+2*browsing_window])

            if mean_0 < mean_next > mean_next_next: # Condition to preselect the peak
                pre_peaks.append(np.argmax(sweep_filt[m:m+browsing_window]) + m)
            
            mean_0 = mean_next
            
        #################################
            
        for pre_peak in pre_peaks:
            sv = int(sampling_rate/500)     # smoothing index
            peaks_list = []
            
            ### Finding the peak using the 2nd derivative
            
            pre_peak_trace = sweep_filt[pre_peak-browsing_window:pre_peak+browsing_window]
            
            pre_peak_diff = ndimage.uniform_filter(np.diff(ndimage.uniform_filter(np.diff(pre_peak_trace),sv)),sv)      # Smooth version of the differential trace

            diff_max_list = [np.argmax(pre_peak_diff[l:l+bw_25])+l for l in range(bw_25,len(pre_peak_diff)-bw_25,bw_25)]
            diff_max_list_refine = [n for n in diff_max_list if np.argmax(pre_peak_diff[n-bw_25:n+bw_25]) + n-bw_25 == n]


            psc_start_list = [j for j in diff_max_list_refine if pre_peak_diff[j] > np.median(pre_peak_diff[j-bw_25:j+bw_25])]
            psc_peak_list = [(np.argmin(pre_peak_diff[h:h+bw_25])+h) for h in psc_start_list]


            ###  Refining the indexes found from the differential using the original values
            
            # Finds the most prominent dip as the refined starting point
            psc_start_list_refine = []
            
            for u in psc_start_list:                
                psc_start_bin = pre_peak_trace[u-int(bw_25/2):u+int(bw_25/2)]
                interpolate = np.linspace(psc_start_bin[0],psc_start_bin[-1],num=len(psc_start_bin))
                diff = interpolate - ndimage.uniform_filter(psc_start_bin, 10)                
                base = np.argmax(diff) + u-int(bw_25/2)
                psc_start_list_refine.append(base)
            
            # Finds the most prominent peak as the refined peak of the event
            psc_peak_list_refine = []
            
            for w in psc_peak_list:                
                psc_start_bin = pre_peak_trace[w-int(bw_25/2):w+int(bw_25/2)]
                new_peak = np.argmax(psc_start_bin) + w-int(bw_25/2)
                
                psc_peak_list_refine.append(new_peak)
            ###
            
            for b in range(len(psc_start_list_refine)):
                height = pre_peak_trace[psc_peak_list_refine[b]]-pre_peak_trace[psc_start_list_refine[b]]
                if height >= parameters[2]:
                    peaks_list.append([psc_start_list_refine[b]+pre_peak-browsing_window,psc_peak_list_refine[b]+pre_peak-browsing_window])

            if len(peaks_list) > 0:
                peaks.append(peaks_list)

        return(peaks)

    """A series of mathematical functions used in the analysis"""
    
    def logifunc(x,A,x0,k,off):
        return A / (1 + np.exp(-k*(x-x0)))+off
    
    def mono_decay(t, a, b):
        return a * np.exp(b * t)
    
    def bi_decay(t, a, b, c, d):
        return a * np.exp(b * t) + c * np.exp(d * t)

    def inverse_log(y,A,x0,k,off):  # x = -(log(A/(y-off)-1)/k + x0
        x = -np.log(A/(y-off)-1)/k + x0
        # check the conditions when the log becomes 
        # if either y = off or when off = 0 and A = y
        return x
    
    def inverse_exp(y,a,b):
        x = np.log(y/a)/b
        return x


    def peak_typer(time_seriesX,time_seriesY,parameters,sampling_rate, raw_trace = False, filename = None, path = None, sweep = None, cursors = None):
        """Peak typer uses the preselected peaks by peak finder and cross-checks them based on a series of given parameters"""        
        
        # Calling on to peak finder
        peaks = signal_processing.peak_finder(time_seriesY,parameters,sampling_rate)

        sweep_filt = ndimage.uniform_filter(time_seriesY, int(sampling_rate/5000)) # the filtered trace to use to extract the parameters
        browsing_window = int(parameters[0]/1000*sampling_rate)
        bw_25 = int(browsing_window*0.25)   # 25% of the browsing window
        bw_50 = int(browsing_window*0.5)    # 50% of the browsing window
        bw_75 = int(browsing_window*0.75)   # 75% of the browsing window

        counter = 0
        
        if raw_trace == True: # If true the average trace will be extracted from the unfiltered trace
            abf_raw = pyabf.ABF(f'{path}/{filename}.abf') # To load the .abf files
            abf_raw.setSweep(sweep)
            time_series_raw = abf_raw.sweepY[int(cursors[0]*abf_raw.dataRate):int(cursors[1]*abf_raw.dataRate)]
        
        
        events_list = []
        
        for peak in peaks: # Looping through each peak to extract parameters
            psc = {} # the dictionary where the parameters will be saved
            psc_start = peak[0][0]
            if psc_start > bw_50 and len(sweep_filt) - psc_start > bw_75:
                new_trace = time_seriesY[psc_start-bw_25:psc_start+bw_75]
                new_trace_x = time_seriesX[psc_start-bw_25:psc_start+bw_75]
                
                if raw_trace == True:
                    new_trace_raw = time_series_raw[psc_start-bw_25:psc_start+bw_75]
                else:
                    new_trace_raw = new_trace
                    
                psc_peak = peak[0][1] - psc_start    # Index value relative to psc_start
                peak_amp = new_trace[bw_25+psc_peak]-new_trace[bw_25]   # Amplitude of the event
                peak_slope = (peak_amp/(time_seriesX[psc_start+psc_peak]-time_seriesX[psc_start]))/1000 # Slope of the event
                peak_rise_time = (time_seriesX[psc_start + psc_peak] - time_seriesX[psc_start])*1000 # Rise time of the event
                
                decay_trace = new_trace[psc_peak+bw_25:] # Decay trace, starting from the index of the max value of the event
                peak_end = np.argmin(abs(decay_trace[parameters[0]:]), axis=0) + parameters[0] + psc_peak # Finding the end of the event
                peak_width = (new_trace_x[peak_end] - new_trace_x[bw_25])*1000 # Total duration of the event

                # selection exerted pased on the values decided in parameters [Index 1 for slope, 2 for amplitude, 3 rise time, 4 event duration ]
                if peak_slope > parameters[1] and peak_amp > parameters[2] and peak_rise_time > parameters[3] and peak_width > parameters[4] and abs(new_trace[bw_25+psc_peak]) > 0:
                    psc['trace'] = new_trace_raw
                    psc['event_index'] = psc_start
                    psc['event_time'] = time_seriesX[psc_start]
                    psc['peak'] = psc_peak
                    psc['peaks_list'] = peak
                    psc['amplitude'] = peak_amp
                    psc['rise_slope'] = peak_slope
                    psc['rise_time'] = peak_rise_time
                                    
                    ### Logistic fit of the rise phase
                    rise_trace = new_trace[bw_25:bw_25+psc_peak]

                    log_p, log_pcov = curve_fit(signal_processing.logifunc,
                                list(range(0, len(rise_trace))),
                                rise_trace, 
                                maxfev=100000,
                                p0=[np.max(rise_trace), 1, 1,np.min(rise_trace)])
                    
                    psc['rise_fit_log'] = [log_p, log_pcov]
                    
                    ###
                    
                    
                    rise_tau = signal_processing.inverse_log(peak_amp*0.63,log_p[0],log_p[1],log_p[2],log_p[3])
                    
                    # Pass from fit by index to time and then from seconds to msec
                    psc['rise_tau'] = rise_tau / sampling_rate * 1000 # 

                    psc['peak_end'] = peak_end
                    psc['width_ms'] = peak_width
                    
                    fit_decay_trace = decay_trace[0:peak_end-psc_peak]
                    
                    ### Mono-exponential decay fit ###
                    
                    a_initial = 10    # usually -200
                    b_initial = 0.01    # usually -0.1
                    popt1, pcov1,infodict1, errmsg1, ier1 = curve_fit(signal_processing.mono_decay,
                                            list(range(0, len(fit_decay_trace))), 
                                            fit_decay_trace,p0=(a_initial,b_initial),
                                            maxfev=100000,full_output=True)
                                        
                    rmsd1 = np.sqrt(np.sum(np.square(infodict1['fvec']))/len(infodict1['fvec']))
                    
                    
                    ### Bi-exponential decay fit ###
                    
                    a_initial = 10    # usually -200
                    b_initial = 0.01    # usually -0.1
                    c_initial = 10    # usually -200
                    d_initial = 0.01    # usually -0.1
                    
                    popt2, pcov2,infodict2, errmsg2, ier2 = curve_fit(signal_processing.bi_decay,
                                            list(range(0, len(fit_decay_trace))), 
                                            fit_decay_trace,p0=(a_initial,b_initial,c_initial,d_initial),
                                            maxfev=100000,full_output=True)
                                        
                    rmsd2 = np.sqrt(np.sum(np.square(infodict2['fvec']))/len(infodict1['fvec']))

                    ###

                    if rmsd1 <= rmsd2:
                        psc['decay_fit'] = [popt1,pcov1,rmsd1,'Mono']
                        decay_tau = abs(1/popt1[1])/sampling_rate*1000

                    elif rmsd1 > rmsd2:
                        psc['decay_fit'] = [popt2,pcov2,rmsd2,'Bi']
                        decay_tau1 = abs(1/popt2[1])/sampling_rate*1000
                        decay_tau2 = abs(1/popt2[3])/sampling_rate*1000
                        
                        if (decay_tau1/decay_tau2 >= 4):
                            decay_tau = [decay_tau2,decay_tau1]
                        elif (decay_tau2/decay_tau1 >= 4):
                            decay_tau = [decay_tau1,decay_tau2]
                        else:
                            decay_tau = np.mean([decay_tau1,decay_tau2])

                    psc['decay_tau'] = decay_tau
                    
                    psc['AUC'] = auc(new_trace_x[bw_25:peak_end],new_trace[bw_25:peak_end])
                                            
                    if counter == 0:
                        psc['ISI'] = np.nan
                        psc['inst_freq'] = np.nan
                        counter += 1
                        isi_prev = time_seriesX[psc_start]
                    else:
                        isi = time_seriesX[psc_start] - isi_prev
                        psc['ISI'] = isi

                        if isi != 0:
                            psc['inst_freq'] = 1/isi
                        else:
                            psc['inst_freq'] = np.nan


                        isi_prev = time_seriesX[psc_start]
                        counter += 1
                        
                        events_list.append(psc)

        return(events_list)


    def sPSC_analysis(path,save_path,filename,abf,sampling_rate,n_sweeps,cursors,direction, raw_trace = False):
        
        ### Parameters for spontaneous events ######################################
        """ You can change these parameters based on the type of events """
        
        parameters_IPSC = [20,      # time window in ms
                            2,      # threshold for slope
                            15,     # threshold for amp
                            3,     # threshold for peak rise time in ms, ideal is 10 ms
                            5,      # Minimal event duration
                            500,     # Max decay tau
                            ]  
        
        # ## Parameters for EPSCs ###
        parameters_EPSC = [20,      # time window in ms
                            5,      # threshold for slope
                            7,      # threshold for amp
                            2,     # threshold for peak rise time in ms, ideal is 10 ms
                            1,      # Minimal event duration
                            60,     # Max decay tau
                            ]  
        
        ############################################################################
      
        
        columns = ['filename',
                   'type',
                   'sweep',
                   'event_index',
                   'event_time',
                    'peak',
                    'peaks_list',
                    'peak_end',
                    'width_ms',
                    'amplitude',
                    'rise_slope',
                    'rise_time',
                    'rise_tau',
                    'rise_fit_log',
                    'decay_tau',
                    'decay_fit',
                    'ISI',
                    'inst_freq',
                    'AUC',
                  ]


        database = pd.DataFrame(columns=columns)

        ############################################################################

        sweepY_list = []
        peaks_traces = []
        freq_list = []
        noise_frequencies_list = []
        
        # Deciding which parameters to use based on the type of event, whether
        # excitatory or inhibitory
        
        for sweep in range(n_sweeps):
            sweepY,sweepX = signal_processing.load_sweep(abf,sampling_rate,sweep,cursors)
            if direction == 'negative':
                sweepY = -sweepY
                parameters = parameters_EPSC
                psc = 'EPSC'
            else:
                parameters = parameters_IPSC
                psc = 'IPSC'
        
            # Filtering the signal, low-pass + notch filters
            sweepY_filt = signal_processing.bessel_filter_combo(sweepY, sampling_rate,cutoff_frequency=500)
            sweepY_filt,noise_frequencies = signal_processing.denoise(sweepY_filt, sweepX, sampling_rate)
            sweepY_list.append(sweepY_filt)
            ###
            
            # Running peak typer to extract the events
            peaks = signal_processing.peak_typer(sweepX, sweepY_filt, parameters, sampling_rate, raw_trace == raw_trace, filename = filename, path = path, sweep = sweep, cursors = cursors)
            noise_frequencies_list.append(noise_frequencies)
            
            for peak in peaks:
                peak_list = [filename,
                             psc,
                             sweep,
                             peak['event_index'],
                             peak['event_time'],
                             peak['peak'],
                             peak['peaks_list'],
                             peak['peak_end'],
                             peak['width_ms'],
                             peak['amplitude'],
                             peak['rise_slope'],
                             peak['rise_time'],
                             peak['rise_tau'],
                             peak['rise_fit_log'],
                             peak['decay_tau'],
                             peak['decay_fit'],
                             peak['ISI'],
                             peak['inst_freq'],
                             peak['AUC']]
                    
                peaks_traces.append(peak['trace'])
                
                database.loc[len(database)] = peak_list
            freq_list.append(len(peaks)/(cursors[1]-cursors[0]))
            
        database.to_excel(save_path + '/' + filename + '_' + "output.xlsx")
        # noise_frequencies_list = pd.DataFrame(noise_frequencies_list)
        # noise_frequencies_list.to_csv(save_path + 'noise_frequencies_list.csv')
        
        n_events = int(len(database))
        avg_amplitude =     database['amplitude'].mean()
        avg_rise_slope =    database['rise_slope'].mean()
        avg_rise_time =     database['rise_time'].mean()
        avg_rise_tau =      database['rise_tau'].mean()
        avg_width =         database['width_ms'].mean()
        avg_decay_tau =     database.decay_tau.loc[[type(database.loc[v,'decay_tau']) == float for v in range(len(database))]].median()
        
        bi_decay_tau =     database.decay_tau.loc[[type(database.loc[v,'decay_tau']) == list for v in range(len(database))]].reset_index(drop=True)
        avg_fast_decay_tau = np.median([bi_decay_tau[l][0] for l in range(len(bi_decay_tau))])
        avg_slow_decay_tau = np.median([bi_decay_tau[l][1] for l in range(len(bi_decay_tau))])
        
        avg_ISI =           database['ISI'].mean()
        avg_inst_freq =     database['inst_freq'].mean()
        avg_auc =           database['AUC'].mean()
        total_auc =         database['AUC'].sum()
        avg_freq =          np.mean(freq_list)
        
        avg_trace_name = filename + '_' + 'Mean_event_trace.npy'
        avg_trace =         np.mean(peaks_traces,axis=0)
        np.save(save_path + '/' + avg_trace_name,avg_trace)

        
        H,X1 = np.histogram(database.amplitude,bins=50,density=True,range=(0,150))
        dx = X1[1] - X1[0]
        cumulative_amplitude = np.cumsum(H)*dx
        cumulative_amplitude_name = filename + '_' + 'cumulative_amplitude.npy'
        np.save(save_path + '/' + cumulative_amplitude_name,cumulative_amplitude)

        # avg_trace_path = save_path + '/' + filename + '_' + 'mean_event_trace.csv'
        # avg_trace.to_file(avg_trace_path,sep=',',)
        
        
        file_df = pd.DataFrame(columns=['n_events',
                                         'amplitude',
                                         'rise_slope',
                                         'rise_time',
                                         'rise_tau',
                                         'width_ms',
                                         'decay_tau',
                                         'decay_tau_fast',
                                         'decay_tau_slow',
                                         'ISI',
                                         'inst_freq',
                                         'mean_AUC',
                                         'total_AUC',
                                         'mean_freq',
                                         'mean_event_trace',
                                         'cumulative_amplitude'
                                         ]
                               )
        
        
        data = [n_events,
                avg_amplitude,
                avg_rise_slope,
                avg_rise_time,
                avg_rise_tau,
                avg_width,
                avg_decay_tau,
                avg_fast_decay_tau,
                avg_slow_decay_tau,
                avg_ISI,
                avg_inst_freq,
                avg_auc,
                total_auc,
                avg_freq,
                avg_trace_name,
                cumulative_amplitude_name,
                ]

        file_df.loc[0] = data    
        
        return(database,file_df)

###############################################################################

class evoked_stimulus_analysis:
    """
    This class contains functions useful to evoked currents from
    electrophysiological recordings.
    """

    def evoked_currents_analysis(time_seriesY,sampling_rate,direction,frequency,n_pulses,cursor,art):
        """
        This function allows to analyse individual evoked currents from a train
        of stimulating pulses.
        """
        
        max_values_pulse = np.empty(n_pulses)           # Creating the array onto which append the peak amplitudes
        tails_pulse = np.empty(n_pulses)                # Creating the array onto which append the tails values
        peak_indexes = np.empty(n_pulses,dtype=int)     # Creating the array onto which append the peak indexes
        decay_tau_pulse = np.empty(n_pulses,dtype=int)     # Creating the array onto which append the decay taus
        
        ### The section of the time series where the stimulation occurs --> Trace sweep
        trace_sweep = ndimage.uniform_filter(time_seriesY[int(cursor*sampling_rate):int(cursor*sampling_rate)+int(1/frequency*n_pulses*sampling_rate)],5)
       
        ### The length of each pulse bin, based on the frequency and the sampling rate
        bin_length = int(1/frequency*sampling_rate)
        
        ## Loops through the sweep, pulse by pulse (bin_int)
        for bin_int in range(0,bin_length*n_pulses,bin_length): 
           
            bin_trace = trace_sweep[bin_int+int(art*sampling_rate):bin_int+bin_length]
            bin_trace_smooth = ndimage.uniform_filter(trace_sweep[bin_int+int(art*sampling_rate):bin_int+bin_length],10)

            # Find the peak
            if direction == 'negative':
                peak_index = np.argmin(bin_trace_smooth[0:int(bin_length/4)])+int(art*sampling_rate)#+pre_max_index-10#+pre_max_index+bin_int+int(art*sampling_rate)
            else:
                peak_index = np.argmax(bin_trace_smooth[0:int(bin_length/4)])+int(art*sampling_rate)#+pre_max_index-10#+pre_max_index+bin_int+int(art*sampling_rate)

            peak_indexes[int(bin_int/bin_length)] = peak_index
            
            amplitude = abs(trace_sweep[peak_index])
            max_values_pulse[int(bin_int/bin_length)] = amplitude
            
            # Find the value at the end of the pulse (tail)
            tail_trace = bin_trace[int(bin_length*0.75):bin_length]
            tail = statistics.median(tail_trace)
            tails_pulse[int(bin_int/bin_length)] = tail

            # Decay fit of the evoked current
            fit_decay_trace = abs(bin_trace_smooth[peak_index:len(bin_trace)-1])
            a_initial = 10    # usually -200
            b_initial = 0.01    # usually -0.1
            popt, pcov = curve_fit(lambda t, a, b: a * np.exp(b * t),
                                    list(range(0, len(fit_decay_trace))), 
                                    fit_decay_trace,p0=(a_initial,b_initial),
                                    maxfev=100000)
            
            decay_tau = abs(1/popt[1])/sampling_rate*1000
            decay_tau_pulse[int(bin_int/bin_length)] = decay_tau

            
        return(trace_sweep,max_values_pulse,peak_indexes,tails_pulse,decay_tau_pulse)
