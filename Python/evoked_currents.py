# -*- coding: utf-8 -*-
"""
Created on Mon May 20 18:12:05 2024

@author: Vincenzo Mastrolia (vincenzo.mastrolia@kcl.ac.uk)
Research associate at the Burrone lab, Centre for Developmental Neurobiology, King's College London

This Script allows to extract parameters from evoked currents from patch clamp recordings in neurons
"""

import pyabf
import numpy as np
import pandas as pd
from signal_processing_tools import signal_processing as sp
from signal_processing_tools import evoked_stimulus_analysis as es

###############################################################################

path = 'please/insert/your/path/here' # Path where file is loc
file = 'insert your .abf file here'
save_path = 'please/insert/your/path/here' # Path where files are saved

frequency = 20 # The frequency in Hz used for the evoked currents
n_pulses = 2 # the number of stimulating pulses used
cursor = 10 # Position of the first oulse within the sweep, in seconds
artifact = [0.004,0.003] # The artifact length in seconds, usually 0.005
direction = 'negative' # The direction of the evoked current, whether inward or
#                       outward. Usually eEPSC are inward, while eIPSC are 
#                       outward.
###############################################################################

# Creating the database onto which save the values

columns=['Date',
        'File',
        'Age',
        'Cell',
        'Frequency',
        'n_pulses',
        'cursor',
        'Direction',
        'Type',
        'n_sweeps',
        'trace',
        'amplitudes',
        'avg_amplitudes',
        'PPR',
        ]

evoked_db = pd.DataFrame(columns=columns)

evoked_db['File'] = file
evoked_db['Date'] = [file.split('_')[0] + '_' + file.split('_')[1] + '_' + file.split('_')[2]]
evoked_db['Frequency'] = frequency
evoked_db['n_pulses'] = n_pulses
evoked_db['cursor'] = cursor
evoked_db['Direction'] = direction

evoked_db.set_index('File',inplace=True)

###############################################################################

abf = pyabf.ABF(path + "/" + file + ".abf") # To load the .abf files
sampling_rate = abf.dataRate # acquisition rate from the .abf file
n_sweeps = abf.sweepCount # number of sweeps in the recording
evoked_db['n_sweeps'] = n_sweeps

# Setting up lists and arrays where to save the values
stim_traces = np.empty([int(1/frequency*n_pulses*sampling_rate),n_sweeps])
sweepY_list = []
max_values = np.empty([n_pulses,n_sweeps])
tails_values = np.empty([n_pulses,n_sweeps])
decay_tau_values = np.empty([n_pulses,n_sweeps])
###

# Selection of the artifact's duration based on the type of current
if direction == 'negative':
    art = artifact[1]
elif direction == 'positive':
    art = artifact[0]
###

for sweep in range(n_sweeps):   # Looping through each sweep
    cursors = [0,abf.sweepLengthSec]
    
    # Loading, baseline subtraction and low-pass filtering of the sweep's signal
    sweepY,sweepX = sp.load_sweep(abf,sampling_rate,sweep,cursors)            
    sweepY_filt = sp.bessel_filter_combo(sweepY, sampling_rate)
    sweepY_list.append(sweepY_filt)
    
    # Extracting the parameters for each evoked response
    trace_sweep,max_values_pulse,peak_indexes,tails_pulse,decay_tau_pulse = es.evoked_currents_analysis(sweepY_filt, sampling_rate, direction, frequency, n_pulses, cursor, art)
    
    for pulse in range(0,n_pulses): # Substracting the relative baseline from the previous event, to avoid cumulative amplitude
        if pulse > 0:
            max_values[pulse,sweep] = abs(max_values_pulse[pulse])-abs(tails_pulse[pulse-1]) # subtracting the tail
        else:
            max_values[pulse,sweep] = abs(max_values_pulse[pulse])
            
    max_values[:,sweep] = max_values_pulse    
    stim_traces[:,sweep] = trace_sweep
    decay_tau_values[:,sweep] = decay_tau_pulse

# Extracting amplitudes from the first pulse
max_values_avg_1st = np.array([np.nanmedian(max_values[0,o:o+5]) for o in range(0,max_values.shape[1],5)])
max_values_std_1st = np.array([np.nanstd(max_values[0,o:o+5]) for o in range(0,max_values.shape[1],5)])
max_values_max_norm_avg_1st = np.array([np.nanmean(max_values[0,o:o+5])/np.max(max_values[0,o:o+5]) for o in range(0,max_values.shape[1],5)])
max_values_norm_avg_1st = np.array([np.nanmean(max_values[0,o:o+5])/np.nanmean(max_values[0,0:5]) for o in range(0,max_values.shape[1],5)])
decay_tau_values_avg_1st = np.array([np.nanmedian(decay_tau_values[0,o:o+5]) for o in range(0,decay_tau_values.shape[1],5)])

# Paired-pulse ratio calculaton
PPR = np.array([np.nanmean(max_values[1,o:o+5])/np.nanmean(max_values[0,o:o+5]) for o in range(0,max_values.shape[1],5)])    
# PPR = np.array([(max_values[1,sweep]/max_values[0,sweep]) for sweep in range(max_values.shape[1])])
PPR_name = file + '_' + 'PPR.npy'
np.save(f'{save_path}/Single_files/{PPR_name}',PPR)
evoked_db.at[file,'PPR'] = PPR_name
####

# Saving values as single numpy arrays and into the database evoked_db
stim_traces_name = file + '_' + 'stim_traces.npy'
np.save(f'{save_path}/Single_files/{stim_traces_name}',stim_traces)
evoked_db.at[file,'trace'] = stim_traces_name

max_values_name = file + '_' + 'max_values.npy'
np.save(f'{save_path}/Single_files/{max_values_name}',max_values)
evoked_db.at[file,'amplitudes'] = max_values_name

max_values_avg_1st_name = file + '_' + 'max_values_avg_1st.npy'
np.save(f'{save_path}/Single_files/{max_values_avg_1st_name}', max_values_avg_1st)
evoked_db.at[file,'avg_amplitudes'] = max_values_avg_1st_name

max_values_norm_avg_1st_name = file + '_' + 'max_values_norm_avg_1st.npy'
np.save(f'{save_path}/Single_files/{max_values_norm_avg_1st_name}', max_values_norm_avg_1st)
evoked_db.at[file,'norm_avg_amplitudes'] = max_values_norm_avg_1st_name

max_values_max_norm_avg_1st_name = file + '_' + 'max_values_max_norm_avg_1st.npy'
np.save(f'{save_path}/Single_files/{max_values_max_norm_avg_1st_name}', max_values_max_norm_avg_1st)
evoked_db.at[file,'max_norm_avg_amplitudes'] = max_values_max_norm_avg_1st_name

max_values_std_1st_name = file + '_' + 'max_values_std_1st.npy'
np.save(f'{save_path}/Single_files/{max_values_std_1st_name}', max_values_std_1st)
evoked_db.at[file,'std_amplitudes'] = max_values_std_1st_name

f_rate = np.array([np.sum((max_values[0,o:o+5] < 15 )*0.2) for o in range(0,max_values.shape[1],5)])
f_rate_name = file + '_' + 'f_rate.npy'
np.save(f'{save_path}/Single_files/{f_rate_name}',f_rate)
evoked_db.at[file,'f_rate'] = f_rate_name

decay_tau_values_avg_1st_name = file + '_' + 'max_values_avg_1st.npy'
np.save(f'{save_path}/Single_files/{decay_tau_values_avg_1st_name}', decay_tau_values_avg_1st)
evoked_db.at[file,'decay_tau'] = max_values_avg_1st_name

###### Saving the database as a .csv file
evoked_db.to_csv(f'{save_path}/Evoked_list.csv',index_label=('File')) 