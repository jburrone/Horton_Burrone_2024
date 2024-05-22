# -*- coding: utf-8 -*-
"""
Created on Mon May 20 18:41:15 2024

@author: Vincenzo Mastrolia (vincenzo.mastrolia@kcl.ac.uk)
Research associate at the Burrone lab, Centre for Developmental Neurobiology, King's College London

This Script allows to extract parameters from spontaneous events from patch clamp recordings in neurons
"""

import pyabf
import numpy as np
from signal_processing_tools import signal_processing as sp


###############################################################################

path = 'please/insert/your/path/here' # Path where file is loc
file = 'insert your .abf file here'
save_path = 'please/insert/your/path/here' # Path where files are saved

cursors = [0,60] # Interval to analyse within the sweep, in seconds
direction = 'negative' # The direction of the spontaneous currents, whether inward or
#                       outward. Usually sEPSC are inward, while sIPSC are 
#                       outward.


abf = pyabf.ABF(path + "/" + file + ".abf") # To load the .abf files
n_sweeps = abf.sweepCount # number of sweeps to analyse
sampling_rate = abf.dataRate # Sampling rate of the acquired trace

# Script to analyse spontaneous events
file_events,file_parameters = sp.sPSC_analysis(path, save_path, file, abf, sampling_rate, n_sweeps, cursors, direction, raw_trace = True)
