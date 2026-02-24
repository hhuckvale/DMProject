import sys
import time
import matplotlib.pyplot as plt # https://matplotlib.org/
import numpy as np # http://www.numpy.org/
import tekwfm
import os
from os.path import exists
from scipy import signal
from scipy import integrate
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from scipy.special import gamma
import matplotlib.cm as cm
from collections import Counter

import gc
import re
import tempfile
from collections import defaultdict, namedtuple
from pathlib import Path
import pandas as pd

doVerbose = False

#timing window HAS to be -20ns to 80ns

#arguments for terminal
folder_in = sys.argv[1]
filename_in = sys.argv[2]
foldername_out = sys.argv[3]


#outputted file
filename_out = f"{filename_in}_processed.csv"
file_path = os.path.join(os.getcwd(), foldername_out, filename_out)


startEvent = 0
nEvents = 10000


volts, tstart, tscale, tfrac, tdatefrac, tdate = tekwfm.read_wfm(folder_in + "/" + filename_in+".wfm")
if doVerbose:
    print('>>>>', volts, tstart, tscale, tfrac, tdatefrac, tdate)

samples = volts.shape

print(len(volts))
print(tscale)
print('samples', samples[0])
tstop = samples[0]*tscale+tstart
sampleTimes = [tstart+x*tscale for x in range(samples[0])]


negVolts = [volts[:, event] for event in range(startEvent, startEvent + nEvents)] # all events, array of array of voltages
#tempVolts = -1 * negVolts # makes pulses positive
tempVolts = [-v for v in negVolts]

"""variables which can be changed"""
baseline_end_frac = 0.15
baseline_start_frac_last = 0.9
pulse_end_frac = 0.99
amplitude_start_frac = 0.1
amplitude_end_frac = 1

timediffbins=1000 #number of bins in the histogram for time between each event
zoomedtimediffbins=10000 #number of bins in the zoomed in histogram for time between each event
zoom=0.0005 #the max x coordinate for zooming in on this histogram (looking at range 0 to 0.0005)



def calcbaseline(tempVolt, samples):
    baseline_voltages = tempVolt[0: int(baseline_end_frac*samples[0])]
    mean_baseline_voltages = tempVolt[0: int(baseline_end_frac*samples[0])]
    mean_b = np.mean(mean_baseline_voltages)
    sigma_b = np.std(baseline_voltages, ddof=1)
    correct_tempVolt = tempVolt - mean_b
    return mean_b, sigma_b, correct_tempVolt



def calcamplitude(corrected_tempVolt, samples):
    #amplitude_voltages = corrected_tempVolt[int(amplitude_start_frac*samples[0]):int(amplitude_end_frac*samples[0])]
    threshold = 3 * sigma_b
    if threshold > 0.003:
        min_height = 3 * sigma_b
    else:
        min_height = 0.003
    min_separation = int(0.05*samples[0])
    prominance = 0.005
    #peaks, _ = find_peaks(amplitude_voltages, height=min_height, distance=min_separation, prominence=prominance)
    peaks, properties = find_peaks(corrected_tempVolt, height=min_height, distance=min_separation, prominence=prominance)
    number_of_peaks = len(peaks)
    if len(peaks) > 0:
        #amplitude = amplitude_voltages[peaks[0]]
        amplitude = corrected_tempVolt[peaks[0]]
        #first_peak_index = peaks[0] + int(amplitude_start_frac*samples[0])
        first_peak_index = peaks[0]
        first_peak_time = sampleTimes[first_peak_index]
    else:
        #amplitude = np.max(amplitude_voltages)
        amplitude = np.max(corrected_tempVolt)
        amp_index = np.argmax(corrected_tempVolt)
        #amp_index = np.argmax(amplitude_voltages)
        #first_peak_index = amp_index + int(amplitude_start_frac*samples[0])
        first_peak_index = amp_index
        first_peak_time = sampleTimes[first_peak_index]
    if amplitude > 0.22:
       number_of_peaks = 1

    return amplitude, first_peak_index, first_peak_time, number_of_peaks

event_index_array=[]
mean_b_array = []
sigma_b_array = []
corrected_tempVolts = []
all_amplitudes = []
all_peak_times = []
all_peak_indices = []

for tempVolt in tempVolts:
    mean_b, sigma_b, corrected_tempVolt = calcbaseline(tempVolt, samples)
    mean_b_array.append(mean_b)
    sigma_b_array.append(sigma_b)
    corrected_tempVolts.append(corrected_tempVolt)

mean_b_array = np.array(mean_b_array)
sigma_b_array = np.array(sigma_b_array)
corrected_tempVolts = np.array(corrected_tempVolts)

for event_index, corrected_tempVolt in enumerate(corrected_tempVolts):
    true_event_number = startEvent + event_index
    amplitude, first_peak_index, first_peak_time, number_of_peaks = calcamplitude(corrected_tempVolt, samples)

    event_index_array.append(true_event_number)
    all_amplitudes.append(amplitude)
    all_peak_times.append(first_peak_time)
    all_peak_indices.append(first_peak_index)

all_amplitudes = np.array(all_amplitudes)
all_peak_times = np.array(all_peak_times)
all_peak_indices = np.array(all_peak_indices)

"""deadtime analysis"""
event_timestamps = np.array(tdatefrac) + np.array(tdate)
event_time_diff = np.diff(event_timestamps)

df = pd.DataFrame({
            'event_id': event_index_array, 
            'baseline': mean_b_array,
            'sd_baseline': sigma_b_array,
            'amplitude': all_amplitudes,
            'peak_time': all_peak_times,
            'event_timestamps' : event_timestamps[:-1]             
            })

df.to_csv(file_path, sep=',', encoding='utf-8-sig', index=False, header=True)
