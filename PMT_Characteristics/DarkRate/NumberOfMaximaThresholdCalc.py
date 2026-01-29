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
import matplotlib.cm as cm
from collections import Counter

import gc
import re
import tempfile
from collections import defaultdict, namedtuple
from pathlib import Path
import pandas as pd

doVerbose = False
filename_base = "PMT1_rawdata/141125_darkrate_1850_10000_2_on"

startEvent = 0
nEvents = 10000

# code here to pass arguments

volts, tstart, tscale, tfrac, tdatefrac, tdate = tekwfm.read_wfm(filename_base+".wfm")
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

baseline_end_frac = 0.15
amplitude_start_frac = 0
amplitude_end_frac = 0.6
pulse_end_frac = 0.99

# ranges for No cleaning cuts
peak_t_range = (-20e-9,80e-9)
first_crossing_range = (-20e-9,80e-9)
v_window_size = (-25e-3,225e-3) #remember flipping of y axis occurs
tolerance = 0 #discard events if they are 1uV away from being cut off by the window

def calcbaseline(tempVolt, samples):
    baseline_voltages = tempVolt[0: int(baseline_end_frac*samples[0])]
    mean_b = np.mean(baseline_voltages)
    sigma_b = np.std(baseline_voltages, ddof=1)
    correct_tempVolt = tempVolt - mean_b
    return mean_b, sigma_b, correct_tempVolt

def calcamplitude(corrected_tempVolt, samples):
    amplitude_voltages = corrected_tempVolt[int(amplitude_start_frac*samples[0]):int(amplitude_end_frac*samples[0])]
    amplitude = np.max(amplitude_voltages)
    amp_index = np.argmax(amplitude_voltages)
    peak_index = amp_index + int(amplitude_start_frac*samples[0])
    peak_time = sampleTimes[peak_index]
    return amplitude, peak_time, peak_index

def correcting_tempVolts(tempVolts, mean_b):
    corrected_tempVolts = []
    for i, tempVolt in enumerate(tempVolts):
        corrected = tempVolt - mean_b[i]
        corrected_tempVolts.append(corrected)
    return corrected_tempVolts

mean_b_array = []
sigma_b_array = []
corrected_tempVolts = []
discarded_event_indices = [] 

for i, tempVolt in enumerate(tempVolts):
    if np.any((tempVolt>=v_window_size[0]-tolerance)&(tempVolt<=v_window_size[0])):
        discarded_event_indices.append(i)
    else: 
        mean_b, sigma_b, corrected_tempVolt = calcbaseline(tempVolt, samples)
        mean_b_array.append(mean_b)
        sigma_b_array.append(sigma_b)
        corrected_tempVolts.append(corrected_tempVolt)

mean_b_array = np.array(mean_b_array)
sigma_b_array = np.array(sigma_b_array)
corrected_tempVolts = np.array(corrected_tempVolts)


event_numbers_by_multiplier = []
multiple_maxima_by_multiplier = []
maxima_counts_by_multiplier = []
maxima_distribution_by_multiplier = [] # stores {0: x, 1: y, 2: z, ...} for each multiplier

#maximum_threshold = 0.03 #threshold for code to consider a peak a maximum in the event, eg. 30mV
min_peak_sep = 0.05*samples[0] #minimum separate between 'peaks' eg. 10% of all samples
multipliers = np.linspace(3, 23, 11)

for multiplier in multipliers:
    number_of_maxima_all = []   
    number_of_multiple_maxima = []
    event_numbers = []

#loop over events
    for event_index, corrected_tempVolt in enumerate(corrected_tempVolts):
        count = 0  #reset counter of maxima
        prev_max_index = -min_peak_sep
        maximum_threshold = multiplier * sigma_b_array[event_index]
        #loop over all voltage points in each event
    
        for i, voltage in enumerate(corrected_tempVolt):
            if voltage > maximum_threshold:
                if i - prev_max_index >= min_peak_sep:
                    count += 1  #count number of times passes below threshold
                    prev_max_index = i
    #add ALL event counts to all maxima array
        final_count = max(1, count)
        number_of_maxima_all.append(final_count)
    

        if final_count > 1: #store events with more than 1 maximum
            event_numbers.append(event_index) 
            number_of_multiple_maxima.append(count) 

    maxima_counts_by_multiplier.append(np.array(number_of_maxima_all))
    event_numbers_by_multiplier.append(np.array(event_numbers))
    multiple_maxima_by_multiplier.append(np.array(number_of_multiple_maxima))
   
    distribution = Counter(number_of_maxima_all)
    maxima_distribution_by_multiplier.append(distribution)

#plot histogram of number of maxima
total_maxima_by_multiplier = []


plt.figure(figsize=(6, 6))
colors = cm.get_cmap('tab20', len(multipliers))

for i, multiplier in enumerate(multipliers):
    distribution = maxima_distribution_by_multiplier[i]
    x = sorted(distribution.keys())  # number of maxima
    y = [distribution[k] for k in x]  # number of events
    errors = np.sqrt(y)
    plt.errorbar(x, y, yerr=errors, markersize=5, fmt='o-', color=colors(i), label=f'{multiplier:.1f}×σ', capsize=3)

plt.xlabel('Number of Maxima in Event')
plt.ylabel('Number of Events')
plt.title('Event Count by Maxima Count Across Thresholds')
plt.legend()
plt.grid(True)
plt.show()

# counts as a function of threshold for one maxima
single_maxima_counts = []
for distribution in maxima_distribution_by_multiplier:
    count = distribution.get(1,0)
    single_maxima_counts.append(count)
single_maxima_counts = np.array(single_maxima_counts)
single_maxima_counts_errors = np.sqrt(single_maxima_counts)

two_maxima_counts = []
for distribution in maxima_distribution_by_multiplier:
    count2 = distribution.get(2,0)
    two_maxima_counts.append(count2)
two_maxima_counts = np.array(two_maxima_counts)
two_maxima_counts_error = np.sqrt(two_maxima_counts)

three_maxima_counts = []
for distribution in maxima_distribution_by_multiplier:
    count3 = distribution.get(3,0)
    three_maxima_counts.append(count3)
three_maxima_counts = np.array(three_maxima_counts)
three_maxima_counts_error = np.sqrt(three_maxima_counts)


plt.figure(figsize=(6,6))
plt.errorbar(multipliers, single_maxima_counts, yerr=single_maxima_counts_errors, fmt='o-', markersize=4, capsize=4, color='darkblue', label='Events with 1 Maximum')
plt.errorbar(multipliers, two_maxima_counts, yerr=two_maxima_counts_error, fmt='o-', markersize=4, capsize=4, color='red', label='Events with 2 Maxima')
plt.errorbar(multipliers, three_maxima_counts, yerr=three_maxima_counts_error, fmt='o-', markersize=4, capsize=4, color='green', label='Events with 3 Maxima')
plt.axvline(15, color='black', linestyle='--', label='Selected threshold = 15 $\sigma$')
plt.xlabel('Threshold Multiplier ($\sigma$ Units)')
plt.ylabel('Number of Events')
plt.ylim(0, 9000)
plt.legend()
plt.show()