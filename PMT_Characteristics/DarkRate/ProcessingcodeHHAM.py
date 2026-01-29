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

#arguments for terminal
filename_in = sys.argv[1]
n_str = sys.argv[2]   
foldername_out = sys.argv[3]

n = int(n_str)
#outputted file
filename_out = f"{filename_in}_processed_threshold{n_str}.csv"
file_path = os.path.join(os.getcwd(), foldername_out, filename_out)


startEvent = 0
nEvents = 10000


volts, tstart, tscale, tfrac, tdatefrac, tdate = tekwfm.read_wfm(filename_in+".wfm")
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
amplitude_start_frac = 0
amplitude_end_frac = 0.6
pulse_end_frac = 0.99

timediffbins=1000 #number of bins in the histogram for time between each event
zoomedtimediffbins=10000 #number of bins in the zoomed in histogram for time between each event
zoom=0.0005 #the max x coordinate for zooming in on this histogram (looking at range 0 to 0.0005)



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



def peak_FWHM(corrected_tempVolt, amplitude, peak_index, sampleTimes):   
    half_maximum = amplitude/2
    
    #left crossing of half maximum
    left = peak_index #begin at index of peak
    #move left until voltage rises to half maximum, stopping if left>0 (beginning of window reached)
    while left>0 and corrected_tempVolt[left]>half_maximum:
        left -= 1 #index moves left
    #interpolate between left and left+1
    if left == 0: #if beginning of window reached, take first time value as left crossing
        left_time = sampleTimes[0]
    else: 
        y1, y2 = corrected_tempVolt[left], corrected_tempVolt[left+1]
        x1, x2 = sampleTimes[left], sampleTimes[left+1]
        left_time = x1 + (half_maximum - y1)/(y2-y1) * (x2-x1)
        
    #same again for right
    right = peak_index
    while right<len(corrected_tempVolt)-1 and corrected_tempVolt[right]>half_maximum:
        right += 1
    #interpolate between right and right-1
    if right == len(corrected_tempVolt)-1: #if end of window reached, take last time value as right crossing
        right_time = sampleTimes[-1]
    else: 
        y1, y2 = corrected_tempVolt[right], corrected_tempVolt[right-1]
        x1, x2 = sampleTimes[right], sampleTimes[right-1]
        right_time = x1 + (half_maximum - y1)/(y2-y1) * (x2-x1)
    
    singleFWHM = right_time - left_time
    return singleFWHM, left_time, right_time



#threshold passed as an argument, n
def time_above_threshold(corrected, times, sigma_b, mean_b, peak_index, n):
    threshold = mean_b + n*sigma_b

    i = peak_index
    while i > 0 and corrected[i] > threshold:
        i -= 1
    # If we hit the boundary and are STILL above threshold → no crossing
    if i == 0 and corrected[i] > threshold:
        return None, None, None
    # Interpolate only if we actually crossed between i and i+1
    y1, y2 = corrected[i], corrected[i+1]
    x1, x2 = times[i], times[i+1]
    # Ensure interpolation makes physical sense
    if (y2 - y1) == 0:
        return None, None, None
    frac = (threshold - y1) / (y2 - y1)
    if not (0 <= frac <= 1):  
        return None, None, None
    first_cross = x1 + frac * (x2 - x1)

    j = peak_index
    while j < len(corrected)-1 and corrected[j] > threshold:
        j += 1
    if j == len(corrected)-1 and corrected[j] > threshold:
        return None, None, None
    y1, y2 = corrected[j-1], corrected[j]
    x1, x2 = times[j-1], times[j]
    if (y2 - y1) == 0:
        return None, None, None
    frac = (threshold - y1) / (y2 - y1)
    if not (0 <= frac <= 1):
        return None, None, None
    second_cross = x1 + frac * (x2 - x1)

    time_above = second_cross - first_cross
    if time_above <= 0:
        return None, None, None

    return time_above, first_cross, second_cross



#use 3sigma for the integral bound 1
def calc_time_above_3sigma(corrected, times, sigma_b, mean_b, peak_index):
    threshold = mean_b + 3*sigma_b

    i = peak_index
    while i > 0 and corrected[i] > threshold:
        i -= 1
    # If we hit the boundary and are STILL above threshold → no crossing
    if i == 0 and corrected[i] > threshold:
        return None, None, None
    # Interpolate only if we actually crossed between i and i+1
    y1, y2 = corrected[i], corrected[i+1]
    x1, x2 = times[i], times[i+1]
    # Ensure interpolation makes physical sense
    if (y2 - y1) == 0:
        return None, None, None
    frac = (threshold - y1) / (y2 - y1)
    if not (0 <= frac <= 1):  
        return None, None, None
    first_3sigma_cross = x1 + frac * (x2 - x1)

    j = peak_index
    while j < len(corrected)-1 and corrected[j] > threshold:
        j += 1
    if j == len(corrected)-1 and corrected[j] > threshold:
        return None, None, None
    y1, y2 = corrected[j-1], corrected[j]
    x1, x2 = times[j-1], times[j]
    if (y2 - y1) == 0:
        return None, None, None
    frac = (threshold - y1) / (y2 - y1)
    if not (0 <= frac <= 1):
        return None, None, None
    second_3sigma_cross = x1 + frac * (x2 - x1)

    time_above_3sigma = second_3sigma_cross - first_3sigma_cross
    if time_above_3sigma <= 0:
        return None, None, None

    return time_above_3sigma, first_3sigma_cross, second_3sigma_cross



def time_index(first_3sigma_cross, sampleTimes):
    if first_3sigma_cross is None:
        return len(sampleTimes) - 1
    else:
        time_difference = np.abs(np.array(sampleTimes) - first_3sigma_cross)
        closest_index = np.argmin(time_difference)
        return closest_index
    


#integrate a single event
def single_charge_integral(corrected_tempVolt, sampleTimes, closest_index):
    # define bound 1 and 2
    bound_1 = closest_index
    bound_2 = int(pulse_end_frac*len(corrected_tempVolt))
    integral_voltages = corrected_tempVolt[bound_1:bound_2]
    time_slice = sampleTimes[bound_1:bound_2]
    single_integral = integrate.trapezoid(integral_voltages, x=time_slice)
    #return a single number for the integral of one event
    return single_integral



def correcting_tempVolts(tempVolts, mean_b):
    corrected_tempVolts = []
    for i, tempVolt in enumerate(tempVolts):
        corrected = tempVolt - mean_b[i]
        corrected_tempVolts.append(corrected)
    return corrected_tempVolts



"""Working on all events"""
mean_b_array = []
sigma_b_array = []
corrected_tempVolts = []
all_amplitudes = []
all_peak_times = []
all_peak_indices = []
all_FWHM = []
all_t1 = []
all_t2 = []
all_time_above_thresh = []
all_first_crossing = []
all_second_crossing = []
all_tabove_3sig = []
all_first_t_3sig = []
all_second_t_3sig = []
all_closest_indices = []
all_integrals = []

for tempVolt in tempVolts:
    mean_b, sigma_b, corrected_tempVolt = calcbaseline(tempVolt, samples)
    mean_b_array.append(mean_b)
    sigma_b_array.append(sigma_b)
    corrected_tempVolts.append(corrected_tempVolt)

mean_b_array = np.array(mean_b_array)
sigma_b_array = np.array(sigma_b_array)
corrected_tempVolts = np.array(corrected_tempVolts)

for event_index, corrected_tempVolt in enumerate(corrected_tempVolts):

    amplitude, peak_time, peak_index = calcamplitude(corrected_tempVolt, samples)
    singleFWHM, left_time, right_time = peak_FWHM(corrected_tempVolt, amplitude, peak_index, sampleTimes)
    time_above, first_cross, second_cross = time_above_threshold(corrected_tempVolt, sampleTimes, sigma_b_array[event_index], mean_b_array[event_index], peak_index, n)
    time_above_3sigma, first_3sigma_cross, second_3sigma_cross = calc_time_above_3sigma(corrected_tempVolt, sampleTimes, sigma_b_array[event_index], mean_b_array[event_index], peak_index)
    closest_index = time_index(first_3sigma_cross, sampleTimes)
    single_integral = single_charge_integral(corrected_tempVolt, sampleTimes, closest_index)

    all_amplitudes.append(amplitude)
    all_peak_times.append(peak_time)
    all_peak_indices.append(peak_index)
    all_FWHM.append(singleFWHM)
    all_t1.append(left_time)
    all_t2.append(right_time)
    all_time_above_thresh.append(time_above)
    all_first_crossing.append(first_cross)
    all_second_crossing.append(second_cross)
    all_tabove_3sig.append(time_above_3sigma)
    all_first_t_3sig.append(first_3sigma_cross)
    all_second_t_3sig.append(second_3sigma_cross) 
    all_closest_indices.append(closest_index)
    all_integrals.append(single_integral)

all_amplitudes = np.array(all_amplitudes)
all_peak_times = np.array(all_peak_times)
all_peak_indices = np.array(all_peak_indices)
all_FWHM = np.array(all_FWHM)
all_t1 = np.array(all_t1)
all_t2 = np.array(all_t2)
all_time_above_thresh = np.array(all_time_above_thresh)
all_first_crossing = np.array(all_first_crossing)
all_second_crossing = np.array(all_second_crossing)
all_tabove_3sig = np.array(all_tabove_3sig)
all_first_t_3sig = np.array(all_first_t_3sig)
all_second_t_3sig = np.array(all_second_t_3sig)
all_closest_indices = np.array(all_closest_indices)
all_integrals = np.array(all_integrals)

print('Analysis done, now printing plots')

"""deadtime analysis"""
event_timestamps = np.array(tdatefrac) + np.array(tdate)
event_time_diff = np.diff(event_timestamps)

event_numbers = []
number_of_multiple_maxima = []
number_of_maxima_all = []

#maximum_threshold = 0.03 #threshold for code to consider a peak a maximum in the event, eg. 30mV
min_peak_sep = 0.1*samples[0] #minimum separate between 'peaks' eg. 10% of all samples

#loop over events
for event_index, corrected_tempVolt in enumerate(corrected_tempVolts):
    count = 0  #reset counter of maxima
    prev_max_index = -min_peak_sep
    maximum_threshold = 15 * sigma_b_array[event_index]
    #loop over all voltage points in each event
    
    for i, voltage in enumerate(corrected_tempVolt):
        if voltage > maximum_threshold:
            if i - prev_max_index >= min_peak_sep:
                count += 1  #count number of times passes below threshold
                prev_max_index = i
    #add ALL event counts to all maxima array
    number_of_maxima_all.append(count)

    if count > 1: #store events with more than 1 maximum
        event_numbers.append(event_index) 
        number_of_multiple_maxima.append(count) 

event_numbers = np.array(event_numbers)
number_of_multiple_maxima = np.array(number_of_multiple_maxima)   
number_of_maxima_all = np.array(number_of_maxima_all)

Resistance = ((1/10000) + (1/50))**(-1)
all_integrals_picocharge = (all_integrals * 1e12) / Resistance

"""export data in a csv"""

df = pd.DataFrame({
            'baseline': mean_b_array,
            'sd_baseline': sigma_b_array,
            'amplitude': all_amplitudes,
            'peak_time': all_peak_times,
            'FWHM': all_FWHM,
            't1_of_FWHM': all_t1,
            't2_of_FWHM': all_t2,
            'time_above_threshold': all_time_above_thresh,
            'first_threshold_crossing': all_first_crossing,
            'second_threshold_crossing': all_second_crossing,
            'integral_pC': all_integrals_picocharge,
            'peaks_over_threshold': number_of_maxima_all,
            'event_timestamps' : event_timestamps[:-1]                
            })


df.to_csv(file_path, sep=',', encoding='utf-8-sig', index=False, header=True)

print('CSV done!')

