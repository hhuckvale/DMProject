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
n_str = sys.argv[3]   
foldername_out = sys.argv[4]

n = int(n_str)
#outputted file
filename_out = f"{filename_in}_processed_threshold{n_str}.csv"
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
pulse_end_frac = 0.99
amplitude_start_frac = 0.1
amplitude_end_frac = 1

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
    #amplitude_voltages = corrected_tempVolt[int(amplitude_start_frac*samples[0]):int(amplitude_end_frac*samples[0])]
    threshold = 5 * sigma_b
    if threshold > 0.003:
        min_height = 5 * sigma_b
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
        number_of_peaks = 1

    return amplitude, first_peak_index, first_peak_time, number_of_peaks


def peak_FWHM(corrected_tempVolt, amplitude, first_peak_index, sampleTimes):   
    half_maximum = amplitude/2
    
    #left crossing of half maximum
    left = first_peak_index #begin at index of peak
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
    right = first_peak_index
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
def time_above_threshold(corrected, times, sigma_b, mean_b, first_peak_index, n):
    threshold = mean_b + n*sigma_b

    i = first_peak_index
    while i > 0 and corrected[i] > threshold:
        i -= 1
    # If we hit the boundary and are STILL above threshold → no crossing
    if i == 0 and corrected[i] > threshold:
        return None, None, None
    if i == len(corrected) - 1:
        # This implies the peak is at the very end, and the left search stopped
        # right at the end of the array, meaning we can't interpolate i to i+1.
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

    j = first_peak_index
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
def calc_time_above_3sigma(corrected, times, sigma_b, mean_b, first_peak_index):
    threshold = mean_b + 3*sigma_b

    i = first_peak_index
    while i > 0 and corrected[i] > threshold:
        i -= 1
    # If we hit the boundary and are STILL above threshold → no crossing
    if i == 0 and corrected[i] > threshold:
        return None, None, None
    if i == len(corrected) - 1:
        # This implies the peak is at the very end, and the left search stopped
        # right at the end of the array, meaning we can't interpolate i to i+1.
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

    j = first_peak_index
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
all_number_of_peaks = []

for tempVolt in tempVolts:
    mean_b, sigma_b, corrected_tempVolt = calcbaseline(tempVolt, samples)
    mean_b_array.append(mean_b)
    sigma_b_array.append(sigma_b)
    corrected_tempVolts.append(corrected_tempVolt)

mean_b_array = np.array(mean_b_array)
sigma_b_array = np.array(sigma_b_array)
corrected_tempVolts = np.array(corrected_tempVolts)

for event_index, corrected_tempVolt in enumerate(corrected_tempVolts):

    amplitude, first_peak_index, first_peak_time, number_of_peaks = calcamplitude(corrected_tempVolt, samples)
    singleFWHM, left_time, right_time = peak_FWHM(corrected_tempVolt, amplitude, first_peak_index, sampleTimes)
    time_above, first_cross, second_cross = time_above_threshold(corrected_tempVolt, sampleTimes, sigma_b_array[event_index], mean_b_array[event_index], first_peak_index, n)
    time_above_3sigma, first_3sigma_cross, second_3sigma_cross = calc_time_above_3sigma(corrected_tempVolt, sampleTimes, sigma_b_array[event_index], mean_b_array[event_index], first_peak_index)
    closest_index = time_index(first_3sigma_cross, sampleTimes)
    single_integral = single_charge_integral(corrected_tempVolt, sampleTimes, closest_index)

    all_amplitudes.append(amplitude)
    all_peak_times.append(first_peak_time)
    all_peak_indices.append(first_peak_index)
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
    all_number_of_peaks.append(number_of_peaks)

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
all_number_of_peaks = np.array(all_number_of_peaks)

print('Analysis done, now printing plots')

"""deadtime analysis"""
event_timestamps = np.array(tdatefrac) + np.array(tdate)
event_time_diff = np.diff(event_timestamps)


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
            'peaks_over_threshold': all_number_of_peaks,
            'event_timestamps' : event_timestamps[:-1]                
            })


df.to_csv(file_path, sep=',', encoding='utf-8-sig', index=False, header=True)

print('CSV done!')

