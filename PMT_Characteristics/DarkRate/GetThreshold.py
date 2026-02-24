import pandas as pd
import numpy as np
import sys

# Load data from the CSV generated in Step A
# Usage: python GetThreshold.py temp_stats.csv
input_csv = sys.argv[1]
df = pd.read_csv(input_csv)

# Configuration
threshold_factors = np.arange(0, 40) # Integer values 0 to 39
sigma_column = df['sd_baseline']
amplitude_column = df['amplitude']
timestamp_column = df['event_timestamps']

timetaken = timestamp_column.iloc[-1] - timestamp_column.iloc[0]

# 1. Calculate Dark Rates for all factors
dark_rates = []
for factor in threshold_factors:
    # count how many pulses exceed (factor * sigma)
    count = np.sum(amplitude_column > (sigma_column * factor))
    dark_rates.append(count / timetaken)

dark_rates = np.array(dark_rates)

# 2. Define the starting point (3 sigma)
# Since threshold_factors is [0, 1, 2, 3...], index 3 is factor 3
start_rate = dark_rates[3]
target_rate = start_rate * 0.98  # 2% drop means 98% of original

## 3. Find the CLOSEST integer sigma (starting from 3 upwards)
# We slice from index 3 to 39
search_range_rates = dark_rates[3:]
search_range_factors = threshold_factors[3:]

# Calculate absolute difference between each rate and the target
differences = np.abs(search_range_rates - target_rate)

# Find the index of the minimum difference
best_index = np.argmin(differences)
n_threshold = search_range_factors[best_index]

# 5. Output for Bash
print(int(n_threshold))