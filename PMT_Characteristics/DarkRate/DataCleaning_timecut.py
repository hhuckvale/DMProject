import pandas as pd
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# arguments for terminal
folder_in = sys.argv[1]       
filename_in = sys.argv[2]     
folder_out = sys.argv[3]      

# outputted file
filename_out = f"{filename_in}_cleaned.csv"
file_path = os.path.join(os.getcwd(), folder_out, filename_out)

# read CSV
csv_path = os.path.join(os.getcwd(), folder_in, f"{filename_in}.csv")
df = pd.read_csv(csv_path)

# make data frame for cleaned data and drop NaNs
cleaned_df = df.dropna(subset=['amplitude', 'time_above_threshold', 'total_time_above', 'event_timestamps', 'value_at_end']).copy()

# --- VOLTAGE EXTRACTION ---
# We use filename_in directly since it contains the voltage string
extracted_volts = pd.Series([filename_in]).str.extract(r"darkrate_(\d+)_")
current_voltage = int(extracted_volts.iloc[0,0])
cleaned_df["voltage"] = current_voltage

'''Deadtime cut'''
cleaned_df = cleaned_df.sort_values('event_timestamps')
cleaned_df['delta_t'] = cleaned_df['event_timestamps'].diff()

# 1. Define the desired bin width
bin_width = 1e-6

# 2. Create bin edges from the minimum to the maximum of your data
# We use np.floor and np.ceil to ensure we cover the whole range
bins_edges = np.arange(np.floor(cleaned_df['delta_t'].min()), 
                       np.ceil(cleaned_df['delta_t'].max()) + bin_width, 
                       bin_width)



counts, bins = np.histogram(cleaned_df['delta_t'].dropna(), bins=bins_edges, range=(0, 0.0002))
mask = bins[:-1] > 1e-8 
max_idx = np.argmax(counts[mask])
bin_centers = 0.5 * (bins[1:] + bins[:-1])
peak_location = bin_centers[mask][max_idx]

afterpulse_threshold = 5 * peak_location
cleaned_df = cleaned_df[cleaned_df['delta_t'] > afterpulse_threshold]

'''Voltage-Specific Integral Cut'''
# Define hard-coded cuts in pC
voltage_cuts = {
    1750: 1.0,
    1800: 1.1,
    1850: 1.5, # Bumped to 2.0 to fix the positive residual/step
    1900: 1.55,
    1950: 1.6,
    2000: 1.9
}

# Get the specific cut for this file
target_cut = voltage_cuts.get(current_voltage, 1.0)
cleaned_df = cleaned_df[cleaned_df['integral_pC'] > target_cut]

'''Other cuts'''
# 15sigma threshold clean (keep events > 0)
cleaned_df = cleaned_df[cleaned_df["time_above_threshold"] > 0]

# Total time above threshold clean (15ns)
cleaned_df = cleaned_df[cleaned_df['total_time_above'] < 2e-8]

# Disregard if value_at_end > 1mV (baseline recovery)
cleaned_df = cleaned_df[cleaned_df['value_at_end'] < 1]

# REMOVED: good_integral = cleaned_df['integral_pC'] > 0.7 
# (This would have overwritten your hard-coded voltage cuts!)

print(f"Processing {current_voltage}V")
print(f"Final Pulse Count: {len(cleaned_df)}")

# export the cleaned data
cleaned_df.to_csv(file_path, sep=',', encoding='utf-8-sig', index=True, header=True)