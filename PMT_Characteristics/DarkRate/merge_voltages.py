import pandas as pd
import glob
import os
import re

# --- CONFIGURATION ---
input_folder = "PMT7_Processed"  # Folder where your individual CSVs are
output_folder = "PMT7_Voltages_combined"      # Folder to save the combined CSVs
file_extension = ".csv"       # The suffix your files have

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Get all CSV files in the input folder
all_files = glob.glob(os.path.join(input_folder, f"*{file_extension}"))

# Group files by voltage using regex
# This looks for the pattern '_darkrate_NUMBER_' in the filename
voltage_groups = {}

for f in all_files:
    match = re.search(r'_darkrate_(\d+)_', os.path.basename(f))
    if match:
        voltage = match.group(1)
        if voltage not in voltage_groups:
            voltage_groups[voltage] = []
        voltage_groups[voltage].append(f)

# Merge and save
for voltage, files in voltage_groups.items():
    print(f"Merging {len(files)} files for voltage {voltage}...")
    
    # List comprehension to read all CSVs for this voltage
    combined_df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    
    # Sort by timestamp if column exists to keep data chronological
    if 'event_timestamps' in combined_df.columns:
        combined_df = combined_df.sort_values(by='event_timestamps')

    output_path = os.path.join(output_folder, f"Combined_Voltage_{voltage}.csv")
    combined_df.to_csv(output_path, index=False)
    print(f"Saved: {output_path}")

print("\nAll merges complete.")