#!/bin/bash

input_folder="PMT7_Processed" 
output_folder="PMT7_Cleaned"

for csv_file in "$input_folder"/*.csv;
do
    base_name=$(basename "$csv_file" .csv)
    echo "Cleaning $base_name.csv"

    python DataCleaning_timecut.py "$input_folder" "$base_name" "$output_folder"

done