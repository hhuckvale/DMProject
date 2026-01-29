#!/bin/bash

input_folder="PMT1_Processed_16" 
output_folder="PMT1_Cleaned_16"

for csv_file in "$input_folder"/*.csv;
do
    base_name=$(basename "$csv_file" .csv)
    echo "Cleaning $base_name.csv"

    python FINALDataCleaning.py "$input_folder" "$base_name" "$output_folder"

done