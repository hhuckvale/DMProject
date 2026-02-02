#!/bin/bash

input_folder="PMT4_Processed_totaltimecolumn" 
output_folder="PMT4_Cleaned_allTgreaterthan0"

for csv_file in "$input_folder"/*.csv;
do
    base_name=$(basename "$csv_file" .csv)
    echo "Cleaning $base_name.csv"

    python DataCleaning_Tabovecutimproved.py "$input_folder" "$base_name" "$output_folder"

done