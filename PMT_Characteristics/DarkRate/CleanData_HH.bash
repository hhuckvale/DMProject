#!/bin/bash

input_folder="PMT5_AllProcessed" 
output_folder="PMT5_Cleaned"

for csv_file in "$input_folder"/*.csv;
do
    base_name=$(basename "$csv_file" .csv)
    echo "Cleaning $base_name.csv"

    python DataCleaning_final.py "$input_folder" "$base_name" "$output_folder"

done