#!/bin/bash

input_folder="PMT7_Processed_valueatend" 
output_folder="PMT7_Clean_valueinlast10"

for csv_file in "$input_folder"/*.csv;
do
    base_name=$(basename "$csv_file" .csv)
    echo "Cleaning $base_name.csv"

    python DataCleaning_valueinlast10.py "$input_folder" "$base_name" "$output_folder"

done