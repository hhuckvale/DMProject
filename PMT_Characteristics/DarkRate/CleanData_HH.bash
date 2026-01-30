#!/bin/bash

input_folder="CORRECTED_INDICES_processed" 
output_folder="CORRECTED_INDICES_AMPCUTS"

for csv_file in "$input_folder"/*.csv;
do
    base_name=$(basename "$csv_file" .csv)
    echo "Cleaning $base_name.csv"

    python FINALDataCleaning_extracutsonamp.py "$input_folder" "$base_name" "$output_folder"

done