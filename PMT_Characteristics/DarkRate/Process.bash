#!/bin/bash

sub_folder="PMT3_rawdata"
voltages=("1750" "1800" "1850" "1900" "1950" "2000")
date=190126
extract_folder="PMT3_Stats_Extraction"

# Arrays of trials
arrayOff0=("1" "2" "3"); arrayOff1=("1" "2" "3"); arrayOff2=("1" "2" "3")
arrayOff3=("1" "2" "3"); arrayOff4=("1" "2" "3"); arrayOff5=("1" "2" "3")

arrayOn0=("1" "2" "3"); arrayOn1=("1" "2" "3"); arrayOn2=("1" "2" "3")
arrayOn3=("1" "2" "3"); arrayOn4=("1" "2" "3"); arrayOn5=("1" "2" "3")

arrayOffNames=(arrayOff0 arrayOff1 arrayOff2 arrayOff3 arrayOff4 arrayOff5)
arrayOnNames=(arrayOn0 arrayOn1 arrayOn2 arrayOn3 arrayOn4 arrayOn5)

# Loop over voltages
for i in "${!voltages[@]}"; do
    voltage=${voltages[$i]}
    fnm="${date}_darkrate_${voltage}_10000"

    eval "offArray=(\"\${${arrayOffNames[$i]}[@]}\")"
    eval "onArray=(\"\${${arrayOnNames[$i]}[@]}\")"

    # Combine on/off loops into a temporary list to avoid repeating code logic
    for state in "off" "on"; do
        eval "currentArray=(\"\${${state}Array[@]}\")"
        
        for run in "${currentArray[@]}"; do
            base_filename="${fnm}_${run}_${state}"
            echo "--------------------------------------------"
            echo "Processing: $base_filename"

            # STEP A: Extract stats to a CSV
            # We create a unique CSV name for this specific run
            python ProcessingCode1.py "$sub_folder" "$base_filename" "$extract_folder"
            temp_csv="${extract_folder}/${base_filename}_processed.csv"

            # STEP B: Calculate the threshold 'n' from that CSV
            # This assumes GetThreshold.py prints ONLY the number to the terminal
            n=$(python GetThreshold.py "$temp_csv")

            echo "Calculated threshold n = $n"

            # STEP C: Run final processing using the calculated 'n'
            # I replaced the hardcoded '16' with '$n'
            python ProcessingcodeHHAM_final.py "$sub_folder" "$base_filename" "$n" PMT3_Processed
            
            # Optional: Clean up the temporary CSV
            # rm "$temp_csv"
        done
    done
done
