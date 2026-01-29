#!/bin/bash

sub_folder="PMT1_rawdata"
#voltages=("1850" "1900" "1950" "2000")
#date=141125
#date=241125
voltages=("1750" "1800")
date=131125

#arrays of trials for each voltage (on and off)
arrayOff0=("1" "2" "3") #for 1850
arrayOff1=("1" "2" "3") #for 1900
arrayOff2=("1" "2" "3") #for 1950
arrayOff3=("1" "2" "3")     #for 2000
#arrayOff0=("1" "2" "3") #for 1750
#arrayOff1=("1" "2" "3") #for 1800

arrayOn0=("1" "2" "3")
arrayOn1=("1" "2" "3")
#arrayOn2=("1" "3")
#arrayOn3=("1" "2" "3")
arrayOn2=("1" "2" "3")
arrayOn3=("1" "2" "3")

#names of sub-arrays
arrayOffNames=(arrayOff0 arrayOff1 arrayOff2 arrayOff3)
arrayOnNames=(arrayOn0 arrayOn1 arrayOn2 arrayOn3)
#arrayOffNames=(arrayOff0 arrayOff1)
#arrayOnNames=(arrayOn0 arrayOn1)

#loop over voltages
for i in "${!voltages[@]}"; 
do
    voltage=${voltages[$i]}
    fnm="${date}_darkrate_${voltage}_10000"

    offArrayName=${arrayOffNames[$i]}
    eval "offArray=(\"\${${offArrayName}[@]}\")"

    onArrayName=${arrayOnNames[$i]}
    eval "onArray=(\"\${${onArrayName}[@]}\")"

    #loop over off trials
    for run in "${offArray[@]}"; 
    do
        filename="${sub_folder}/${fnm}_${run}_off"
        echo "Doing this file now: $filename"
        python ProcessingcodeHHAM_findpeaks_PMT1.py "$sub_folder" "${fnm}_${run}_off" 16 PMT1_Processed_16
    done

    #loop over on trials
    for run in "${onArray[@]}"; 
    do
        filename="${sub_folder}/${fnm}_${run}_on"
        echo "Doing this file now: $filename"
        python ProcessingcodeHHAM_findpeaks_PMT1.py "$sub_folder" "${fnm}_${run}_on" 16 PMT1_Processed_16
    done
done


#type ./ProcessData.bash in terminal to run