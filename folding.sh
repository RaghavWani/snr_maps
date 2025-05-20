#!/bin/bash

n_beams=800 #number of beams
par_file="B0329+54.par" #par file for pulsar
fil_dir="SPLT_1Jan_B0329_B3" #path of directory containing filterbnk files

# Check if required arguments are provided
if [[ -z $n_beams || -z $par_file || -z $fil_dir ]]; then
    echo "Usage: $0 <n_beams> <par_file> <fil_dir>"
    exit 1
fi

#parent_dir=$(dirname "$fil_dir")
out_dir="${fil_dir}/FOLDING_OUTPUTS" #output directory
mkdir -p "$out_dir"

for i in $(seq 0 $((n_beams-1))); do #Loop over beams
        prepfold -timing "$par_file" "${fil_dir}/BM${i}.fil" -noxwin -o "${out_dir}/L${i}aa"
        echo "-----------------Beam number ${i} done--------------------------------"
done