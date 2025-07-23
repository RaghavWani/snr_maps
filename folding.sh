#!/bin/bash

n_beams=800 #number of beams
par_file="B0329+54.par" #parameter file for pulsar
fil_dir="SPLT_1Jan_B0329_B3" #path of directory containing filterbank files

out_dir="${fil_dir}/FOLDING_OUTPUTS" #output directory
mkdir -p "$out_dir"

# Beam wise folding
for i in $(seq 0 $((n_beams-1))); do #Loop over beams
        prepfold -timing "$par_file" "${fil_dir}/BM${i}.fil" -noxwin -o "${out_dir}/L${i}aa"
        echo "------------ Folding for Beam number ${i} done --------------------------------"
done