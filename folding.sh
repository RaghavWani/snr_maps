#!/bin/bash

# Script to fold filterbank files for a pulsar using prepfold
# Source relevant env before running this script (/lustre_archive/apps/tdsoft/env.sh)
# Usage: ./folding.sh
# Last updated: 23 July 2025; Author: Raghav Wani

n_beams=800 #number of beams
par_file="B0329+54.par" #parameter file for pulsar from ATNF catalogue
fil_dir="SPLT_1Jan_B0329_B3" #path of directory containing filterbank files

out_dir="${fil_dir}/FOLDING_OUTPUTS" #output directory
mkdir -p "$out_dir"

# Beam wise folding
for i in $(seq 0 $((n_beams-1))); do #Loop over beams
        prepfold -timing "$par_file" "${fil_dir}/BM${i}.fil" -noxwin -o "${out_dir}/L${i}aa"
        echo "------------ Folding for Beam number ${i} done --------------------------------"
done