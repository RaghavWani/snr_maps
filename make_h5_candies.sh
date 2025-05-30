#!/bin/bash

# path to the candydates.csv file
CSV_PATH="$1"

# path to the directory containing the .fil files
FIL_DIR="$2"

if [[ -z "$CSV_PATH" || -z "$FIL_DIR" ]]; then
    echo "Run the script with: $0 <csv_path> <fil_dir>!"
    exit 1
fi


# loop through each .fil file in the directory
for filfile in "$FIL_DIR"/*.fil; do
    # getting the beam name for the folder
    beam_name=$(basename "$filfile" .fil)
    mkdir -p "$beam_name"
    echo "Running candies for $filfile in $beam_name/"
    
    # running candies inside the created beam folder
    (cd "$beam_name" && candies make "$CSV_PATH" --fil "$filfile" && candies plot *.h5)
done

