#!/bin/bash

make_beam_file=$1
sel_beam=$2
nbeams=$3
radec_offset=$4

if [ -z "$make_beam_file" ] || [ -z "$sel_beam" ] || [ -z "$nbeams" ] || [ -z "$radec_offset" ]; then
    echo "Usage: $0 <make_beam_file> <sel_beam> <nbeams> <radec_offset>"
    echo "Example: $0 make_beam_HHMMSS.pkl g0 160 0,0"
    exit 1
fi

if [ ! -f "$make_beam_file" ]; then
    echo "Error: File '$make_beam_file' not found!"
    exit 1
fi

if ! [[ "$nbeams" =~ ^[0-9]+$ ]]; then
    echo "Error: nbeams must be a positive integer!"
    exit 1
fi

if ! [[ "$radec_offset" =~ ^-?[0-9]+,-?[0-9]+$ ]]; then
    echo "Error: radec_offset must be in the format 'RA_offset,DEC_offset' (e.g., '0,0' or '-1.5,2.5')!"
    exit 1
fi

echo "Simulating beam with the following parameters:"
echo "  Make beam file: $make_beam_file"
echo "  Selected beam: $sel_beam"
echo "  Number of beams: $nbeams"
echo "  RA/DEC offset: $radec_offset"

/lustre_archive/apps/correlator/code/conda/bin/python3.1 /lustre_archive/apps/correlator/code/BMSTEERING-NEW/at/make_tile.py -F ${make_beam_file} -${sel_beam} -n ${nbeams} -C -R ${radec_offset}

