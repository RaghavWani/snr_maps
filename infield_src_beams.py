# Python script to plot beam pattern from ahdr files
# Last updated: 09th January 2026; Author: Raghav Wani

import numpy as np
import matplotlib.pyplot as plt
import argparse
import pandas as pd
import os, pytz, glob, sys
from astropy import units as u
from astropy.coordinates import SkyCoord, TETE
from astropy.time import Time
from datetime import datetime

def getmjd(t: datetime):
    """
    Function to convert an IST Date/Time to MJD
    """

    localtz = pytz.timezone("Asia/Kolkata")
    localdt = localtz.localize(t, is_dst=None)
    utcdt = localdt.astimezone(pytz.utc)
    mjd = Time(utcdt).mjd
    return mjd

def get_precessed_coords(src_ra, src_dec, mjd):
    """
    Function to precess the source coordinates to the specified epoch.
    Input: src_ra (Right Ascension in degrees), src_dec (Declination in degrees), epoch (e.g., 'J2000').
    Returns: Precessed SkyCoord object.
    """
    src_ra = src_ra.split(":")
    src_ra = f"{src_ra[0]}h{src_ra[1]}m{src_ra[2]}s"

    src_dec = src_dec.split(":")
    src_dec = f"{src_dec[0]}d{src_dec[1]}m{src_dec[2]}s"

    obstime = Time(mjd, format="mjd")
    coords = SkyCoord(src_ra, src_dec, frame="icrs")
    precessed_coord = coords.transform_to(TETE(obstime=obstime))
  
    return precessed_coord

def get_ahdr_files(directory_path):
    ahdr_files = glob.glob(directory_path)
    return sorted(glob.glob(os.path.join(directory_path, "*.ahdr")))

def ra_dec_from_ahdr(ahdr_files):
    """
    Extracts RA, DEC, BM-Idx, and BM-SubIdx values from ahdr files in a directory.
    Input: Directory containing ahdr files.
    Returns: A DataFrame containing the extracted data: RA, DEC, BM-Idx, BM-SubIdx.
    """
    data_columns = ["RA", "DEC", "BM-Idx", "BM-SubIdx"]
    ahdr_data = pd.DataFrame(columns=data_columns)

    for file in ahdr_files:
        if os.path.exists(file):
            with open(file, "r") as infile:
                lines = infile.readlines()

                src_name = lines[7].strip().split()[-1]
                pc_ra_line = lines[8]
                pc_dec_line = lines[9]
                nbeams_line = lines[25]
                beams_per_host_line = lines[26]
                freq = float(lines[12].strip().split()[-1]) #Frequency in Hz

                isttime_line = lines[-1]
                date_line = lines[-2]
                ist = " ".join([date_line.strip().split()[-1], isttime_line.strip().split()[-1][:-3]])
                istdatetime = datetime.strptime(ist, "%d/%m/%Y %H:%M:%S.%f")
                mjd = getmjd(istdatetime)

                if int(freq) == 550000000:
                    band = 3
                elif int(freq) == 750000000:
                    band = 4
                elif int(freq) == 1460000000:
                    band = 5
                else:
                    raise ValueError(f"Unsupported frequency {freq/1e6} MHz. Supported frequencies are 550, 750, and 1460 MHz.")

                pc_ra = float(pc_ra_line.strip().split()[-1])
                pc_dec = float(pc_dec_line.strip().split()[-1])
                beam_per_host = int(beams_per_host_line.strip().split()[-1])
                nbeams = int(nbeams_line.strip().split()[-1])

                ra_dec_lines = lines[28:28+beam_per_host]
                temp_df = pd.DataFrame([line.strip().split() for line in ra_dec_lines], columns=data_columns)
                ahdr_data = pd.concat([ahdr_data, temp_df], ignore_index=True)
        else:
            print(f"File {file} not found!")

    ahdr_data = ahdr_data.apply(lambda col: pd.to_numeric(col, errors='coerce'))
    return ahdr_data, src_name, pc_ra, pc_dec, nbeams, band, mjd

def src_data_from_text(filename):
    """
    Read alternative source list and precess coordinates to MJD. This file has ra dec in HHhMMmSSs and DDdMMmSSs format.
    """
    cols = ['Src_name', 'RA', 'DEC', 'Epoch', 'Ang_dist']
    df = pd.read_csv(filename, sep=' ', names=cols, skipinitialspace=True)

    coords = SkyCoord(df['RA'].values,
                      df['DEC'].values,
                      frame='icrs')

    df['RA']  = coords.ra.to_value(u.rad)
    df['DEC'] = coords.dec.to_value(u.rad)
    return df

def beam_pattern_plot(ahdr_data, src_name, band, src_ra_dec, pc_ra, pc_dec, output_dir, df):
    '''
    Plots beam pattern from above merged dataframe. Each beam is annotated with its beam number.
    Saves the plot to output directory.
    '''
    ra = ahdr_data["RA"]
    dec = ahdr_data["DEC"]
    beam_numbers = ahdr_data["BM-Idx"]

    fig, ax = plt.subplots(figsize=(10,8), constrained_layout=True)
    scatter = ax.plot(ra, dec, 'o', markersize=5, alpha=1)
    ax.plot(pc_ra, pc_dec, '*', markersize=10, label="Phase Centre", color='red')
    if src_ra_dec is not None:
        ax.plot(src_ra_dec.ra.rad, src_ra_dec.dec.rad, '+', markersize=15, label='Precessed Target Coord', color='forestgreen', alpha=1)
    ax.scatter(df['RA'], df['DEC'], s=80, marker='*', color='green', label='Infield Sources')   

    # Annotate each point with its beam number
    for i in range(len(ra)):
        ax.annotate((beam_numbers[i]), (ra[i], dec[i]),
            textcoords="offset points", xytext=(0, 4), ha='right', fontsize=6, color='blue', weight='bold')

    ax.set_xlabel('Right Ascension')
    ax.set_ylabel('Declination')
    ax.set_title('Beam Pattern plot')
    plt.legend()
    output_path = os.path.join(output_dir, f"Infield_BeamPattern_{src_name}_B{band}.png")
    plt.savefig(output_path)
    plt.show()

def infield_src_plot(ahdr_data, src_name, band, src_ra_dec, pc_ra, pc_dec, output_dir, df):
    '''
    Plots beam pattern from above merged dataframe. Each beam is annotated with its beam number.
    Saves the plot to output directory.
    '''
    ra = ahdr_data["RA"]
    dec = ahdr_data["DEC"]
    beam_numbers = ahdr_data["BM-Idx"]

    fig, ax = plt.subplots(figsize=(10,8), constrained_layout=True)
    ax.scatter(df['RA'], df['DEC'], s=80, marker='*', color='green', label='Infield Sources')
    ax.plot(ra, dec, 'o', markersize=5, color='k')
    ax.plot(pc_ra, pc_dec, '*', markersize=10, label="Phase Centre", color='red')
    if src_ra_dec is not None:
        ax.plot(src_ra_dec.ra.rad, src_ra_dec.dec.rad, '+', markersize=15, label='Precessed Target Coord', color='forestgreen', alpha=1)
    
    for _, row in df.iterrows():
        ax.text(
            row['RA'],
            row['DEC'],
            row['Src_name'],
            fontsize=8,
            ha='left',
            va='bottom',
            color='purple'
        )    

    # Annotate each point with its beam number
    for i in range(len(ra)):
        ax.annotate((beam_numbers[i]), (ra[i], dec[i]),
            textcoords="offset points", xytext=(0, 4), ha='right', fontsize=10, color='blue', weight='bold')
     
    max_dist_ra = max(df['RA']) - min(df['RA'])
    max_dist_dec = max(df['DEC']) - min(df['DEC'])

    ax.set_xlabel('Right Ascension')
    ax.set_ylabel('Declination')
    ax.set_title('Infield sources in FOV')
    ax.set_xlim(min(df['RA']) - max_dist_ra, max(df['RA']) + max_dist_ra)
    ax.set_ylim(min(df['DEC']) - max_dist_dec, max(df['DEC']) + max_dist_dec)
    plt.legend()
    output_path = os.path.join(output_dir, f"BeamPattern_{src_name}_B{band}.png")
    plt.savefig(output_path)
    plt.show()

def main():
    parser = argparse.ArgumentParser(description="Plot beam pattern from AHDR files")
    parser.add_argument("--ahdr_files", required=True, nargs="+", help="List of AHDR files")
    parser.add_argument("--src_list", required=True, help="Path to source list text file")
    parser.add_argument("--src_ra", help="Source Right Ascension in HH:MM:SS") #optional
    parser.add_argument("--src_dec", help="Source Declination in DD:MM:SS") #optional
    args = parser.parse_args()

    infield_data = src_data_from_text(args.src_list)
    ahdr_data, src_name, pc_ra, pc_dec, nbeams, band, mjd = ra_dec_from_ahdr(args.ahdr_files)
    if args.src_ra is None or args.src_dec is None:
        print("Target source RA/DEC not provided. Skipping plotting target coordinates...")
        src_ra_dec = None
    else:
        src_ra_dec = get_precessed_coords(args.src_ra, args.src_dec, mjd)
        print("precessed_src_ra_dec:", src_ra_dec)


    print(f"\nList of Infield Sources in the FOV: \n{infield_data}")
    ahdr_dir_path = os.path.dirname(args.ahdr_files[0])
    beam_pattern_plot(ahdr_data, src_name, band, src_ra_dec, pc_ra, pc_dec, ahdr_dir_path, infield_data)
    infield_src_plot(ahdr_data, src_name, band, src_ra_dec, pc_ra, pc_dec, ahdr_dir_path, infield_data)
    print(f"\nSaved beam pattern and infield sources plots at {ahdr_dir_path}")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"Error: {e}")