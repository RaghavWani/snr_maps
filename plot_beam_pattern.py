# Python script to plot beam pattern from ahdr files
# Last updated: 30th July 2025; Author: Raghav Wani

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os, pytz
from astropy import units as u
from astropy.coordinates import SkyCoord, TETE
from astropy.time import Time
from datetime import datetime

ahdr_dir_path = "/lustre_data/spotlight/data/TEST_RT_20250813_203023/BeamData" 
src_ra_dec = "05:34:31.93357 +22:00:52.1927"

def getmjd(t: datetime):
    """
    Function to convert an IST Date/Time to MJD
    """

    localtz = pytz.timezone("Asia/Kolkata")
    localdt = localtz.localize(t, is_dst=None)
    utcdt = localdt.astimezone(pytz.utc)
    mjd = Time(utcdt).mjd
    return mjd

def get_precessed_coords(src_ra_dec, mjd):
    """
    Function to precess the source coordinates to the specified epoch.
    Input: src_ra (Right Ascension in degrees), src_dec (Declination in degrees), epoch (e.g., 'J2000').
    Returns: Precessed SkyCoord object.
    """
    src_ra_dec = src_ra_dec.split(" ")
    src_ra = src_ra_dec[0].split(":")
    src_ra = f"{src_ra[0]}h{src_ra[1]}m{src_ra[2]}s"

    src_dec = src_ra_dec[1].split(":")
    src_dec = f"{src_dec[0]}d{src_dec[1]}m{src_dec[2]}s"

    obstime = Time(mjd, format="mjd")
    coords = SkyCoord(src_ra, src_dec, frame="icrs")
    precessed_coord = coords.transform_to(TETE(obstime=obstime))

    return precessed_coord

def ra_dec_from_ahdr(directory_path):
    """
    Extracts RA, DEC, BM-Idx, and BM-SubIdx values from ahdr files in a directory.
    Input: Directory containing ahdr files.
    Returns: A DataFrame containing the extracted data: RA, DEC, BM-Idx, BM-SubIdx.
    """
    data_columns = ["RA", "DEC", "BM-Idx", "BM-SubIdx"]
    ahdr_data = pd.DataFrame(columns=data_columns)

    # List all .ahdr files in the directory
    ahdr_files = [os.path.join(directory_path, file) for file in os.listdir(directory_path) if file.endswith(".ahdr")]
    if not ahdr_files:
        print(f"No header files found in directory: {directory_path}")
        return None

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

                if int(freq) == 500000000:
                    band = 3
                elif int(freq) == 550000000:
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

def beam_pattern_plot(ahdr_data, src_name, band, src_ra_dec, pc_ra, pc_dec, output_dir):
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
    ax.plot(src_ra_dec.ra.rad, src_ra_dec.dec.rad, '+', markersize=15, label='Precessed Pulsar Coord', color='forestgreen', alpha=1)
    # Annotate each point with its beam number
    for i in range(len(ra)):
        ax.annotate((beam_numbers[i]), (ra[i], dec[i]),
            textcoords="offset points", xytext=(0, 4), ha='center', fontsize=6, color='blue', weight='bold')

    ax.set_xlabel('Right Ascension')
    ax.set_ylabel('Declination')
    ax.set_title('Beam Pattern plot')
    plt.legend()
    output_path = os.path.join(output_dir, f"BeamPattern_{src_name}_B{band}.png")
    plt.savefig(output_path)
    plt.show()

ahdr_data, src_name, pc_ra, pc_dec, nbeams, band, mjd = ra_dec_from_ahdr(ahdr_dir_path)
src_ra_dec = get_precessed_coords(src_ra_dec, mjd)
beam_pattern_plot(ahdr_data, src_name, band, src_ra_dec, pc_ra, pc_dec, ahdr_dir_path)
print(f"Saved beam pattern plot at {ahdr_dir_path}")
