# Essential utilities for the functioning of the SNR map plotting code in the SPOTLIGHT pipeline

import os
import yaml
import numpy as np
import pandas as pd

from astropy import units as u
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
from astropy.time import Time
from datetime import datetime
import pytz

import h5py
import ultraplot as uplt

def load_config(config_path="config.yaml"):
    """
    Function to load the configuration file
    """

    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file '{config_path}' not found. Please create one.")
    
    with open(config_path, "r") as file:
        try:
            config = yaml.safe_load(file)
        except yaml.YAMLError as e:
            raise ValueError(f"Error parsing the YAML file: {e}")
    
    return config

def getmjd(t: datetime):
    """
    Function to convert an IST Date/Time to MJD
    """

    localtz = pytz.timezone("Asia/Kolkata")
    localdt = localtz.localize(t, is_dst=None)
    utcdt = localdt.astimezone(pytz.utc)
    mjd = Time(utcdt).mjd
    return mjd

def ra_dec_from_ahdr(observation_path, scan_id):
    """
    Extracts RA, DEC, BM-Idx, and BM-SubIdx values from ahdr files
    Inputs: Observation file directory
    Returns: A DataFrame containing the extracted data: RA, DEC, BM-Idx, BM-SubIdx from the .ahdr files and the target source's RA, DEC
    """

    # taking a default 10 beams per compute node and default 160 total beams (global var - will get overwritten from the beams per node from the .ahdr files)
    beams_per_node = 10
    num_beams = 160    

    # obtaining the ra, dec of the source (used for finding central beam later)
    source_ra = 0
    source_dec = 0
    
    data_columns = ["RA", "DEC", "BM-Idx", "BM-SubIdx"]
    stacked_data = pd.DataFrame(columns=data_columns)

    # list all .ahdr files in the given directory
    ahdr_files = [os.path.join(observation_path, "BeamData", file) for file in os.listdir(os.path.join(observation_path, "BeamData")) if (file.endswith(".ahdr") and file.startswith(scan_id))]
    if not ahdr_files:
        print(f"No .ahdr files found in directory: {ahdr_dir_path}")
        return None

    for file in ahdr_files:
        if os.path.exists(file):
            with open(file, "r") as infile:
                lines = infile.readlines()
                
                # reading some additional info other than beam ra, dec
                beams_per_node_line = lines[26]
                num_beams_line = lines[25]
                source_ra_line = lines[8]
                source_dec_line = lines[9]

                # reading the observation time in MJD
                isttime_line = lines[-1]
                date_line = lines[-2]
                ist = " ".join([date_line.strip().split()[-1], isttime_line.strip().split()[-1][:-3]])
                istdatetime = datetime.strptime(ist, "%d/%m/%Y %H:%M:%S.%f")
                mjd = snr_plt_utils.getmjd(istdatetime)
                
                source_ra = float(source_ra_line.strip().split()[-1])
                source_dec = float(source_dec_line.strip().split()[-1])
                beams_per_node = int(beams_per_node_line.strip().split()[-1])
                num_beams = int(num_beams_line.strip().split()[-1])
                
                # reading the required dataframe info from the file
                info_selected_lines = lines[28:28+beams_per_node]

                temp_df = pd.DataFrame([line.strip().split() for line in info_selected_lines], columns=data_columns)
                stacked_data = pd.concat([stacked_data, temp_df], ignore_index=True)
        else:
            print(f"File {file} not read!")

    stacked_data = stacked_data.apply(lambda col: pd.to_numeric(col, errors='coerce'))
    return stacked_data, source_ra, source_dec, beams_per_node, num_beams, mjd

def hdms_to_rad(hhmmss, coord): 
    """
    Converts hhmmss to radians
    Inputs: hhmmss: int, coord: str ('RA' or 'DEC')
    Returns: Angle in radians
    """
    
    try:
        hh = int(hhmmss/10000)
        mm = int(hhmmss/100 - hh*100)
        ss = float(hhmmss - hh*10000 - mm*100)
        
        if coord == "RA":
            angle = Angle(str(hh)+"h"+str(mm)+"m"+str(ss)+"s")
        elif coord == "DEC":
            angle = Angle(str(hh)+"d"+str(mm)+"m"+str(ss)+"s")
        else:
            raise ValueError("Invalid coordinate name entered. Use 'RA' or 'DEC'.")
        return angle.to(u.radian).value
    
    except Exception as e:
        log.error(f"Error converting hhmmss to radians: {e}")
        return None

def h5_to_dataframe(observation_path, scan_id, header_dataframe):
    """
    Converts information in all the h5 files present in FRBPipeData/BM* to a dataframe
    Input: directory path of observation, scan identifier (name of the folder inside /FRBPipeData), header dataframe (returned by 'ra_dec_from_ahdr()')
    Returns: DataFrame with columns RA, DEC, BM_Idx, DM, Time, SNR
    """

    rows = []
    h5files = []
    
    # navigating to the folder of that particualar scan in the observation
    scan_path = os.path.join(observation_path, "FRBPipeData", scan_id)
    # print(scan_path)
    if not os.path.isdir(scan_path):
        raise FileNotFoundError(f"H5 file directory '{scan_path}' does not exist. Please check the observation path in the .config file.")

    # going through the individual beam folders to extract the h5 files and appending them to rows
    for beam_folder in os.listdir(scan_path):
        if beam_folder.startswith("BM") and os.path.isdir(os.path.join(scan_path, beam_folder)):
            for filename in os.listdir(os.path.join(scan_path, beam_folder)):
                file_path = os.path.join(scan_path, beam_folder, filename)
                if os.path.isfile(file_path) and filename.endswith('.h5'):
                    h5files.append(file_path)
                    with h5py.File(file_path, "r") as h5_file:
                        try:
                            beam_no = int(beam_folder[2:])
                            ra = float(header_dataframe[header_dataframe['BM-Idx'] == beam_no]['RA'].iloc[0])
                            dec = float(header_dataframe[header_dataframe['BM-Idx'] == beam_no]['DEC'].iloc[0])
                            snr = h5_file.attrs.get('snr', None)
                            time = h5_file.attrs.get('t0', None)
                            dm = h5_file.attrs.get('dm', None)
                            rows.append({'RA': ra, 'DEC': dec, 'BM_Idx': beam_no, 'DM': dm, 'Time': time, 'SNR': snr})
                        except Exception as e:
                            log.error(f"Error processing file {filename}: {e}")

    return pd.DataFrame(rows)

def classification_to_df(observation_path, scan_id, header_dataframe):
    """
    Reads all the positively classified candidates' h5 files and converts it into a dataframe
    Inputs: Directory paths and header dataframe
    Returns: Dataframe with all positively classified candidates
    """

    rows = []

    # reading classification.csv to a df
    classification_path = os.path.join(observation_path, "FRBPipeData", scan_id, "classification.csv")
    if not (os.path.exists(classification_path)):
        print("classification.csv not found in the scan directory!")
        return
    classification_df = pd.read_csv(classification_path)

    # extracting only the label 1 candidates
    classification_df = classification_df[classification_df['labels'] == 1.0]

    # reading the .h5 files of the label 1 candidates
    for i, row in classification_df.iterrows():
        with h5py.File(row['candidate'], "r") as h5_file:
            try:
                beam_no = int(str(row['candidate']).strip().split("/")[-2][2:])
                ra = float(header_dataframe[header_dataframe['BM-Idx'] == beam_no]['RA'].iloc[0])
                dec = float(header_dataframe[header_dataframe['BM-Idx'] == beam_no]['DEC'].iloc[0])
                snr = h5_file.attrs.get('snr', None)
                time = h5_file.attrs.get('t0', None)
                dm = h5_file.attrs.get('dm', None)
                rows.append({'RA': ra, 'DEC': dec, 'BM_Idx': beam_no, 'DM': dm, 'Time': time, 'SNR': snr})
            except Exception as e:
                log.error(f"Error processing file {filename}: {e}")

    return pd.DataFrame(rows)

def candidates_to_df(observation_path, scan_id, header_dataframe):
    """
    Converts information in all the candidates.csv files present in FRBPipeData/BM* to a dataframe
    Input: directory path of observation, scan identifier (name of the folder inside /FRBPipeData), header dataframe (returned by 'ra_dec_from_ahdr()')
    Returns: DataFrame with columns RA, DEC, BM_Idx, DM, Time, SNR
    """

    rows = []
    csvfiles = []
    
    # navigating to the folder of that particualar scan in the observation
    scan_path = os.path.join(observation_path, "FRBPipeData", scan_id)
    # print(scan_path)
    if not os.path.isdir(scan_path):
        raise FileNotFoundError(f"Scan file directory '{scan_path}' does not exist. Please check the paths in the .config file.")

    # going through the individual beam folders to extract the candidates.csv files and appending them to rows
    for beam_folder in os.listdir(scan_path):
        if beam_folder.startswith("BM") and os.path.isdir(os.path.join(scan_path, beam_folder)):
            for filename in os.listdir(os.path.join(scan_path, beam_folder)):
                file_path = os.path.join(scan_path, beam_folder, filename)
                if os.path.isfile(file_path) and (filename == 'candidates.csv'):
                    csvfiles.append(file_path)
                    csv_df = pd.read_csv(file_path)
                    for i, row in csv_df.iterrows():
                        try:
                            beam_no = int(beam_folder[2:])
                            ra = float(header_dataframe[header_dataframe['BM-Idx'] == beam_no]['RA'].iloc[0])
                            dec = float(header_dataframe[header_dataframe['BM-Idx'] == beam_no]['DEC'].iloc[0])
                            snr = row['snr']
                            time = row['time']
                            dm = row['dm']
                            rows.append({'RA': ra, 'DEC': dec, 'BM_Idx': beam_no, 'DM': dm, 'Time': time, 'SNR': snr})
                        except Exception as e:
                            log.error(f"Error processing file {filename}: {e}")
            print(f"{beam_folder} done!")
                            
    return pd.DataFrame(rows)

def grpd_candidates_to_df(observation_path, scan_id, header_dataframe, toa_ver, dm_ver, time_thresh, dm_thresh):
    """
    Converts information in all the candidates.csv files present in FRBPipeData/BM* to a dataframe
    Input: directory path of observation, scan identifier (name of the folder inside /FRBPipeData), header dataframe (returned by 'ra_dec_from_ahdr()')
    Returns: DataFrame with columns RA, DEC, BM_Idx, DM, Time, SNR

    NOTE : this function does the DM-Time thresholding before appending to the candidates dataframe
    """

    rows = []
    csvfiles = []

    # navigating to the folder of that particualar scan in the observation
    scan_path = os.path.join(observation_path, "FRBPipeData", scan_id)
    # print(scan_path)
    if not os.path.isdir(scan_path):
        raise FileNotFoundError(f"Scan file directory '{scan_path}' does not exist. Please check the paths in the .config file.")

    # going through the individual beam folders to extract the candidates.csv files and appending them to rows
    for beam_folder in os.listdir(scan_path):
        if beam_folder.startswith("BM") and os.path.isdir(os.path.join(scan_path, beam_folder)):
            for filename in os.listdir(os.path.join(scan_path, beam_folder)):
                file_path = os.path.join(scan_path, beam_folder, filename)
                if os.path.isfile(file_path) and (filename == 'candidates.csv'):
                    csvfiles.append(file_path)
                    csv_df = pd.read_csv(file_path)
                    for i, row in csv_df.iterrows():
                        try:
                            beam_no = int(beam_folder[2:])
                            ra = float(header_dataframe[header_dataframe['BM-Idx'] == beam_no]['RA'].iloc[0])
                            dec = float(header_dataframe[header_dataframe['BM-Idx'] == beam_no]['DEC'].iloc[0])
                            snr = row['snr']
                            time = row['time']
                            dm = row['dm']
                            if (time <= toa_ver + time_thresh) and (time >= toa_ver - time_thresh) and (dm <= dm_ver + dm_thresh) and (dm >= dm_ver - dm_thresh):
                                rows.append({'RA': ra, 'DEC': dec, 'BM_Idx': beam_no, 'DM': dm, 'Time': time, 'SNR': snr})
                        except Exception as e:
                            log.error(f"Error processing file {filename}: {e}")
            print(f"{beam_folder} done!")

    return pd.DataFrame(rows)

def clustering_to_df(observation_path, scan_id, header_dataframe):
    """
    Converts information in all the filtered_candidates.csv files present in FRBPipeData/BM* to a dataframe
    Input: directory path of observation, scan identifier (name of the folder inside /FRBPipeData), header dataframe (returned by 'ra_dec_from_ahdr()')
    Returns: DataFrame with columns RA, DEC, BM_Idx, DM, Time, SNR
    """

    rows = []
    csvfiles = []
    
    # navigating to the folder of that particualar scan in the observation
    scan_path = os.path.join(observation_path, "FRBPipeData", scan_id)
    # print(scan_path)
    if not os.path.isdir(scan_path):
        raise FileNotFoundError(f"Scan file directory '{scan_path}' does not exist. Please check the paths in the .config file.")

    # going through the individual beam folders to extract the candidates.csv files and appending them to rows
    for beam_folder in os.listdir(scan_path):
        if beam_folder.startswith("BM") and os.path.isdir(os.path.join(scan_path, beam_folder)):
            for filename in os.listdir(os.path.join(scan_path, beam_folder)):
                file_path = os.path.join(scan_path, beam_folder, filename)
                if os.path.isfile(file_path) and (filename == 'filtered_candidates.csv'):
                    csvfiles.append(file_path)
                    csv_df = pd.read_csv(file_path)
                    for i, row in csv_df.iterrows():
                        try:
                            beam_no = int(beam_folder[2:])
                            ra = float(header_dataframe[header_dataframe['BM-Idx'] == beam_no]['RA'].iloc[0])
                            dec = float(header_dataframe[header_dataframe['BM-Idx'] == beam_no]['DEC'].iloc[0])
                            snr = row['snr']
                            time = row['stime']
                            dm = row['dm']
                            rows.append({'RA': ra, 'DEC': dec, 'BM_Idx': beam_no, 'DM': dm, 'Time': time, 'SNR': snr})
                        except Exception as e:
                            log.error(f"Error processing file {filename}: {e}")
                            
    return pd.DataFrame(rows)

def beamfilter_to_df(observation_path, scan_id, header_dataframe):
    """
    Reads all the beam-filtered candidates' h5 files and converts it into a dataframe
    Inputs: Directory paths and header dataframe
    Returns: Dataframe with all beam filtered candidates
    """

    rows = []

    # reading classification.csv to a df
    beamfilter_path = os.path.join(observation_path, "FRBPipeData", scan_id, "sorted_files_filtered.csv")
    if not (os.path.exists(beamfilter_path)):
        print("sorted_files_filtered.csv not found in the scan directory!")
        return
    bf_df = pd.read_csv(beamfilter_path, header=None)

    # reading the .h5 files
    for i, row in bf_df.iterrows():
        beam_str = row[0]
        h5filename = row[1]
        h5filepath = os.path.join(observation_path, "FRBPipeData", scan_id, beam_str, h5filename)
        with h5py.File(h5filepath, "r") as h5_file:
            try:
                beam_no = int(beam_str[2:])
                ra = float(header_dataframe[header_dataframe['BM-Idx'] == beam_no]['RA'].iloc[0])
                dec = float(header_dataframe[header_dataframe['BM-Idx'] == beam_no]['DEC'].iloc[0])
                snr = h5_file.attrs.get('snr', None)
                time = h5_file.attrs.get('t0', None)
                dm = h5_file.attrs.get('dm', None)
                rows.append({'RA': ra, 'DEC': dec, 'BM_Idx': beam_no, 'DM': dm, 'Time': time, 'SNR': snr})
            except Exception as e:
                log.error(f"Error processing file {filename}: {e}")

    return pd.DataFrame(rows)

def get_central_beam_no(header_df, source_ra, source_dec):
    """
    Gets the RA, DEC of the beam that is closest to the RA, DEC of the source observed
    Input: Header DataFrame
    Output: Beam Index of the central beam
    """

    # computing the angular distance using spherical coordinates
    coords = SkyCoord(ra=header_df["RA"].values*u.rad, dec=header_df["DEC"].values*u.rad)
    center = SkyCoord(ra=source_ra*u.rad, dec=source_dec*u.rad)
    separations = coords.separation(center)
    
    # computing the angular distance to the geometric center of the FoV (Euclidean Approximation)
    # ang_dist = (df["RA"] - ra_center)**2 + (df["DEC"] - dec_center)**2 
    # central_df_index = ang_dist.idxmin()

    # return df.index[separations.argmin()]
    return int(header_df['BM-Idx'][header_df.index[separations.argmin()]])

def get_burst_dm_toa(cand_h5file_path):
    """
    Gets the DM and ToA of the verified burst
    Input: .h5 file path of the verified candidate (in FRBPipeData)
    Outputs: DM and ToA of the verified candidate
    """
    
    cand_h5file_path.strip().split('/')[-1].split('_')
    ToA_ver = float(cand_h5file_path.strip().split('/')[-1].split('_')[1][1:])
    DM_ver = float(cand_h5file_path.strip().split('/')[-1].split('_')[2][2:])
    return DM_ver, ToA_ver

def plot_beam_pattern(header_df, source_ra, source_dec, output_dir):
    """
    Function to plot the observation-specific beam pattern
    """

    central_beam_index = get_central_beam_no(header_df, source_ra, source_dec)
    
    fig = uplt.figure(width=7.5, height=5)
    ax = fig.subplot()
    for i, row in header_df.iterrows():
        ax.text(row["RA"] - 1.5e-5, row["DEC"] - 1.5e-5, int(row["BM-Idx"]), transform="data", size=5)
    ax.scatter(header_df["RA"], header_df["DEC"], facecolor="none", edgecolor="black", markersize=100)
    ax.scatter(header_df[header_df['BM-Idx'] == central_beam_index]["RA"], header_df[header_df['BM-Idx'] == central_beam_index]["DEC"], c="red", marker="*", markersize=50, label="Phase center")
    fig.savefig(os.path.join(output_dir, "SNR_BeamPattern.png"), dpi=150)
