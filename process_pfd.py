import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import logging
from astropy import units as u
from astropy.coordinates import Angle
import h5py
import os
from pathlib import Path
import warnings

warnings.filterwarnings('ignore')

logging.basicConfig(level="INFO", datefmt="[%X]", format="%(message)s")
log = logging.getLogger("snr_plot")

# Function to get ra dec from ahdr files 
# [This function WILL NOT BE USED Once Ra Dec in filterbank issue is fixed]
def ra_dec_from_ahdr(directory_path, beam_per_host):
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
        print(f"No .ahdr files found in directory: {directory_path}")
        return None
    
    for file in ahdr_files:
        if os.path.exists(file):
            with open(file, "r") as infile:
                lines = infile.readlines()
            
                selected_lines = lines[28:28+beam_per_host]
                
                temp_df = pd.DataFrame([line.strip().split() for line in selected_lines], columns=data_columns)
                ahdr_data = pd.concat([ahdr_data, temp_df], ignore_index=True)
        else:
            print(f"File {file} not found!")

    ahdr_data = ahdr_data.apply(lambda col: pd.to_numeric(col, errors='coerce'))
    return ahdr_data

# Snr data from pfd files
def extract_snr(pfd_dir,ahdr_data,nbeams):
    '''
    Extracts beam index and snr from folding outputs *.pfd.bestprof files 
    Returns a new dataframe RA, DEC, BM-Idx, BM-SubIdx, SNR
    '''
    data_columns = ["BM-Idx", "SNR"]
    pfd_data = pd.DataFrame(columns=data_columns)

    for i in range(nbeams):
        filename = f"L{i}aa_PSR_0332+5434.pfd.bestprof"
        file_path = os.path.join(pfd_dir, filename)
        
        if os.path.isfile(file_path) and filename.endswith('.pfd.bestprof'):
            with open(file_path, "r") as f:
                lines = f.readlines()
            
                selected_line = lines[13]
                
                temp_df = pd.DataFrame({
                    'BM-Idx': [i],
                    'SNR': [np.float64(selected_line.split("(")[2].split("~")[1].split(" ")[0])]                      
                })
                pfd_data = pd.concat([pfd_data, temp_df], ignore_index=True)
        else:
            raise ValueError(f"No pfd file found for Beam {i}!")
        
    merged_df = pd.merge(ahdr_data, pfd_data,on='BM-Idx', how='inner')
    return merged_df

# SNR plotting
def snr_plot(merged_df):
    '''
    Plots SNR scatter plot for merged dataframe
    '''
        
    fig, ax = plt.subplots(figsize=(5, 4))
    scatter = ax.scatter(merged_df['RA'], merged_df['DEC'], c=merged_df['SNR']/max(merged_df['SNR']), s=50, cmap='plasma', edgecolors='black', alpha=1)
    fig.colorbar(scatter, ax=ax, label='SNR')
    
    ax.set_xlabel('Right Ascension (rad)')
    ax.set_ylabel('Declination (rad)')
    ax.set_title("SNR Scatter Plot from folding")

    # beamnum = merged_df.loc[merged_df['SNR'].idxmax(), 'BM_Idx']
    # ax.annotate((beamnum).astype(int), (merged_df.loc[merged_df['SNR'].idxmax(), 'RA'], merged_df.loc[merged_df['SNR'].idxmax(), 'DEC']), 
    #         textcoords="offset points", xytext=(0, 4), ha='center', fontsize=8, color='blue', weight='bold')

    #plt.savefig(f"SNR_Scatter_Plot_{idx + 1}.png")
    plt.show()
    
def main():
    '''
    Main function to process h5 files
    Input: directory path of h5 files, time tolerance, dm tolerance
    Output: SNR scatter plot for each group
    '''
    parser = argparse.ArgumentParser(prog=__file__)
    parser.add_argument("-D1", "--ahdr_dir_path", type=Path, required=True)
    parser.add_argument("-D2", "--pfd_dir_path", type=Path, required=True)
    parser.add_argument("-bph", "--beam_per_host", type=int)
    parser.add_argument("-n", "--nbeams", type=int)
    args = parser.parse_args()

    log.info(f"Plotting SNR Map for pfd files in directory: {args.pfd_dir_path}")
    ahdr_data = ra_dec_from_ahdr(args.ahdr_dir_path,args.beam_per_host)
    new_data = extract_snr(args.pfd_dir_path,ahdr_data,args.nbeams)
    snr_plot(new_data)

if __name__ == "__main__":
    main() 