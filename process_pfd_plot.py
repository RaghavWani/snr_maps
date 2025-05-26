import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import logging
from matplotlib.ticker import FormatStrFormatter
from astropy import units as u
from astropy.coordinates import Angle
from tabulate import tabulate
import os
import re
import warnings
import yaml

warnings.filterwarnings('ignore')
logging.basicConfig(level="INFO", datefmt="[%X]", format="%(message)s")
log = logging.getLogger("snr_plot")

# Function to load YAML config file
def load_config(config_path="config.yaml"):
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file '{config_path}' not found. Please create one.")
    
    with open(config_path, "r") as file:
        try:
            config = yaml.safe_load(file)
        except yaml.YAMLError as e:
            raise ValueError(f"Error parsing the YAML file: {e}")
    
    return config

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
        print(f"No header files found in directory: {directory_path}")
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

# Function to get PSR code from source name
def get_psr_code(src_name):
    match = re.match(r'^[A-Z](\d{4}[+-]\d{4})$', src_name)
    if match:
        return f"PSR_{match.group(1)}"
    else:
        raise ValueError("Invalid code format")
    
# Snr data from pfd files
def extract_snr(pfd_dir, ahdr_data, nbeams, src_name):
    '''
    Extracts beam index and snr from folding outputs *.pfd.bestprof files 
    Returns a new dataframe RA, DEC, BM-Idx, BM-SubIdx, SNR
    '''
    data_columns = ["BM-Idx", "SNR"]
    pfd_data = pd.DataFrame(columns=data_columns)

    for i in range(nbeams):
        filename = f"L{i}aa_{get_psr_code(src_name)}.pfd.bestprof"
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

# Beam Pattern plotting
def beam_pattern_plot(merged_df, src_name, band, output_dir):
    '''
    Plots beam pattern from above merged dataframe. Each beam is annotated with its beam number.
    Saves the plot to output directory.
    '''
    ra = merged_df['RA']
    dec = merged_df['DEC']
    beam_numbers = merged_df['BM-Idx']
    
    fig, ax = plt.subplots(figsize=(10,8), constrained_layout=True)
    scatter = ax.plot(ra, dec, 'o',markersize=5, alpha=1) 
    # Annotate each point with its beam number
    for i in range(len(ra)):
        ax.annotate((beam_numbers[i]), (ra[i], dec[i]),
                    textcoords="offset points", xytext=(0, 4), ha='center', fontsize=6, color='blue', weight='bold')
    
    ax.set_xlabel('Right Ascension')
    ax.set_ylabel('Declination')
    ax.set_title('Beam Pattern plot')
    output_path = os.path.join(output_dir, f"BeamPattern_{src_name}_B{band}.png")
    plt.savefig(output_path)
    plt.show()

# Folded SNR plotting
def fold_snr_plot(merged_df, src_name, band, pc_ra,pc_dec, output_dir, normal, saveif):
    '''
    Plots SNR scatter plot for merged dataframe. Saves the data and folded snr map to output directory.
    '''
    merged_df.sort_values(by=['BM-Idx'],ascending=True,inplace=True)
    merged_df.reset_index(drop=True, inplace=True)
    
    if int(normal) == 1:
       #Shift and renormalization
       merged_df['SNR'] -= min(merged_df['SNR']) #Note negative sign
       merged_df['SNR'] /= max(merged_df['SNR'])

    #Highest SNR RA DEC used as source coordinates
    src_ra = merged_df.loc[merged_df['SNR'].idxmax(),'RA']
    src_dec = merged_df.loc[merged_df['SNR'].idxmax(),'DEC']
    src_bm = merged_df.loc[merged_df['SNR'].idxmax(),'BM-Idx']

    #SNR Map from data:
    fig, ax = plt.subplots(figsize=(7.5,5), constrained_layout=True)
    scatter = ax.scatter(merged_df['RA'], merged_df['DEC'], c=merged_df['SNR'], s=150, cmap='viridis', edgecolors='black', alpha=1)
    fig.colorbar(scatter, ax=ax, label='SNR')
    ax.plot(pc_ra, pc_dec, '*',markersize=8, label="Phase Centre", color='red')
    ax.plot(src_ra, src_dec, 'o',markersize=2, label=f"Max SNR Beam {src_bm}", color='k')
    ax.set_xlabel('Right Ascension (rad)')
    ax.set_ylabel('Declination (rad)')
    ax.set_title(f'Folded SNR Map: {src_name}, Band {band}')
    plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%.4f'))
    plt.legend()
    output_path = os.path.join(output_dir, f"SNRMapFOLD_{src_name}_B{band}.png")
    plt.savefig(output_path)
    plt.show()

    if int(saveif) == 1:
       csv_path = os.path.join(output_dir, f"SNRMapFOLD_{src_name}_B{band}.csv")
       merged_df.to_csv(csv_path, index=False)
    return merged_df, src_ra, src_dec, src_bm

# Simulated SNR plotting
def sim_snr_plot(sim_file, src_name, src_ra, src_dec, band, output_dir, normal):
    '''
    Plots SNR scatter plot for simulation data. Saves the simulated snr map to output directory.
    '''

    # Simulation SNR map data:
    sim_sm = pd.read_csv(sim_file, delim_whitespace=True, header=None) 
    sim_sm.columns = ['RA', 'DEC', 'SNR', 'Beam-Idx']

    if int(normal) == 1:
       #Shift and renormalization
       sim_sm['SNR'] -= min(sim_sm['SNR']) #Note negative sign
       sim_sm['SNR'] /= max(sim_sm['SNR'])
    
    #arcsec to radian [REQUIRED IF SIMULATION DATA IS IN ARCSEC]
    # sim_sm['RA'] = sim_sm['RA'].apply(lambda x: (x * u.arcsec).to(u.rad).value)
    # sim_sm['DEC'] = sim_sm['DEC'].apply(lambda x: (x * u.arcsec).to(u.rad).value)

    #Linear transformation (ra dec shift): [REQUIRED IF SIMULATION DATA IS NOT CENTERED AROUND PHASE CENTER]
    # sim_sm['RA'] += pc_ra 
    # sim_sm['DEC'] += pc_dec

    #SNR Map from simulation:
    fig, ax = plt.subplots(figsize=(7.5,5), constrained_layout=True)
    scatter = ax.scatter(sim_sm['DEC'], sim_sm['RA'], c=sim_sm['SNR'], s=150, cmap='viridis', edgecolors='black', alpha=1)
    fig.colorbar(scatter, ax=ax, label='SNR')
    ax.plot(src_ra, src_dec, 'o',markersize=2, label='Pulsar Coord', color='red')
    ax.set_xlabel('Right Ascension (rad)')
    ax.set_ylabel('Declination (rad)')
    ax.set_title(f'Simulation SNR Map: {src_name}, Band {band}')
    plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%.4f'))
    plt.legend()

    output_path = os.path.join(output_dir, f"SNRMapSIM_{src_name}_B{band}.png")
    plt.savefig(output_path)
    plt.show()

    return sim_sm

# Residual SNR plotting
def residual_plot(fold_sm, sim_sm, src_name, src_ra, src_dec, src_bm, band, output_dir):
    '''
    Plots residual scatter plot for simulation and folded data.
    '''
    residual = fold_sm.copy()
    residual['SNR'] = abs(fold_sm['SNR'] - sim_sm['SNR']) # Residual calculation

    # Residual Details:
    res_details = [
    ["Source Name", src_name],
    ["Source RA (rad)", "{:.7f}".format(src_ra.iloc[0])],
    ["Source DEC (rad)", "{:.7f}".format(src_dec.iloc[0])],
    ["Source BM Idx", src_bm.iloc[0]],
    ["Max Residual SNR", f"{max(residual['SNR']):.3f} (at BM {fold_sm['BM-Idx'][residual['SNR'].idxmax()]})"],
    ["Min Residual SNR", f"{min(residual['SNR']):.3f} (at BM {fold_sm['BM-Idx'][residual['SNR'].idxmin()]})"],
    ["Residual at Source Coord", "{:.3f}".format(residual[residual['BM-Idx'] == src_bm.iloc[0]]['SNR'].iloc[0])],
    ]

    print("Residaul Details:")
    print(tabulate(res_details, tablefmt="plane"))
    
    #Residual Plot:
    fig, ax = plt.subplots(figsize=(7.5,5), constrained_layout=True)
    scatter = ax.scatter(sim_sm['DEC'], sim_sm['RA'], c=residual['SNR'], s=150, cmap='viridis', edgecolors='black', alpha=1)
    fig.colorbar(scatter, ax=ax, label='SNR')
    ax.plot(src_ra, src_dec, 'o',markersize=2, label='Pulsar Coord', color='red')
    ax.set_xlabel('Right Ascension (rad)')
    ax.set_ylabel('Declination (rad)')
    ax.set_title(f'SNR Residual Plot: {src_name}, Band {band}')
    plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%.4f'))
    plt.legend()

    output_path = os.path.join(output_dir, f"SNRMapRES_{src_name}_B{band}.png")
    plt.savefig(output_path)
    plt.show()

# Main function
def main():
    try:
        config = load_config()
        
        # Details from config file
        src_name = config.get("src_name")
        pfd_code = config.get("pfd_code")
        band = config.get("band")
        pc_ra = config.get("pc_ra")
        pc_dec = config.get("pc_dec")
        normal = config.get("normal")
        saveif = config.get("saveif")

        header_dir_path = config.get("header_dir_path")
        pfd_dir_path = config.get("pfd_dir_path")
        sim_file_path = config.get("sim_file_path")
        output_dir_path = config.get("output_dir_path")
        bph = config.get("beam_per_host")
        nbeams = config.get("nbeams")

        log.info(f"Processing PFD files from {pfd_dir_path} with {nbeams} beams...")

        #Main script
        ahdr_data = ra_dec_from_ahdr(header_dir_path,bph)
        new_data = extract_snr(pfd_dir_path,ahdr_data,nbeams, pfd_code)

        log.info(f"Plotting Beam Pattern...")
        beam_pattern_plot(new_data,src_name, band,output_dir_path)

        log.info(f"Plotting Folded SNR Map...")
        fold_sm, src_ra, src_dec, src_bm = fold_snr_plot(new_data, src_name, band, pc_ra, pc_dec, output_dir_path, normal, saveif)
        log.info(f"Plotting Simulated SNR Map...")
        sim_sm = sim_snr_plot(sim_file_path, src_name, src_ra, src_dec, band, output_dir_path, normal)
        
        log.info(f"Plotting Residual SNR Map...")
        residual_plot(fold_sm, sim_sm, src_name, src_ra, src_dec ,src_bm, band, output_dir_path)

        print(f"\nAll SNR Map plots and folded data saved to {output_dir_path}")
    
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main() 