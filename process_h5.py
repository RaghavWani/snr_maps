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
def ra_dec_from_ahdr(directory_path,beam_per_host):
    """
    Extracts RA, DEC, BM-Idx, and BM-SubIdx values from ahdr files in a directory.
    Input: Directory containing ahdr files.
    Returns: A DataFrame containing the extracted data: RA, DEC, BM-Idx, BM-SubIdx.
    """
    data_columns = ["RA", "DEC", "BM-Idx", "BM-SubIdx"]
    stacked_data = pd.DataFrame(columns=data_columns)

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
                stacked_data = pd.concat([stacked_data, temp_df], ignore_index=True)
        else:
            print(f"File {file} not found!")

    stacked_data = stacked_data.apply(lambda col: pd.to_numeric(col, errors='coerce'))
    return stacked_data

# Function to explore h5 file structure
def explore_h5(group, structure):
    for key, item in group.items():
        if isinstance(item, h5py.Group):
            structure[key] = {"type": "Group", "attributes": dict(item.attrs), "children": {}}
            explore_h5(item, structure[key]["children"])
        elif isinstance(item, h5py.Dataset):
            structure[key] = {
                "type": "Dataset",
                "attributes": dict(item.attrs),
                "shape": item.shape,
                "dtype": str(item.dtype),
            }

#RA DEC in radians:
# [This function WILL BE USED Once Ra Dec in filterbank issue is fixed]
def hdms_to_rad(hhmmss, coord): 
    '''
    Converts hhmmss to radians
    hhmmss: int, hhmmss format
    coord: str, 'RA' or 'DEC'
    Returns: Angle in radians
        
    '''
    try:
        hh = int(hhmmss/10000)
        mm = int(hhmmss/100 - hh*100)
        ss = float(hhmmss - hh*10000 - mm*100)
        
        if coord == "RA":
            angle = Angle(str(hh)+"h"+str(mm)+"m"+str(ss)+"s")
        elif coord == "DEC":
            angle = Angle(str(hh)+"d"+str(mm)+"m"+str(ss)+"s")
        else :
            raise ValueError("Invalid Coordinate name entered. Use 'RA' or 'DEC'.")
        return angle.to(u.radian).value
    except Exception as e:
        log.error(f"Error converting hhmmss to radians: {e}")
        return None

# Function to convert h5 files to DataFrame
def h5_to_dataframe(directory,header_dataframe):
    '''
    Converts h5 files information to dataframe
    Input: directory path of h5 files
    Output: DataFrame with columns RA, DEC, BM_Idx, DM, Time, SNR
    '''
    rows = []
    for filename in os.listdir(directory):
        file_structure = {}
        file_path = os.path.join(directory, filename)
        if os.path.isfile(file_path) and filename.endswith('.h5'):
            with h5py.File(file_path, "r") as h5_file:
                try:
                    beam_no = int(filename.split("_")[0].split("BM")[1])  # Adjust filename parsing
                    ra = float(header_dataframe['RA'][beam_no])
                    dec = float(header_dataframe['DEC'][beam_no])
                
                    # explore_h5(h5_file, file_structure)
                    # ra = hdms_to_rad(file_structure['extras']['attributes']['src_raj'], "RA")
                    # dec = hdms_to_rad(file_structure['extras']['attributes']['src_dej'], "DEC")
                    # print(beam_no, file_structure['extras']['attributes']['src_raj'], file_structure['extras']['attributes']['src_dej'])
                    
                    snr = h5_file.attrs.get('snr', None)
                    time = h5_file.attrs.get('t0', None)
                    dm = h5_file.attrs.get('dm', None)
                    rows.append({'RA': ra, 'DEC': dec, 'BM_Idx': beam_no, 'DM': dm, 'Time': time, 'SNR': snr})
                except Exception as e:
                    log.error(f"Error processing file {filename}: {e}")
    
    return pd.DataFrame(rows)

def time_filter(df,tolerance):
    '''
    Groups the dataframe based on time of arrival of burst (within tolerance)
    Input: DataFrame, tolerance in seconds
    Output: List of grouped DataFrames 
    '''
    
    df = df.sort_values('Time') # Sort by time of arrival

    groups = []
    current_group = [df.iloc[0]]
    for i in range(1, len(df)):
        if (df.iloc[i]['Time'] - current_group[-1]['Time']) <= tolerance:
            current_group.append(df.iloc[i])
        else:
            groups.append(pd.DataFrame(current_group))
            current_group = [df.iloc[i]]
    groups.append(pd.DataFrame(current_group))

    return groups

def dm_filter(df,tolerance):
    '''
    Checks if all DM values in time filtered group are approximately equal within the tolerance
    Input: DataFrame, tolerance
    Output: Boolean, True if all DM values are approximately equal, False otherwise
    '''

    # Check if all values are equal within the tolerance
    are_values_equal = np.abs(df['DM'].max() - df['DM'].min()) <= tolerance

    if are_values_equal:
        print(f"All values in 'DM' are approximately equal within the tolerance of {tolerance}.")
        return True
    else:
        print(f"Values in 'DM' vary beyond the tolerance of {tolerance}.")
        return False


def snr_plot(df,groups, dm_tol):
    '''
    Plots SNR scatter plot for each group
    Input: List of grouped DataFrames, DM tolerance
    Output: SNR scatter plot for each group
    '''
    i=0
    for idx, group in enumerate(groups):
        print(f"Group {idx + 1}:")
        
        if not dm_filter(group, dm_tol):
            continue  # Skip this group if DM filtering fails
        i+=1
        fig, ax = plt.subplots(figsize=(5, 4))
        scatter = ax.scatter(group['RA'], group['DEC'], c=group['SNR']/max(group['SNR']), s=50, cmap='plasma', edgecolors='black', alpha=1)
        fig.colorbar(scatter, ax=ax, label='SNR')
        
        ax.set_xlabel('Right Ascension')
        ax.set_ylabel('Declination')
        ax.set_title(f"SNR Scatter Plot, T={group['Time'].iloc[0]}")
        ax.set_aspect('auto')

        # beamnum = group.loc[group['SNR'].idxmax(), 'BM_Idx']
        # ax.annotate((beamnum).astype(int), (group.loc[group['SNR'].idxmax(), 'RA'], group.loc[group['SNR'].idxmax(), 'DEC']), 
        #         textcoords="offset points", xytext=(0, 4), ha='center', fontsize=8, color='blue', weight='bold')

        plt.text(0.5, 0.9, f"Number of candidates {len(group['SNR'])}", fontsize=12, ha='center', va='center', transform=ax.transAxes)
        plt.xlim(min(df['RA'])-0.0005, max(df['RA'])+0.0005)
        plt.ylim(min(df['DEC'])-0.0005, max(df['DEC'])+0.0005)
        plt.grid()
        #plt.savefig(f"SNR_Scatter_Plot_{idx + 1}.png")
        plt.show()
    print("All groups plotted.\n")
    print(f"Total {i} ToA groups passed DM filter.")
    
def main():
    '''
    Main function to process h5 files
    Input: directory path of h5 files, time tolerance, dm tolerance
    Output: SNR scatter plot for each group
    '''
    parser = argparse.ArgumentParser(prog=__file__)
    parser.add_argument("-D1", "--ahdr_dir_path", type=Path, required=True)
    parser.add_argument("-D2", "--h5_dir_path", type=Path, required=True)
    parser.add_argument("-bph", "--beam_per_host", type=int)
    parser.add_argument("-Tt", "--time_tol", type=float)
    parser.add_argument("-DMt", "--dm_tol", type=float, default=0.1)
    args = parser.parse_args()

    log.info(f"Plotting SNR Map for h5 files in directory: {args.h5_dir_path}")
    header_df = ra_dec_from_ahdr(args.ahdr_dir_path,args.beam_per_host)
    df = h5_to_dataframe(args.h5_dir_path,header_df)
    groups = time_filter(df,args.time_tol)
    print(f" Total {len(groups)} different ToA found.")
    print()
    snr_plot(df,groups, args.dm_tol)

if __name__ == "__main__":
    main() 
