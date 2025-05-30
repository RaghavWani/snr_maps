import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from astropy import units as u
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
from datetime import datetime

import ultraplot as uplt

import os
from pathlib import Path
import warnings
import shutil

import snr_plt_utils

warnings.filterwarnings('ignore')

def fetch_highest_snr(header_df, df, source_ra, source_dec, time_thresh, toa_ver):
    """
    Gets the ToA of the burst and return dataframe with details of beams having Time close to toA.
    Inputs: Header DataFrame, Candidates DF, Phase Center Coords, Threshold time
    Returns: New DataFrame with details of beams having Time close to ToA, ToA of the burst
    """

    central_beam_no = snr_plt_utils.get_central_beam_no(header_df, source_ra, source_dec)
    # print("The central beam index is: ", central_beam_no)
    central_beam_df = df[df['BM_Idx'] == central_beam_no]
    # toA = central_beam_df[central_beam_df['SNR'] == max(central_beam_df['SNR'])]['Time']
    # toA = df[df['SNR'] == max(df['SNR'])]['Time']
    
    # Setting SNR of beams that dont temporally coincide (time threshold) with the central burst to 0
    new_df = df.copy()
    # new_df.loc[new_df['Time'] != toA.iloc[0], 'SNR'] = 0
    mask = (new_df['Time'] < toa_ver - time_thresh) | (new_df['Time'] > toa_ver + time_thresh)
    new_df.loc[mask, 'SNR'] = 0
    # note: this dataframe contains all the beams which have candidates wherein atleast one candidate was detected within the range (not necessarily the burst)

    # Creating dataframe that contains all the beams where that burst was detected
    mask = (df['Time'] >= toa_ver - time_thresh) & (df['Time'] <= toa_ver + time_thresh)
    new_df_2 = df[mask]

    new_df_unique = new_df.sort_values('SNR')
    new_df_unique = new_df_unique.drop_duplicates(subset='BM_Idx', keep='last') # keeping the highest SNR cand in that beam if multiple cands present

    new_df_2_u = new_df_2.sort_values('SNR')
    new_df_2_u = new_df_2_u.drop_duplicates(subset='BM_Idx', keep='last')
    
    Bool_val = new_df_unique['BM_Idx'].duplicated().any() # Check if there are duplicate BM_Idx values
    
    if Bool_val:
        print(new_df_unique)
        raise ValueError("There are rows in the DataFrame with the same value of BM_Idx.")
    return new_df_unique, new_df_2_u, toa_ver

def dm_filter(df, tolerance, dm_ver):
    """
    Checks if all DM values in time filtered group are approximately equal within the DM tolerance
    Input: DataFrame, tolerance
    Output: Boolean, True if all DM values are approximately equal, False otherwise
    """
    # print(np.var(df['DM']))
    mask = (df['DM'] >= dm_ver - tolerance) & (df['DM'] <= dm_ver + tolerance)
    dm_fil_df = df[mask]
    print("DM filtering done!")
    return dm_fil_df
    

def snr_plot(header_df, df, new_df, dm_tol, toA, dm_ver, source_ra, source_dec, output_dir):
    """
    Plots SNR scatter plot for each group
    Input: List of grouped DataFrames, DM tolerance
    Output: SNR scatter plot for each group
    """
    new_df = dm_filter(new_df, dm_tol, dm_ver)

    central_beam_index = snr_plt_utils.get_central_beam_no(header_df, source_ra, source_dec)
    snrmaxsp = new_df[new_df["SNR"] == new_df["SNR"].max()]
    print(snrmaxsp)
    fig = uplt.figure(width=7.5, height=5)
    ax = fig.subplot()
    #scatter_header_df = ax.scatter(header_df['RA'], header_df['DEC'], vmin=0.0, c=pd.Series(0, index=header_df['RA'].index), cmap='viridis', edgecolor="black", label="__nolegend__", markersize=150)
    scatter = ax.scatter(new_df['RA'], new_df['DEC'], vmin=0.0, c=new_df['SNR'], cmap='viridis', edgecolor='black', label="Beams", markersize=150)
    ax.scatter(snrmaxsp["RA"], snrmaxsp["DEC"], c="red", markersize=50, label=f"Beam with maximum SNR")
    ax.scatter(header_df[header_df['BM-Idx'] == central_beam_index]["RA"], header_df[header_df['BM-Idx'] == central_beam_index]["DEC"], c="red", marker="*", markersize=50, label="Phase center")
    #ax.scatter(hdrs[0]["ra"], hdrs[0]["dec"], c="red", marker="*", markersize=300, label="Phase center")
    ax.colorbar(scatter, label='SNR')
    ax.legend(loc="top")
    ax.format(suptitle=f"Single Pulse SNR Map (Post-Classify), T={toA}")   
    ax.set_xlabel('Right Ascension')
    ax.set_ylabel('Declination')
    ax.set_aspect('auto')

    # beamnum = group.loc[group['SNR'].idxmax(), 'BM_Idx']
    # ax.annotate((beamnum).astype(int), (group.loc[group['SNR'].idxmax(), 'RA'], group.loc[group['SNR'].idxmax(), 'DEC']), 
    #         textcoords="offset points", xytext=(0, 4), ha='center', fontsize=8, color='blue', weight='bold')

    # plt.text(0.5, 0.9, f"Number of candidates {len(new_df['SNR'])}", fontsize=12, ha='center', va='top', transform=ax.transAxes)
    #plt.savefig(f"SNR_Scatter_Plot_{idx + 1}.png")
    fig.savefig(os.path.join(output_dir, "SNR_postclassify.png"), dpi=150)

def main():
    """
    Main function to plot the beamwise-SNR maps from the data files
    Inputs: directory path of h5 files, time tolerance, dm tolerance
    Returns: SNR scatter plot for each group
    """

    # reading the config file
    config = snr_plt_utils.load_config("config-h5.yaml")
    scan_id = config.get("scan_id")
    observation_path = config.get("observation_path")
    dm_thresh = config.get("dm_thresh")
    time_thresh = config.get("time_thresh")
    output_dir = config.get("output_dir")
    cand_h5file_path = config.get("cand_h5file_path")

    # copying the config file to the output directory
    shutil.copy("config-h5.yaml", os.path.join(traced_h5files_path, "run_config.yaml"))

    dm_ver, toa_ver = snr_plt_utils.get_burst_dm_toa(cand_h5file_path)

    # computing the dataframe for plotting the SNR maps
    header_df, source_ra, source_dec, beams_per_node, num_beams = snr_plt_utils.ra_dec_from_ahdr(observation_path, scan_id)
    df = snr_plt_utils.classification_to_df(observation_path, scan_id, header_df)
    new_df, new_df_2, toA = fetch_highest_snr(header_df, df, source_ra, source_dec, time_thresh, toa_ver)

    # plotting the SNR map
    snr_plot(header_df, df, new_df_2, dm_thresh, toA, dm_ver, source_ra, source_dec, output_dir)

if __name__ == "__main__":
    main()
