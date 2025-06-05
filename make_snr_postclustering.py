import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from astropy import units as u
from astropy.coordinates import Angle, SkyCoord, TETE
from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel
from scipy.ndimage import gaussian_filter
from astropy.time import Time
from datetime import datetime

import ultraplot as uplt

import os
from pathlib import Path
import logging
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
 

def snr_plot(header_df, df, new_df, dm_tol, toA, dm_ver, source_ra, source_dec, output_dir, mjd, cand_ra_dec):
    """
    Plots SNR scatter plot for each group
    Input: List of grouped DataFrames, DM tolerance
    Output: SNR scatter plot for each group
    """
    new_df = dm_filter(new_df, dm_tol, dm_ver)

    central_beam_index = snr_plt_utils.get_central_beam_no(header_df, source_ra, source_dec)
    snrmaxsp = new_df[new_df["SNR"] == new_df["SNR"].max()]

    fig = uplt.figure(width=7.5, height=5)
    ax = fig.subplot()
    #scatter_header_df = ax.scatter(header_df['RA'], header_df['DEC'], vmin=0.0, c=pd.Series(0, index=header_df['RA'].index), cmap='viridis', edgecolor="black", label="__nolegend__", markersize=150)
    scatter = ax.scatter(new_df['RA'], new_df['DEC'], vmin=0.0, c=new_df['SNR'], cmap='viridis', edgecolor='black', label="__nolegend__", markersize=150)
    
    if cand_ra_dec:
        # plotting the precessed coords
        obstime = Time(mjd, format="mjd")
        coords = SkyCoord(cand_ra_dec, frame="icrs")
        tc = coords.transform_to(TETE(obstime=obstime))
        ax.scatter(tc.ra.rad, tc.dec.rad, c="grey", marker=".", markersize=50, label="Precessed coords")
    else:
        print("Not plotting precessed coords~!")

    ax.scatter(snrmaxsp["RA"], snrmaxsp["DEC"], c="red", markersize=50, label=f"Beam with maximum SNR")
    ax.scatter(header_df[header_df['BM-Idx'] == central_beam_index]["RA"], header_df[header_df['BM-Idx'] == central_beam_index]["DEC"], c="red", marker="*", markersize=50, label="Phase center")
    ax.colorbar(scatter, label='SNR')
    ax.legend(loc="top")
    ax.format(suptitle=f"Single Pulse SNR Map (Post-Clustering), T={toA}")   
    ax.set_xlabel('Right Ascension')
    ax.set_ylabel('Declination')
    ax.set_aspect('auto')

    fig.savefig(os.path.join(output_dir, "SNR_postclustering.png"), dpi=150)

def spatial_snr_plot(header_df, df, new_df, dm_tol, toA, dm_ver, source_ra, source_dec, output_dir, mjd, cand_ra_dec):
    new_df = dm_filter(new_df, dm_tol, dm_ver)
    central_beam_index = snr_plt_utils.get_central_beam_no(header_df, source_ra, source_dec)
    
    n_pix = 256

    ra_min, ra_max = new_df["RA"].min(), new_df["RA"].max()
    dec_min, dec_max = new_df["DEC"].min(), new_df["DEC"].max()
    
    # adding artificial padding if necessary
    eps = 1e-6  # ~0.2 arcsec in radians

    if np.isclose(ra_max, ra_min):
        ra_min -= eps / 2
        ra_max += eps / 2

    if np.isclose(dec_max, dec_min):
        dec_min -= eps / 2
        dec_max += eps / 2

    ra_min_deg, ra_max_deg = np.degrees([ra_min, ra_max])
    dec_min_deg, dec_max_deg = np.degrees([dec_min, dec_max])
    
    ra_grid = np.linspace(ra_min, ra_max, n_pix)
    dec_grid = np.linspace(dec_min, dec_max, n_pix)

    snr_image = np.zeros((n_pix, n_pix))

    w = WCS(naxis=2)
    w.wcs.crpix = [n_pix // 2, n_pix // 2]
    w.wcs.cdelt = [(ra_max_deg - ra_min_deg) / n_pix, (dec_max_deg - dec_min_deg) / n_pix]
    w.wcs.crval = [(ra_max_deg + ra_min_deg) / 2, (dec_max_deg + dec_min_deg) / 2]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    kernel_sigma_pix = 8.0  # scale with SNR if desired

    for _, row in new_df.iterrows():
        coord = SkyCoord(ra=row["RA"] * u.rad, dec=row["DEC"] * u.rad)
        x_pix, y_pix = skycoord_to_pixel(coord, w)
    
        x_pix = int(np.round(x_pix))
        y_pix = int(np.round(y_pix))
    
        if 0 <= x_pix < n_pix and 0 <= y_pix < n_pix:
            temp_image = np.zeros_like(snr_image)
            temp_image[y_pix, x_pix] = row["SNR"]
    
            smoothed = gaussian_filter(temp_image, sigma=kernel_sigma_pix)
            snr_image += smoothed
    
    fig = uplt.figure(width=7.5, height=5)
    ax = fig.subplot()
    skysnr = ax.imshow(snr_image, origin="lower", cmap="viridis", extent=[ra_min, ra_max, dec_min, dec_max])
    ax.scatter(header_df[header_df['BM-Idx'] == central_beam_index]["RA"], header_df[header_df['BM-Idx'] == central_beam_index]["DEC"], c="red", marker="*", markersize=25, label="Phase center")
    
    if cand_ra_dec:
        # plotting the precessed coords
        obstime = Time(mjd, format="mjd")
        coords = SkyCoord(cand_ra_dec, frame="icrs")
        tc = coords.transform_to(TETE(obstime=obstime))
        ax.scatter(tc.ra.rad, tc.dec.rad, c="grey", marker=".", markersize=50, label="Precessed coords")
    else:
        print("Not plotting precessed coords~!")
    
    ax.colorbar(skysnr, label='SNR')
    ax.legend(loc="top")
    ax.format(suptitle=f"Single Pulse Sky-SNR Map (Post-Clustering), T={toA}")
    ax.set_xlabel('Right Ascension')
    ax.set_ylabel('Declination')
    ax.set_aspect('auto')

    fig.savefig(os.path.join(output_dir, "SNR_postclustering_skymap.png"), dpi=150)

def main():
    """
    Main function to plot the beamwise-SNR maps from the data files
    Inputs: directory path of h5 files, time tolerance, dm tolerance
    Returns: SNR scatter plot for each group
    """

    # reading the config file
    config = snr_plt_utils.load_config("config.yaml")
    scan_id = config.get("scan_id")
    observation_path = config.get("observation_path")
    dm_thresh = config.get("dm_thresh")
    time_thresh = config.get("time_thresh")
    output_dir = config.get("output_dir")
    cand_h5file_path = config.get("cand_h5file_path")
    cand_ra_dec = config.get("cand_ra_dec")

    # copying the config file to the output directory
    shutil.copy("config.yaml", os.path.join(output_dir, "run_config.yaml"))

    dm_ver, toa_ver = snr_plt_utils.get_burst_dm_toa(cand_h5file_path)

    # computing the dataframe for plotting the SNR maps
    header_df, source_ra, source_dec, beams_per_node, num_beams, mjd = snr_plt_utils.ra_dec_from_ahdr(observation_path, scan_id)
    df = snr_plt_utils.clustering_to_df(observation_path, scan_id, header_df)
    new_df, new_df_2, toA = fetch_highest_snr(header_df, df, source_ra, source_dec, time_thresh, toa_ver)

    # plotting the SNR map
    snr_plot(header_df, df, new_df_2, dm_thresh, toA, dm_ver, source_ra, source_dec, output_dir, mjd, cand_ra_dec)
    spatial_snr_plot(header_df, df, new_df_2, dm_thresh, toA, dm_ver, source_ra, source_dec, output_dir, mjd, cand_ra_dec)

if __name__ == "__main__":
    main()
