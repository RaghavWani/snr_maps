# This program plots the SNR map for a verified source (burst) by tracing the candidate from the output of the pipeline in the raw (filterbank files)
# Look at README for more info on the script, inputs, outputs and the configuration file

# AUTHOR: Vyaas

# Essential imports
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from astropy import units as u
from astropy.coordinates import SkyCoord, TETE, Angle
from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel
from scipy.ndimage import gaussian_filter
from astropy.time import Time
from datetime import datetime

from candies.base import Candidate
from spyden import snratio, Template
import ultraplot as uplt

import h5py
import yaml

import os
from pathlib import Path
import subprocess
import warnings
import shutil

# SNR plotting utility module
import snr_plt_utils

warnings.filterwarnings('ignore')

def make_candysv(traced_h5files_path, cand_h5file_path, observation_path, scan_id):
    """
    This function creates the .csv file that candies takes as input to perform feature extraction
    """

    # this is the BeamID of the beam wherein the verified burst is in
    beam_number = int(cand_h5file_path.split("/")[-2][2:])
    beam_str = cand_h5file_path.split("/")[-2]
    burst_fil_name = (beam_str + ".fil")
    target_snr = float(cand_h5file_path.split("/")[-1].split("_")[-1][3:-3])

    # this is the filterbank file containing the verified burst
    burst_fil_path = os.path.join(observation_path, "FilData", scan_id, burst_fil_name)

    # this is the filtered_candidates.csv file for the beam containing the verified burst
    filtered_cands_csv_path = os.path.join(observation_path, "FRBPipeData", scan_id, beam_str, "filtered_candidates.csv")
    
    # making a similar .csv for candies to use to make the .h5 files
    fcands = pd.read_csv(filtered_cands_csv_path)
    matching_rows = fcands[np.isclose(fcands['snr'], target_snr, atol=1e-5)]
    candysv = pd.DataFrame([matching_rows.iloc[0]]).reset_index(drop=True)
    candysv = candysv.drop(candysv.columns[0], axis=1)
    
    # removing the redundant "BMx" part in the directory as .fil files are directly stored in the /FilData
    parts = candysv["file"].values[0].strip().split("/")
    if beam_str in parts:
        parts.remove(beam_str)
    new_path = "/".join(parts)
    candysv.at[0, "file"] = new_path

    if os.path.exists(os.path.join(traced_h5files_path, "candydates.csv")):
        print("A .csv file already exists in the output directory - using the existing file!")
    else:
        # writing the candydates.csv (input to candies make --fil) to the output dir
        print("Writing a new .csv to the output directory!")
        candysv.to_csv(f"{traced_h5files_path}/candydates.csv", index=False)


def load_cands_from_h5(traced_h5files_path, cand_h5file_path, numBeams):
    """
    This function loads the candidates from all the h5 files after running candies
    """

    cands = []
    cand_h5file_name = cand_h5file_path.split("/")[-1]
    for i in range(numBeams):
        cands.append([i, Candidate.load(f"{traced_h5files_path}/BM{i}/{cand_h5file_name}")])
    cands = pd.DataFrame(cands, columns=["BM-Idx", "Candidate"])
    cands.set_index("BM-Idx", drop=True, inplace=True)
    return cands


def snr_plot(header_df, cands, traced_h5files_path, source_ra, source_dec, mjd, cand_ra_dec):
    """
    This function plots the SNR map from the feature extracted .h5 files
    """

    df = pd.merge(header_df, cands, on="BM-Idx")
    # df = df.sort_values("BM-Idx", inplace=True)
    
    candsnrs = []
    candamps = []
    
    for i, row in df.iterrows():
        c = row["Candidate"]
    
        # Bandpass normalisation.
        dynspec = c.dedispersed.data
        mu = np.median(dynspec, axis=1)
        sigma = np.std(dynspec, axis=1)
        sigma[sigma == 0] = 1.0
        normalised = (dynspec - mu.reshape(-1, 1)) / sigma.reshape(-1, 1)
    
        # Get profile.
        prof = dynspec.sum(axis=0)
    
        # Calculate SNR.
        template = Template.boxcar(int(c.wbin))
        snrs, _, _, models = snratio(prof, template)
        model = models[0, :]
        snrs = snrs[0, 0, :]
        snr = snrs[snrs.argmax()]
    
        # Save SNR.
        candsnrs.append(snr)
        candamps.append(prof.max())
        
    df["SNRsp"] = candsnrs
    df["AMPsp"] = candamps
    snrmaxsp = df[df["SNRsp"] == df["SNRsp"].max()]

    # normalizing SNR
    df["SNRsp"] = (df["SNRsp"] - df["SNRsp"].min()) / (df["SNRsp"].max() - df["SNRsp"].min())

    # generating consistent x and y bounds
    window_rad = 0.00005  # adjust as needed
    ra_min, ra_max = header_df["RA"].min() - window_rad, header_df["RA"].max() + window_rad
    dec_min, dec_max = header_df["DEC"].min() - window_rad, header_df["DEC"].max() + window_rad

    fig = uplt.figure(width=7.5, height=5)
    ax = fig.subplot()
    sm = ax.scatter(df["RA"], df["DEC"], vmin=0.0, c=df["SNRsp"], cmap="viridis", edgecolor="black", label="__nolegend__", markersize=150)
    
    if cand_ra_dec:
        # plotting the precessed coords
        obstime = Time(mjd, format="mjd")
        coords = SkyCoord(cand_ra_dec, frame="icrs")
        tc = coords.transform_to(TETE(obstime=obstime))
        ax.scatter(tc.ra.rad, tc.dec.rad, c="grey", marker=".", markersize=50, label="Precessed coords")
    else:
        print("Not plotting precessed coords~!")

    ax.scatter(snrmaxsp["RA"], snrmaxsp["DEC"], c="red", markersize=50, label=f"Beam with maximum SNR")
    ax.scatter(source_ra, source_dec, c="red", marker="*", markersize=50, label="Phase center")
    ax.colorbar(sm)
    ax.legend(loc="top")
    ax.format(xlim=(ra_min, ra_max), ylim=(dec_min, dec_max), suptitle="Single Pulse SNR Map (Raw Data)")
    fig.savefig(os.path.join(traced_h5files_path, "rSNR.png"), dpi=150)
    uplt.show()

def spatial_snr_plot(header_df, cands, traced_h5files_path, source_ra, source_dec, mjd, cand_ra_dec):
    """
    This function plots the SNR map from the feature extracted .h5 files
    """

    df = pd.merge(header_df, cands, on="BM-Idx")
    # df = df.sort_values("BM-Idx", inplace=True)
    
    candsnrs = []
    candamps = []
    
    for i, row in df.iterrows():
        c = row["Candidate"]
    
        # Bandpass normalisation.
        dynspec = c.dedispersed.data
        mu = np.median(dynspec, axis=1)
        sigma = np.std(dynspec, axis=1)
        sigma[sigma == 0] = 1.0
        normalised = (dynspec - mu.reshape(-1, 1)) / sigma.reshape(-1, 1)
    
        # Get profile.
        prof = dynspec.sum(axis=0)
    
        # Calculate SNR.
        template = Template.boxcar(int(c.wbin))
        snrs, _, _, models = snratio(prof, template)
        model = models[0, :]
        snrs = snrs[0, 0, :]
        snr = snrs[snrs.argmax()]
    
        # Save SNR.
        candsnrs.append(snr)
        candamps.append(prof.max())
        
    df["SNRsp"] = candsnrs
    df["AMPsp"] = candamps
    snrmaxsp = df[df["SNRsp"] == df["SNRsp"].max()]

    # central_beam_index = snr_plt_utils.get_central_beam_no(header_df, source_ra, source_dec)
    
    n_pix = 256

    #ra_min, ra_max = df["RA"].min(), df["RA"].max()
    #dec_min, dec_max = df["DEC"].min(), df["DEC"].max()

    # generating consistent x and y bounds
    window_rad = 0.00005  # adjust as needed
    ra_min, ra_max = header_df["RA"].min() - window_rad, header_df["RA"].max() + window_rad
    dec_min, dec_max = header_df["DEC"].min() - window_rad, header_df["DEC"].max() + window_rad

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

    for _, row in df.iterrows():
        coord = SkyCoord(ra=row["RA"] * u.rad, dec=row["DEC"] * u.rad)
        x_pix, y_pix = skycoord_to_pixel(coord, w)
    
        x_pix = int(np.round(x_pix))
        y_pix = int(np.round(y_pix))
    
        if 0 <= x_pix < n_pix and 0 <= y_pix < n_pix:
            temp_image = np.zeros_like(snr_image)
            temp_image[y_pix, x_pix] = row["SNRsp"]
    
            smoothed = gaussian_filter(temp_image, sigma=kernel_sigma_pix)
            snr_image += smoothed

    # normalizing SNR
    snr_image = (snr_image - snr_image.min()) / (snr_image.max() - snr_image.min())

    # saving spatial SNR map as numpy array to load later
    np.save(os.path.join(output_dir, "rSNR_skymap.npy"), snr_image)

    fig = uplt.figure(width=7.5, height=5)
    ax = fig.subplot()
    skysnr = ax.imshow(snr_image, origin="lower", cmap="viridis", extent=[ra_min, ra_max, dec_min, dec_max])
    ax.scatter(source_ra, source_dec, c="red", marker="*", markersize=25, label="Phase center")
    
    if cand_ra_dec:
        # plotting the precessed coords
        obstime = Time(mjd, format="mjd")
        coords = SkyCoord(cand_ra_dec, frame="icrs")
        tc = coords.transform_to(TETE(obstime=obstime))
        ax.scatter(tc.ra.rad, tc.dec.rad, c="grey", marker=".", markersize=25, label="Precessed coords")
    else:
        print("Not plotting precessed coords~!")
    
    ax.colorbar(skysnr, label='SNR')
    ax.legend(loc="top")
    ax.format(xlim=(ra_min, ra_max), ylim=(dec_min, dec_max), suptitle=f"Single Pulse Sky-SNR Map (Raw Data)") 
    ax.set_xlabel('Right Ascension')
    ax.set_ylabel('Declination')
    ax.set_aspect('auto')
    
    fig.savefig(os.path.join(traced_h5files_path, "rSNR_skymap.png"), dpi=150)


def main():
    """
    Main function to plot the beamwise-SNR maps from the RAW data files
    """

    # reading the config file
    config = snr_plt_utils.load_config("config.yaml")
    scan_id = config.get("scan_id")
    observation_path = config.get("observation_path")
    dm_thresh = config.get("dm_thresh")
    time_thresh = config.get("time_thresh")
    cand_h5file_path = config.get("cand_h5file_path")
    traced_h5files_path = config.get("output_dir")
    cand_ra_dec = config.get("cand_ra_dec")

    # copying the config file to the output directory
    shutil.copy("config.yaml", os.path.join(traced_h5files_path, "run_config.yaml"))
    
    # computing the dataframe for plotting the SNR maps
    header_df, source_ra, source_dec, beams_per_node, num_beams, mjd = snr_plt_utils.ra_dec_from_ahdr(observation_path, scan_id)
    header_df = header_df.sort_values('BM-Idx')
    header_df = header_df.set_index('BM-Idx')

    # making the candydates.csv as an input to candies make
    make_candysv(traced_h5files_path, cand_h5file_path, observation_path, scan_id)
    
    h5_files = list(Path(traced_h5files_path).rglob("*.h5"))
    if h5_files:
        raise RuntimeError(f"Found {len(h5_files)} .h5 file(s) in {traced_h5files_path}. Aborting script ~ CLEAR O/P DIRECTORY BEFORE RUNNING CANDIES MAKE")

    # running candies make on all the filterbank files for the scan provided using a bash script
    bash_script_path = str((Path.cwd() / "make_h5_candies.sh").resolve())
    try:
        subprocess.run(
            [bash_script_path, os.path.join(traced_h5files_path, "candydates.csv"), os.path.join(observation_path, "FilData", scan_id)],
            check=True,  # raises error if script fails
            cwd=os.path.join(traced_h5files_path)
        )
    except subprocess.CalledProcessError as e:
        print("..make_h5_candies.. script failed with error:", e)
    
    cands = load_cands_from_h5(traced_h5files_path, cand_h5file_path, num_beams)
    snr_plot(header_df, cands, traced_h5files_path, source_ra, source_dec, mjd, cand_ra_dec)
    spatial_snr_plot(header_df, cands, traced_h5files_path, source_ra, source_dec, mjd, cand_ra_dec)

if __name__ == "__main__":
    main()
