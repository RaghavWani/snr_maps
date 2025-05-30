# This program plots the SNR map for a verified source (burst) by tracing the candidate from the output of the pipeline in the raw (filterbank files)
# Look at README for more info on the script, inputs, outputs and the configuration file

# AUTHOR: Vyaas

# Essential imports
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from astropy import units as u
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
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


def snr_plot(header_df, cands, traced_h5files_path, source_ra, source_dec):
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
    
    fig = uplt.figure(width=7.5, height=5)
    ax = fig.subplot()
    sm = ax.scatter(df["RA"], df["DEC"], vmin=0.0, c=df["SNRsp"], cmap="viridis", edgecolor="black", label="Beams", markersize=150)
    ax.scatter(snrmaxsp["RA"], snrmaxsp["DEC"], c="red", markersize=50, label=f"Beam with maximum SNR")
    ax.scatter(source_ra, source_dec, c="red", marker="*", markersize=50, label="Phase center")
    ax.colorbar(sm)
    ax.legend(loc="top")
    ax.format(suptitle="Single Pulse SNR Map (Raw Data)")
    fig.savefig(os.path.join(traced_h5files_path, "rSNR.png"), dpi=150)
    uplt.show()


def main():
    """
    Main function to plot the beamwise-SNR maps from the RAW data files
    """

    # reading the config file
    config = snr_plt_utils.load_config("config-h5.yaml")
    scan_id = config.get("scan_id")
    observation_path = config.get("observation_path")
    dm_thresh = config.get("dm_thresh")
    time_thresh = config.get("time_thresh")
    cand_h5file_path = config.get("cand_h5file_path")
    traced_h5files_path = config.get("output_dir")

    # copying the config file to the output directory
    shutil.copy("config-h5.yaml", os.path.join(traced_h5files_path, "run_config.yaml"))
    
    # computing the dataframe for plotting the SNR maps
    header_df, source_ra, source_dec, beams_per_node, num_beams, mjd = snr_plt_utils.ra_dec_from_ahdr(observation_path, scan_id)
    header_df = header_df.sort_values('BM-Idx')
    header_df = header_df.set_index('BM-Idx')

    # making the candydates.csv as an input to candies make
    make_candysv(traced_h5files_path, cand_h5file_path, observation_path, scan_id)

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
    snr_plot(header_df, cands, traced_h5files_path, source_ra, source_dec)

if __name__ == "__main__":
    main()
