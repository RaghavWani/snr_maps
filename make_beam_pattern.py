# Simple code to generate the beam pattern from the header files

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from astropy import units as u
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
from datetime import datetime

from candies.base import Candidate
from spyden import snratio, Template
from priwo import *
import ultraplot as uplt

import h5py
import argparse
import yaml

import os
from pathlib import Path
import logging
import warnings

import snr_plt_utils

def main():
    
    # reading the config file
    config = snr_plt_utils.load_config("config.yaml")
    scan_id = config.get("scan_id")
    observation_path = config.get("observation_path")
    dm_thresh = config.get("dm_thresh")
    time_thresh = config.get("time_thresh")
    output_dir = config.get("output_dir")
    cand_h5file_path = config.get("cand_h5file_path")
    dm_ver, toa_ver = snr_plt_utils.get_burst_dm_toa(cand_h5file_path)

    header_df, source_ra, source_dec, beams_per_node, num_beams = snr_plt_utils.ra_dec_from_ahdr(observation_path, scan_id)

    snr_plt_utils.plot_beam_pattern(header_df, source_ra, source_dec, output_dir)

if __name__ == "__main__":
    main()
