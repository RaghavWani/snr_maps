# Python script to convert GMRT tile pickle file to RA, Dec, SNR, and BM index format
# Author: Jayaram Chengalur, Mekhala Muley

import sys
import getopt
import numpy as np
from   astropy import units as u
from   astropy.coordinates import SkyCoord
import pickle
import matplotlib.pyplot as plt
import os
 
script_dir=os.path.abspath('/lustre_archive/apps/correlator/code/BMSTEERING-NEW/at/')
sys.path.append(script_dir)
from make_tile import GMRTTile,TileType

def do_tile2radec():
    """
    Function converts the output format from the make_tile.py pickle format to the 
    text format and saves it in file (DEC_rad, RA_rad, unnormalized_SNR, BM_Indx ).
    
    Usage:
     python tile2radec.py -f gmrt_tile.pkl
    """    
    fname="gmrt_tile.pkl"
    #print(fname)
    try:
        opts,args=getopt.getopt(sys.argv[1:],"f:h",["file=","help"])
    except getopt.GetoptError as err:
        print(err)
        help(do_tile2radec)
        sys.exit(1)
    for opt, arg in opts:
        if opt in ("-f","--file"):
            fname=arg
        elif opt in ("-h","--help"):
            help(do_tile2radec)
            sys.exit(0)

    ofname="sim_snr.txt"
    of=open(ofname,"w")
    i=0
    
    f=open(fname,"rb");
    gtile=pickle.load(f);
    f.close()
    nb=gtile.n_beam
    crd=SkyCoord(gtile.beam_crd[:nb])
    amp=gtile.grid[:nb,2]   
    
    for c in crd:
        line=str(c.dec.radian)+"\t"+str(c.ra.radian)+"\t"+str(amp[i])+"\t"+str(i)+"\n"
        i=i+1
        of.write(line)
    of.close()

do_tile2radec()    
