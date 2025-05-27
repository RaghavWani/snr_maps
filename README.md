# SNR Map
Repository containing beamwise Signal-to-Noise Ratio (SNR) Map plotting code that validates the true candidates and their SNR distribution over sky coordinates by visuallizing the SNR across multiple beams. This code is a part of the SPOTLIGHT FRB detection pipeline.

## Usage
To run the code, use the following commands:

* For Folding - `folding.sh`:

Modify the paths for Pulsar parameter file and the directory for beam wise filterbank files. Run the shell script as `./folding.sh`

* For Plotting - `process_pfd_plot.py`:

Fill all the observational parameters correctly in the configuration file `config.yaml` and then execute the code as `python process_pfd_plot.py`. Make sure to source the correct Python environment that includes all the required libraries.

## Description
*  `process_pfd_plot.py` \
This script extracts SNR (Signal-to-Noise Ratio) information from pulsar folding results for all observed beams and plots the SNR distribution across sky coordinates. It is specifically designed for `.pfd` files obtained from _prepfold_, a tool from the _PRESTO_ ([Pulsar REsearch Software Toolkit](http://www.cv.nrao.edu/~sransom/presto/)) software package, which is widely used for processing and analyzing pulsar data. Below is a brief workflow overview: 

1. Extarcts SNR values for individual beams from the `.pfd.bestprof` file and stores them in a pandas DataFrame containing RA, DEC, Beam Index, and SNR (_Note: The RA DEC values are read from associated header files, as `.pfd` file doesn't include this information_) </li>
2. Reads observation details and plotting settings from the configuration file `config.yaml`</li>
3. Requires simulated SNR Map data (`.txt`) generated from the beam simulation code to enable comparison between observed and simulated maps. </li>
4. Following plots are created sequentially: 
    <ul>
        <li> Beam formation Map </li>
        <li> Folded Profile SNR Map </li>
        <li> Folded Profile Data (optional)</li>
        <li> Simulated SNR Map </li>
        <li> Residual results and Residual SNR Map </li>
    </ul>
    </li>
5. Ensure that `config.yaml` is in the same directory as `process_pfd_plot.py`</li>
</ol>

