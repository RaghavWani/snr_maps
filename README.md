# SNR Map
Repository containing beamwise Signal-to-Noise Ratio (SNR) Map plotting code that validates the true candidates and their SNR distribution over sky coordinates by visuallizing the SNR across multiple beams. This code is a part of the SPOTLIGHT FRB detection pipeline.

## Usage
To run the code, use the following commands:

* For `process_pfd_plot.py`:

Fill all the observational parameters correctly in the configuration file `config.yaml` and then execute the code as follows:
`python process_pfd_plot.py`

* For `folding.sh`:

Modify the paths for Pulsar parameter file and directory for beam wise filterbank files

## Description
*  `process_pfd_plot.py` \
This code is similar to `process_h5.py`, with the significant difference being that this code extracts SNR information from Pulsar folding. As a result, it works only for pfd files obtained from _prepfold_, a tool part of the _PRESTO_ software package. _PRESTO_ ([Pulsar REsearch Software Toolkit](http://www.cv.nrao.edu/~sransom/presto/)) is widely used for processing and analyzing pulsar data. Below is a brief outline of the code: 

1. Extarcts SNR information for individual beams from the `.pfd.bestprof` file and stores it as pandas data frame - RA, DEC, Beam Index, SNR (NOTE: The RA DEC values are taken from header files as pfd file doesn't have this information) </li>
2. Extracts other observation related details from the configuration file `config.yaml`</li>
3. Also, simulated SNR Map data (`.txt`) generated from the beam simulation code must be provided in order to plot Simulated and Residual SNR maps. </li>
4. Following plots are generated sequentially: 
    <ol>
        <li> The Beam formation Map </li>
        <li> The Folded SNR Map </li>
        <li> The Simulated SNR Map </li>
        <li> Residual results and Residual SNR Map </li>
    </ol>
    </li>
5. Note that the configuration file should be in the same directory as `process_pfd_plot.py`</li>
</ol>

* `process_h5.py` \
This code is similar to `centralized_process_h5.py` but with additional changes to DM-Time filtering instead of only DM filtering. As this code uses HDF5 files as input, it works only for data processed using the FRB detection pipeline. Below is a brief outline of the code: 
1. Converts HDF5 files into a pandas data frame by extracting information like - RA, DEC, Beam Index, SNR, Time of arrival, and DM. (NOTE: Due to some issue with RA DEC values in filterbank files, the code takes these values from ahdr files. This will be modified once the issue with the filterbank file is resolved)
2. Given the time tolerance, this data frame is divided into several "groups", each having approximately the same time of arrival.
3. The code then does DM filtering and processes only bursts that have their DM variation across all beams within the expected range of DM distribution of the observing source.
4. The SNR Map is then plotted for these DM-Time sorted bursts.


