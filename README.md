# coincidence_filter
Repository containing coincidence - anti-coincidence burst filtering code, that validate classified burst candidates across multiple beams. This filter is a part of SPOTLIGHT FRB detection pipeline.

### centralized_process_h5.py
This code is one of the initial development in building coincidence filter. The results are meaningful only when observing known source (such as pulsar) and beams synthesized in a way that the source lie in the central beam. The data can be processed using many ways - real-time FRB detection pipeline, folding, Pulsar search pipeline, or single pulse search - but as this code takes h5 files as input, the data should be processed only using FRB detection pipeline. For other methods, the development is under progress. Below is brief outline of code: 
1. Converts h5 files into pandas dataframe, by extracting information like - RA, DEC, Beam Index, SNR, Time of arrival, and DM. (NOTE: Due to some issue with RA DEC values in filterbank files, the code takes these values from ahdr files. This will be modified once issue with filterbank file is resolved)
2. Given the central beam index, the code looks for highest SNR burst from central beam and gets its time of arrival.
3. The code then fetches the bursts that had this time of arrival across all beams.
4. The SNR Map is then plot for these sorted burst, but first ensuring it has DM variation within expected range of DM distribution. 
