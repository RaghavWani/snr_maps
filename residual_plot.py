import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import FormatStrFormatter
from astropy import units as u
import warnings

warnings.filterwarnings('ignore')

# Observation details:
src_name = 'B2154+40'#Source name
src_ra = 5.749967  	#Source RA in radians
src_dec = 0.707463 #Source DEC in radians

pc_ra =  5.752066 #Pointing center RA in radians
pc_dec = 0.706426 #Pointing center DEC in radians
band = 4

# Folded SNR map data:
data1 = f"7Jan25/SnrMap_{src_name}_band{band}_data.csv" #folding snr map
fold_sm = pd.read_csv(data1,delimiter=',')
fold_sm.sort_values(by=['BM-Idx'],ascending=True,inplace=True)
fold_sm.drop(columns=['I'],inplace=True)
fold_sm.reset_index(drop=True, inplace=True)

print('Max:',fold_sm['SNR'].max(),', Min:',fold_sm['SNR'].min())

fold_sm['SNR'] -= min(fold_sm['SNR']) #Shift and renormalization
fold_sm['SNR'] /= max(fold_sm['SNR'])
print('Renormalizing...')
print('Max:',fold_sm['SNR'].max(),', Min:',fold_sm['SNR'].min()) #should be 1 and 0

# print(fold_sm[fold_sm['SNR'] == 1])


# Simulation SNR map data:
data2 = f"7Jan25/ra-dec_centers_snr_{src_name}_b{band}_offset.txt" #folding snr map
sim_sm = pd.read_csv(data2, delim_whitespace=True, header=None) 
sim_sm.columns = ['RA', 'DEC', 'SNR']

print('Max:',sim_sm['SNR'].max(),', Min:',sim_sm['SNR'].min())

sim_sm['SNR'] -= min(sim_sm['SNR']) #Shift and renormalization
sim_sm['SNR'] /= max(sim_sm['SNR'])
print('Renormalizing...')
print('Max:',sim_sm['SNR'].max(),', Min:',sim_sm['SNR'].min()) #should be 1 and 0

#arcsec to radian 
sim_sm['RA'] = sim_sm['RA'].apply(lambda x: (x * u.arcsec).to(u.rad).value)
sim_sm['DEC'] = sim_sm['DEC'].apply(lambda x: (x * u.arcsec).to(u.rad).value)

#Linear transformation (ra dec shift): replace the RA and DEC values with the phase center ra dec values
sim_sm['RA'] += pc_ra 
sim_sm['DEC'] += pc_dec

# Residual calculation
residual = abs(fold_sm['SNR'] - sim_sm['SNR'])

#SNR Map from data:
fig,ax = plt.subplots(figsize=(5,4))
scatter = ax.scatter(fold_sm['RA'], fold_sm['DEC'], c=fold_sm['SNR'], s=50, cmap='viridis', edgecolors='black', alpha=1)
fig.colorbar(scatter, ax=ax, label='SNR')
ax.plot(src_ra, src_dec, 'o',markersize=2, label='Pulsar Coord', color='red')
ax.set_xlabel('Right Ascension (rad)')
ax.set_ylabel('Declination (rad)')
ax.set_title(f'Folded SNR Map: {src_name}, Band {band}')
plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%.4f'))
plt.legend()
plt.show()

#SNR Map from simulation:
fig,ax = plt.subplots(figsize=(5,4))
scatter = ax.scatter(sim_sm['RA'], sim_sm['DEC'], c=sim_sm['SNR'], s=50, cmap='viridis', edgecolors='black', alpha=1)
fig.colorbar(scatter, ax=ax, label='SNR')
ax.plot(src_ra, src_dec, 'o',markersize=2, label='Pulsar Coord', color='red')
ax.set_xlabel('Right Ascension (rad)')
ax.set_ylabel('Declination (rad)')
ax.set_title(f'Simulation SNR Map: {src_name}, Band {band}')
plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%.4f'))
plt.legend()
plt.show()

#Residual Plot:
print(f"Max residual:\n{round(max(residual),2)}")
fig,ax = plt.subplots(figsize=(5,4))
scatter = ax.scatter(sim_sm['RA'], sim_sm['DEC'], c=residual, s=50, cmap='viridis', edgecolors='black', alpha=1)
fig.colorbar(scatter, ax=ax, label='SNR')
ax.plot(src_ra, src_dec, 'o',markersize=2, label='Pulsar Coord', color='red')
ax.set_xlabel('Right Ascension (rad)')
ax.set_ylabel('Declination (rad)')
ax.set_title(f'SNR Residual Plot: {src_name}, Band {band}')
plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%.4f'))
plt.legend()
plt.show()