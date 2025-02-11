import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import FormatStrFormatter
import warnings
warnings.filterwarnings('ignore')

# Folded SNR map data:
band = 3
src_name = 'B2154+40' #Source name

#Pointing center coordinates from header file
pc_ra =  5.751332  
pc_dec = 0.705683 

data1 = f"SnrMap_{src_name}_band{band}_data.csv" #folding snr map
fold_sm = pd.read_csv(data1,delimiter=',')
fold_sm.sort_values(by=['BM-Idx'],ascending=True,inplace=True)
fold_sm.drop(columns=['I'],inplace=True)
fold_sm.reset_index(drop=True, inplace=True)

#Highest SNR RA DEC used as source coordinates
src_ra = fold_sm[fold_sm['SNR'] == 1]['RA']
src_dec = fold_sm[fold_sm['SNR'] == 1]['DEC']

print('Max SNR:',fold_sm['SNR'].max(),', Min SNR:',fold_sm['SNR'].min())
#Shift and renormalization
fold_sm['SNR'] -= min(fold_sm['SNR']) #Note negative sign
fold_sm['SNR'] /= max(fold_sm['SNR'])
print('Renormalizing...')
print('Max SNR:',fold_sm['SNR'].max(),', Min SNR:',fold_sm['SNR'].min())

#SNR Map from data:
fig,ax = plt.subplots(figsize=(5,4))
scatter = ax.scatter(fold_sm['RA'], fold_sm['DEC'], c=fold_sm['SNR'], s=50, cmap='viridis', edgecolors='black', alpha=1)
fig.colorbar(scatter, ax=ax, label='SNR')
ax.plot(src_ra, src_dec, 'o',markersize=2, label='Pulsar Coord', color='red')
ax.set_xlabel('Right Ascension (rad)')
ax.set_ylabel('Declination (rad)')
ax.set_title(f'Folded SNR Map: {src_name}, Band {band}')
plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%.4f'))
# plt.grid()
plt.legend()
plt.show()
