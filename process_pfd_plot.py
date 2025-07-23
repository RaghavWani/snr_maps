import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import yaml, os, re, pytz
from matplotlib.ticker import FormatStrFormatter
from matplotlib.collections import PatchCollection
from astropy import units as u
from astropy.coordinates import SkyCoord, TETE
from astropy.time import Time
from datetime import datetime
from matplotlib.patches import Ellipse
from tabulate import tabulate
import warnings, logging


warnings.filterwarnings('ignore')
logging.basicConfig(level="INFO", datefmt="[%X]", format="%(message)s")
log = logging.getLogger("snr_plot")

# Function to load YAML config file
def load_config(config_path="config.yaml"):
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file '{config_path}' not found. Please create one.")
    
    with open(config_path, "r") as file:
        try:
            config = yaml.safe_load(file)
        except yaml.YAMLError as e:
            raise ValueError(f"Error parsing the YAML file: {e}")
    
    return config

def getmjd(t: datetime):
    """
    Function to convert an IST Date/Time to MJD
    """

    localtz = pytz.timezone("Asia/Kolkata")
    localdt = localtz.localize(t, is_dst=None)
    utcdt = localdt.astimezone(pytz.utc)
    mjd = Time(utcdt).mjd
    return mjd

def get_precessed_coords(src_ra_dec, mjd):
    """
    Function to precess the source coordinates to the specified epoch.
    Input: src_ra (Right Ascension in degrees), src_dec (Declination in degrees), epoch (e.g., 'J2000').
    Returns: Precessed SkyCoord object.
    """
    src_ra_dec = src_ra_dec.split(" ")
    src_ra = src_ra_dec[0].split(":")
    src_ra = f"{src_ra[0]}h{src_ra[1]}m{src_ra[2]}s"

    src_dec = src_ra_dec[1].split(":")
    src_dec = f"{src_dec[0]}d{src_dec[1]}m{src_dec[2]}s"

    obstime = Time(mjd, format="mjd")
    coords = SkyCoord(src_ra, src_dec, frame="icrs")
    precessed_coord = coords.transform_to(TETE(obstime=obstime))

    return precessed_coord

# Function to get ra dec from ahdr files 
# [This function WILL NOT BE USED Once Ra Dec in filterbank issue is fixed]
def ra_dec_from_ahdr(directory_path):
    """
    Extracts RA, DEC, BM-Idx, and BM-SubIdx values from ahdr files in a directory.
    Input: Directory containing ahdr files.
    Returns: A DataFrame containing the extracted data: RA, DEC, BM-Idx, BM-SubIdx.
    """
    data_columns = ["RA", "DEC", "BM-Idx", "BM-SubIdx"]
    ahdr_data = pd.DataFrame(columns=data_columns)

    # List all .ahdr files in the directory
    ahdr_files = [os.path.join(directory_path, file) for file in os.listdir(directory_path) if file.endswith(".ahdr")]
    if not ahdr_files:
        print(f"No header files found in directory: {directory_path}")
        return None
    
    for file in ahdr_files:
        if os.path.exists(file):
            with open(file, "r") as infile:
                lines = infile.readlines()

                src_name = lines[7].strip().split()[-1]
                pc_ra_line = lines[8]
                pc_dec_line = lines[9]
                nbeams_line = lines[25]
                beams_per_host_line = lines[26]
                freq = float(lines[12].strip().split()[-1]) #Frequency in Hz

                isttime_line = lines[-1]
                date_line = lines[-2]
                ist = " ".join([date_line.strip().split()[-1], isttime_line.strip().split()[-1][:-3]])
                istdatetime = datetime.strptime(ist, "%d/%m/%Y %H:%M:%S.%f")
                mjd = getmjd(istdatetime)

                if int(freq) == 500000000:
                    band = 3
                elif int(freq) == 550000000:
                    band = 4
                elif int(freq) == 1460000000:
                    band = 5
                else:
                    raise ValueError(f"Unsupported frequency {freq} MHz. Supported frequencies are 550, 750, and 1460 MHz.")
                
                pc_ra = float(pc_ra_line.strip().split()[-1])
                pc_dec = float(pc_dec_line.strip().split()[-1])
                beam_per_host = int(beams_per_host_line.strip().split()[-1])
                nbeams = int(nbeams_line.strip().split()[-1])

                ra_dec_lines = lines[28:28+beam_per_host]
                temp_df = pd.DataFrame([line.strip().split() for line in ra_dec_lines], columns=data_columns)
                ahdr_data = pd.concat([ahdr_data, temp_df], ignore_index=True)
        else:
            print(f"File {file} not found!")

    ahdr_data = ahdr_data.apply(lambda col: pd.to_numeric(col, errors='coerce'))
    return ahdr_data, src_name, pc_ra, pc_dec, nbeams, band, mjd

# Function to get PSR code from source name
def get_psr_code(src_name):
    match = re.match(r'^[A-Z](\d{4}[+-]\d{4})$', src_name)
    if match:
        return f"PSR_{match.group(1)}"
    else:
        raise ValueError("Invalid code format")
    
# Snr data from pfd files
def extract_snr(pfd_dir, ahdr_data, nbeams, src_name):
    '''
    Extracts beam index and snr from folding outputs *.pfd.bestprof files 
    Returns a new dataframe RA, DEC, BM-Idx, BM-SubIdx, SNR
    '''
    data_columns = ["BM-Idx", "SNR"]
    pfd_data = pd.DataFrame(columns=data_columns)

    for i in range(nbeams):
        filename = f"L{i}aa_{get_psr_code(src_name)}.pfd.bestprof"
        file_path = os.path.join(pfd_dir, filename)
        
        if os.path.isfile(file_path) and filename.endswith('.pfd.bestprof'):
            with open(file_path, "r") as f:
                lines = f.readlines()
            
                selected_line = lines[13]
                
                temp_df = pd.DataFrame({
                    'BM-Idx': [i],
                    'SNR': [np.float64(selected_line.split("(")[2].split("~")[1].split(" ")[0])]                      
                })
                pfd_data = pd.concat([pfd_data, temp_df], ignore_index=True)
        else:
            raise ValueError(f"No pfd file found for Beam {i}!")
        
    merged_df = pd.merge(ahdr_data, pfd_data,on='BM-Idx', how='inner')
    return merged_df

# Beam Pattern plotting
def beam_pattern_plot(merged_df, src_name, band, output_dir):
    '''
    Plots beam pattern from above merged dataframe. Each beam is annotated with its beam number.
    Saves the plot to output directory.
    '''
    ra = merged_df['RA']
    dec = merged_df['DEC']
    beam_numbers = merged_df['BM-Idx']
    
    fig, ax = plt.subplots(figsize=(10,8), constrained_layout=True)
    scatter = ax.plot(ra, dec, 'o',markersize=5, alpha=1) 
    # Annotate each point with its beam number
    for i in range(len(ra)):
        ax.annotate((beam_numbers[i]), (ra[i], dec[i]),
                    textcoords="offset points", xytext=(0, 4), ha='center', fontsize=6, color='blue', weight='bold')
    
    ax.set_xlabel('Right Ascension')
    ax.set_ylabel('Declination')
    ax.set_title('Beam Pattern plot')
    output_path = os.path.join(output_dir, f"BeamPattern_{src_name}_B{band}.png")
    plt.savefig(output_path)
    plt.show()

# Folded SNR plotting
def fold_snr_plot(merged_df, src_name, pc_ra, pc_dec, band, a, b, angle, output_dir, normal, saveif, src_ra_dec=None):
    '''
    Plots SNR scatter plot for merged dataframe. Saves the data and folded snr map to output directory.
    '''
    merged_df.sort_values(by=['BM-Idx'],ascending=True,inplace=True)
    merged_df.reset_index(drop=True, inplace=True)
    
    if normal:
       #Normalization to [0,1] range
       merged_df['SNR'] -= min(merged_df['SNR']) #Note negative sign
       merged_df['SNR'] /= max(merged_df['SNR'])

    #Highest SNR beam-
    high_snr_ra = merged_df.loc[merged_df['SNR'].idxmax(),'RA']
    high_snr_dec = merged_df.loc[merged_df['SNR'].idxmax(),'DEC']
    high_snr_bm = merged_df.loc[merged_df['SNR'].idxmax(),'BM-Idx']

    # Source closest beam:
    df_coords = SkyCoord(merged_df['RA'].values*u.rad, merged_df['DEC'].values *u.rad, frame='icrs')
    src_coord = SkyCoord(src_ra_dec.ra.rad*u.rad, src_ra_dec.dec.rad*u.rad, frame='icrs')
    closest_idx = df_coords.separation(src_coord).argmin()
    src_closest_bm = merged_df.iloc[closest_idx]['BM-Idx']

    #SNR Map from data:
    fig, ax = plt.subplots(figsize=(7.5,5), constrained_layout=True)
    ellipses=[Ellipse(xy=(x,y), width=a, height=b, angle = angle, edgecolor='blue', fill=True)
            for (x,y) in zip(merged_df['RA'], merged_df['DEC'])]
    e_col=PatchCollection(ellipses)
    e_col.set(array=merged_df['SNR'],cmap='viridis')
    ax.add_collection(e_col)
    plt.colorbar(e_col)
    
    if src_ra_dec:
        ax.plot(src_ra_dec.ra.rad, src_ra_dec.dec.rad, '+', markersize=15, label='Precessed Pulsar Coord', color='forestgreen', alpha=1)
    
    # scatter = ax.scatter(merged_df['RA'], merged_df['DEC'], c=merged_df['SNR'], s=150, cmap='viridis', edgecolors='black', alpha=1)
    # fig.colorbar(scatter, ax=ax, label='SNR')
    ax.plot(pc_ra, pc_dec, '*', markersize=10, label="Phase Centre", color='red')
    ax.plot(high_snr_ra, high_snr_dec, 'o', markersize=3, label=f"Max SNR Beam {high_snr_bm}", color='k')
    
    ax.set_xlabel('Right Ascension (rad)')
    ax.set_ylabel('Declination (rad)')
    ax.set_title(f'Folded SNR Map: {src_name}, Band {band}')
    plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%.4f'))
    plt.legend()

    output_path = os.path.join(output_dir, f"SNRMapFOLD_{src_name}_B{band}.png")
    plt.savefig(output_path)
    plt.show()

    if saveif:
        csv_path = os.path.join(output_dir, f"SNRMapFOLD_{src_name}_B{band}.csv")
        merged_df.to_csv(csv_path, index=False)
    return merged_df, src_closest_bm

# Simulated SNR plotting
def sim_snr_plot(sim_file, src_name, pc_ra, pc_dec, band, a, b, angle, output_dir, normal, src_ra_dec):
    '''
    Plots SNR scatter plot for simulation data. Saves the simulated snr map to output directory.
    '''

    # Simulation SNR map data:
    sim_sm = pd.read_csv(sim_file, delim_whitespace=True, header=None) 
    sim_sm.columns = ['DEC', 'RA', 'SNR', 'BM-Idx']

    sim_sm.sort_values(by=['BM-Idx'],ascending=True,inplace=True)
    sim_sm.reset_index(drop=True, inplace=True)
    
    if normal:
       #Shift and renormalization to [0,1] range
       sim_sm['SNR'] -= min(sim_sm['SNR']) #Note negative sign
       sim_sm['SNR'] /= max(sim_sm['SNR'])
    
    #arcsec to radian [REQUIRED IF SIMULATION DATA IS IN ARCSEC]
    # sim_sm['RA'] = sim_sm['RA'].apply(lambda x: (x * u.arcsec).to(u.rad).value)
    # sim_sm['DEC'] = sim_sm['DEC'].apply(lambda x: (x * u.arcsec).to(u.rad).value)

    #Linear transformation (ra dec shift): [REQUIRED IF SIMULATION DATA IS NOT CENTERED AROUND PHASE CENTER]
    # sim_sm['RA'] += pc_ra 
    # sim_sm['DEC'] += pc_dec

    #Highest SNR beam:
    high_snr_ra = sim_sm.loc[sim_sm['SNR'].idxmax(),'RA']
    high_snr_dec = sim_sm.loc[sim_sm['SNR'].idxmax(),'DEC']
    high_snr_bm = sim_sm.loc[sim_sm['SNR'].idxmax(),'BM-Idx']

    #SNR Map from simulation:
    fig, ax = plt.subplots(figsize=(7.5,5), constrained_layout=True)
    ellipses=[Ellipse(xy=(x,y), width=a, height=b, angle = angle, edgecolor='blue', fill=True)
            for (x,y) in zip(sim_sm['RA'], sim_sm['DEC'])]
    e_col=PatchCollection(ellipses)
    e_col.set(array=sim_sm['SNR'],cmap='viridis')
    ax.add_collection(e_col)
    plt.colorbar(e_col)
    
    if src_ra_dec:
        ax.plot(src_ra_dec.ra.rad, src_ra_dec.dec.rad, '+', markersize=15, label='Precessed Pulsar Coord', color='forestgreen', alpha=1)

    # scatter = ax.scatter(sim_sm['DEC'], sim_sm['RA'], c=sim_sm['SNR'], s=150, cmap='viridis', edgecolors='black', alpha=1)
    # fig.colorbar(scatter, ax=ax, label='SNR')
    ax.plot(pc_ra, pc_dec, '*', markersize=10, label="Phase Centre", color='red')
    ax.plot(high_snr_ra, high_snr_dec, 'o', markersize=3, label=f"Max SNR Beam {high_snr_bm}", color='k')

    ax.set_xlabel('Right Ascension (rad)')
    ax.set_ylabel('Declination (rad)')
    ax.set_title(f'Simulation SNR Map: {src_name}, Band {band}')
    plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%.4f'))
    plt.legend()

    output_path = os.path.join(output_dir, f"SNRMapSIM_{src_name}_B{band}.png")
    plt.savefig(output_path)
    plt.show()

    return sim_sm

# Residual SNR plotting
def residual_plot(fold_sm, sim_sm, src_name, src_bm, pc_ra, pc_dec, band, a, b, angle, output_dir, src_ra_dec):
    '''
    Plots residual scatter plot for simulation and folded data.
    '''

    residual = fold_sm.copy()
    residual['SNR'] = abs(fold_sm['SNR'] - sim_sm['SNR']) # Residual calculation

    # Residual Details:
    res_details = [
    ["Source Name", src_name],
    ["Source RA (rad)", "{:.7f}".format(src_ra_dec.ra.rad)],
    ["Source DEC (rad)", "{:.7f}".format(src_ra_dec.dec.rad)],
    ["Phase Centre RA (rad)", "{:.7f}".format(pc_ra)],
    ["Phase Centre DEC (rad)", "{:.7f}".format(pc_dec)],
    ["Source closest BM Idx", src_bm],
    ["Max Residual SNR", f"{max(residual['SNR']):.3f} (at BM {fold_sm['BM-Idx'][residual['SNR'].idxmax()]})"],
    ["Min Residual SNR", f"{min(residual['SNR']):.3f} (at BM {fold_sm['BM-Idx'][residual['SNR'].idxmin()]})"],
    ["Residual at Source Coord", "{:.3f}".format(residual[residual['BM-Idx'] == src_bm]['SNR'].iloc[0])],
    ]

    print("\nResidual Details:")
    print(tabulate(res_details, tablefmt="plane"))
    
    #Residual Plot:
    fig, ax = plt.subplots(figsize=(7.5,5), constrained_layout=True)
    ellipses=[Ellipse(xy=(x,y), width=a, height=b, angle = angle, edgecolor='blue', fill=True)
            for (x,y) in zip(sim_sm['RA'], sim_sm['DEC'])]
    e_col=PatchCollection(ellipses)
    e_col.set(array=residual['SNR'],cmap='viridis')
    ax.add_collection(e_col)
    plt.colorbar(e_col)

    if src_ra_dec:
        ax.plot(src_ra_dec.ra.rad, src_ra_dec.dec.rad, '+', markersize=15, label='Precessed Pulsar Coord', color='forestgreen', alpha=1)

    # scatter = ax.scatter(sim_sm['DEC'], sim_sm['RA'], c=residual['SNR'], s=150, cmap='viridis', edgecolors='black')
    # fig.colorbar(scatter, ax=ax, label='SNR')
    ax.plot(pc_ra, pc_dec, '*', markersize=8, label="Phase Centre", color='red')
    ax.set_xlabel('Right Ascension (rad)')
    ax.set_ylabel('Declination (rad)')
    ax.set_title(f'SNR Residual Plot: {src_name}, Band {band}')
    plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%.4f'))
    plt.legend()

    output_path = os.path.join(output_dir, f"SNRMapRES_{src_name}_B{band}.png")
    plt.savefig(output_path)
    plt.show()

# Main function
def main():
    try:
        config = load_config()
        
        # Details from config file
        src_ra_dec = config.get("src_ra_dec")  # Expected format: "RA DEC" (e.g., "07:42:48.9513465 -28:22:44.30356")
        normal = config.get("normal")
        saveif = config.get("saveif")
        a = config.get("ellipse_a")
        b = config.get("ellipse_b")
        angle = config.get("ellipse_angle")
        
        header_dir_path = config.get("header_dir_path")
        pfd_dir_path = config.get("pfd_dir_path")
        sim_file_path = config.get("sim_file_path")
        output_dir_path = config.get("output_dir_path")
        
        a = (a * u.arcsec).to(u.rad).value
        b = (b * u.arcsec).to(u.rad).value

        # Main script
        if not isinstance(normal, bool):
            raise ValueError("Config error: 'normal' must be either True or False (boolean).")

        if not isinstance(saveif, bool):
                    raise ValueError("Config error: 'saveif' must be either True or False (boolean).")

        ahdr_data, src_name, pc_ra, pc_dec,  nbeams, band, mjd = ra_dec_from_ahdr(header_dir_path)
        log.info(f"Processing PFD files from {pfd_dir_path} with {nbeams} beams...")
        if src_ra_dec:
            src_ra_dec = get_precessed_coords(src_ra_dec, mjd)
        else:
            raise ValueError("Please provide pulsar coordinates from ATNF PSR Catalogue in the format 'RA DEC' (e.g., '07:42:48.9513465 -28:22:44.30356')")

        new_data = extract_snr(pfd_dir_path, ahdr_data, nbeams, src_name)

        log.info(f"Plotting Beam Pattern...")
        beam_pattern_plot(new_data,src_name, band,output_dir_path)

        log.info(f"Plotting Folded SNR Map...")
        fold_sm, src_bm = fold_snr_plot(new_data, src_name, pc_ra, pc_dec, band, a, b, angle, output_dir_path, normal, saveif, src_ra_dec)

        log.info(f"Plotting Simulated SNR Map...")
        sim_sm = sim_snr_plot(sim_file_path, src_name, pc_ra, pc_dec, band, a, b, angle, output_dir_path, normal, src_ra_dec)

        if normal != True:
            print("Residual plot requires normalized data. Set 'normal' to True in config file.")
        else:
            log.info(f"Plotting Residual SNR Map...")
            residual_plot(fold_sm, sim_sm, src_name, src_bm, pc_ra, pc_dec, band,  a, b, angle, output_dir_path, src_ra_dec)

        print(f"\nAll SNR Map plots and folded data saved to {output_dir_path}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main() 