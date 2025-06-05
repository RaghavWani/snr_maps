# Simple code to visualize the beam pattern from the header files from the SPOTLIGHT Pipeline

import snr_plt_utils

def main():
    
    # reading the config file
    config = snr_plt_utils.load_config("config.yaml")
    scan_id = config.get("scan_id")
    observation_path = config.get("observation_path")
    output_dir = config.get("output_dir")

    header_df, source_ra, source_dec, beams_per_node, num_beams = snr_plt_utils.ra_dec_from_ahdr(observation_path, scan_id)

    snr_plt_utils.plot_beam_pattern(header_df, source_ra, source_dec, output_dir)

if __name__ == "__main__":
    main()
