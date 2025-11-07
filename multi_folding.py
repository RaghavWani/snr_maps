####################################################################
# Python script to fold pulsar data using PRESTO's folding tool 
# namely "prepfold", which is developed by Scott Ransom. The script  
# parallelizes the folding task across CPU threads. 
# 
# Usage: 
#       python multi_folding.py -D /path/to/filterbanks -s "source_name"
# Output:
#	Folding outputs saved at /path/to/filterbanks/FOLDING_OUTPUTS
#
# NOTE: Use appropriate env before running this script:
# $ source /lustre_archive/apps/tdsoft/env.sh
#
#  Last Update: 07th Nov 2025; ~ Raghav Wani
####################################################################

import shutil
import os
import glob
import argparse
import subprocess
import time
import fcntl
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
from rich.console import Console
from rich_argparse import RichHelpFormatter

console = Console()

LOG_FILE_NAME = "folding_log.txt"

def safe_log_write(log_path, message):
    with open(log_path, "a") as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        f.write(message + "\n")
        fcntl.flock(f, fcntl.LOCK_UN)

def run_folding(fil_path, par_file, fold_dir, log_path=None):
    """Run prepfold on a single .fil file."""
    fil_path = Path(fil_path)
    cmd = [
        "prepfold",
        "-timing", str(par_file),
        str(fil_path),
        "-noxwin",
        "-o", str(fold_dir / fil_path.stem)
    ]
    subprocess.run(cmd, check=True)

    # Log completion
    if log_path:
        safe_log_write(log_path, f"{fil_path.name} Processed")

def main():
    parser = argparse.ArgumentParser(
        prog='multi_folding',
        description="[bold cyan] Pulsar folding using PRESTO's prepfold[/bold cyan]",
        epilog="[gold3]For any issues, feel free to contact me (Raghav)[/gold3]",
        formatter_class=RichHelpFormatter
    )
    parser.add_argument("-D", "--input_dir", type=Path, required=True, help="Absolute Path to directory containing filterbanks (REQUIRED)")
    parser.add_argument("-s", "--src_name", type=str, required=True, help="Pulsar name (REQUIRED)")
    parser.add_argument("-w", "--workers", type=int, default=4, help="Max parallel workers (DEFAULT=4)")
    args = parser.parse_args()
    
    fil_files = sorted(glob.glob(str(args.input_dir / "*.fil")))

    #fil_files = [f for f in fil_files]

    print("""
      _______  _______  _        ______  _________ _        _______ 
     (  ____ \(  ___  )( \      (  __  \ \__   __/( (    /|(  ____ \ 
     | (    \/| (   ) || (      | (  \  )   ) (   |  \  ( || (    \/
     | (__    | |   | || |      | |   ) |   | |   |   \ | || |      
     |  __)   | |   | || |      | |   | |   | |   | (\ \) || | ____ 
     | (      | |   | || |      | |   ) |   | |   | | \   || | \_  )
     | )      | (___) || (____/\| (__/  )___) (___| )  \  || (___) |
     |/       (_______)(_______/(______/ \_______/|/    )_)(_______)
                                                             
    """)

    if not fil_files:
        console.print("[red]No .fil files found![/red]")
        return

    par_file= Path(f"/lustre_archive/spotlight/data/{args.src_name}.par") #parameter file for pulsar from ATNF catalogue
    if not par_file.exists():
        raise FileNotFoundError(f"Parameter file not found for {args.src_name}. Please create one at '/lustre_archive/spotlight/data/'")

    log_path = args.input_dir / LOG_FILE_NAME
    open(log_path, "w").close()

    console.print(f"[cyan] Input directory: {args.input_dir}[/cyan]")
    console.print(f"[green]Found {len(fil_files)} .fil files[/green]")
    console.rule("[bold blue] Beginning pulsar folding using prepfold")

    fold_dir = args.input_dir / "FOLDING_OUTPUTS"
    os.makedirs(fold_dir, exist_ok=True)

    start_time = time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(run_folding, f, par_file, fold_dir, log_path): f for f in fil_files}
        for fut in tqdm(as_completed(futures), total=len(fil_files), desc="Folding", unit="file"):
            try:
                fut.result()
            except Exception as e:
                console.print(f"[red]Error processing {futures[fut]}: {e}[/red]")

    total_time = time.time() - start_time
    h, rem = divmod(total_time, 3600)
    m, s = divmod(rem, 60)
    formatted_time = f"{int(h):02d}:{int(m):02d}:{s:05.2f}"

    # Prepend final summary
    with open(log_path, "r+") as f:
        old_content = f.read()
        f.seek(0)
        f.write(f"Folding done for all the files\nTime taken: {formatted_time}\n{old_content}")

    console.print(f"\n [bold green]Done with pulsar folding for {args.src_name}[/bold green]")
    console.print(f"\n📝 Log saved at: [cyan]{log_path}[/cyan]")

if __name__ == "__main__":
    main()

