import os
import subprocess
import sys
from pathlib import Path
import time

import benchmark_config as cfg

# Helper function to run shell comands.
def run(cmd):
    """
    Run a shell command
    
    :param cmd: The command to run
    """
    print(f"[ENV]: {cmd}")
    result = subprocess.run(cmd, shell = True)
    if result.returncode != 0:
        print(f"[ERROR] Command failed: {cmd}")
        sys.exit(1)

def conda_env_exists(env):
    """
    Check if a conda environment exists
    
    :param env: The environment to check
    """
    result = subprocess.run(
        f"conda env list | grep '^{env}", # Check if the environment exists
        shell=True, # Run in shell
        stdout=subprocess.PIPE, # Capture stdout
        stderr=subprocess.PIPE, # Capture stderr
    )
    return result.returncode == 0

# ProteinGenerator Installation
def install_proteingenerator():
    """
    Installs ProteinGenerator
    """
    # Clone tool from GitHub
    print(f"[TOOLS] Installing ProteinGenerator")
    print(f"[TOOLS] Downloading ProteinGenerator from GitHub")
    if os.path.exists(cfg.PG_DIR):
        print("[TOOLS] ProteinGenerator already installed, skipping.")
    else:
        run(f"git clone https://github.com/RosettaCommons/protein_generator {cfg.PG_DIR}")
        print(f"[TOOLS] ProteinGenerator downloaded successfully.")
    
    print(f"[TOOLS] Downloading model checkpoints")
    
    # Download model checkpoints
    model_dir = os.path.join(cfg.PG_DIR, "model")
    os.makedirs(model_dir, exist_ok = True)
    ckpts = [
        "http://files.ipd.uw.edu/pub/sequence_diffusion/checkpoints/SEQDIFF_221219_equalTASKS_nostrSELFCOND_mod30.pt",
        "http://files.ipd.uw.edu/pub/sequence_diffusion/checkpoints/SEQDIFF_230205_dssp_hotspots_25mask_EQtasks_mod30.pt",
    ]
    
    for url in ckpts:
        file = os.path.join(model_dir, url.split("/")[-1])
        if not os.path.exists(file):
            run(f"wget {url} -O {file}")
        else:
            print(f"[TOOLS] Checkpoint already exists: {file}")

    print(f"[TOOLS] ProteinGenerator fully installed.")

def install_rfaa():
    """
    Installs RoseTTAFold-All-Atom
    """
    print(f"[TOOLS] Installing RoseTTAFold-All-Atom")
    if os.path.exists(cfg.RF_DIR):
        print("[TOOLS] RoseTTAFold-All-Atom already installed, skipping.")
    else:
        run(f"git clone https://github.com/baker-laboratory/RoseTTAFold-All-Atom.git {cfg.RF_DIR}")
        print(f"[TOOLS] RoseTTAFold-All-Atom downloaded successfully.")
    rfaa_dependencies()
    print(f"[TOOLS] RoseTTAFold-All-Atom fully installed.")

def rfaa_dependencies():
    """
    Loads all RFAA dependencies.
    """
    # Running setup scripts according to RFAA documentation.
    run(f"conda run {cfg.RF_DIR}/rf2aa/SE3Transformer/ pip3 install --no-cache-dir -r requirements.txt")
    run(f"python3 {cfg.RF_DIR}/rf2aa/setup.py install")

    # Prompting user to download signalp6.
    signalp_dir = os.path.join(cfg.TOOLS_DIR, "signalp6")
    os.makedirs(signalp_dir, exist_ok = True)

    tarball_path = os.path.join(signalp_dir, "signalp-6.0.tar.gz")

    print(f"[TOOLS] For installation of RFAA, a licensed copy of signalp6 is required.")
    print(f"[TOOLS] This unforunately must be downloaded manually.")
    print(f"[TOOLS] Download it from:")
    print(f"    https://services.healthtech.dtu.dk/services/SignalP-6.0/")
    print(f"[TOOLS] Move the downloaded tarball to:\n {tarball_path}\n")
    input(f"[TOOLS] Press ENTER once you have placed signalp-6.0h.fast.tar.gz...\n")

    if not os.path.exists(tarball_path):
        print(f"[ERROR] File does not exist: {tarball_path}")
        print(f"[ERROR] You must download and place it manually. Aborting.\n")
        return

    print("[TOOLS] Running signalp6-register...")
    run(f"conda run -n {cfg.SIGNALP_ENV} signalp6-register {tarball_path}")

    # Locating installed weights inside the conda environment
    conda_prefix = subprocess.check_output(
        f"conda run -n {cfg.SIGNALP_ENV} bash -c 'echo $CONDA_PREFIX'",
        shell=True
    ).decode().strip()

    # Locating weights
    weights_dir = os.path.join(
        conda_prefix,
        "lib",
        "python3.10",
        "site-packages",
        "signalp",
        "model_weights"
    )

    # Setting paths
    distilled = os.path.join(weights_dir, "distilled_model_signalp6.pt")
    ensemble = os.path.join(weights_dir, "ensemble_model_signalp6.pt")

    # Checking for distilled weights path
    if not os.path.exists(distilled):
        print(f"[ERROR] Expected distilled weight not found:\n   {distilled}")
        print(f"[ERROR] Registration may have failed.")
        return
    
    # Rename model files (instructed by documentation)
    os.rename(distilled, ensemble)

    print("\n[TOOLS] Success - SignalP-6.0 fully installed and registered.")
    print(f"[TOOLS] Weights installed in:\n   {weights_dir}\n")

    # Running install dependencies script.
    run(f"bash {cfg.RF_DIR}/install_dependencies.sh")

    # Downloading paper weights.
    if not os.path.exists(f"{cfg.RF_DIR}/RFAA_paper_weights.pt"):
        print("[TOOLS] Downloading RFAA paper weights...")
        run(f"wget http://files.ipd.uw.edu/pub/RF-All-Atom/weights/RFAA_paper_weights.pt -P {cfg.RF_DIR}")

    print("[TOOLS] The following files are VERY large and need to be in the tool directory.")
    input("[TOOLS] Please ensure you have enough space. Press ENTER to continue...\n")

    # Downloading UniRef30
    print("[TOOLS] Downloading UniRef30...")
    if not os.path.exists(f"{cfg.RF_DIR}/UniRef30_2020_06"):
        run(f"mkdir -p {cfg.RF_DIR}/UniRef30_2020_06")
        run(f"aria2c -x 16 -s 16 {cfg.RF_DIR} UniRef30_2020_06_hhsuite.tar.gz http://wwwuser.gwdg.de/~compbiol/uniclust/2020_06/UniRef30_2020_06_hhsuite.tar.gz")
        run(f"tar -xvzf UniRef30_2020_06_hhsuite.tar.gz -C {cfg.RF_DIR}/UniRef30_2020_06")
        run(f"rm {cfg.RF_DIR}/UniRef30_2020_06_hhsuite.tar.gz")
    else:
        print("[TOOLS] UniRef30 already downloaded.")
    
    # Downloading BFD
    print("[TOOLS] Downloading BFD...")
    if not os.path.exists(f"{cfg.RF_DIR}/bfd"):
        run(f"mkdir -p {cfg.RF_DIR}/bfd")
        run(f"aria2c -x 16 -s 16 {cfg.RF_DIR} bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt.tar.gz https://bfd.mmseqs.com/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt.tar.gz")
        run(f"tar -xvzf bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt.tar.gz -C {cfg.RF_DIR}/bfd")
        run(f"rm {cfg.RF_DIR}/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt.tar.gz")
    else:
        print("[TOOLS] BFD already downloaded.")

    # Downloading Structure Templates
    if not os.path.exists(f"{cfg.RF_DIR}/pdb100_2021Mar03"):
        run(f"aria2c -x 16 -s 16 {cfg.RF_DIR}/pdb100_2021Mar03.tar.gz https://files.ipd.uw.edu/pub/RoseTTAFold/pdb100_2021Mar03.tar.gz")
        run(f"tar -xvzf pdb100_2021Mar03.tar.gz -C {cfg.RF_DIR}")
        run(f"rm {cfg.RF_DIR}/pdb100_2021Mar03.tar.gz")
    else:
        print("[TOOLS] Structure Templates already downloaded.")

    # Downloading BLAST
    if not os.path.exists(f"{cfg.RF_DIR}/blast-2.2.26"):
        run(f"mkdir -p {cfg.RF_DIR}/blast-2.2.26")
        run(f"aria2c -x 16 -s 16 {cfg.RF_DIR}/blast-2.2.26-x64-linux.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.26/blast-2.2.26-x64-linux.tar.gz")
        run(f"tar -xvzf blast-2.2.26-x64-linux.tar.gz -C {cfg.RF_DIR}/blast-2.2.26")
        run(f"cp -r {cfg.RF_DIR}blast-2.26/blast-2.26/ {cfg.RF_DIR}/blast-2.26_bk")
        run(f"rm -r {cfg.RF_DIR}/blast-2.2.26")
        run(f"mv {cfg.RF_DIR}/blastr-2.26_bk {cfg.RF_DIR}/blast-2.2.26")
        run(f"rm {cfg.RF_DIR}/blast-2.2.26-x64-linux.tar.gz")
    else:
        print("[TOOLS] BLAST already downloaded.")

    print("[TOOLS] Installation complete.")