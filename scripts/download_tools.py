import os
import sys
import subprocess
import benchmark_config as cfg

# Checking appropriate environment used.
REQUIRED_ENV = cfg.TOOLS_ENV
current_env = os.environ.get("CONDA_DEFAULT_ENV")
if current_env != REQUIRED_ENV:
    print(f"[DOWNLOAD] Relaunching download_dataset.py in {REQUIRED_ENV} environment.")
    cmd = ["conda", "run", "-n", REQUIRED_ENV, "python"] + sys.argv
    subprocess.run(cmd)
    sys.exit()


from pathlib import Path
import time
import yaml

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

def rfaa_databases():
    """
    Configures downloaded RFAA databases. 
    """
    ur30 = os.path.join(cfg.RF_DIR, "UniRef30_2020_06")
    bfd = os.path.join(cfg.RF_DIR, "bfd")
    blastmat = os.path.join(cfg.RF_DIR, "blast-2.2.26")
    pdb100 = os.path.join(cfg.RF_DIR, "pdb100_2021Mar03")
    
    for p in [ur30, bfd, blastmat, pdb100]:
        if not p.exists():
            print(f"[WARNING] Expected DB directory missing: {p}")

        # Export environment variables permanently with modification to bashrc - recommended by documentation.
        bashrc = Path.home() / ".bashrc"

        def append_export(name, value):
            """
            Appends export line to bashrc

            :param name: name of environment variable
            :param value: value of environment variable
            """
            line = f"export {name}=\"{value}\"" # line to append
            existing = bashrc.read_text() if bashrc.exists() else "" # existing content

            if line not in existing: # if not already in bashrc
                with open(bashrc, "a") as f: # append to bashrc
                    f.write("\n" + line + "\n")
                print("[TOOLS] Appended export line to bashrc.")
            else:
                print("[TOOLS] bashrc already contains {name}")
        append_export("DB_UR30", str(ur30)) # Append export line for ur30
        append_export("DB_BFD", str(bfd)) # Append export line for bfd
        append_export("DB_PDB", str(blastmat)) # Append export line for blast
        append_export("DB_BLASTMAT", str(pdb100)) # Append export line for pdb100

def rfaa_dependencies():
    """
    Loads all RFAA dependencies.
    """
    # Running setup scripts according to RFAA documentation.
    run(f"conda run pip3 install --no-cache-dir -r {cfg.RF_DIR}/rf2aa/SE3Transformer/requirements.txt")
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
    
    # Configuring Databases
    print(f"[TOOLS] Configuring databases...")
    rfaa_databases()

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

def prepare_diffdock_weights():
    """
    Unzips DiffDock weights
    """
    # This is a file that only the author of this pipeline has - it automatically downloads the correct zip file.
    if os.path.exists(cfg.PROJECT_ROOT / "scripts" / "diffdock_weights.py"):
        run(f"run python {cfg.PROJECT_ROOT}/scripts/diffdock_weights.py", shell = True)
    # Unzips weights.
    MODEL_DIR = os.path.join(cfg.DD_DIR, "models")
    ZIP_DIR = os.path.join(cfg.DD_DIR, "diffdock_weights.zip")
    run(f"unzip {ZIP_DIR} 'diffdock-main/workdir/'-d {cfg.DD_DIR}/tmp_dd", shell = True)
    run(f"mv {cfg.DD_DIR}/tmp_dd/diffdock-main/workdir/* {MODEL_DIR}/", shell = True)
    run("rm -rf tmp_dd diffdock_weights.zip", shell = True)

def has_diffdock_weights():
    """
    Checks if DiffDock weights are installed.
    """
    MODEL_DIR = os.path.join(cfg.DD_DIR, "models")
    REQUIRED_FILES = ["diffdock_model.pt",
                    "score_model.pt",
                    "noise_scheduler.pt",
                    "inference_args.json"]
    os.makedirs(MODEL_DIR, exist_ok = True)
    all_exist = all(os.path.exists(os.path.join(MODEL_DIR, f)) for f in REQUIRED_FILES)
    return all_exist
    
def install_diffdock():
    """
    Installs DiffDock
    """
    print(f"[TOOLS] Installing DiffDock")
    print(f"[TOOLS] Downloading DiffDock from GitHub")
    if os.path.exists(cfg.DD_DIR):
        if os.path.exists(cfg.DD_DIR):
            print("[TOOLS] DiffDock already installed, skipping.")
        else:
            run(f"git clone https://github.com/gcorso/DiffDock.git {cfg.DD_DIR}")
            print(f"[TOOLS] DiffDock downloaded successfully.")
        if has_diffdock_weights():
            print("[TOOLS] DiffDock weights already installed, skipping.")
        else:
            print(f"[TOOLS] The next step is to obtain the model weights.")
            print(f"[TOOLS] However, due to copyright infringement issues, the model weights are not publicly available.")
            print(f"[TOOLS] This pipeline was successfully run despite this, because there may or may not be another copy... somewhere...")
            print(f"[TOOLS] Please download the .zip file and put it in the root of DiffDock.")
            print(f"[TOOLS] Alternatively, reach out to Emil Cacayan for help - he has a script that automatically does this.")
            print(f"[TOOLS] Name it diffdock_weights.zip.")
            print(f"[TOOLS] Recommend saving the zip to another place, as this tool will delete it upon unpacking.")
            input(f"[TOOLS] Press ENTER to continue.")
            try: 
                prepare_diffdock_weights()
            except:
                print("[TOOLS] Please download the .zip file and put it in the root of DiffDock.")
                print("[TOOLS] Alternatively, reach out to Emil Cacayan for help - he has a script that automatically does this.")
                print("[TOOLS] Name it diffdock_weights.zip.")
                sys.exit(1)
    print(f"[TOOLS] DiffDock fully installed.")

def install_vina():
    """
    Installs Vina
    """
    if os.path.exists(cfg.VINA_BIN):
        print("[TOOLS] Vina already installed, skipping.")
        return
    os.makedirs(cfg.VINA_DIR, exist_ok = True)
    url = "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_linux_x86_64"
    out = os.path.join(cfg.VINA_DIR, "vina")
    run(f"wget {url} -O {out}")
    run(f"chmod +x {out}")
    print(f"[TOOLS] AutoDock Vina downloaded successfully.")
    
def main():
    """
    Main function
    """
    print("[TOOLS] Downloading tools...")
    install_proteingenerator()
    install_rfaa()
    install_diffdock()
    install_vina()
    print("[TOOLS] Tools downloaded successfully.")

if __name__ == "__main__":
    main()





            