import os
import subprocess
import sys
from pathlib import Path

import benchmark_config as cfg

# Tools and their corresponding .yml files.
ENV_DIR = Path(cfg.PROJECT_ROOT) / "envs"
ENV_FILES = {
    cfg.PG_ENV: Path(ENV_DIR) / "run_proteingenerator.yml",
    cfg.RF_ENV: Path(ENV_DIR) / "run_rosettafold.yml",
    cfg.DD_ENV: Path(ENV_DIR) / "run_diffdock.yml",
    cfg.VINA_ENV: Path(ENV_DIR) / "run_vina.yml",
    cfg.DATA_ENV: Path(ENV_DIR) / "download_dataset.yml",
    cfg.SIGNALP_ENV: Path(ENV_DIR) / "signalp.yml",
    cfg.TOOLS_ENV: Path(ENV_DIR) / "download_tools.yml",
    cfg.MUTATE_ENV: Path(ENV_DIR) / "mutate_dataset.yml"
}

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

# Function which resets and rebuilds all environments
def recreate_env(env_name, yaml_path):
    """
    Deletes old environment and populates a new conda environment.
    
    :param env_name: The name of the environment to create
    :param yaml_path: The path to the .yml file
    """
    print(f"[ENV] Resetting environment: {env_name}")
    print(f"[ENV] YAML file: {yaml_path}")

    # Remove old environment if exists.
    print(f"[ENV] Removing old environment: {env_name}")
    try:
        run(f"conda env remove -n {env_name} --yes")
    except:
        print(f"[ENV] Environment {env_name} does not exist, creating now...")
        pass

    # Create new environment.
    if yaml_path.exists():
        print(f"[ENV] Creating new environment: {env_name} from {yaml_path}.")
        run(f"conda env create -n {env_name} -f \"{yaml_path}\"")
    else:
        print(f"[WARNING] YAML file not found for {env_name}")
        print("[WARNING] Creating a minimal environment instead (Python 3.10)")
        run(f"conda create -y -n {env_name} python=3.10")

    print(f"[ENV] Environment {env_name} created successfully.")

def main():
    print(f"[ENV] Restting environments.")
    for env_name, yaml_path in ENV_FILES.items():
        recreate_env(env_name, yaml_path)
    print(f"[ENV] Finished resetting environments.")

if __name__ == "__main__":
    main()
