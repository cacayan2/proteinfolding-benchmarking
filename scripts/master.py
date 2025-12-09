# master.py

import os

import sys

import subprocess

from pathlib import Path

import benchmark_config as cfg

ROOT = Path(__file__).resolve().parent

SCRIPTS = ROOT  # folder containing run_*.py


def run_step(env: str, script_relpath: str):
   script_path = str(SCRIPTS / script_relpath)
   print("\n" + "="*90)
   print(f"[MASTER] Running {script_relpath} (env: {env})")
   print("="*90 + "\n")
   # Non-interactive conda activation: THIS FIXES THE STOP ISSUE
   cmd = (
       f"bash -c '"
       f"source ~/.bashrc >/dev/null 2>&1 ; "
       f"conda activate {env} ; "
       f"python -u {script_path}"
       f"'"
   )
   process = subprocess.Popen(
       cmd,
       shell=True,
       stdout=subprocess.PIPE,
       stderr=subprocess.STDOUT,
       bufsize=1,
       text=True
   )
   for line in process.stdout:
       print(line, end="", flush=True)
   rc = process.wait()
   if rc != 0:
       raise RuntimeError(f"[MASTER] FAILED: {script_relpath}")
   print(f"\n[MASTER] COMPLETED: {script_relpath}\n")


def main():

    print("\n==============================")

    print("     BENCHMARK PIPELINE       ")

    print("==============================\n")

    # All steps in order

    STEPS = [

        (cfg.TOOLS_ENV,        "download_tools.py"),

        (cfg.DATA_ENV,         "download_dataset.py"),

        (cfg.MUTATE_ENV,       "mutate_dataset.py"),

        (cfg.PG_ENV,           "run_proteingenerator.py"),

        (cfg.RF_ENV,           "run_rosettafold.py"),

        (cfg.DD_ENV,           "run_diffdock.py"),

        (cfg.VINA_ENV,         "run_vina.py"),

    ]

    for env, script in STEPS:

        run_step(env, script)

    print("\n==============================")

    print("       PIPELINE COMPLETE      ")

    print("==============================\n")

    print(f"Results:  {cfg.RESULTS_DIR}")

    print(f"Plots:    {cfg.PLOTS_DIR}\n")


if __name__ == "__main__":

    main()
 