# Contains central config for pipeline
from pathlib import Path
import os

print("[CONFIG] Setting Benchmark Config Variables...")

# Project Root (auto-detected)
PROJECT_ROOT = Path(__file__).resolve().parent.parent
PROJECT_ROOT = str(PROJECT_ROOT)
print(f"[CONFIG] Project Root: {PROJECT_ROOT}")

# Data Directories
DATA_DIR = os.path.join(PROJECT_ROOT, "data")
PDB_LIST = ["1HVR", "1STP", "1OHR", "2VTA", "1BMA"]
DATA_ENV = "download_dataset"
print(f"[CONFIG] Data Directory: {DATA_DIR}")

# Results Directories
RESULTS_DIR = os.path.join(PROJECT_ROOT, "results")
PG_OUT = os.path.join(RESULTS_DIR, "proteingenerator")
RF_OUT = os.path.join(RESULTS_DIR, "rosettafold")
DD_OUT = os.path.join(RESULTS_DIR, "diffdock")
VINA_OUT = os.path.join(RESULTS_DIR, "vina")
print(f"[CONFIG] ProteinGenerator Output Directory: {PG_OUT}")
print(f"[CONFIG] RosettaFold Output Directory: {RF_OUT}")
print(f"[CONFIG] DiffDock Output Directory: {DD_OUT}")
print(f"[CONFIG] Vina Output Directory: {VINA_OUT}")

# Plots Directory
PLOTS_DIR = os.path.join(PROJECT_ROOT, "plots")
print(f"[CONFIG] Plots Directory: {PLOTS_DIR}")

# Tools Directories
TOOLS_DIR = os.path.join(PROJECT_ROOT, "tools")
TOOLS_ENV = "download_tools"
MUTATE_ENV = "mutate_dataset"
print(f"[CONFIG] Tools Directory: {TOOLS_DIR}")

## ProteinGenerator
PG_DIR = os.path.join(TOOLS_DIR, "ProteinGenerator")
PG_ENV = "proteingenerator"
PG_INFERENCE = os.path.join(TOOLS_DIR, "ProteinGenerator/inference.py")
print(f"[CONFIG] ProteinGenerator Directory: {PG_DIR}")
print(f"[CONFIG] ProteinGenerator Environment: {PG_ENV}")
print(f"[CONFIG] ProteinGenerator Inference Script: {PG_INFERENCE}")

## RoseTTAFold-All-Atom
RF_DIR = os.path.join(TOOLS_DIR, "RoseTTAFold-All-Atom")
RF_ENV = "RFAA"
RF_SCRIPT = os.path.join(RF_DIR, "rf2aa/run_inference.py")
SIGNALP_ENV = "signalp6"
print(f"[CONFIG] RoseTTAFold-All-Atom Directory: {RF_DIR}")
print(f"[CONFIG] RoseTTAFold-All-Atom Environment: {RF_ENV}")
print(f"[CONFIG] RoseTTAFold-All-Atom Script: {RF_SCRIPT}")

## DiffDock
DD_DIR = os.path.join(TOOLS_DIR, "DiffDock")
DD_ENV = "diffdock"
print(f"[CONFIG] DiffDock Directory: {DD_DIR}")
print(f"[CONFIG] DiffDock Environment: {DD_ENV}")

## Vina
VINA_DIR = os.path.join(TOOLS_DIR, "vina")
VINA_BIN = os.path.join(VINA_DIR, "vina")
VINA_ENV = "vina"
print(f"[CONFIG] Vina Directory: {VINA_DIR}")
print(f"[CONFIG] Vina Environment: {VINA_ENV}")

# Pipeline Settings
## ProteinGenerator
## We need to enforce ligand trimming & chain selection,
## otherwise ProteinGenerator throws an error.
PG_MAX_LIGAND_HEAVY_ATOMS = 40 # This avoids tensor mismatch errors.
PG_SELECT_FIRST_CHAIN = True # Ensure we use the first chain available (usually chain A).
print(f"[CONFIG] ProteinGenerator Max Ligand Heavy Atoms: {PG_MAX_LIGAND_HEAVY_ATOMS}")
print(f"[CONFIG] ProteinGenerator Select First Chain: {PG_SELECT_FIRST_CHAIN}")

## RFAA
RFAA_MODE = "full" # Alllows RFAA to run with all options available (MSA + templating).
RFAA_CONFIG_RELATIVE = "config/inference/protein.yaml" # Sets the path to the config file.
LOG_MSA_DEPTH = True
MSA_DEPTH_FILE = os.path.join(RESULTS_DIR, "msa_depths.csv")
print(f"[CONFIG] RoseTTAFold-All-Atom Mode: {RFAA_MODE}")
print(f"[CONFIG] RoseTTAFold-All-Atom Config Relative Path: {RFAA_CONFIG_RELATIVE}")
print(f"[CONFIG] RoseTTAFold-All-Atom Log MSA Depth: {LOG_MSA_DEPTH}")
print(f"[CONFIG] RoseTTAFold-All-Atom MSA Depth File: {MSA_DEPTH_FILE}")

## DiffDock
DD_NUM_SAMPLES = 20
print(f"[CONFIG] DiffDock Number of Samples: {DD_NUM_SAMPLES}")

## Vina
VINA_BOX_PADDING = 10.0 # Å, padding around a ligand when generating auto-sized box

## Performance Logging
LOG_TIME_AND_MEMORY = True
LOG_CPU_TIME = True
LOG_MEMORY_TRACE = True
MEMORY_SAMPLING_INTERVAL = 0.5
LOG_IO_STATS = True
LOG_GPU_STATS = True
GPU_SAMPLING_INTERVAL = 0.5
print(f"[CONFIG] Log Time and Memory: {LOG_TIME_AND_MEMORY}")
print(f"[CONFIG] Log CPU Time: {LOG_CPU_TIME}")
print(f"[CONFIG] Log Memory Trace: {LOG_MEMORY_TRACE}")
print(f"[CONFIG] Memory Sampling Interval: {MEMORY_SAMPLING_INTERVAL}")
print(f"[CONFIG] Log IO Stats: {LOG_IO_STATS}")
print(f"[CONFIG] Log GPU Stats: {LOG_GPU_STATS}")
print(f"[CONFIG] GPU Sampling Interval: {GPU_SAMPLING_INTERVAL}")

## Performance Logging Results Directories
COMPLEXITY_LOG = os.path.join(RESULTS_DIR, "complexity_log.csv")
MEMORY_TRACE_DIR = os.path.join(RESULTS_DIR, "memory_traces")
GPU_TRACE_DIR = os.path.join(RESULTS_DIR, "gpu_traces")
IO_TRACE_DIR = os.path.join(RESULTS_DIR, "io_traces")
print(f"[CONFIG] Complexity Log: {COMPLEXITY_LOG}")
print(f"[CONFIG] Memory Trace Directory: {MEMORY_TRACE_DIR}")
print(f"[CONFIG] GPU Trace Directory: {GPU_TRACE_DIR}")
print(f"[CONFIG] IO Trace Directory: {IO_TRACE_DIR}")


