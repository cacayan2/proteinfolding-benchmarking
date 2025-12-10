# run_rosettafold.py

# Runs RoseTTAFold-All-Atom (RFAA) inference with full benchmarking instrumentation

import os

import sys

import time

import json

import subprocess

import threading

from pathlib import Path

import psutil

# Load pipeline config

sys.path.append(str(Path(__file__).parent))

import benchmark_config as cfg


###############################################################################

# Utility: Directory Setup

###############################################################################

def ensure_dirs():

    """Ensure benchmark output directories exist."""

    os.makedirs(cfg.RF_OUT, exist_ok=True)

    if cfg.LOG_MEMORY_TRACE:

        os.makedirs(cfg.MEMORY_TRACE_DIR, exist_ok=True)

    if cfg.LOG_GPU_STATS:

        os.makedirs(cfg.GPU_TRACE_DIR, exist_ok=True)

    if cfg.LOG_IO_STATS:

        os.makedirs(cfg.IO_TRACE_DIR, exist_ok=True)


###############################################################################

# Profilers

###############################################################################

def monitor_gpu(pdb_id, stop_event):

    """Poll nvidia-smi and log GPU usage."""

    if not cfg.LOG_GPU_STATS:

        return

    out_path = os.path.join(cfg.GPU_TRACE_DIR, f"{pdb_id}_rfaa_gpu.csv")

    with open(out_path, "w") as f:

        f.write("timestamp_ms,utilization.gpu [%],memory.used [MiB]\n")

        while not stop_event.is_set():

            try:

                smi = subprocess.check_output(

                    ["nvidia-smi", "--query-gpu=utilization.gpu,memory.used",

                     "--format=csv,noheader,nounits"],

                    encoding="utf-8"

                ).strip()

                util, mem = smi.split(", ")

                ts = int(time.time() * 1000)

                f.write(f"{ts},{util},{mem}\n")

            except Exception:

                pass

            time.sleep(cfg.GPU_SAMPLING_INTERVAL)


def monitor_memory(pdb_id, pid, stop_event):

    """Log memory usage of the RFAA process."""

    if not cfg.LOG_MEMORY_TRACE:

        return

    out_path = os.path.join(cfg.MEMORY_TRACE_DIR, f"{pdb_id}_rfaa_mem.csv")

    with open(out_path, "w") as f:

        f.write("timestamp_ms,rss_mb,vms_mb\n")

        proc = psutil.Process(pid)

        while not stop_event.is_set():

            try:

                mem = proc.memory_info()

                ts = int(time.time() * 1000)

                f.write(f"{ts},{mem.rss/1e6:.3f},{mem.vms/1e6:.3f}\n")

            except psutil.NoSuchProcess:

                break

            time.sleep(cfg.MEMORY_SAMPLING_INTERVAL)


def monitor_io(pdb_id, pid, stop_event):

    """Log I/O read/write bytes."""

    if not cfg.LOG_IO_STATS:

        return

    out_path = os.path.join(cfg.IO_TRACE_DIR, f"{pdb_id}_rfaa_io.csv")

    with open(out_path, "w") as f:

        f.write("timestamp_ms,read_bytes,write_bytes\n")

        proc = psutil.Process(pid)

        while not stop_event.is_set():

            try:

                io = proc.io_counters()

                ts = int(time.time() * 1000)

                f.write(f"{ts},{io.read_bytes},{io.write_bytes}\n")

            except psutil.NoSuchProcess:

                break

            time.sleep(cfg.MEMORY_SAMPLING_INTERVAL)


###############################################################################

# Complexity Logging

###############################################################################

def log_complexity(pdb_id, start, end, peak_mem_mb):

    """Write summary performance metrics."""

    if not cfg.LOG_TIME_AND_MEMORY:

        return

    header_needed = not os.path.exists(cfg.COMPLEXITY_LOG)

    with open(cfg.COMPLEXITY_LOG, "a") as f:

        if header_needed:

            f.write("tool,pdb_id,time_sec,mem_mb\n")

        f.write(f"RFAA,{pdb_id},{end-start:.2f},{peak_mem_mb:.2f}\n")


###############################################################################

# RFAA Execution

###############################################################################

def run_rfaa_for(pdb_id):

    print(f"\n==== Running RFAA for {pdb_id} ====")

    pdb_dir = os.path.join(cfg.DATA_DIR, pdb_id)

    fasta_path = os.path.join(pdb_dir, "protein.fasta")

    if not os.path.exists(fasta_path):

        print(f"[ERROR] Missing FASTA for {pdb_id}. Run download_dataset.py.")

        return False

    out_dir = os.path.join(cfg.RF_OUT, pdb_id)

    os.makedirs(out_dir, exist_ok=True)

    python_cmd = (

        f"export PYTHONPATH={cfg.RF_DIR}:$PYTHONPATH && "

        f"python {os.path.join(cfg.RF_DIR, 'rf2aa/run_inference.py')} "

        f"--config-name=protein.yaml "

        f"--config-path={os.path.join(cfg.RF_DIR, 'rf2aa/config/inference')} "

        f"job_name={pdb_id} "

        f"output_path={out_dir} "

        f"protein_inputs.A.fasta_file={fasta_path}"

    )

    cmd = ["conda", "run", "-n", cfg.RF_ENV, "bash", "-c", python_cmd]

    print("[RFAA CMD]", " ".join(cmd))

    start = time.time()

    # Launch RFAA process

    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, cwd=cfg.RF_DIR)

    # Start profilers

    stop_event = threading.Event()

    gpu_thread = threading.Thread(target=monitor_gpu, args=(pdb_id, stop_event))

    mem_thread = threading.Thread(target=monitor_memory, args=(pdb_id, process.pid, stop_event))

    io_thread = threading.Thread(target=monitor_io, args=(pdb_id, process.pid, stop_event))

    if cfg.LOG_GPU_STATS:

        gpu_thread.start()

    if cfg.LOG_MEMORY_TRACE:

        mem_thread.start()

    if cfg.LOG_IO_STATS:

        io_thread.start()

    # Wait for RFAA to finish

    stdout, stderr = process.communicate()

    stop_event.set()

    end = time.time()

    print(stdout)

    print(stderr)

    # Measure peak memory via psutil

    try:

        peak_mem = psutil.Process(process.pid).memory_info().rss / 1e6

    except psutil.NoSuchProcess:

        peak_mem = 0.0

    log_complexity(pdb_id, start, end, peak_mem)

    if process.returncode != 0:

        print(f"[ERROR] RFAA failed for {pdb_id}")

        return False

    print(f"[RFAA] Completed {pdb_id}")

    return True


###############################################################################

# Entrypoint

###############################################################################

def main():

    ensure_dirs()

    pdb_ids = sorted([

        d for d in os.listdir(cfg.DATA_DIR)

        if os.path.isdir(os.path.join(cfg.DATA_DIR, d))

    ])

    print("[INFO] Running RFAA on:", pdb_ids)

    for pdb_id in pdb_ids:

        run_rfaa_for(pdb_id)


if __name__ == "__main__":

    main()
 