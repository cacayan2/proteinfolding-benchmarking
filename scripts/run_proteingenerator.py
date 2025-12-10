# run_proteingenerator.py

"""

Runs ProteinGenerator on each PDB, logs GPU + memory usage,

renames outputs, extracts sequence from .trb, and prepares PG output

for downstream DiffDock / AutoDock steps.

NO RDKit, NO Meeko, NO pdbqt generation here.

"""

import os

import sys

import json

import time

import subprocess

import threading

from pathlib import Path

# Load config

sys.path.append(str(Path(__file__).parent))

import benchmark_config as cfg


# ======================================================================

# GPU / Memory Tracking

# ======================================================================

def track_gpu(pid, outfile, stop_event):

    if not cfg.LOG_GPU_STATS:

        return

    os.makedirs(cfg.GPU_TRACE_DIR, exist_ok=True)

    with open(outfile, "w") as f:

        f.write("timestamp_sec,gpu_util,memory_used_mb\n")

    while not stop_event.is_set():

        try:

            out = subprocess.check_output(

                "nvidia-smi --query-gpu=utilization.gpu,memory.used "

                "--format=csv,noheader,nounits",

                shell=True

            ).decode().strip()

            util, mem = out.split(", ")

            with open(outfile, "a") as f:

                f.write(f"{time.time()},{util},{mem}\n")

        except:

            pass

        time.sleep(cfg.GPU_SAMPLING_INTERVAL)


def track_memory(pid, outfile, stop_event):

    if not cfg.LOG_MEMORY_TRACE:

        return

    os.makedirs(cfg.MEMORY_TRACE_DIR, exist_ok=True)

    with open(outfile, "w") as f:

        f.write("timestamp_sec,vmrss_kb,vmhwm_kb\n")

    while not stop_event.is_set():

        try:

            with open(f"/proc/{pid}/status") as s:

                rss = hwm = None

                for line in s:

                    if line.startswith("VmRSS:"):

                        rss = line.split()[1]

                    elif line.startswith("VmHWM:"):

                        hwm = line.split()[1]

                if rss and hwm:

                    with open(outfile, "a") as f:

                        f.write(f"{time.time()},{rss},{hwm}\n")

        except:

            pass

        time.sleep(cfg.MEMORY_SAMPLING_INTERVAL)


# ======================================================================

# Complexity Logging

# ======================================================================

def log_complexity(tool, pdb_id, start, end, mem_mb):

    if not cfg.LOG_TIME_AND_MEMORY:

        return

    need_header = not os.path.exists(cfg.COMPLEXITY_LOG)

    with open(cfg.COMPLEXITY_LOG, "a") as f:

        if need_header:

            f.write("tool,pdb_id,time_sec,mem_mb\n")

        f.write(f"{tool},{pdb_id},{end - start:.2f},{mem_mb:.2f}\n")


# ======================================================================

# PG Inference

# ======================================================================

def run_proteingenerator_once(pdb_id, protein_pdb, out_dir):

    """

    Runs PG inside its conda env, tracks GPU & memory, and renames outputs.

    """

    meta_path = os.path.join(cfg.DATA_DIR, pdb_id, "metadata.json")

    with open(meta_path) as f:

        meta = json.load(f)

    chain = meta["chain_used"]

    length = meta["residues"]

    contigs = f"{chain}1-{length}"

    os.makedirs(out_dir, exist_ok=True)

    pg_cmd = [

        "conda", "run",

        "--cwd", out_dir,

        "-n", cfg.PG_ENV,

        "python", cfg.PG_INFERENCE,

        "--pdb", protein_pdb,

        "--num_designs", "1",

        "--contigs", contigs,

        "--dump_pdb",

        "--out", "."

    ]

    wrapped_cmd = f"/usr/bin/time -v {' '.join(pg_cmd)}"

    print(f"[PG CMD] {wrapped_cmd}")

    start = time.time()

    proc = subprocess.Popen(

        wrapped_cmd,

        shell=True,

        cwd=out_dir,

        stdout=subprocess.PIPE,

        stderr=subprocess.PIPE,

        text=True

    )

    stop_event = threading.Event()

    # tracking

    gpu_file = os.path.join(cfg.GPU_TRACE_DIR, f"PG_{pdb_id}.csv")

    mem_file = os.path.join(cfg.MEMORY_TRACE_DIR, f"PG_{pdb_id}.csv")

    if cfg.LOG_GPU_STATS:

        gpu_thread = threading.Thread(target=track_gpu, args=(proc.pid, gpu_file, stop_event))

        gpu_thread.start()

    if cfg.LOG_MEMORY_TRACE:

        mem_thread = threading.Thread(target=track_memory, args=(proc.pid, mem_file, stop_event))

        mem_thread.start()

    stdout, stderr = proc.communicate()

    stop_event.set()

    # extract memory from time output

    mem_mb = 0.0

    for line in stderr.splitlines():

        if "Maximum resident set size" in line:

            mem_kb = float(line.split()[-1])

            mem_mb = mem_kb / 1024.0

    end = time.time()

    log_complexity("ProteinGenerator", pdb_id, start, end, mem_mb)

    if proc.returncode != 0:

        print("[ERROR] PG failed:", stderr)

        return False

    rename_pg_outputs(out_dir, pdb_id)

    return True


# ======================================================================

# Output Handling

# ======================================================================

def rename_pg_outputs(out_dir, pdb_id):

    std_pdb = os.path.join(out_dir, "._000000.pdb")

    std_trb = os.path.join(out_dir, "._000000.trb")

    new_pdb = os.path.join(out_dir, f"{pdb_id}.pdb")

    new_trb = os.path.join(out_dir, f"{pdb_id}.trb")

    if os.path.exists(std_pdb):

        os.replace(std_pdb, new_pdb)

    if os.path.exists(std_trb):

        os.replace(std_trb, new_trb)

    extract_pg_sequence(new_trb, os.path.join(out_dir, "sequence.fasta"))


def extract_pg_sequence(trb_path, fasta_out):

    if not os.path.exists(trb_path):

        return

    try:

        import pickle

        with open(trb_path, "rb") as f:

            trb = pickle.load(f)

        seq = trb.get("best_sequence") or trb.get("sequence")

        if not seq:

            return

        with open(fasta_out, "w") as f:

            f.write(">PG_design\n" + seq + "\n")

    except Exception as e:

        print(f"[WARN] Could not extract PG sequence: {e}")


# ======================================================================

# Pipeline Loop

# ======================================================================

def process_pdb(pdb_id):

    pdb_dir = os.path.join(cfg.DATA_DIR, pdb_id)

    protein_pdb = os.path.join(pdb_dir, "protein.pdb")

    if not os.path.exists(protein_pdb):

        print(f"[SKIP] Missing protein.pdb for {pdb_id}")

        return

    out_dir = os.path.join(cfg.PG_OUT, pdb_id)

    print(f"\n=== ProteinGenerator: {pdb_id} ===")

    ok = run_proteingenerator_once(pdb_id, protein_pdb, out_dir)

    if not ok:

        print(f"[ERROR] PG failed for {pdb_id}")


def main():

    pdb_ids = sorted(

        d for d in os.listdir(cfg.DATA_DIR)

        if os.path.isdir(os.path.join(cfg.DATA_DIR, d))

    )

    print("[INFO] Running ProteinGenerator on:", pdb_ids)

    for pdb_id in pdb_ids:

        process_pdb(pdb_id)


if __name__ == "__main__":

    main()
 