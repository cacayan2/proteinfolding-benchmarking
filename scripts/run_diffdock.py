# run_diffdock.py

import os

import subprocess

import time

import resource

import psutil

import threading

import benchmark_config as cfg


# =====================================================================

# Helper: Run a command with error handling

# =====================================================================

def run_cmd(cmd):

    print("[CMD]", " ".join(cmd))

    result = subprocess.run(cmd, shell=False)

    if result.returncode != 0:

        raise RuntimeError(f"[ERROR] Command failed: {cmd}")


# =====================================================================

# Complexity Log

# =====================================================================

def _log_complexity(tool, pdb_id, elapsed):

    if not cfg.LOG_TIME_AND_MEMORY:

        return

    mem_mb = psutil.Process().memory_info().rss / (1024**2)

    header = "tool,pdb_id,time_sec,mem_mb\n"

    row = f"{tool},{pdb_id},{elapsed:.2f},{mem_mb:.2f}\n"

    if not os.path.exists(cfg.COMPLEXITY_LOG):

        with open(cfg.COMPLEXITY_LOG, "w") as f:

            f.write(header)

    with open(cfg.COMPLEXITY_LOG, "a") as f:

        f.write(row)


# =====================================================================

# Auxiliary Tracing Utilities (Memory / IO / GPU)

# =====================================================================

def start_memory_trace(tag):

    if not cfg.LOG_MEMORY_TRACE:

        return None

    os.makedirs(cfg.MEMORY_TRACE_DIR, exist_ok=True)

    path = os.path.join(cfg.MEMORY_TRACE_DIR, f"{tag}_diffdock_memory.csv")

    f = open(path, "w")

    f.write("time_sec,memory_mb\n")

    start = time.time()

    def trace():

        while True:

            mem = psutil.Process().memory_info().rss / (1024**2)

            f.write(f"{time.time() - start:.2f},{mem:.2f}\n")

            f.flush()

            time.sleep(cfg.MEMORY_SAMPLING_INTERVAL)

    t = threading.Thread(target=trace, daemon=True)

    t.start()

    return f


def start_io_trace(tag):

    if not cfg.LOG_IO_STATS:

        return None

    os.makedirs(cfg.IO_TRACE_DIR, exist_ok=True)

    path = os.path.join(cfg.IO_TRACE_DIR, f"{tag}_diffdock_io.csv")

    f = open(path, "w")

    f.write("time_sec,read_MB,write_MB\n")

    proc = psutil.Process()

    start = time.time()

    def trace():

        while True:

            io = proc.io_counters()

            f.write(f"{time.time() - start:.2f},{io.read_bytes/1e6:.2f},{io.write_bytes/1e6:.2f}\n")

            f.flush()

            time.sleep(cfg.MEMORY_SAMPLING_INTERVAL)

    threading.Thread(target=trace, daemon=True).start()

    return f


def start_gpu_trace(tag):

    if not cfg.LOG_GPU_STATS:

        return None

    os.makedirs(cfg.GPU_TRACE_DIR, exist_ok=True)

    path = os.path.join(cfg.GPU_TRACE_DIR, f"{tag}_diffdock_gpu.csv")

    f = open(path, "w")

    f.write("time_sec,gpu_mem,gpu_util\n")

    start = time.time()

    # CPU-only: GPU=0 always

    def trace():

        while True:

            f.write(f"{time.time() - start:.2f},0,0\n")

            f.flush()

            time.sleep(cfg.GPU_SAMPLING_INTERVAL)

    threading.Thread(target=trace, daemon=True).start()

    return f


# =====================================================================

# Build the DiffDock command

# =====================================================================

def build_diffdock_cmd(protein, ligand, out_dir):

    """Returns full shell command string for conda+CPU DiffDock run."""

    model_dir = os.path.join(cfg.DD_DIR, "models")

    # Disable GPU:

    env_prefix = "CUDA_VISIBLE_DEVICES=\"\" "

    cmd = (

        f"cd {cfg.DD_DIR} && "

        f"{env_prefix}"

        f"OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 NUMEXPR_NUM_THREADS=8 "

        f"python inference.py "

        f"--protein_path {protein} "

        f"--ligand_description {ligand} "

        f"--out_dir {out_dir} "

        f"--samples_per_complex {cfg.DD_NUM_SAMPLES} "

        f"--model_dir {model_dir}"

    )

    # conda execution wrapper:

    return ["conda", "run", "-n", cfg.DD_ENV, "bash", "-c", cmd]


# =====================================================================

# Docking Functions

# =====================================================================

def run_diffdock_native(pdb_id):

    print(f"\n=== Running DiffDock (Native) for {pdb_id} ===")

    protein = os.path.join(cfg.DATA_DIR, pdb_id, "protein.pdb")

    ligand  = os.path.join(cfg.DATA_DIR, pdb_id, "ligand.sdf")

    if not os.path.isfile(protein) or not os.path.isfile(ligand):

        print(f"[SKIP] Missing native protein/ligand for {pdb_id}")

        return

    out_dir = os.path.join(cfg.RESULTS_DIR, "diffdock_native", pdb_id)

    os.makedirs(out_dir, exist_ok=True)

    # tracing

    mem_f = start_memory_trace(f"{pdb_id}_native")

    io_f  = start_io_trace(f"{pdb_id}_native")

    gpu_f = start_gpu_trace(f"{pdb_id}_native")

    start = time.time()

    cmd = build_diffdock_cmd(protein, ligand, out_dir)

    run_cmd(cmd)

    elapsed = time.time() - start

    _log_complexity("diffdock_native", pdb_id, elapsed)

    if mem_f: mem_f.close()

    if io_f: io_f.close()

    if gpu_f: gpu_f.close()

    print(f"[SUCCESS] DiffDock native finished for {pdb_id}")


def run_diffdock_generated(pdb_id):

    print(f"\n=== Running DiffDock (Generated Structure) for {pdb_id} ===")

    pg_pdb = os.path.join(cfg.PG_OUT, pdb_id, f"{pdb_id}.pdb")

    ligand = os.path.join(cfg.DATA_DIR, pdb_id, "ligand.sdf")

    if not os.path.isfile(pg_pdb):

        print(f"[SKIP] Missing ProteinGenerator output for {pdb_id}")

        return

    if not os.path.isfile(ligand):

        print(f"[SKIP] Missing ligand.sdf for {pdb_id}")

        return

    out_dir = os.path.join(cfg.RESULTS_DIR, "diffdock_generated", pdb_id)

    os.makedirs(out_dir, exist_ok=True)

    mem_f = start_memory_trace(f"{pdb_id}_generated")

    io_f  = start_io_trace(f"{pdb_id}_generated")

    gpu_f = start_gpu_trace(f"{pdb_id}_generated")

    start = time.time()

    cmd = build_diffdock_cmd(pg_pdb, ligand, out_dir)

    run_cmd(cmd)

    elapsed = time.time() - start

    _log_complexity("diffdock_generated", pdb_id, elapsed)

    if mem_f: mem_f.close()

    if io_f: io_f.close()

    if gpu_f: gpu_f.close()

    print(f"[SUCCESS] DiffDock generated finished for {pdb_id}")


# =====================================================================

# Run for All PDB IDs

# =====================================================================

def run_diffdock_all():

    print("[INFO] Running DiffDock for:", cfg.PDB_LIST)

    for pdb_id in cfg.PDB_LIST:

        run_diffdock_native(pdb_id)

        run_diffdock_generated(pdb_id)

    print("\n[INFO] All DiffDock runs complete.")


if __name__ == "__main__":

    run_diffdock_all()
 