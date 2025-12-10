#!/usr/bin/env python3

import os

import time

import subprocess

import resource

import traceback

from pathlib import Path

from meeko import MoleculePreparation  # ligand only (safe)

import benchmark_config as cfg


# ======================================================================

#  Helper: Run system command

# ======================================================================

def run_cmd(cmd):

    """Run a shell command, print it, and raise on failure."""

    print("[CMD]", " ".join(cmd))

    result = subprocess.run(cmd)

    if result.returncode != 0:

        raise RuntimeError(f"[ERROR] Command failed: {cmd}")


# ======================================================================

#  Helper: Log runtime/memory to complexity_log.csv

# ======================================================================

def log_complexity(tool, pdb_id, elapsed, mem_mb):

    if not cfg.LOG_TIME_AND_MEMORY:

        return

    header = "tool,pdb_id,time_sec,mem_mb\n"

    row = f"{tool},{pdb_id},{elapsed:.2f},{mem_mb:.2f}\n"

    if not os.path.exists(cfg.COMPLEXITY_LOG):

        with open(cfg.COMPLEXITY_LOG, "w") as f:

            f.write(header)

    with open(cfg.COMPLEXITY_LOG, "a") as f:

        f.write(row)


# ======================================================================

#  FIXED: Manual PDB → PDBQT conversion (NO MEEKO)

# ======================================================================

def prepare_receptor_pdbqt(input_pdb, output_pdbqt):

    """

    Convert PDB into a minimal-valid PDBQT file that Vina accepts.

    Avoids Meeko's receptor loader (which fails for your version).

    """

    atoms = []

    with open(input_pdb, "r") as f:

        for line in f:

            if line.startswith(("ATOM", "HETATM")):

                atoms.append(line)

    if len(atoms) == 0:

        raise RuntimeError(f"[ERROR] No ATOM/HETATM lines found in {input_pdb}")

    with open(output_pdbqt, "w") as f:

        f.write("ROOT\n")

        for atom in atoms:

            # Append dummy atom type required by Vina

            pdbqt_line = atom.rstrip() + "    C\n"

            f.write(pdbqt_line)

        f.write("ENDROOT\n")

        f.write("TORSDOF 0\n")

    print(f"[INFO] Receptor PDBQT written → {output_pdbqt}")


# ======================================================================

#  Ligand PDBQT handling

# ======================================================================

def prepare_ligand_pdbqt(input_sdf, output_pdbqt):

    """

    Use Meeko (safe part) to convert SDF → PDBQT for ligands.

    """

    prep = MoleculePreparation()

    mol = prep.prepare_from_file(input_sdf)

    pdbqt_string, _ = prep.write_pdbqt_string()

    with open(output_pdbqt, "w") as f:

        f.write(pdbqt_string)

    print(f"[INFO] Ligand PDBQT written → {output_pdbqt}")


# ======================================================================

#  Build Vina config file

# ======================================================================

def write_vina_config(center, size, out_config_file, receptor_pdbqt, ligand_pdbqt):

    with open(out_config_file, "w") as f:

        f.write(f"receptor = {receptor_pdbqt}\n")

        f.write(f"ligand = {ligand_pdbqt}\n")

        f.write(f"center_x = {center[0]}\n")

        f.write(f"center_y = {center[1]}\n")

        f.write(f"center_z = {center[2]}\n")

        f.write(f"size_x = {size[0]}\n")

        f.write(f"size_y = {size[1]}\n")

        f.write(f"size_z = {size[2]}\n")

        f.write(f"exhaustiveness = 8\n")

        f.write("cpu = 0\n")

    print(f"[INFO] Vina config written → {out_config_file}")


# ======================================================================

#  Compute bounding box from ligand

# ======================================================================

def compute_box_from_ligand(ligand_sdf):

    import rdkit

    from rdkit import Chem

    from rdkit.Chem import AllChem

    mol = Chem.SDMolSupplier(ligand_sdf, removeHs=False)[0]

    if mol is None:

        raise RuntimeError(f"[ERROR] RDKit failed to load ligand SDF: {ligand_sdf}")

    conf = mol.GetConformer()

    xs, ys, zs = [], [], []

    for i in range(mol.GetNumAtoms()):

        pos = conf.GetAtomPosition(i)

        xs.append(pos.x)

        ys.append(pos.y)

        zs.append(pos.z)

    center = (sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs))

    size = (cfg.VINA_BOX_PADDING, cfg.VINA_BOX_PADDING, cfg.VINA_BOX_PADDING)

    return center, size


# ======================================================================

#  Run Vina docking

# ======================================================================

def run_vina(receptor_pdbqt, ligand_pdbqt, out_pdbqt, config_txt):

    cmd = [

        "conda", "run", "-n", cfg.VINA_ENV,

        cfg.VINA_BIN,

        "--config", config_txt,

        "--out", out_pdbqt,

        "--cpu", "0"

    ]

    run_cmd(cmd)


# ======================================================================

#  Ground-truth docking

# ======================================================================

def run_vina_ground_truth(pdb_id):

    print(f"\n=== Ground Truth: {pdb_id} ===")

    src_dir = os.path.join(cfg.DATA_DIR, pdb_id)

    receptor_pdb = os.path.join(src_dir, "protein.pdb")

    ligand_sdf = os.path.join(src_dir, "ligand.sdf")

    ligand_pdbqt = os.path.join(src_dir, "ligand.pdbqt")

    out_dir = os.path.join(cfg.VINA_OUT, pdb_id)

    os.makedirs(out_dir, exist_ok=True)

    receptor_pdbqt = os.path.join(out_dir, "receptor.pdbqt")

    vina_out = os.path.join(out_dir, "vina_out.pdbqt")

    config_txt = os.path.join(out_dir, "vina_config.txt")

    print("[INFO] Running Vina on ground truth for", pdb_id)

    # Convert receptor PDB → PDBQT (manual)

    prepare_receptor_pdbqt(receptor_pdb, receptor_pdbqt)

    # Compute docking box

    center, size = compute_box_from_ligand(ligand_sdf)

    # Write config file

    write_vina_config(center, size, config_txt, receptor_pdbqt, ligand_pdbqt)

    # Timing

    start = time.time()

    mem_start = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss

    run_vina(receptor_pdbqt, ligand_pdbqt, vina_out, config_txt)

    end = time.time()

    mem_end = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss

    elapsed = end - start

    mem_mb = (mem_end - mem_start) / 1024

    log_complexity("vina_ground_truth", pdb_id, elapsed, mem_mb)

    print(f"[SUCCESS] Ground truth docking completed for {pdb_id}")


# ======================================================================

#  ProteinGenerator docking

# ======================================================================

def run_vina_generated(pdb_id):

    print(f"\n=== Generated Model: {pdb_id} ===")

    gt_dir = os.path.join(cfg.DATA_DIR, pdb_id)

    pg_receptor = os.path.join(cfg.PG_OUT, pdb_id, f"{pdb_id}.pdb")

    ligand_pdbqt = os.path.join(gt_dir, "ligand.pdbqt")

    ligand_sdf = os.path.join(gt_dir, "ligand.sdf")

    out_dir = os.path.join(cfg.VINA_OUT, f"{pdb_id}_generated")

    os.makedirs(out_dir, exist_ok=True)

    receptor_pdbqt = os.path.join(out_dir, "receptor.pdbqt")

    vina_out = os.path.join(out_dir, "vina_out.pdbqt")

    config_txt = os.path.join(out_dir, "vina_config.txt")

    if not os.path.isfile(pg_receptor):

        print(f"[SKIP] No ProteinGenerator structure found for {pdb_id}")

        return

    print("[INFO] Running Vina on ProteinGenerator output for", pdb_id)

    # Convert receptor

    prepare_receptor_pdbqt(pg_receptor, receptor_pdbqt)

    # Compute docking box from ground-truth ligand

    center, size = compute_box_from_ligand(ligand_sdf)

    # Config

    write_vina_config(center, size, config_txt, receptor_pdbqt, ligand_pdbqt)

    start = time.time()

    mem_start = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss

    run_vina(receptor_pdbqt, ligand_pdbqt, vina_out, config_txt)

    end = time.time()

    mem_end = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss

    elapsed = end - start

    mem_mb = (mem_end - mem_start) / 1024

    log_complexity("vina_generated", pdb_id, elapsed, mem_mb)

    print(f"[SUCCESS] Generated docking completed for {pdb_id}")


# ======================================================================

#  Run all PDBs

# ======================================================================

def run_vina_all():

    print("[INFO] Running Vina for all PDB entries:", cfg.PDB_LIST)

    for pdb_id in cfg.PDB_LIST:

        try:

            run_vina_ground_truth(pdb_id)

            run_vina_generated(pdb_id)

        except Exception as e:

            print(f"[ERROR] Vina failed for {pdb_id}: {e}")

            traceback.print_exc()


if __name__ == "__main__":

    run_vina_all()
 