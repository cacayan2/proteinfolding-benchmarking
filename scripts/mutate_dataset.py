# A dataset mutation and regeneration pipeline for the benchmarking project.
import os
import sys
import subprocess
import random
import json
import warnings
from pathlib import Path
import benchmark_config as cfg
from Bio import SeqIO

from Bio.PDB import PDBParser, PPBuilder
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from rdkit import Chem
from meeko import MoleculePreparation, PDBQTWriterLegacy
import pdbfixer
from openmm.app import PDBFile

warnings.simplefilter("ignore", PDBConstructionWarning)
sys.setrecursionlimit(10000)

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

def mutate_sequence(seq, rate=0.02):
    """
    Returns a mutated sequence and a list of mutation dictionaries.
    Uses the simplest approach: random point mutations with probability `rate`.
    """
    seq = list(seq)
    mutations = []

    for i, aa in enumerate(seq): # iterate over sequence
        if random.random() <= rate and aa in AMINO_ACIDS: # mutate
            new_aa = random.choice(AMINO_ACIDS) # random mutation
            while new_aa == aa: # make sure it's not the same
                new_aa = random.choice(AMINO_ACIDS) # random mutation
            mutations.append({ # add mutation to list
                "position": i + 1,
                "original_residue": aa,
                "mutated_residue": new_aa,
                "type": "point"
            })
            seq[i] = new_aa # mutate

    return "".join(seq), mutations

def refold_sequence_with_rosettafold(fasta_path, outdir):
   """
   Runs RFAA (factory repo) with MSA + templating enabled.
   fasta_path: path to FASTA file
   outdir: output directory where model_001.pdb will be written
   """
   if os.path.exists(outdir):
       print(f"[RFAA] Skipping {fasta_path}, output directory already exists...")
       return
   print(outdir)
   os.makedirs(outdir, exist_ok=True)
   python_cmd = (
       f"export PYTHONPATH={cfg.RF_DIR}:$PYTHONPATH && "
       f"python {os.path.join(cfg.RF_DIR, 'rf2aa/run_inference.py')} "
       f"--config-name=protein.yaml "  
       f"--config-path={os.path.join(cfg.RF_DIR, 'rf2aa/config/inference')} "
       f"job_name=mutant "
       f"output_path={outdir} "
       f"protein_inputs.A.fasta_file={fasta_path}"
   )
   cmd = ["conda", "run", "-n", cfg.RF_ENV, "bash", "-c", python_cmd]
   print("\n[RFAA CMD]", " ".join(cmd))
   subprocess.run(cmd, check=True, cwd=cfg.RF_DIR)

def extract_fasta_from_pdb(pdb_path, fasta_path):
    """
    Extracts a FASTA file from a PDB file.
    
    :param pdb_path: The path to the PDB file
    :param fasta_path: The path to save the FASTA file to
    """
    print(f"[FASTA] Extracting chain A from {fasta_path}")
    parser = PDBParser(QUIET = True) # Suppress warnings
    structure = parser.get_structure("X", pdb_path) # Parse the PDB file
    seq = None # Initialize the sequence
    ppb = PPBuilder() # Initialize the PPBuilder
    for model in structure:
        for chain in model:
            chain_id = chain.id.strip() # Get the chain ID
            if cfg.PG_SELECT_FIRST_CHAIN: # If we want to use the first chain
                if seq is None: # If we haven't found a sequence yet
                    peptides = ppb.build_peptides(chain) # Build the peptides
                    if peptides:
                        seq = peptides[0].get_sequence() # Get the sequence
                        break
            else:
                if chain_id == "A": # If the chain ID is A
                    peptides = ppb.build_peptides(chain) # Build the peptides
                    if peptides:
                        seq =  "".join(str(pp.get_sequence())) # Get the sequence
                        break
            if seq:
                break # If we found a sequence, break out of the loop
    if not seq:
        print(f"[WARNING] Could not extract sequence - writing an empty FASTA file.")
        seq = ""
    with open(fasta_path, "w") as f:
        f.write(f">chainA\n{seq}\n")
    return seq

def parse_heavy_atoms(line):
    """
    Extracts heavy atoms from a PDB file.
    Lots of trial and error with this one to figure out
    which line contained the heavy atom.
    
    :param line: The line to parse
    """
    element = line[76:78].strip() # Get the element
    if not element:
        element = line[12:14].strip()
    return 0 if element.upper() == "H" else 1 # Return 0 if hydrogen, 1 otherwise

# Preprocessing steps
def preprocess_pdb(raw_pdb, clean_protein_out, original_directory, new_directory, metadata_out, fasta_out):
    """
    Preprocesses a PDB file.
    
    :param raw_pdb: The path to the raw PDB file
    :param clean_protein_out: The path to save the cleaned protein to
    :param ligand_out: The path to save the ligands to
    :param metadata_out: The path to save the metadata to
    :param fasta_out: The path to save the FASTA file to
    """
    sequence = extract_fasta_from_pdb(raw_pdb, fasta_out) # Extract the sequence
    target_chain = None # Initialize the target chain
    protein_lines = [] # Initialize the protein lines
    ligands = {} # Initialize the ligands
    with open(raw_pdb) as f:
        for line in f:
            if line.startswith("ATOM"):
                ch = line[21].strip() # Get the chain
                if cfg.PG_SELECT_FIRST_CHAIN:
                    if target_chain is None: # If we haven't found the target chain
                        target_chain = ch # Set the target chain
                    if ch != target_chain: # If the chain is not the target chain
                        continue # Skip the line
                else:
                    if ch != "A": # If the chain is not A
                        continue # Skip the line
                line = list(line) # Convert the line to a list
                if line[16] != " ": # If the residue name is not a space
                    line[16] = " " # Replace the residue name with a space
                line = "". join(line) # Convert the list back to a string
                protein_lines.append(line) # Add the line to the protein lines
            elif line.startswith("HETATM"): # If the line starts with HETATM
                resn = line[17:20].strip() # Get the residue name
                if resn in ("HOH", "WAT"): # If the residue name is HOH or WAT i.e. water
                    continue # Skip the line
                ligands.setdefault(resn, []).append(line) # Add the line to the ligands
    renumbered = [] # Initialize the renumbered lines
    prev_res = None # Initialize the previous residue
    new_idx = 0 # Initialize the new index
    for line in protein_lines: # For each line in the protein lines
        old = int(line[22:26]) # Get the old index
        if old != prev_res: # If the old index is different from the previous residue
            new_idx += 1 # Increment the new index
            prev_res = old # Set the previous residue
        newline = line[:22] + f"{new_idx:4d}" + line[26:] # Add the new index
        renumbered.append(newline) # Add the line to the renumbered lines
    with open(clean_protein_out, "w") as f: # Write the cleaned protein
        f.writelines(renumbered) # Write the renumbered lines
    receptor_pdbqt = clean_protein_out.replace(".pdb", ".pdbqt") # Get the receptor PDBQT file
    # Converting receptor to PDBQT using RDKit and Meeko
    print(f"[RECEPTOR] Repairing receptor and adding hydrogens with PDBFixer")
    fixer = pdbfixer.PDBFixer(filename=clean_protein_out) # Initialize the PDBFixer
    fixer.findMissingResidues() # Find missing residues
    fixer.findMissingAtoms() # Find missing atoms
    fixer.addMissingAtoms() # Add missing atoms
    fixer.addMissingHydrogens(pH = 7.4) # Add missing hydrogens
    fixed_pdb = clean_protein_out.replace(".pdb", "_fixed.pdb") # Get the fixed PDB file
    with open(fixed_pdb, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f) # Write the fixed PDB file
    print(f"[RECEPTOR] Preparing rigid receptor PDBQT (Meeko) to write {receptor_pdbqt}")
    # Load receptor PDB with RDKit
    protein_mol = Chem.MolFromPDBFile(fixed_pdb, sanitize = False, removeHs = False) # Read the PDB file
    if protein_mol is None:
        raise RuntimeError(f"[ERROR] RDKit failed to load the PDB: {fixed_pdb}")
    # Make all hydrogens explicit
    protein_mol = Chem.AddHs(protein_mol) # Add hydrogen atoms
    # Force receptor into one fragment (required for Meeko)
    frags = Chem.GetMolFrags(protein_mol, asMols = True, sanitizeFrags = False) # Get the fragments
    if len(frags) > 1:
        print(f"[RECEPTOR] RDKit detected {len(frags)} fragments, stitching fragments...")
        # Combine all fragments into a single molecule container
        combo = frags[0]
        for frag in frags[1:]: 
            combo = Chem.CombineMols(combo, frag) # Combine the fragments
        # Add dummy single-bonds to link fragments so Meeko allows them
        em = Chem.EditableMol(combo)
        # We need to track any offsets for atom indices
        frag_starts = []
        offset = 0
        for frag in frags:
            frag_starts.append(offset) # Get the start of the fragment
            offset += frag.GetNumAtoms() # Get the number of atoms in the fragment
        # Connect every fragment to the first
        anchor = 0 # Anchor atom index
        for i in range(1, len(frags)):
            em.AddBond(anchor, frag_starts[i], Chem.rdchem.BondType.SINGLE) # Add a single bond
        stitched = em.GetMol() # Get the stitched molecule
        # Sanitize the stitch molecule without valence errors
        Chem.SanitizeMol(stitched, sanitizeOps=Chem.SANITIZE_NONE)
        protein_mol = stitched
    # Meeko Receptor PDBQT Writing
    prep = MoleculePreparation() # Initialize the MoleculePreparation class
    setup_list = prep.prepare(protein_mol) # Prepare the protein, returns MoleculeSetup
    setup = setup_list[0] # Get the first setup
    setup.is_rigid = True # Set the setup to be rigid
    setup.flexibility_model = None # Set the flexibility model to None
    setup.bonds_to_break = [] # Set the bonds to break to an empty list
    # Write the PDBQT file safely
    writer = PDBQTWriterLegacy() # Initialize the PDBQTWriterLegacy class
    txt = writer.write_string(setup) # Write the PDBQT file
    if isinstance(txt, tuple): # Meeko 0.5 returns (text, error)
        txt = txt[0] # Get the text
    with open(receptor_pdbqt, "w") as f:
        f.write(txt) # Write the PDBQT file
    print(f"[RECEPTOR] Receptor PDBQT successfully written.")
    # Generating metadata

    meta = dict(
        sequence_length = len(sequence), # Sequence length
        chain_used = target_chain, # Chain used
        residues = new_idx, # Residues
        ligand_file = os.path.join(new_directory, "ligand.pdb"),
        ligand_sdf = os.path.join(new_directory, "ligand.sdf"),
        ligand_pdbqt = os.path.join(new_directory, "ligand.pdbqt"),
        receptor_pdbqt = receptor_pdbqt, # Receptor PDBQT file
        fasta_file = fasta_out, # FASTA file
        raw_pdb = raw_pdb, # Raw PDB file
    )

    with open(metadata_out, "w") as f: # Write the metadata
        json.dump(meta, f, indent = 2)
    return meta


def main():
    print("[MUTATE] Generating mutated dataset")

    # Create mutated copies for each PDB folder inside DATA_DIR
    for pdb_id in os.listdir(cfg.DATA_DIR): # Iterate over the PDB folders
        wt_dir = os.path.join(cfg.DATA_DIR, pdb_id) # WT directory
        
        if not os.path.isdir(wt_dir): # If it's not a directory
            continue

        mut_id = pdb_id + "m" # Mutated ID
        if 'm' in pdb_id: # If it already exists
            print(f"[MUTATE] Skipping {pdb_id}, mutated directory already exists...")
            continue
        mut_dir = os.path.join(cfg.DATA_DIR, mut_id) # Mutated directory
        os.makedirs(mut_dir, exist_ok=True) # Create the directory

        print(f"[MUTATE] Processing {pdb_id}: {mut_id}") # Print

        fasta_path = os.path.join(wt_dir, "protein.fasta") # WT FASTA
        mut_fasta = os.path.join(mut_dir, "protein.fasta") # Mutated FASTA
        if not os.path.exists(fasta_path): # If it doesn't exist
            print(f"[WARNING] No protein.fasta for {pdb_id}, skipping.") 
            continue

        mutated_seq = None # Mutated sequence
        mutations = None # Mutations
        if not os.path.exists(mut_fasta) and not os.path.exists(os.path.join(mut_dir, "mutations.json")): # If it doesn't exist
            with open(fasta_path) as f: # Open
                lines = f.read().splitlines() # Read
                seq = "" if len(lines) < 2 else lines[1] # Get the sequence

            #  Mutate sequence 
            mutated_seq, mutations = mutate_sequence(seq) # Mutate
            print(f"[MUTATE] {len(mutations)} mutations applied.") # Print

            #  Save mutated FASTA 
            mut_fasta = os.path.join(mut_dir, "protein.fasta") # Mutated FASTA

            with open(mut_fasta, "w") as f:
                f.write(f">chainA_mutated\n{mutated_seq}\n") # Write
        else:
            print(f"[MUTATE] Using existing mutated FASTA: {mut_fasta}")
            with open(mut_fasta) as f:
                all_records = list(SeqIO.parse(f, "fasta"))
                seq = all_records[0].seq
                mutated_seq = str(seq)
            with open(os.path.join(mut_dir, "mutations.json"), "r") as f:
                mutations = json.load(f)
        #  Refold mutated structure

        refold_outdir = os.path.join(mut_dir, "rfaa_out")  # RFAA output directory

        refold_sequence_with_rosettafold(mut_fasta, refold_outdir)

        # RFAA (factory default) writes 'mutant.pdb' at rfaa_out root

        # but we keep this flexible in case configs change.

        mut_pdb = "mutant.pdb"

        print(f"[MUTATE] Using RFAA PDB: {mut_pdb}")

        #  Preprocess PDB just like download_dataset
        mut_raw_path = os.path.join(mut_dir, "rfaa_out/mutant.pdb")
        mut_clean_path = os.path.join(mut_dir, "protein.pdb")
        mut_fasta_path = os.path.join(mut_dir, "protein.fasta")
        mut_meta_path = os.path.join(mut_dir, "metadata.json")

        print("[MUTATE] Preprocessing mutated PDB...")
        meta = preprocess_pdb(mut_raw_path, mut_clean_path, wt_dir, mut_dir, mut_meta_path, mut_fasta_path)
        #  Copy ligand files / raw PDB exactly 
        for fname in ["ligand.pdb", "ligand.pdbqt", "ligand.sdf"]: 
            src = os.path.join(wt_dir, fname) # Source
            dst = os.path.join(mut_dir, fname) # Destination
            if os.path.exists(src): # If it exists
                print(f"[COPY] {fname}") # Print
                with open(src, "r") as fin, open(dst, "w") as fout: # Open
                    fout.write(fin.read()) # Write

        #  Update metadata 
        meta_in = os.path.join(wt_dir, "metadata.json") # WT metadata
        meta_out = os.path.join(mut_dir, "metadata.json") # Mutated metadata

        if os.path.exists(meta_in): # If it exists
            with open(meta_in) as f: # Open
                meta = json.load(f) # Load
        else:
            meta = {} # Empty

        meta["mutated_from"] = pdb_id # Mutated from
        meta["mutated_id"] = mut_id # Mutated ID
        meta["mutations"] = mutations # Mutations
        meta["sequence_length"] = len(mutated_seq) # Sequence length

        with open(meta_out, "w") as f: # Open
            json.dump(meta, f, indent=2) # Write

        #  Write mutations.json 
        mut_info_path = os.path.join(mut_dir, "mutations.json") # Mutations info
        with open(mut_info_path, "w") as f: # Open
            json.dump({
                "original_id": pdb_id, # Original
                "mutated_id": mut_id, # Mutated
                "mutations": mutations # Mutations
            }, f, indent=2)

        print(f"[MUTATE] Finished {pdb_id}: {mut_id}")
    print("[MUTATE] All mutations complete!")

if __name__ == "__main__":
    main()
 