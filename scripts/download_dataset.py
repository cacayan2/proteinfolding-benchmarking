# A dataset downloading and preprocessing pipeline for the benchmarking project.
import os
import sys
import subprocess
import benchmark_config as cfg
import json
import requests
import warnings
from pathlib import Path
from Bio.PDB import PDBParser, PPBuilder
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from rdkit import Chem
from meeko import MoleculePreparation, PDBQTWriterLegacy
from openbabel import openbabel as ob
import pdbfixer
from openmm.app import PDBFile

sys.setrecursionlimit(10000) # Increase the recursion limit for meeko

# Loading the config file.
sys.path.append(str(Path(__file__).parent))
warnings.simplefilter("ignore", PDBConstructionWarning)

# Some helper functions.
def fetch_pdb(pdb_id, out_path):
    """
    Obtains a pdb file from rcsb given a PDB ID and output path.
    
    :param pdb_id: The PDB ID to download
    :param out_path: The path to save the PDB file to
    """
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    print(f"[DOWNLOAD] Fetching {pdb_id} from {url}")
    r = requests.get(url, timeout = 30)
    if r.status_code != 200: # Check if the request was successful for debugging.
        raise RuntimeError(f"Failed to download PDB {pdb_id}")
    with open(out_path, "w") as f:
        f.write(r.text)

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
def preprocess_pdb(raw_pdb, clean_protein_out, ligand_out, metadata_out, fasta_out):
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
    organic_candidates = [] # Initialize the organic candidates
    for resn, atoms in ligands.items(): # For each ligand
        heavy = sum(parse_heavy_atoms(a) for a in atoms) # Count the heavy atoms
        if heavy < 5: # If the ligand has less than 5 heavy atoms
            continue # Skip the ligand
        has_carbon = any(" C" in a[12:20] or "C" == a[76:78].strip() for a in atoms) # Check if the ligand has carbon
        if not has_carbon: # If the ligand does not have carbon
            continue # Skip the ligand
        organic_candidates.append((heavy, resn, atoms)) # Add the ligand to the organic candidates
    if not organic_candidates: # If there are no organic candidates
        print(f"[WARNING] No organic ligand found - receptor only.")
        ligand_found = False # Set ligand_found to False
        selected_atoms = None # Set selected_atoms to None
    else:
        organic_candidates.sort(reverse = True) # Sort the organic candidates
        heavy, resn, selected_atoms = organic_candidates[0] # Get the first organic candidate
        print(f"[LIGAND] Selected {resn} with {heavy} heavy atoms for downstream docking.")
        ligand_found = True # Set ligand_found to True
    ligand_sdf = ligand_out.replace(".pdb", ".sdf") # Get the ligand SDF file
    ligand_pdbqt = ligand_out.replace(".pdb", ".pdbqt") # Get the ligand PDBQT file
    receptor_pdbqt = clean_protein_out.replace(".pdb", ".pdbqt") # Get the receptor PDBQT file
    if ligand_found: # If a ligand was found
        print(f"[LIGAND] Ligand found for this receptor, writing {ligand_out}")
        with open(ligand_out, "w") as f:
            f.writelines(selected_atoms) # Write the selected atoms
        raw_mol = Chem.MolFromPDBFile(ligand_out, sanitize = False, removeHs = False) # Read the ligand PDB file
        if raw_mol is None:
            raise RuntimeError(f"RDKit failed to load the ligand PDB: {ligand_out}")
        frags = Chem.GetMolFrags(raw_mol, asMols = True, sanitizeFrags = True) # Get the ligand fragments
        organic_frags = [m for m in frags if any (a.GetSymbol() == "C" for a in m.GetAtoms())] # Get the organic fragments
        largest = max(organic_frags, key = lambda m: m.GetNumAtoms()) # Get the largest organic fragment
        largest = Chem.AddHs(largest) # Add hydrogen atoms
        Chem.SanitizeMol(largest) # Sanitize the ligand
    # Write the SDF for DiffDock
        print(f"[LIGAND] Writing {ligand_sdf}")
        Chem.MolToMolFile(largest, ligand_sdf) # Write the ligand SDF
    # Convert the ligand using meeko
        mp = MoleculePreparation() # Initialize the MoleculePreparation class
        setup_list = mp.prepare(largest) # Prepare the ligand, returns MoleculeSetup
        setup = setup_list[0] # Get the first setup
        writer = PDBQTWriterLegacy() # Initialize the PDBQTWriterLegacy class
        txt = writer.write_string(setup) # Write the PDBQT file
        if isinstance(txt, tuple): # Meeko 0.5 returns (text, error)
            txt = txt[0] # Get the text
        with open(ligand_pdbqt, "w") as f:
            f.write(txt) # Write the PDBQT file
    else:
        print("[WARNING] No ligand found, skipping PDBQT generation.")
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
        ligand_found = ligand_found, # Ligand found
        ligand_file = ligand_out if ligand_found else None, # Ligand file
        ligand_sdf = ligand_sdf if ligand_found else None, # Ligand SDF file
        ligand_pdbqt = ligand_pdbqt if ligand_found else None, # Ligand PDBQT file
        receptor_pdbqt = receptor_pdbqt, # Receptor PDBQT file
        fasta_file = fasta_out, # FASTA file
        raw_pdb = raw_pdb, # Raw PDB file
    )

    with open(metadata_out, "w") as f: # Write the metadata
        json.dump(meta, f, indent = 2)
    return meta

def main():
    """
    Main function.
    """
    print(f"[DOWNLOAD] Downloading dataset")
    for pdb_id in cfg.PDB_LIST:
        pdb_dir = os.path.join(cfg.DATA_DIR, pdb_id)
        os.makedirs(pdb_dir, exist_ok = True)
        raw_path = os.path.join(pdb_dir, f"{pdb_id}.pdb")
        clean_path = os.path.join(pdb_dir, f"protein.pdb")
        ligand_path = os.path.join(pdb_dir, f"ligand.pdb")
        fasta_path = os.path.join(pdb_dir, f"protein.fasta")
        meta_path = os.path.join(pdb_dir, f"metadata.json")

        print(f"[DOWNLOAD] Processing {pdb_id}...")
        if not os.path.exists(raw_path):
            fetch_pdb(pdb_id, raw_path)
            meta = preprocess_pdb(raw_path, clean_path, ligand_path, meta_path, fasta_path)
            print(f"[PREPROCESS] Finished processing {pdb_id}!")
        else:
            print(f"[DOWNLOAD] Skipping {pdb_id}, directory already exists...")
    print(f"[DOWNLOAD] Finished downloading datasets!")

if __name__ == "__main__":
    main()
