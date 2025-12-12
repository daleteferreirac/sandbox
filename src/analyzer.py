# Main module of the structural PDB analyzer project.
from utils import load_pdb, extract_residues, count_residue_types

atoms = load_pdb("../data/1CRN.pdb")
print(len(atoms))

residues = extract_residues(atoms)
count_residues = count_residue_types(residues)

print(f"Total residues: {len(residues)}")
print("First five residues:", list(residues.items())[:5])
print("Count of residues: ", count_residues)
