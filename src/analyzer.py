# Main module of the structural PDB analyzer project.
from utils import load_pdb, extract_residues

atoms = load_pdb("../data/1CRN.pdb")
print(len(atoms))

residues = extract_residues(atoms)

print(f"Total residues: {len(residues)}")
print("First five residues:", list(residues.items())[:5])