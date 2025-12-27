# Main module of the structural PDB analyzer project.
from utils import load_pdb, extract_residues, classify_residues, multiple_protein_analyzer

# Load atoms from PDB file
atoms = load_pdb("../data/1A3N.pdb")
print(f"Number of atoms: {len(atoms)}")

# Extract amino acid residues
residues, chain_info = extract_residues(atoms)


print(f"Total residues: {len(residues)}")
print("First five residues: ", list(residues.items())[:5])
print(f"chain information: {chain_info}")



classes, counts, chain_counts = classify_residues(residues) # tuple unpacking
print("Total residue classes:", classes)
print("Total residue counts:", counts)
print("% per chain:", chain_counts)

protein1 = "../data/1A3N.pdb"
protein2 = "../data/1CRN.pdb"
multiple_protein_analyzer(protein1, protein2)