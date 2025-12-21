# Main module of the structural PDB analyzer project.
from utils import load_pdb, extract_residues, count_residue_types, classify_residues

# Load atoms from PDB file
atoms = load_pdb("../data/1CRN.pdb")
print(f"Number of atoms: {len(atoms)}")

# Extract amino acid residues
residues, chain_info = extract_residues(atoms)
count_residues = count_residue_types(residues)

print(f"Total residues: {len(residues)}")
print("First five residues: ", list(residues.items())[:5])
print(f"chain information: {chain_info}")

print("Residue composition:")
for res, count in count_residues.items():
    print(res, count)

classes, counts = classify_residues(residues) # tuple unpacking
print("Residue classes:", classes)
print("Residue counts:", counts)

