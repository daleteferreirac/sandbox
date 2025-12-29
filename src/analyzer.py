# Main module of the structural PDB analyzer project.
from utils import analyze_proteins

protein1 = "../data/1A3N.pdb"
protein2 = "../data/1CRN.pdb"

results = analyze_proteins(protein1, protein2)
for protein, data in results.items():
    print(f"Protein: {protein}")
    for key, value in data.items():
        print(f"\t{key}: {value}")
