# # Compares proteins and generates a CSV file.
from utils import analyze_proteins
import csv

protein1 = "../data/1A3N.pdb"
protein2 = "../data/1CRN.pdb"
protein3 = "../data/4AG8.pdb"

results = analyze_proteins(protein1, protein2, protein3)

with open("protein-summary.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile) # create a CSV writer object

    #head
    writer.writerow(["Protein ID", "atoms", "residues", "chains", "waters", "ligands"])

    #rows
    for protein, data in results.items():
        writer.writerow([
            protein,
            data["atoms"],
            data["residues"],
            len(data["residue-chain"]),
            data["waters"],
            len(data["ligands"])
        ])

# the results is generation of a .csv file with the information