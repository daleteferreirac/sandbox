# Main module of the structural PDB analyzer project.
from utils import analyze_proteins
import csv
import matplotlib.pyplot as plt

protein1 = "../data/1A3N.pdb"
protein2 = "../data/1CRN.pdb"
protein3 = "../data/4AG8.pdb"

results = analyze_proteins(protein1, protein2, protein3)
for protein, data in results.items():
    print(f"Protein: {protein}")
    for key, value in data.items():
        print(f"\t{key}: {value}")

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

for protein, data in results.items():
    # Extract residue class counts
    class_data = data["residues-classes-counts"]

    # Prepare data for plotting
    classes = list(class_data.keys())
    counts = [class_data[cls]["count"] for cls in classes]

    plt.figure()     # Create bar plot
    plt.bar(classes, counts)

    plt.xlabel("Residue class")
    plt.ylabel("Number of residues")
    plt.title(f"Residue class composition – {protein}")
    plt.show()
