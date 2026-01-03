# Analyzes a single protein and returns information in dictionary and graph format.
from utils import analyze_proteins
import csv
import matplotlib.pyplot as plt


protein = "../data/4AG8.pdb"

results = analyze_proteins(protein)
for protein, data in results.items():
    print(f"Protein: {protein}")
    for key, value in data.items():
        print(f"\t{key}: {value}")

# create the plot: Residue class composition – {protein}
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
