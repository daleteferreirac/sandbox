# Analyzes a single protein and returns information in dictionary and graph format.
from utils import parse_pdb, analyze_proteins, detect_residue_contacts
import csv
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

protein = "../data/1A6M.pdb"

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
#    plt.show()

atoms = parse_pdb("../data/1A6M.pdb")
contacts, residues = detect_residue_contacts(atoms, cutoff=4.5)

print(f"Number of residue contacts: {len(contacts)}")
print("10 contacts:")
for c in list(contacts)[:10]: # pair of contact: (('chain', 'residue number'), ('chain', 'residue number'))
    print(c)

