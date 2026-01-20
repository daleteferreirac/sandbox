# Analyzes a single protein and returns information in dictionary and graph format.
from utils import analyze_proteins, load_pdb, extract_atom_coordinates, detect_residue_contacts
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
#    plt.show()

atoms_lines = load_pdb("../data/4AG8.pdb")
coords_data = extract_atom_coordinates(atoms_lines)

contacts, residues = detect_residue_contacts(coords_data, cutoff=4.5)
print(f"Number of residue contacts: {len(contacts)}")
print("10 contacts:")
for c in list(contacts)[:10]: # pair of contact: (('chain', 'residue number'), ('chain', 'residue number'))
    print(c)

# matrix of contacts
res = len(residues)
matrix = []
for i in range(res):
    row = []
    for j in range(res):
        row.append(0)
    matrix.append(row)

index_residues = {} # key: i, value: res number, {'900': 0, '1073': 1, '1129': 2, '1147': 3,...}
for i, res in enumerate(residues):
    index_residues[res] = i

for item in contacts:
    res1 = item[0][1] #
    res2 = item[1][1]

    i = index_residues[res1]
    j = index_residues[res2]

    matrix[i][j] = 1
    matrix[j][i] = 1
for line in matrix:
    print(line)
