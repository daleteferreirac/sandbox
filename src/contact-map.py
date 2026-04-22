from utils import parse_pdb, detect_residue_contacts
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import numpy as np
import os

pdb_path = "../data/4AG8.pdb"
atoms = parse_pdb(pdb_path)
pdb_name = os.path.basename(pdb_path).replace(".pdb", "")

contacts, residues = detect_residue_contacts(atoms, cutoff=4.5)
print(f"Number of residue contacts: {len(contacts)}")
print("10 contacts:")
for c in list(contacts)[:10]: # pair of contact: (('chain', 'residue number'), ('chain', 'residue number'))
    print(c)


# matrix of contacts
res = len(residues)
matrix = np.zeros((res, res))

index_residues = {} # key: i, value: res number, {'900': 0, '1073': 1, '1129': 2, '1147': 3,...}
for i, res in enumerate(residues):
    index_residues[res] = i

for item in contacts:
    res1 = item[0][1] #
    res2 = item[1][1]

    i = index_residues[res1]
    j = index_residues[res2]

    matrix[i, j] = 1
    matrix[j, i] = 1
for line in matrix:
    print(line)

# colors: 0 = black, 1 = yellow
cmap = ListedColormap(["white", "black"])
norm = BoundaryNorm([0, 0.5, 1], cmap.N)

plt.figure(figsize=(6, 6))
plt.imshow(matrix, cmap=cmap, norm=norm, origin="lower")
plt.xlabel("Residue index")
plt.ylabel("Residue index")
plt.title(f"Residue Contact Map ({pdb_name})")
plt.colorbar(ticks=[0, 1], label="Contact (0 = no, 1 = yes)")

plt.savefig(f"../outputs/{pdb_name}_contact_map.png", dpi=300, bbox_inches="tight")
plt.close()

