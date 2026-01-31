# Helper functions for the structural analyzer

import os
import math

def load_pdb(filepath):
    """
    Load ATOM and HETATM records from a PDB file.
    :param : str
        path to PDB file
    Returns
    list of str
        Lines corresponding to ATOM and HETATM records. (the entire line of the pdb file)
    """
    atoms_lines = []
    with open(filepath, "r") as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"): # atoms or ligands
                atoms_lines.append(line.strip()) # removes \n

    return atoms_lines

def extract_residues(atom_lines):
    """
    Extracts all unique residues from a list of ATOM/HETATM lines.
    :param : list of str
        ATOM and HETATM lines from a PDB file.

    Returns a dictionary:
        { (chain, residue_number) : residue_name }
        the key is a tuple of chain and residue number.
        the value is a string of the residue name.
    chain_info dictionary:
        {chain: number of residues}
    """
    residues = {} # ('A', '1'): 'THR', (('A', '2'), 'LEU'), (('A', '3'), 'SER'),...
    waters = {}
    ligands = {}
    chain_info = {} # {'HYDROPHOBICS': {'count': 19, 'percentage': 45.24}, ...}

    for line in atom_lines:
        if line.startswith("ATOM"):
            res_name = line[17:20].strip()   # VAL, ALA, ARG...
            chain = line[21].strip()         # A, B, C...
            res_num = int(line[22:26].strip())    # 1, 2, 3, 8, 54...

            key = (chain, res_num) # this key(a tuple) >
            if key not in residues:
                residues[key] = res_name # will be assigned to residues dict and receive this value : res_name

                if chain not in chain_info:
                    chain_info[chain] = 1
                else:
                    chain_info[chain] += 1
        elif line.startswith("HETATM"):
            name_hetatm = line [17:20].strip()
            chain_hetatm = line[21].strip()
            num_hetatm = line[22:26].strip()
            if name_hetatm == "HOH":
                key = (chain_hetatm, num_hetatm)
                if key not in waters:
                    waters[key] = name_hetatm

            else:
                key = (chain_hetatm, num_hetatm)
                if key not in ligands:
                    ligands[key] = name_hetatm


    return residues, chain_info, waters, ligands

# def count_residue_types(residues):
#     """
#     Counts how many times each amino acid appears.
#     residues: dict {(chain, res_num): res_name}
#     returns: dict {res_name: count}
#     """
#     counts = {}
#
#     for res_name in residues.values():
#         if res_name not in counts:
#             counts[res_name] = 1
#         else:
#             counts[res_name] += 1
#
#     return counts

def classify_residues(residues):
    """
    Classify residues based on hydrophobic, polar and charged.
        use the dictionary residues: {(chain, res_num): res_name}

    :param : dict
        {(chain, residue_number): residue_name}
    :return: dict
    classes : dict
        {class_name: set(residue_names)}

    counts : dict
        {class_name: {'count': int, 'percentage': float}}

    chain_counts : dict
        {chain: {class_name: {'count': int, 'percentage': float}}}
    """
    hydrophobics = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"}
    polars = {"SER", "THR", "ASN", "GLN", "TYR", "CYS"}
    charged = {"ARG", "LYS", "HIS", "ASP", "GLU"}

    counts = {
        "HYDROPHOBICS": 0,
        "POLAR": 0,
        "CHARGED": 0
    }

    chain_counts = {}

    for (chain, _),res_name in residues.items(): # ignore the residue number
        if chain not in chain_counts:
            chain_counts[chain] = {
                "HYDROPHOBICS": 0,
                "POLAR": 0,
                "CHARGED": 0
            }

        if res_name in hydrophobics:
            counts["HYDROPHOBICS"] += 1
            chain_counts[chain]["HYDROPHOBICS"] += 1
        elif res_name in polars:
            counts["POLAR"] += 1
            chain_counts[chain]["POLAR"] += 1
        elif res_name in charged:
            counts["CHARGED"] += 1
            chain_counts[chain]["CHARGED"] += 1

    total = sum(counts.values())
    for key in counts:
        counts[key] = {"count": counts[key], "percentage": round((counts[key] / total)*100, 2) if total > 0 else 0}

    for chain, class_counts in chain_counts.items():
        total_chain = sum(class_counts.values())
        for cls in class_counts:
            class_counts[cls] = {
                "count": class_counts[cls],
                "percentage": round((class_counts[cls]/total_chain) * 100, 2)
                if total_chain > 0 else 0
            }

    return counts, chain_counts

def analyze_proteins(*filepaths):
    """
    Receives multiple PDB file paths and returns the information organized.
    :param filepaths: tuple of str
        PDB file paths
    :return: dict
    """
    results = {}

    for path in filepaths:
        protein_name = os.path.splitext(os.path.basename(path))[0]

        atoms = load_pdb(path)
        # extract amino acid residues and hetatm
        residues, chain_info, waters, ligands = extract_residues(atoms)
        counts, chain_counts = classify_residues(residues)  # tuple unpacking
        results[protein_name] = {
            "atoms": len(atoms), "residues": len(residues), "First five residues": (list(residues.items())[:5]),
            "residue-chain": chain_info, "residues-classes-counts": counts, "classes-chain": chain_counts, "ligands": ligands, "waters": len(waters)
        }

    return results

def extract_atom_coordinates(atom_lines):
    """
    Returns a list of atoms with coordinates.
    Each atom: dict with keys: chain, res_num, res_name, atom_name, x, y, z
    [{'chain': 'A', 'res_num': '1', 'res_name': 'THR', 'atom_name': 'N', 'x': 17.047, 'y': 14.099, 'z': 3.625} ....
    """
    atoms = []

    for line in atom_lines:
        if line.startswith("ATOM"):
            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            chain = line[21].strip()
            res_num = int(line[22:26].strip())

            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())

            atoms.append({
                "chain": chain,
                "res_num": res_num,
                "res_name": res_name,
                "atom_name": atom_name,
                "x": x,
                "y": y,
                "z": z
            })

    return atoms

def distance(atom1, atom2):
    """

    :param atom1:
    :param atom2:
    :return:
    """
    dx = atom1["x"] - atom2["x"]
    dy = atom1["y"] - atom2["y"]
    dz = atom1["z"] - atom2["z"]
    return math.sqrt(dx * dx + dy * dy + dz * dz)

def atomic_distance(atoms, cutoff=4.5): # contacts atom-atom
    """
        Returns list of tuples:
        (atom1, atom2, distance)
        """
    contacts = []

    n = len(atoms)
    for i in range(n):
        for j in range(i + 1, n):
            d = distance(atoms[i], atoms[j])
            if d <= cutoff:
                contacts.append((atoms[i], atoms[j], d))

    return contacts

def detect_residue_contacts(atoms, cutoff=4.5):
    """
    Returns set of residue pairs in contact and a liste of the residues numbers
    ((chain1, res1), (chain2, res2))
    """
    contacts = set() # avoid duplicates

    n = len(atoms)
    for i in range(n):
        for j in range(i+1, n):
            d = distance(atoms[i], atoms[j])
            if d <= cutoff:
                res1 = (atoms[i]["chain"], atoms[i]["res_num"])
                res2 = (atoms[j]["chain"], atoms[j]["res_num"])

                # remove trivial contacts
                if res1 == res2:
                    continue # a residue is always in contact with itself
                if abs(res1[1] - res2[1]) <= 1: # immediate neighbors
                    continue # skips contacts between adjacent residues covalently linked (not folding)

                contacts.add((res1, res2))

            residues = sorted({res[1] for res_pair in contacts for res in res_pair}) # {...} →(set), result: ['845', '846', '872', '910', ...]

    return contacts, residues # made a set of residues to construct the matrix of contacts
