import os
# Helper functions for the structural analyzer
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
            res_num = line[22:26].strip()    # 1, 2, 3, 8, 54...

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

    for (chain, _),res_name in residues.items():
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

