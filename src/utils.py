# Helper functions for the structural analyzer

def load_pdb(filepath):
    """
    Loads a PDB file.
    filepath: path to PDB file
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
    Returns a dictionary:
        { (chain, residue_number) : residue_name }
        the key is a tuple of chain and residue number.
        the value is a string of the residue name.
    chain_info dictionary:
        {chain: number of residues}
    """
    residues = {}
    chain_info = {}

    for line in atom_lines:
        res_name = line[17:20].strip()   # VAL, ALA, ARG...
        chain = line[21].strip()         # A, B, C...
        res_num = line[22:26].strip()    # 1, 2, 3, 8, 54...

        key = (chain, res_num)
        if key not in residues:
            residues[key] = res_name # ('A', '1'): 'THR', (('A', '2'), 'LEU'), (('A', '3'), 'SER'),...

            if chain not in chain_info:
                chain_info[chain] = 1
            else:
                chain_info[chain] += 1

    return residues, chain_info

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
    returns: dict {class_residue: residues}, dict {class_residue: count}
    """
    hydrophobics = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"}
    polars = {"SER", "THR", "ASN", "GLN", "TYR", "CYS"}
    charged = {"ARG", "LYS", "HIS", "ASP", "GLU"}

    classes = {
        "HYDROPHOBICS": set(),
        "POLAR": set(),
        "CHARGED": set()
    }

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
            classes["HYDROPHOBICS"].add(res_name)
            counts["HYDROPHOBICS"] += 1
            chain_counts[chain]["HYDROPHOBICS"] += 1
        if res_name in polars:
            classes["POLAR"].add(res_name)
            counts["POLAR"] += 1
            chain_counts[chain]["POLAR"] += 1
        if res_name in charged:
            classes["CHARGED"].add(res_name)
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

    return classes, counts, chain_counts

def multiple_protein_analyzer(*filepaths):
    """
    Receives multiple PDB file paths and returns the information organized.
    :param filepaths: tuple of PDB file paths
    :return:
    """
    for path in filepaths:
        protein_name = path.split("/")[-1]
        print(f"Protein: {protein_name}")
        atoms = load_pdb(path)
        print(f"    Number of atoms: {len(atoms)}")


    return

