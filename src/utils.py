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
    """
    residues = {}

    for line in atom_lines:
        res_name = line[17:20].strip()   # VAL, ALA, ARG...
        chain = line[21].strip()         # A, B, C...
        res_num = line[22:26].strip()    # 1, 2, 3, 8, 54...

        residues[(chain, res_num)] = res_name

    return residues

def count_residue_types(residues):
    """
    Counts how many times each amino acid appears.
    residues: dict {(chain, res_num): res_name}
    returns: dict {res_name: count}
    """
    counts = {}

    for res_name in residues.values():
        if res_name not in counts:
            counts[res_name] = 1
        else:
            counts[res_name] += 1

    return counts

def classify_residues(residues):
    """
    Classify residues based on hydrophobic, polar and charged.
    use the dictionary residues: {(chain, res_num): res_name}
    returns: dict {class_residue: count}
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

    for res_name in residues.values():
        if res_name in hydrophobics:
            classes["HYDROPHOBICS"].add(res_name)
            counts["HYDROPHOBICS"] += 1
        if res_name in polars:
            classes["POLAR"].add(res_name)
            counts["POLAR"] += 1
        if res_name in charged:
            classes["CHARGED"].add(res_name)
            counts["CHARGED"] += 1

    return classes, counts

def add_count_chain():
    pass