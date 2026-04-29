# Helper functions for the structural analyzer
import os
import math
from Bio.PDB import PDBParser, NeighborSearch

def parse_pdb(filepath):
    """
    Parse ATOM and HETATM records into structured dictionaries.
    list of dict:
        Each atom contains:
        - record (ATOM or HETATM)
        - atom_name
        - res_name
        - chain
        - res_num (int)
        - x, y, z (float coordinates)
    example: [{'record': 'ATOM', 'atom_name': 'N', 'res_name': 'VAL', 'chain': 'A', 'res_num': 1, 'x': -4.004, 'y': 15.224, 'z': 13.636}, ...]
    """

    atoms = []

    with open(filepath, "r") as f:
        for line in f:

            record = line[:6].strip()

            # only consider ATOM and HETATM records, including ligands
            if record not in {"ATOM", "HETATM"}:
                continue

            atom = {
                "record": record,
                "atom_name": line[12:16].strip(),
                "res_name": line[17:20].strip(),
                "chain": line[21].strip(),
                "res_num": int(line[22:26].strip()),
                "x": float(line[30:38].strip()),
                "y": float(line[38:46].strip()),
                "z": float(line[46:54].strip()),
            }

            atoms.append(atom)

    return atoms

def extract_residues(atoms):
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

    for atom in atoms:

        key = (atom["chain"], atom["res_num"])

        if atom["record"] == "ATOM":
            if key not in residues:
                residues[key] = atom["res_name"] # {('A', 1): 'VAL', ('A', 2): 'LEU', ('A', 3): 'SER', ('A', 4): 'GLU', ...}
                chain_info[atom["chain"]] = chain_info.get(atom["chain"], 0) + 1 # {'A': 151}

        elif atom["record"] == "HETATM":
            if atom["res_name"] == "HOH":
                waters[key] = atom["res_name"] # {('A', 1002): 'HOH', ('A', 1003): 'HOH', ('A', 1004): 'HOH',...}
            else:
                ligands[key] = atom["res_name"] # example of 4AG8: {('A', 2000): 'AXI'}

    return residues, chain_info, waters, ligands

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

        atoms = parse_pdb(path)
        # extract amino acid residues and hetatm
        residues, chain_info, waters, ligands = extract_residues(atoms)
        counts, chain_counts = classify_residues(residues)  # tuple unpacking
        results[protein_name] = {
            "atoms": len(atoms), "residues": len(residues), "First five residues": (list(residues.items())[:5]),
            "residue-chain": chain_info, "residues-classes-counts": counts, "classes-chain": chain_counts, "ligands": ligands, "waters": len(waters)
        }

    return results

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


def detect_residue_contacts(atoms, cutoff=4.5):
    """
    Returns set of residue pairs in contact and a liste of the residues numbers
    ((chain1, res1), (chain2, res2))
    """
    contacts = set() # avoid duplicates

    n = len(atoms)

    # use j = i + 1 to avoid comparing an atom with itself and duplicate comparisons (i,j) and (j,i)
    # double loop over all atom pairs.
    for i in range(n):
        for j in range(i+1, n):

            # Ignore non-protein atoms
            if atoms[i]["record"] != "ATOM" or atoms[j]["record"] != "ATOM":
                continue

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

def detect_residues_contacts_fast(pdb_file, cuttoff=4.5):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("prot", pdb_file)    # create objects Atom with .coord, .get_parent(), etc.

    all_atoms = list(structure.get_atoms())  # objects from Atom into a list

    ns = NeighborSearch(all_atoms) # kd tree with coodernates

    contacts = set()

    for atom in all_atoms: # for each atom, search neighbors whithin the cutoff

        if atom.get_parent().id[0] != " ": # ignores hetreroatoms
            continue

        neighbors = ns.search(atom.coord, cutoff, level='A')  # .search() retorn list of atoms whitin the radius at the level 'A'atom

        res1 = atom.get_parent() # atom.get_parent() retorn objet Residue from their atom

        for neighbor in neighbors:

            if neighbor.get_parent().id[0] != " ": # ignores hetatom from the neighbors
                continue

            res2 = neighbor.get_parent()

            # ignores contact itself
            if res1 == res2:
                continue

            num1 = res1.id[1]   #sequence of each atom: res1.id = (' ', 47, ' '), [1] is the number
            num2 = res2.id[1]

            # immediate neighbors
            if abs(num1 - num2) <= 1:
                continue

            # chain ID: 'A', 'B', etc.
            chain1 = res1.get_parent().id
            chain2 = res2.get_parent().id

            # tuple(sorted(...)) garantees that (A,4)-(A,79) and (A,79)-(A,4) is the same key
            key = tuple(sorted([(chain1, num1), (chain2, num2)]))
            contacts.add(key)

    residues = sorted({res[1] for pair in contacts for res in pair})

    return contacts,  residues
