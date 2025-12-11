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
    """
    residues = {}

    for line in atom_lines:
        res_name = line[17:20].strip()   # VAL, ALA, ARG...
        chain = line[21].strip()         # A, B, C...
        res_num = line[22:26].strip()    # 1, 2, 3, 8, 54...

        residues[(chain, res_num)] = res_name

    return residues