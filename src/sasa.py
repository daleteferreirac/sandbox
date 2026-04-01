from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

MAX_SASA = {
    "ALA": 129, "VAL": 174, "LEU": 201, "ILE": 197,
    "MET": 224, "PHE": 240, "TRP": 285, "PRO": 159,
    "SER": 155, "THR": 172, "ASN": 195, "GLN": 223,
    "TYR": 263, "CYS": 167,
    "ARG": 274, "LYS": 236, "HIS": 224, "ASP": 193, "GLU": 223
}

def compute_sasa_per_residue(pdb_file):
    """
    Computes SASA per residue using Shrake–Rupley algorithm.

    Returns:
        total_sasa (float)
        sasa_per_residue (dict): {(chain, res_num): sasa}
    """
    parser = PDBParser(QUIET=True) # quiet suppress format warnings
    structure = parser.get_structure("protein", pdb_file)

    sr = ShrakeRupley() # Initialize Shrake–Rupley SASA calculator
    sr.compute(structure, level="R") # R = residue level (stored in residue.sasa)

    sasa_per_residue = {}
    total_sasa = 0.0

    for residue in structure.get_residues():
        # Skip hetero residues (ligands, ions, water)
        if residue.id[0] != " ":
            continue

        chain = residue.get_parent().id # Chain identifier (e.g. 'A')
        res_num = residue.id[1]  # already int

        res_name = residue.get_resname()
        sasa = residue.sasa

        if res_name in MAX_SASA:
            relative = sasa / MAX_SASA[res_name]
        else:
            relative = None  # fallback for unknown residues

        if relative is not None:
            if relative < 0.2:
                exposure = "buried"
            elif relative > 0.5:
                exposure = "exposed"
            else:
                exposure = "intermediate"
        else:
            exposure = None

        sasa_per_residue[(chain, res_num)] = {
            "absolute": sasa,
            "relative": relative,
            "exposure": exposure
        }
        total_sasa += sasa

    return total_sasa, sasa_per_residue

total_sasa, sasa_res = compute_sasa_per_residue("../data/1A6M.pdb")

print(f"Total SASA: {total_sasa:.2f} Å²")

for key, value in list(sasa_res.items())[:10]:
    print(
        key,
        f"abs={value['absolute']:.2f}",
        f"rel={value['relative']:.2f}" if value["relative"] else "rel=None",
        f"class={value['exposure']}"
    )