from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley


def compute_sasa_per_residue(pdb_file):
    """
    Computes SASA per residue using Shrake–Rupley algorithm.

    Returns:
        total_sasa (float)
        sasa_per_residue (dict): {(chain, res_num): sasa}
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    sr = ShrakeRupley()
    sr.compute(structure, level="R")

    sasa_per_residue = {}
    total_sasa = 0.0

    for residue in structure.get_residues():
        # Skip hetero residues if desired
        if residue.id[0] != " ":
            continue

        chain = residue.get_parent().id
        res_num = residue.id[1]  # already int

        sasa = residue.sasa
        sasa_per_residue[(chain, res_num)] = sasa
        total_sasa += sasa

    return total_sasa, sasa_per_residue


total_sasa, sasa_res = compute_sasa_per_residue("../data/1A6M.pdb")
print(f"Total SASA: {total_sasa:.2f} Å²")
for key, value in list(sasa_res.items())[:10]:
    print(key, f"{value:.2f}")
