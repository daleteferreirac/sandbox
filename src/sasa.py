from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

def compute_total_sasa(pdb_file):
    """
    Computes total SASA of a protein structure using Shrake-Rupley algorithm.

    Returns:
        total_sasa (float): total solvent accessible surface area (Å²)
    """
    parser = PDBParser(QUIET=True) # QUIET suppresses format warnings
    structure = parser.get_structure("protein", pdb_file)

    sr = ShrakeRupley() # create Shrake–Rupley SASA calculator
    sr.compute(structure, level="A")  # A = atom level, compute SASA at atom level and store values

    total_sasa = 0.0
    for atom in structure.get_atoms():
        total_sasa += atom.sasa

    return total_sasa

sasa = compute_total_sasa("../data/1A6M.pdb")
print(f"Total SASA: {sasa:.2f} Å²")
