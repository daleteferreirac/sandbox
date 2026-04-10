from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from utils import parse_pdb,detect_residue_contacts

atoms = parse_pdb("../data/1A6M.pdb")
contacts, _ = detect_residue_contacts(atoms) # e.g: (('A', 4), ('A', 79))

contact_count = {} # count contacts per residue

for res1, res2 in contacts:
    # .get(key, 0) avoids KeyError and initializes missing keys with 0
    contact_count[res1] = contact_count.get(res1, 0) + 1 # key-> res1, value-> contacts that each res have
    contact_count[res2] = contact_count.get(res2, 0) + 1

MAX_SASA = {
    "ALA": 129, "VAL": 174, "LEU": 201, "ILE": 197,
    "MET": 224, "PHE": 240, "TRP": 285, "PRO": 159,
    "SER": 155, "THR": 172, "ASN": 195, "GLN": 223,
    "TYR": 263, "CYS": 167,
    "ARG": 274, "LYS": 236, "HIS": 224, "ASP": 193, "GLU": 223
} # # maxi SASA values for each amino acid (used for normalization)

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

        # extract residue identifiers
        chain = residue.get_parent().id # Chain identifier (e.g. 'A')
        res_num = residue.id[1]  # already int

        res_name = residue.get_resname()
        sasa = residue.sasa

        # classify residue exposure based on relative sasa
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
# count how many residues fall into each exposure classe
buried = sum(1 for v in sasa_res.values() if v["exposure"] == "buried")
exposed = sum(1 for v in sasa_res.values() if v["exposure"] == "exposed")
intermediate = sum(1 for v in sasa_res.values() if v["exposure"] == "intermediate")

total = len(sasa_res)

print("\nExposure summary:")
print(f"Buried: {buried} ({buried/total*100:.1f}%)")
print(f"Exposed: {exposed} ({exposed/total*100:.1f}%)")
print(f"Intermediate: {intermediate} ({intermediate/total*100:.1f}%)\n")

for key, value in list(sasa_res.items())[:10]:
    print(
        key,
        f"abs={value['absolute']:.2f}",
        f"rel={value['relative']:.2f}" if value["relative"] is not None else "rel=None",
        f"class={value['exposure']}"
    )

# Combine SASA with contact information
for key in list(sasa_res.keys())[:10]:
    rel = sasa_res[key]["relative"]
    contacts = contact_count.get(key, 0)
# format relative sasa safely (avoid none formatting error)
    if rel is not None:
        rel_str = f"{rel:.2f}"
    else:
        rel_str = "None"

    print(key, f"rel={rel_str}", f"contacts={contacts}")