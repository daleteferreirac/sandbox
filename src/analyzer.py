# Main module of the structural PDB analyzer project.
from utils import analyze_proteins

protein1 = "../data/1A3N.pdb"
protein2 = "../data/1CRN.pdb"
analyze_proteins(protein1, protein2)