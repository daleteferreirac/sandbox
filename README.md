# Structural PDB Analyzer (learning project)

This project is a Python tool to analyze basic structural information from PDB files.

---

## What this project does (so far)

Given one or more PDB files, the code extracts and organizes:

- Total number of atoms
- Number of protein residues
- Residue count per chain
- Classification of residues into:
  - Hydrophobic
  - Polar
  - Charged
- Percentage of each residue class (global and per chain)
- Identification and separation of:
  - Protein residues (ATOM)
  - Water molecules (HOH)
  - Small-molecule ligands (HETATM, excluding HOH)

Water molecules and ligands are excluded from residue-based structural calculations.
All results are returned as Python dictionaries.

---

## Notes

- The focus is educational. 