# rdkit_playground
First meeting with RDKit

## Plan

### :black_square_button: Block 1 Parsing and canonicalization
Three different SMILES for benzoic acid collapse to one canonical string.

### :black_square_button: Block 2 Standardization.
Spend the most time here. `Cleanup` → `LargestFragmentChooser` → `Uncharger` → tautomer canonicalization is the sequence that makes a compound from ChEMBL match the same compound from a vendor catalog or a patent.

### :black_square_button: Block 3 Descriptors and fingerprints.
MW, cLogP, TPSA, HBD/HBA, rotatable bonds, QED, Lipinski.

### :black_square_button: Block 4SMARTS, scaffolds, alerts
SMARTS is regex for molecules.

### :black_square_button: Block 5 3D conformers
ETKDGv3 embedding, MMFF optimization, SDF output.

### :black_square_button: Block 6 Assemble the pipeline
Raw SMILES → standardize → dedup on InChIKey → descriptors → fingerprints → scaffold-grouped split, with failure counting throughout. Eight inputs in, six unique compounds out, one parse failure, one salt collapsed onto its parent. That's a small annotation pipeline, and it's the concrete thing you point at when someone asks what you've actually done with cheminformatics.