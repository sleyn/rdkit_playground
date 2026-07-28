# rdkit_playground
First meeting with RDKit

## Plan

### :black_square_button: Block 1 Parsing and canonicalization
Three different SMILES for benzoic acid collapse to one canonical string. This is your entry point: canonical SMILES and InChIKey are the join keys, and they're the direct analogue of normalized CHROM:POS:REF:ALT. The gotcha to internalize is that MolFromSmiles returns None on bad input rather than raising — every production loop needs a null check, and RDKit logs parse errors to stderr independently of your exception handling.

### :black_square_button: Block 2 Standardization.
Spend the most time here. Cleanup → LargestFragmentChooser → Uncharger → tautomer canonicalization is the sequence that makes a compound from ChEMBL match the same compound from a vendor catalog or a patent. Skip it and your dedup silently fails, your similarity search misses, and your ML training set has leakage. This is the block that gives you your best interview line, because it's structurally identical to left-alignment and multi-allelic splitting.

### :black_square_button: Block 3 Descriptors and fingerprints.
MW, cLogP, TPSA, HBD/HBA, rotatable bonds, QED, Lipinski. Then Morgan fingerprints via rdFingerprintGenerator — note that AllChem.GetMorganFingerprintAsBitVect still works but is the legacy path, and that radius 2 means ECFP4, because the ECFP number is the diameter. Get comfortable with BulkTanimotoSimilarity; count fingerprints usually beat bit vectors as ML features.

### :black_square_button: Block 4SMARTS, scaffolds, alerts
SMARTS is regex for molecules — you'll pick it up fast. The important concept is Murcko scaffolds, because scaffold splits are the leakage-safe way to build train/test sets. Random splits inflate ADMET model performance the same way splitting variants without holding out whole genomic regions inflates a classifier. If you make one point about ML methodology in the technical round, make it this one.

### :black_square_button: Block 5 3D conformers
ETKDGv3 embedding, MMFF optimization, SDF output — this is literally what you feed a docking program, so it's your handoff into the Tier 2 Smina exercise. Two gotchas that trip up everyone: add explicit hydrogens before embedding, and always set randomSeed, since embedding is stochastic and unseeded runs aren't reproducible.

### :black_square_button: Block 6 Assemble the pipeline
Raw SMILES → standardize → dedup on InChIKey → descriptors → fingerprints → scaffold-grouped split, with failure counting throughout. Eight inputs in, six unique compounds out, one parse failure, one salt collapsed onto its parent. That's a small annotation pipeline, and it's the concrete thing you point at when someone asks what you've actually done with cheminformatics.