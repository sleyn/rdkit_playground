# rdkit_playground

The goal for this repo was to get comfortable with the core RDKit workflows.

Everything lives in a single [marimo](https://marimo.io/) notebook, reactive and runnable end-to-end, walking through the basics one topic at a time.

## Stack

- Python 3.13
- [RDKit](https://www.rdkit.org/) — cheminformatics
- [marimo](https://marimo.io/) — reactive notebook
- [pandas](https://pandas.pydata.org/) — tabular data
- [uv](https://docs.astral.sh/uv/) — dependency management

## Getting started

```bash
uv sync
uv run marimo edit notebooks/rdkit_notebook.py   # interactive, editable
# or
uv run marimo run notebooks/rdkit_notebook.py    # read-only app view
```

## What's inside

The notebook ([notebooks/rdkit_notebook.py](notebooks/rdkit_notebook.py)) covers:

1. **Parsing and canonicalization** — three different SMILES for benzoic acid collapse to one canonical string; handling invalid SMILES; InChIKey generation; inspecting a molecule's atom graph.
2. **Standardization** — `Cleanup` → `LargestFragmentChooser` → `Uncharger` → tautomer canonicalization, the sequence that makes a compound from ChEMBL match the same compound from a vendor catalog or a patent.
3. **Descriptors and fingerprints** — MW, cLogP, TPSA, HBD/HBA, rotatable bonds, QED, and Lipinski's Rule of Five; Morgan (ECFP) fingerprints and Tanimoto similarity between drugs.
4. **SMARTS, scaffolds, and alerts** — substructure matching with SMARTS patterns, Bemis-Murcko scaffold extraction, and PAINS (Pan-Assay Interference) filtering.
5. **3D conformers** — ETKDGv3 embedding, MMFF94 energy minimization, RMSD between conformers, and SDF export.
6. **Pipeline** — stitches the above into a small annotation pipeline: raw SMILES → standardize → dedupe on InChIKey → descriptors/fingerprints → scaffold-grouped split, with failure counting throughout.

## Layout

- `notebooks/rdkit_notebook.py` — the marimo notebook (source of truth; all work happens here)
- `data/` — artifacts produced by the notebook (e.g. `ibuprofen_best.sdf`, the lowest-energy 3D conformer)

## License

Apache 2.0 — see [LICENSE](LICENSE).
