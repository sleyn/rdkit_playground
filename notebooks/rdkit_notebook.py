import marimo

__generated_with = "0.23.14"
app = marimo.App()


@app.cell
def _(mo):
    mo.md(r"""
    # RDKit exploration

    As I see a lot of jobs with cheminformatics requirements I decided to make a small repo to get familiarized with one of the most used cheminformatics software [RDKit](https://www.rdkit.org/).

    I've aseked the Claude to give me some mini tasks to solve that would teach me RDKit basics.

    ## 1. Parsing and canonicalization

    Three different SMILES for benzoic acid collapse to one canonical string.
    """)
    return


@app.cell
def _():
    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem import (
        Descriptors, Draw, QED, rdFingerprintGenerator,
        rdMolDescriptors, Crippen, rdMolAlign, Mol
    )
    from rdkit.Chem.MolStandardize import rdMolStandardize
    from rdkit.Chem.Scaffolds import MurckoScaffold
    from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
    import pandas as pd
    from pandas import DataFrame

    return (
        Chem,
        Crippen,
        DataFrame,
        DataStructs,
        Descriptors,
        FilterCatalog,
        FilterCatalogParams,
        Mol,
        MurckoScaffold,
        QED,
        RDLogger,
        pd,
        rdFingerprintGenerator,
        rdMolAlign,
        rdMolDescriptors,
        rdMolStandardize,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Make a cannonical SMILES

    Three different SMILES formulas represent the same benzoic acid.
    """)
    return


@app.cell
def _(Chem):

    variants: list[str] = ["c1ccccc1C(=O)O", "OC(=O)c1ccccc1", "C1=CC=C(C=C1)C(O)=O"]
    canonnical_forms: set[str] = {Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in variants}
    print(f'Cannonical forms: {",".join(canonnical_forms)}\nNumber of cannonical forms: {len(canonnical_forms)}')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Work with an incorrect output
    """)
    return


@app.cell
def _(Chem, Mol):
    # Bad input gives None
    bad: Mol = Chem.MolFromSmiles("C1CC")
    print(f'Bad input processing gives "{bad}"')
    return


@app.cell
def _(Chem, RDLogger):
    # Silence the parse-error logging
    # Try a SMILES with a wrong valency
    RDLogger.DisableLog("rdApp.error")
    assert Chem.MolFromSmiles("n1cccc1") is None
    return


@app.cell
def _(Chem, RDLogger):
    RDLogger.EnableLog("rdApp.error")
    assert Chem.MolFromSmiles("n1cccc1") is None
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### InChIKey

    InChIKey is a hash-based canonical identifier of a SMILES structure.
    """)
    return


@app.cell
def _(Chem, Mol):
    molecule: Mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")   # aspirin
    print(f'InChIKey: {Chem.MolToInchiKey(molecule)}')
    return (molecule,)


@app.cell
def _(molecule: "Mol"):
    print(f'atoms: {molecule.GetNumAtoms()} | heavy: {molecule.GetNumHeavyAtoms()}')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Atom-level overview

    `Mol.GetNumAtoms` gives number of nodes in the molecule graph and `GetNumHeavyAtoms` gives number of nodes that are not hydrogen. As SMILES does not contain hydrogen RDKit builds implicit-hydrogen graph.

    To build explicit-hydrogen graph add hydrogens using `Mol.AddHs`
    """)
    return


@app.cell
def _(Chem, Mol, molecule: "Mol"):
    # Add explicit Hydrogen nodes to the graph
    molecule_with_h: Mol = Chem.AddHs(molecule)
    print(f'atoms: {molecule_with_h.GetNumAtoms()} | heavy: {molecule_with_h.GetNumHeavyAtoms()}')

    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Similar, RDKit calculates degree of an atom without counting hydrogens if hydrogens were not explicitly added.
    """)
    return


@app.cell
def _(molecule: "Mol"):
    for atom in list(molecule.GetAtoms()):
        print(f"idx={atom.GetIdx()} {atom.GetSymbol()} "
              f"arom={atom.GetIsAromatic()} degree={atom.GetDegree()} "
              f"implicitH={atom.GetTotalNumHs()}")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 2. Standardization

    Tautomer canonicalization is the sequence that makes a compound from ChEMBL match the same compound from a vendor catalog or a patent

    Use several preprocessing modules:

    `rdMolStandardize.Cleanup()`:
    1. Disconnect Organometallics `DisconnectOrganometallics`

       Disconnects metal atoms bonded to non-metals (e.g., organometallic complexes or salts). Leaves the metal and ligand as separate fragments in the molecule structure. RDKit treats metals in SMILES as covalently bonded by default as thir representation does not show a type of bond.

    2. Normalize Functional Groups `Normalize`

       Applies a series of SMARTS-based transformations (normalizations) to convert non-standard or ambiguous functional group representations into standard, canonical forms.

    3. Re-ionize / Balance Charges `Reionize`

       Adjusts charges across functional groups to place charges on the most acidic/basic centers.

    `rdMolStandardize.LargestFragmentChooser` is used to select and extract the primary organic or active fragment from a multi-component molecular structure.

    `rdMolStandardize.Uncharger` is designed to neutralize charged functional groups in a molecule whenever possible by adding or removing protons ($H^+$).
    """)
    return


@app.cell
def _(Chem, Mol, rdMolStandardize):
    dirty: list[str] = [
            "CC(=O)Oc1ccccc1C(=O)[O-].[Na+]",   # aspirin sodium salt
            "CC(=O)Oc1ccccc1C(=O)O",            # aspirin free acid
            "CN1CCC[C@H]1c1cccnc1.Cl",          # nicotine HCl
    ]

    lfc = rdMolStandardize.LargestFragmentChooser()
    uncharger = rdMolStandardize.Uncharger()

    def standardize(smiles):
        m: Mol = Chem.MolFromSmiles(smiles)
        if m is None:
            return None
        m = rdMolStandardize.Cleanup(m)     # sanitize, normalize, disconnect metals
        m = lfc.choose(m)                   # drop counterions / solvent
        m = uncharger.uncharge(m)           # neutralize where chemically sensible

        return Chem.MolToSmiles(m)

    return dirty, lfc, standardize


@app.cell
def _(dirty: list[str], standardize):
    for s in dirty:
            print(f"  {s:38s} -> {standardize(s)}")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Tautomers

    `rdMolStandardize.TautomerEnumerator()` Returns the canonical tautomer for a molecule. Uses scoring from M. Sitzmann et al., “Tautomerism in Large Databases.”, JCAMD 24:521 (2010) [https://doi.org/10.1007/s10822-010-9346-4](https://doi.org/10.1007/s10822-010-9346-4) to select the default one.
    """)
    return


@app.cell
def _(Chem, rdMolStandardize):
    taut_a, taut_b = "O=C1CCCCC1", "OC1=CCCCC1"   # cyclohexanone / its enol
    enumerator = rdMolStandardize.TautomerEnumerator()
    canon_a: str = Chem.MolToSmiles(enumerator.Canonicalize(Chem.MolFromSmiles(taut_a)))
    canon_b: str = Chem.MolToSmiles(enumerator.Canonicalize(Chem.MolFromSmiles(taut_b)))
    print(f"  tautomers -> {canon_a} | {canon_b} | match={canon_a == canon_b}")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 3. Descriptors and fingerprints

    * **MW** - Molecular weght
    * **cLogP** - Calculated Octanol-Water Partition Coefficient. The logarithm of the partition coefficient ($P$) between n-octanol (a hydrophobic organic solvent) and water. It quantifies lipophilicity (fat-solubility vs. water-solubility).
    * **TPSA** - Topological Polar Surface Area. The surface area ($\text{Å}^2$) occupied by polar atoms (predominantly oxygen, nitrogen, and their attached hydrogen atoms). Polar regions interact strongly with water molecules. High TPSA impedes a molecule's ability to passively penetrate cell membranes.
    * **HBD/HBA** - Hydrogen Bond Donors & Acceptors. **HBD (Donors):** The count of hydrogen atoms attached to electronegative heteroatoms (typically $-\text{OH}$ and $-\text{NH}$ groups). **HBA (Acceptors):** The count of electronegative heteroatoms with lone pairs that can accept hydrogen bonds (typically $\text{O}$ and $\text{N}$ atoms). Hydrogen bonding with water molecules must be broken for a compound to move from an aqueous environment into a lipid membrane. Excess donors or acceptors make membrane permeation energetically unfavorable.
    * **rotatable bonds (RotB)** - The number of non-ring, single bonds connected to non-hydrogen heavy atoms (excluding terminal bonds like $-\text{CH}_3$ or amide $\text{C-N}$ bonds due to partial double-bond character). Measures molecular flexibility.
    * **QED** - Quantitative Estimate of Drug-likeness. A continuous metric (ranging from $0$ to $1$) proposed by Bickerton et al. that quantifies overall drug-likeness.
      * $1.0$: Ideal, highly drug-like property profile.
      * $0.0$: Very poor drug-likeness profile.
    * **Lipinski** - Lipinski's Rule of Five. A set of empirical guidelines formulated by Christopher Lipinski to evaluate whether a chemical compound has chemical and physical properties that make it likely to be an orally active drug in humans. The rules are *multiples of 5*, not *5 rules*.
      * MW $\le 500 \text{ Da}$
      * cLogP $\le 5$
      * HBD $\le 5$
      * HBA $\le 10$
    """)
    return


@app.cell
def _():
    DRUGS: dict[str, str] = {
        "aspirin":     "CC(=O)Oc1ccccc1C(=O)O",
        "ibuprofen":   "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        "naproxen":    "COc1ccc2cc(ccc2c1)C(C)C(=O)O",
        "paracetamol": "CC(=O)Nc1ccc(O)cc1",
        "caffeine":    "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
        "imatinib":    "Cc1ccc(cc1Nc1nccc(n1)-c1cccnc1)NC(=O)c1ccc(CN2CCN(C)CC2)cc1",
    }
    return (DRUGS,)


@app.cell
def _(
    Chem,
    Crippen,
    DRUGS: dict[str, str],
    DataFrame,
    Descriptors,
    Mol,
    QED,
    pd,
    rdMolDescriptors,
):
    rows = []
    for name, smi in DRUGS.items():
        m: Mol = Chem.MolFromSmiles(smi)
        rows.append({
            "name": name,
            "MW": round(Descriptors.MolWt(m), 1),
            "cLogP": round(Crippen.MolLogP(m), 2),
            "TPSA": round(rdMolDescriptors.CalcTPSA(m), 1),
            "HBD": rdMolDescriptors.CalcNumHBD(m),
            "HBA": rdMolDescriptors.CalcNumHBA(m),
            "RotB": rdMolDescriptors.CalcNumRotatableBonds(m),
            "QED": round(QED.qed(m), 3),
        })
    
    df: DataFrame = pd.DataFrame(rows)
    df
    return (df,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Add Lipinski Ro5 as an explicit filter (the canonical early triage rule).
    """)
    return


@app.cell
def _(df: "DataFrame"):
    df["Ro5_violations"] = (
        (df.MW > 500).astype(int) + (df.cLogP > 5).astype(int)
        + (df.HBD > 5).astype(int) + (df.HBA > 10).astype(int)
    )
    print("Ro5 violations:\n")
    df[["name", "Ro5_violations"]]
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Fingerprints

    In cheminformatics, molecular fingerprints are digital representations of chemical structures encoded as fixed-length bit vectors (arrays of 0s and 1s) or integer vectors.

    **RDKit Fingerprint** are analyze linear or branched paths of atoms and bonds in the molecular graph up to a certain distance (e.g., 7 bonds), hash these paths, and set bits in a vector.

    Once molecules are converted into binary fingerprints, their similarity is most commonly quantified using the Tanimoto Coefficient (Jaccard index):$$T(A, B) = \frac{\vert{}A \cap B\vert{}}{\vert{}A \cup B\vert{}} = \frac{c}{a + b - c}$$Where:$c$ = number of bits set to 1 in both molecules $A$ and $B$.$a$ = number of bits set to 1 in molecule $A$.$b$ = number of bits set to 1 in molecule $B$.A Tanimoto score ranges from 0.0 (no shared features) to 1.0 (identical bit patterns, though not necessarily identical molecules due to hash collisions or lack of stereochemistry awareness). In drug discovery, an ECFP4 Tanimoto score $\ge 0.85$ between two molecules often suggests they may share similar biological activity.
    """)
    return


@app.cell
def _(Chem, DRUGS: dict[str, str], DataStructs, Mol, rdFingerprintGenerator):
    from rdkit.Chem.AllChem import FingerprintGenerator64
    from rdkit.DataStructs.cDataStructs import ExplicitBitVect, UIntSparseIntVect

    gen: FingerprintGenerator64 = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    mols: dict[str, Mol] = {n: Chem.MolFromSmiles(s) for n, s in DRUGS.items()}
    fps: dict[str, ExplicitBitVect] = {n: gen.GetFingerprint(m) for n, m in mols.items()}

    names: list[str] = list(fps)
    ref: ExplicitBitVect = fps["ibuprofen"]
    sims: list[float] = DataStructs.BulkTanimotoSimilarity(ref, [fps[n] for n in names])
    print("\nTanimoto vs ibuprofen:")
    for n, sim in sorted(zip(names, sims), key=lambda t: -t[1]):
        print(f"  {n:12s} {sim:.3f}")
    return FingerprintGenerator64, UIntSparseIntVect, gen, mols


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Count-based fingerprints are an alternative for bit vectors. Some claim that they work better for ML models.
    """)
    return


@app.cell
def _(
    UIntSparseIntVect,
    gen: "FingerprintGenerator64",
    mols: "dict[str, Mol]",
):
    cfp: UIntSparseIntVect = gen.GetCountFingerprint(mols["imatinib"])
    print(f'nonzero count-fp elements: {len(cfp.GetNonzeroElements())}')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 4. Substructure, SMARTS, scaffolds, alert filters

    ### Substructures

    `Mol.HasSubstructMatch()` checks whether a molecule contains a specific chemical pattern or substructure defined by another molecule or SMARTS string.
    """)
    return


@app.cell
def _(Chem, Mol, mols: "dict[str, Mol]"):
    carboxylic_acid: Mol = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    hits: list[str] = [n for n, m in mols.items() if m.HasSubstructMatch(carboxylic_acid)]
    print("Molecules that have carboxylic acid:", hits)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Count aromatic nitrogens.
    """)
    return


@app.cell
def _(Chem, Mol, mols: "dict[str, Mol]"):
    aromatic_n: Mol = Chem.MolFromSmarts("[n]")
    print(f'aromatic N counts: { {n: len(m.GetSubstructMatches(aromatic_n)) for n, m in mols.items()} }')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Scaffolds

    `MurckoScaffold.MurckoScaffoldSmiles()` generates the Bemis-Murcko Scaffold for each molecule by stripping away all non-ring side chains, leaving only the ring systems and the linkers connecting them.

    In drug discovery and machine learning, Murcko scaffolds are essential for **scaffold splitting**. This ensuring structural family held-outs in train/test splits to evaluate how models generalize to novel chemical series, avoiding data leakage. Scaffold splitting groups data points by a shared core structural feature or graph topology so that entire structural groups go to either the train or test set, but never both.
    """)
    return


@app.cell
def _(MurckoScaffold, mols: "dict[str, Mol]"):
    print("scaffolds:")
    for n2, m2 in mols.items():
        print(f"  {n2:12s} {MurckoScaffold.MurckoScaffoldSmiles(mol=m2)}")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Pan-Assay Interference \(PAINS\) Alert Screening

    Pan-assay interference compounds (PAINS) are chemical molecules that frequently yield false positives in high-throughput drug screening.

    In early-stage drug discovery, researchers screen large libraries of chemicals to find ones that selectively interact with a biological target. PAINS contain disruptive functional groups—such as catechols, enones, rhodanines, and toxoflavins—that interfere with the testing process itself rather than the disease pathway.

    Common mechanisms of interference include:
    * **Protein Reactivity:** Covalently binding to thiols (like cysteine residues) or amines on multiple proteins.
    * **Redox Activity:** Generating reactive oxygen species (ROS) or hydrogen peroxide.
    * **Aggregation:** Forming microscopic particles that nonspecifically trap and inhibit assay proteins.
    * **Chelation & Fluorescence:** Binding to trace metals in the assay buffer or interfering with typical light-absorption or fluorescence readouts.

    `FilterCatalogParams()` creates the parameters container.

    `params.AddCatalog()` loads the PAINS rule sets. PAINS filters were originally defined by Baell & Holloway (2010) to identify chemical substructures that frequently yield false-positive hits in high-throughput biological screening assays.

    `FilterCatalog(params)` instantiates the compiled filter catalog.

    `catalog.HasMatch()` checks whether the molecule m contains any substructure that matches the loaded PAINS catalog.
    """)
    return


@app.cell
def _(FilterCatalog, FilterCatalogParams, mols: "dict[str, Mol]"):
    params = FilterCatalogParams() 
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog(params)
    flagged = {n: catalog.HasMatch(m) for n, m in mols.items()}
    print("\nPAINS flagged:", {k: v for k, v in flagged.items() if v} or "none")

    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 3D conformers

    ETKDGv3 Algorithm (`ETKDGv3()`): Uses Experimental Torsion-Angle Preferences with Knowledge-based Distance Geometry (v3). The standard algorithm in RDKit for 3D conformation sampling based on experimental crystal structure data.

    Random seed `params_ibuprofen.randomSeed` is needed as 3D coordinate generation uses stochastic distance geometry. Setting a random seed guarantees exact reproducibility across runs.

    `EmbedMultipleConfs()` generates an ensemble of 20 distinct 3D spatial arrangements (conformers - conformational isomers) for the molecule and returns a list of their Conformer IDs (cids).
    """)
    return


@app.cell
def _(Chem, DRUGS: dict[str, str], Mol):
    from rdkit.Chem.AllChem import EmbedParameters
    from typing import Sequence
    from rdkit.Chem.rdDistGeom import ETKDGv3, EmbedMultipleConfs
    ibuprofen: Mol = Chem.MolFromSmiles(DRUGS["ibuprofen"])

    # !!!IMPORTANT!!! Add explicit Hs BEFORE embedding or your geometry is wrong.
    ibuprofen_explicit_hydro: Mol = Chem.AddHs(ibuprofen)
    params_ibuprofen: EmbedParameters = ETKDGv3()
    params_ibuprofen.randomSeed = 0xC0FFEE          # embedding is stochastic — always seed
    cids: Sequence[int] = EmbedMultipleConfs(
            mol=ibuprofen_explicit_hydro,
            numConfs=20,
            params=params_ibuprofen
        )
    print(f'conformers generated: {len(cids)}')
    return cids, ibuprofen_explicit_hydro


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    `MMFFOptimizeMoleculeConfs()` runs energy minimization on each generated conformer using the MMFF94 (Merck Molecular Force Field) to relax steric clashes and optimize bond angles/lengths.

    Then code pairs each conformer ID with its resulting potential energy (in $\text{kcal/mol}$) and sorts them so the lowest-energy conformer comes first.
    """)
    return


@app.cell
def _(cids: "Sequence[int]", ibuprofen_explicit_hydro: "Mol"):
    from rdkit.Chem.rdForceFieldHelpers import MMFFOptimizeMoleculeConfs

    ibuprofen_energy_results = MMFFOptimizeMoleculeConfs(ibuprofen_explicit_hydro, maxIters=2000)
    energies: list[tuple[int, float]] = [(cid, e) for cid, (_, e) in zip(cids, ibuprofen_energy_results)]
    energies.sort(key=lambda t: t[1])
    print(f'lowest 3 MMFF energies: {[(c, round(e, 2)) for c, e in energies[:3]]}')
    return (energies,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    `GetBestRMS()` calculates the Root Mean Square Deviation (in Ångströms) between the lowest-energy conformer (`refId`) and the second-lowest (`prbId`).
    """)
    return


@app.cell
def _(
    energies: list[tuple[int, float]],
    ibuprofen_explicit_hydro: "Mol",
    rdMolAlign,
):
    rms: float = rdMolAlign.GetBestRMS(
        prbMol=ibuprofen_explicit_hydro, 
        refMol=ibuprofen_explicit_hydro, 
        prbId=energies[1][0], 
        refId=energies[0][0]    # best
        )
    print(f'RMSD between two lowest conformers: {rms:.2f} A')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Exporting to SDF Format. Structure-Data File (SDF) format stores atomic 3D coordinates, bonding topology, and associated metadata
    """)
    return


@app.cell
def _(
    Chem,
    energies: list[tuple[int, float]],
    ibuprofen_explicit_hydro: "Mol",
):
    with Chem.SDWriter("../data/ibuprofen_best.sdf") as w:
            ibuprofen_explicit_hydro.SetProp("_Name", "ibuprofen")
            ibuprofen_explicit_hydro.SetProp("MMFF_energy", f"{energies[0][1]:.3f}")
            w.write(ibuprofen_explicit_hydro, confId=energies[0][0])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # 6. Pipeline

    Constructs an input list containing the valid SMILES from the DRUGS dictionary, plus two  edge cases:
    * a sodium salt version of aspirin (a duplicate in disguise)
    * an invalid SMILES string (`C1CC` has an unclosed ring).
    """)
    return


@app.cell
def _(DRUGS: dict[str, str]):
    raw_pipedata = list(DRUGS.values()) + [
            "CC(=O)Oc1ccccc1C(=O)[O-].[Na+]",    # duplicate of aspirin, as a salt
            "C1CC",                              # unparseable SMILES
        ]
    return (raw_pipedata,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Processing & Standardizing Each Compound

    Attempts to parse each SMILES string. If it fails (like "C1CC"), it increments the failures counter and skips to the next entry.

    Standardization Cleanup:
    1. `rdMolStandardize.Cleanup()` normalizes functional groups and disconnects metals.
    2. `lfc.choose()` strips counterions/salts.
    3. `unch_pipedata.uncharge` neutralizes charges where chemically sensible.
    """)
    return


@app.cell
def _(
    Chem,
    Descriptors,
    FingerprintGenerator64,
    Mol,
    MurckoScaffold,
    QED,
    RDLogger,
    df: "DataFrame",
    gen: "FingerprintGenerator64",
    lfc,
    pd,
    raw_pipedata,
    rdFingerprintGenerator,
    rdMolStandardize,
):
    lfc_pipedata = rdMolStandardize.LargestFragmentChooser()
    unch_pipedata = rdMolStandardize.Uncharger()
    gen_pipedata: FingerprintGenerator64 = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

    RDLogger.DisableLog("rdApp.error")
    records: list  = []
    failures: int = 0

    for smi_pipe in raw_pipedata:
        mol_pipe: Mol = Chem.MolFromSmiles(smi_pipe)
        if mol_pipe is None:
            failures += 1
            continue
        mol_pipe = unch_pipedata.uncharge(lfc.choose(rdMolStandardize.Cleanup(mol_pipe)))
        records.append({
            "canonical_smiles": Chem.MolToSmiles(mol_pipe),
            "inchikey": Chem.MolToInchiKey(mol_pipe),
            "scaffold": MurckoScaffold.MurckoScaffoldSmiles(mol=mol_pipe),
            "MW": Descriptors.MolWt(mol_pipe),
            "QED": QED.qed(mol_pipe),
            "_fp": gen.GetFingerprint(mol_pipe),
        })
    RDLogger.EnableLog("rdApp.error")

    df_pipe = pd.DataFrame(records)
    print(f'parsed {len(df)} / {len(raw_pipedata)}  ({failures} failed)')

    return (df_pipe,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Deduplication on InChIKey

    Uses the inchikey column to identify duplicate molecules. Because the salt form of aspirin was stripped down to free aspirin during standardization, its InChIKey now matches regular aspirin.
    """)
    return


@app.cell
def _(DataFrame, df_pipe):
    before: int = len(df_pipe)
    df_filtered: DataFrame = df_pipe.drop_duplicates(subset="inchikey").reset_index(drop=True)
    print(f"deduplicated on InChIKey: {before} -> {len(df_filtered)}")
    return (df_filtered,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Scaffold Grouping for Train/Test Splits
    """)
    return


@app.cell
def _(df_filtered: "DataFrame"):
    groups = df_filtered.groupby("scaffold").size().sort_values(ascending=False)
    print(f"distinct scaffolds: {len(groups)}")

    print("final table:")
    df_pipe_final = df_filtered.drop(columns="_fp").round(2)
    df_pipe_final

    return


if __name__ == "__main__":
    app.run()
