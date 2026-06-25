"""RDKit structural annotation for the comparative flavonoid pipeline.

This module reads the occurrence-level CSV exported by Part II, annotates each
unique InChIKey once, and writes both compound-level and occurrence-level
outputs. It preserves invalid structures for audit, detects conflicting SMILES
assigned to the same InChIKey, calculates RDKit descriptors, and generates
Bemis-Murcko scaffolds.

Default execution
-----------------
From the project root:

    python scripts/3_annotate_flavonoids_rdkit_auto.py

Explicit execution (recommended for exact reproducibility):

    python scripts/3_annotate_flavonoids_rdkit_auto.py \
        --input outputs/<run_tag>/PartII_ALL/<base_tag>__flavonoids_for_rdkit.csv \
        --output-dir outputs/<run_tag>/RDKit_ALL

Main outputs
------------
- lotus_flavonoids_rdkit_compounds.csv
    One row per InChIKey. This is the authoritative structural table.
- lotus_flavonoids_rdkit_annotations.csv
    Original occurrence rows joined to the compound-level annotations.
- lotus_flavonoids_rdkit_invalid_structures.csv
    Missing or unparsable SMILES and compounds without a valid representative.
- lotus_flavonoids_rdkit_structure_conflicts.csv
    InChIKeys represented by more than one distinct standardized structure.
- lotus_flavonoids_rdkit_structural_qc.csv
    Structures requiring review, with transparent QC flags and taxonomy context.
- lotus_flavonoids_rdkit_run_metadata.json
    Parameters, software versions, file hashes, and row counts.

Scientific scope
----------------
The sugar detector identifies probable oxygen-rich sugar-like rings. It does
not prove glycosidic attachment and should not be described as definitive
"glycosylation" without additional validation. The prenyl detector identifies
prenyl-like substructures and is likewise intended as an interpretable motif
flag rather than a complete biosynthetic annotation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import shutil
import sys
from collections import Counter
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Optional

import pandas as pd
from rdkit import Chem, rdBase
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold


DEFAULT_INPUT_NAME = "lotus_flavonoids_for_rdkit.csv"
DEFAULT_OCCURRENCE_OUTPUT = "lotus_flavonoids_rdkit_annotations.csv"
DEFAULT_COMPOUND_OUTPUT = "lotus_flavonoids_rdkit_compounds.csv"
DEFAULT_INVALID_OUTPUT = "lotus_flavonoids_rdkit_invalid_structures.csv"
DEFAULT_CONFLICT_OUTPUT = "lotus_flavonoids_rdkit_structure_conflicts.csv"
DEFAULT_METADATA_OUTPUT = "lotus_flavonoids_rdkit_run_metadata.json"
DEFAULT_QC_OUTPUT = "lotus_flavonoids_rdkit_structural_qc.csv"

# Independent motif definitions. Counts and presence/absence are both exported.
SMARTS_PATTERNS: dict[str, str] = {
    "phenolic_OH": "[OX2H][cX3]",
    "methoxy_aryl": "[cX3][OX2][CH3]",
    # Broad isoprenyl fragment. Reported explicitly as prenyl-like, not as a
    # definitive biosynthetic prenylation assignment.
    "prenyl_like": "C=C(C)C",
    # More precise name than the previous generic 'conjugated carbonyl'.
    "alpha_beta_unsaturated_carbonyl": "[CX3](=O)[CX3]=[CX3]",
}

DECORATION_FLAGS = (
    "has_methoxy_aryl",
    "has_prenyl_like",
    "has_probable_sugar_ring",
)


@dataclass(frozen=True)
class ParsedStructure:
    input_smiles: str
    input_canonical_smiles: Optional[str]
    standardized_smiles: Optional[str]
    mol: Optional[Chem.Mol]
    parse_status: str
    error_message: Optional[str]
    rdkit_inchikey: Optional[str]
    input_fragment_count: Optional[int]
    stereo_fallback_used: bool


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def clean_text(value: Any) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def normalize_inchikey(value: Any) -> str:
    return clean_text(value).upper()


def copy_without_stereochemistry(mol: Chem.Mol) -> Chem.Mol:
    """Return a molecule copy with all stereochemical annotations removed."""
    copied = Chem.Mol(mol)
    Chem.RemoveStereochemistry(copied)
    return copied


def has_unsafe_double_bond_stereo(mol: Chem.Mol) -> bool:
    """Detect ambiguous or incomplete double-bond stereo annotations."""
    stereo_none = Chem.rdchem.BondStereo.STEREONONE
    stereo_any = Chem.rdchem.BondStereo.STEREOANY
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.rdchem.BondType.DOUBLE:
            continue
        stereo = bond.GetStereo()
        if stereo == stereo_any:
            return True
        if stereo != stereo_none and len(bond.GetStereoAtoms()) != 2:
            return True
    return False


def safe_mol_to_smiles(
    mol: Chem.Mol, *, isomeric: bool
) -> tuple[str, bool]:
    """Generate canonical SMILES, removing malformed stereo when needed.

    Some public-database structures contain incomplete double-bond stereo
    annotations. Recent RDKit versions correctly reject those annotations
    during canonicalization. Since the present analyses do not use
    stereochemistry, the deterministic fallback removes stereo rather than
    discarding the compound.
    """
    stereo_removed = False
    working = mol

    if isomeric and has_unsafe_double_bond_stereo(working):
        working = copy_without_stereochemistry(working)
        stereo_removed = True
        isomeric = False

    try:
        return (
            Chem.MolToSmiles(
                working, canonical=True, isomericSmiles=isomeric
            ),
            stereo_removed,
        )
    except Exception:
        fallback = copy_without_stereochemistry(mol)
        return (
            Chem.MolToSmiles(
                fallback, canonical=True, isomericSmiles=False
            ),
            True,
        )


def safe_inchikey(mol: Chem.Mol) -> Optional[str]:
    try:
        value = Chem.MolToInchiKey(mol)
        return value if value else None
    except Exception:
        try:
            value = Chem.MolToInchiKey(copy_without_stereochemistry(mol))
            return value if value else None
        except Exception:
            return None


def standardize_molecule(mol: Chem.Mol) -> tuple[Chem.Mol, bool]:
    """Clean a molecule and retain its largest fragment parent.

    Tautomer canonicalization and charge neutralization are intentionally not
    applied because they can modify motif assignments. The fragment-parent
    step removes disconnected counterions while preserving the principal
    natural-product structure.
    """
    stereo_removed = False
    working = mol
    if has_unsafe_double_bond_stereo(working):
        working = copy_without_stereochemistry(working)
        stereo_removed = True

    try:
        cleaned = rdMolStandardize.Cleanup(working)
        parent = rdMolStandardize.FragmentParent(cleaned)
        Chem.SanitizeMol(parent)
        return parent, stereo_removed
    except Exception:
        fallback = copy_without_stereochemistry(mol)
        cleaned = rdMolStandardize.Cleanup(fallback)
        parent = rdMolStandardize.FragmentParent(cleaned)
        Chem.SanitizeMol(parent)
        return parent, True


def parse_structure(smiles: Any) -> ParsedStructure:
    text = clean_text(smiles)
    if not text:
        return ParsedStructure(
            input_smiles="",
            input_canonical_smiles=None,
            standardized_smiles=None,
            mol=None,
            parse_status="missing_smiles",
            error_message="SMILES is missing or empty.",
            rdkit_inchikey=None,
            input_fragment_count=None,
            stereo_fallback_used=False,
        )

    try:
        raw_mol = Chem.MolFromSmiles(text)
        if raw_mol is None:
            raise ValueError("RDKit MolFromSmiles returned None.")

        input_canonical, input_stereo_fallback = safe_mol_to_smiles(
            raw_mol, isomeric=True
        )
        parent, standardization_stereo_fallback = standardize_molecule(raw_mol)
        standardized, output_stereo_fallback = safe_mol_to_smiles(
            parent, isomeric=True
        )
        stereo_fallback_used = bool(
            input_stereo_fallback
            or standardization_stereo_fallback
            or output_stereo_fallback
        )
        return ParsedStructure(
            input_smiles=text,
            input_canonical_smiles=input_canonical,
            standardized_smiles=standardized,
            mol=parent,
            parse_status="valid",
            error_message=None,
            rdkit_inchikey=safe_inchikey(raw_mol),
            input_fragment_count=len(Chem.GetMolFrags(raw_mol)),
            stereo_fallback_used=stereo_fallback_used,
        )
    except Exception as exc:
        return ParsedStructure(
            input_smiles=text,
            input_canonical_smiles=None,
            standardized_smiles=None,
            mol=None,
            parse_status="invalid_smiles",
            error_message=str(exc),
            rdkit_inchikey=None,
            input_fragment_count=None,
            stereo_fallback_used=False,
        )


def compile_smarts_patterns(patterns: dict[str, str]) -> dict[str, Chem.Mol]:
    compiled: dict[str, Chem.Mol] = {}
    for name, smarts in patterns.items():
        query = Chem.MolFromSmarts(smarts)
        if query is None:
            raise ValueError(f"Invalid SMARTS pattern for '{name}': {smarts}")
        compiled[name] = query
    return compiled


def count_unique_substructure_matches(mol: Chem.Mol, query: Chem.Mol) -> int:
    """Count unique atom sets to reduce duplicate symmetric SMARTS matches."""
    matches = mol.GetSubstructMatches(query, uniquify=True)
    unique_sets = {tuple(sorted(match)) for match in matches}
    return len(unique_sets)


def probable_sugar_rings(
    mol: Chem.Mol,
    min_ring_oxygens: int = 3,
) -> list[tuple[int, ...]]:
    """Return oxygen-rich, non-aromatic 5/6-membered sugar-like rings.

    This intentionally conservative heuristic detects probable sugar-like
    rings. It does not distinguish O- from C-glycosides and does not prove a
    glycosidic bond to the flavonoid core.
    """
    accepted: list[tuple[int, ...]] = []

    for ring in mol.GetRingInfo().AtomRings():
        ring_tuple = tuple(sorted(ring))
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]

        if len(ring_atoms) not in (5, 6):
            continue
        if any(atom.GetSymbol() not in {"C", "O"} for atom in ring_atoms):
            continue
        if any(atom.GetIsAromatic() for atom in ring_atoms):
            continue

        sp3_count = sum(
            atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3
            for atom in ring_atoms
        )
        if sp3_count < len(ring_atoms) - 1:
            continue

        oxygen_ids: set[int] = set()
        ring_ids = set(ring)

        for atom in ring_atoms:
            if atom.GetSymbol() == "O":
                oxygen_ids.add(atom.GetIdx())

            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() != "O":
                    continue
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                if bond is None or bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                    continue

                # Exclude carbonyl oxygen atoms.
                is_carbonyl_oxygen = any(
                    other_bond.GetBondType() == Chem.rdchem.BondType.DOUBLE
                    and other_bond.GetOtherAtom(neighbor).GetSymbol() == "C"
                    for other_bond in neighbor.GetBonds()
                )
                if not is_carbonyl_oxygen:
                    oxygen_ids.add(neighbor.GetIdx())

        exocyclic_oxygen_count = sum(idx not in ring_ids for idx in oxygen_ids)
        if len(oxygen_ids) >= min_ring_oxygens and exocyclic_oxygen_count >= 1:
            accepted.append(ring_tuple)

    return sorted(set(accepted))


def bemis_murcko_scaffold(mol: Chem.Mol) -> tuple[str, str, int]:
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold is None or scaffold.GetNumAtoms() == 0:
        return "", "", 0

    scaffold_smiles, _ = safe_mol_to_smiles(scaffold, isomeric=False)
    generic = MurckoScaffold.MakeScaffoldGeneric(scaffold)
    generic_smiles, _ = safe_mol_to_smiles(generic, isomeric=False)
    return scaffold_smiles, generic_smiles, scaffold.GetNumAtoms()


def motif_combination(record: dict[str, Any], flags: Iterable[str]) -> str:
    present = [flag.removeprefix("has_") for flag in flags if bool(record.get(flag))]
    return "+".join(present) if present else "none_of_selected_motifs"


def annotate_molecule(
    mol: Chem.Mol,
    compiled_patterns: dict[str, Chem.Mol],
) -> dict[str, Any]:
    record: dict[str, Any] = {}

    for name, query in compiled_patterns.items():
        count = count_unique_substructure_matches(mol, query)
        record[f"n_{name}"] = count
        record[f"has_{name}"] = count > 0

    sugar_rings = probable_sugar_rings(mol)
    record["n_probable_sugar_rings"] = len(sugar_rings)
    record["has_probable_sugar_ring"] = len(sugar_rings) > 0
    # Backward-compatible alias used by the previous Part III.
    record["has_probable_sugar"] = record["has_probable_sugar_ring"]

    scaffold_smiles, generic_scaffold_smiles, scaffold_atoms = (
        bemis_murcko_scaffold(mol)
    )
    record["bemis_murcko_scaffold_smiles"] = scaffold_smiles
    record["bemis_murcko_generic_scaffold_smiles"] = generic_scaffold_smiles
    record["bemis_murcko_scaffold_atom_count"] = scaffold_atoms
    record["has_bemis_murcko_scaffold"] = scaffold_atoms > 0
    # Compatibility alias for LOTUS/previous R code. This column is generated
    # by RDKit here and should replace InChIKey core14 in scaffold analyses.
    record["murcko_framework"] = scaffold_smiles
    record["murko_framework"] = scaffold_smiles

    record["mw"] = float(Descriptors.MolWt(mol))
    record["exact_mw"] = float(Descriptors.ExactMolWt(mol))
    record["logp"] = float(Crippen.MolLogP(mol))
    record["tpsa"] = float(rdMolDescriptors.CalcTPSA(mol))
    record["fsp3"] = float(rdMolDescriptors.CalcFractionCSP3(mol))
    record["num_HBD"] = int(rdMolDescriptors.CalcNumHBD(mol))
    record["num_HBA"] = int(rdMolDescriptors.CalcNumHBA(mol))
    record["num_rings"] = int(rdMolDescriptors.CalcNumRings(mol))
    record["num_aromatic_rings"] = int(
        rdMolDescriptors.CalcNumAromaticRings(mol)
    )
    record["num_rotatable_bonds"] = int(
        rdMolDescriptors.CalcNumRotatableBonds(mol)
    )
    record["heavy_atom_number"] = int(mol.GetNumHeavyAtoms())
    record["formal_charge"] = int(Chem.GetFormalCharge(mol))
    record["fragment_count_after_standardization"] = len(
        Chem.GetMolFrags(mol)
    )

    record["selected_motif_count"] = sum(
        bool(record.get(flag)) for flag in DECORATION_FLAGS
    )
    record["is_simple_selected_motifs"] = record["selected_motif_count"] == 0
    record["decoration_combination"] = motif_combination(
        record, DECORATION_FLAGS
    )
    all_flags = (
        "has_phenolic_OH",
        "has_methoxy_aryl",
        "has_prenyl_like",
        "has_alpha_beta_unsaturated_carbonyl",
        "has_probable_sugar_ring",
    )
    record["all_motif_combination"] = motif_combination(record, all_flags)
    return record



def collapse_unique_text(values: pd.Series) -> str:
    cleaned = sorted(
        {
            clean_text(value)
            for value in values
            if clean_text(value)
        }
    )
    return ";".join(cleaned)


def structure_qc_record(row: pd.Series) -> dict[str, Any]:
    """Calculate conservative structural QC flags.

    These flags do not replace NPClassifier. They identify structures that
    require exclusion from the primary structural analysis or manual review.
    The minimum carbon threshold is 14 because 2-arylbenzofurans, one of the
    selected subclasses, may contain a C14 core.
    """
    status = clean_text(row.get("annotation_status"))
    standardized_smiles = clean_text(row.get("standardized_smiles"))
    input_smiles = clean_text(row.get("representative_input_smiles"))

    result: dict[str, Any] = {
        "standardized_rdkit_inchikey": "",
        "standardized_inchikey_exact_match": False,
        "standardized_inchikey_connectivity_match": False,
        "num_carbons": pd.NA,
        "num_oxygens": pd.NA,
        "num_nitrogens": pd.NA,
        "num_sulfurs": pd.NA,
        "num_halogens": pd.NA,
        "largest_input_fragment_heavy_atoms": pd.NA,
        "second_input_fragment_heavy_atoms": pd.NA,
        "substantial_multicomponent_input": False,
        "qc_zero_oxygen": False,
        "qc_carbon_count_below_14": False,
        "qc_contains_nitrogen": False,
        "qc_contains_sulfur": False,
        "qc_contains_halogen": False,
        "primary_structure_eligible": False,
        "structural_qc_status": "exclude_primary",
        "structural_qc_reasons": "annotation_not_valid",
    }

    if not status.startswith("valid") or not standardized_smiles:
        return result

    mol = Chem.MolFromSmiles(standardized_smiles)
    if mol is None:
        result["structural_qc_reasons"] = "standardized_smiles_unparsable"
        return result

    element_counts = Counter(atom.GetSymbol() for atom in mol.GetAtoms())
    result["num_carbons"] = int(element_counts.get("C", 0))
    result["num_oxygens"] = int(element_counts.get("O", 0))
    result["num_nitrogens"] = int(element_counts.get("N", 0))
    result["num_sulfurs"] = int(element_counts.get("S", 0))
    result["num_halogens"] = int(
        sum(element_counts.get(symbol, 0) for symbol in ("F", "Cl", "Br", "I"))
    )

    standardized_key = safe_inchikey(mol) or ""
    input_key = normalize_inchikey(row.get("inchikey"))
    result["standardized_rdkit_inchikey"] = standardized_key
    result["standardized_inchikey_exact_match"] = bool(
        standardized_key and standardized_key == input_key
    )
    result["standardized_inchikey_connectivity_match"] = bool(
        standardized_key and input_key and standardized_key[:14] == input_key[:14]
    )

    if input_smiles:
        input_mol = Chem.MolFromSmiles(input_smiles)
        if input_mol is not None:
            fragment_sizes = sorted(
                (
                    int(fragment.GetNumHeavyAtoms())
                    for fragment in Chem.GetMolFrags(
                        input_mol, asMols=True, sanitizeFrags=True
                    )
                ),
                reverse=True,
            )
            if fragment_sizes:
                result["largest_input_fragment_heavy_atoms"] = fragment_sizes[0]
            if len(fragment_sizes) >= 2:
                result["second_input_fragment_heavy_atoms"] = fragment_sizes[1]
                result["substantial_multicomponent_input"] = bool(
                    fragment_sizes[1] >= 10
                    and fragment_sizes[1] / max(fragment_sizes[0], 1) >= 0.25
                )

    result["qc_zero_oxygen"] = result["num_oxygens"] == 0
    result["qc_carbon_count_below_14"] = result["num_carbons"] < 14
    result["qc_contains_nitrogen"] = result["num_nitrogens"] > 0
    result["qc_contains_sulfur"] = result["num_sulfurs"] > 0
    result["qc_contains_halogen"] = result["num_halogens"] > 0

    exclusion_reasons: list[str] = []
    review_reasons: list[str] = []

    if result["qc_zero_oxygen"]:
        exclusion_reasons.append("zero_oxygen")
    if result["qc_carbon_count_below_14"]:
        exclusion_reasons.append("carbon_count_below_14")
    if result["substantial_multicomponent_input"]:
        exclusion_reasons.append("substantial_multicomponent_input")

    if bool(row.get("stereo_fallback_used", False)):
        review_reasons.append("stereo_removed_for_canonicalization")
    if result["qc_contains_nitrogen"]:
        review_reasons.append("contains_nitrogen")
    if result["qc_contains_sulfur"]:
        review_reasons.append("contains_sulfur")
    if result["qc_contains_halogen"]:
        review_reasons.append("contains_halogen")
    if not result["standardized_inchikey_connectivity_match"]:
        review_reasons.append("standardized_inchikey_connectivity_changed")

    result["primary_structure_eligible"] = not exclusion_reasons
    if exclusion_reasons:
        result["structural_qc_status"] = "exclude_primary"
        result["structural_qc_reasons"] = ";".join(
            exclusion_reasons + review_reasons
        )
    elif review_reasons:
        result["structural_qc_status"] = "review"
        result["structural_qc_reasons"] = ";".join(review_reasons)
    else:
        result["structural_qc_status"] = "pass"
        result["structural_qc_reasons"] = ""

    return result


def add_structural_qc(
    compounds: pd.DataFrame,
    input_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    qc_rows = [
        structure_qc_record(row)
        for _, row in compounds.iterrows()
    ]
    qc_df = pd.DataFrame(qc_rows)
    compounds_qc = pd.concat(
        [compounds.reset_index(drop=True), qc_df.reset_index(drop=True)],
        axis=1,
    )

    context_columns = [
        column
        for column in ("class_np", "family", "genus", "species", "ref_id")
        if column in input_df.columns
    ]
    if context_columns:
        aggregation = {
            column: collapse_unique_text
            for column in context_columns
        }
        context = (
            input_df.groupby("inchikey", dropna=False)
            .agg(aggregation)
            .reset_index()
        )
        compounds_qc = compounds_qc.merge(
            context,
            on="inchikey",
            how="left",
            validate="one_to_one",
        )

    review = compounds_qc.loc[
        compounds_qc["structural_qc_status"] != "pass"
    ].copy()
    preferred_columns = [
        "inchikey",
        "annotation_status",
        "structural_qc_status",
        "structural_qc_reasons",
        "primary_structure_eligible",
        "class_np",
        "family",
        "genus",
        "species",
        "representative_input_smiles",
        "standardized_smiles",
        "num_carbons",
        "num_oxygens",
        "num_nitrogens",
        "num_sulfurs",
        "num_halogens",
        "input_fragment_count",
        "largest_input_fragment_heavy_atoms",
        "second_input_fragment_heavy_atoms",
        "substantial_multicomponent_input",
        "stereo_fallback_used",
        "standardized_rdkit_inchikey",
        "standardized_inchikey_exact_match",
        "standardized_inchikey_connectivity_match",
        "mw",
        "logp",
        "tpsa",
        "fsp3",
        "bemis_murcko_scaffold_smiles",
        "ref_id",
    ]
    review = review[
        [column for column in preferred_columns if column in review.columns]
    ].sort_values(
        by=["structural_qc_status", "structural_qc_reasons", "inchikey"],
        kind="stable",
    )
    return compounds_qc, review


def input_candidates(project_root: Path) -> list[Path]:
    candidates: list[Path] = []

    for explicit in (
        project_root / DEFAULT_INPUT_NAME,
        Path(__file__).resolve().parent / DEFAULT_INPUT_NAME,
    ):
        if explicit.exists() and explicit.is_file():
            candidates.append(explicit)

    patterns = (
        "outputs/**/PartII_ALL/*__flavonoids_for_rdkit.csv",
        "**/PartII_ALL/*__flavonoids_for_rdkit.csv",
        "**/*__flavonoids_for_rdkit.csv",
    )
    for pattern in patterns:
        candidates.extend(
            path for path in project_root.glob(pattern) if path.is_file()
        )

    unique: dict[Path, Path] = {}
    for path in candidates:
        unique[path.resolve(strict=False)] = path

    return sorted(
        unique.values(),
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )


def resolve_input_csv(explicit_input: Optional[str], project_root: Path) -> Path:
    if explicit_input:
        path = Path(explicit_input).expanduser().resolve(strict=False)
        if not path.exists() or not path.is_file():
            raise FileNotFoundError(f"Input CSV was not found: {path}")
        return path

    candidates = input_candidates(project_root)
    if not candidates:
        raise FileNotFoundError(
            "No RDKit input CSV was found. Run Part II first or provide "
            "--input <path-to-__flavonoids_for_rdkit.csv>."
        )

    selected = candidates[0]
    print(f"[AUTO-DETECT] Selected newest input: {selected}")
    if len(candidates) > 1:
        print(
            f"[AUTO-DETECT] {len(candidates) - 1} additional candidate(s) "
            "were found. Use --input for exact reproducibility."
        )
    return selected


def default_output_dir(input_csv: Path) -> Path:
    if input_csv.parent.name == "PartII_ALL":
        return input_csv.parent.parent / "RDKit_ALL"
    return input_csv.parent / "RDKit_ALL"


def representative_structure(
    group: pd.DataFrame,
    input_inchikey: str,
) -> tuple[pd.Series, bool, list[str]]:
    valid = group.loc[group["parse_status"] == "valid"].copy()
    if valid.empty:
        return group.iloc[0], False, []

    valid["exact_match"] = (
        valid["rdkit_inchikey"].fillna("").str.upper() == input_inchikey
    )
    valid["connectivity_match"] = (
        valid["rdkit_inchikey"].fillna("").str[:14]
        == input_inchikey[:14]
    )

    frequencies = Counter(valid["standardized_smiles"].dropna())
    valid["structure_frequency"] = valid["standardized_smiles"].map(frequencies)
    valid = valid.sort_values(
        by=[
            "exact_match",
            "connectivity_match",
            "structure_frequency",
            "standardized_smiles",
        ],
        ascending=[False, False, False, True],
        kind="stable",
    )

    structures = sorted(valid["standardized_smiles"].dropna().unique().tolist())
    return valid.iloc[0], len(structures) > 1, structures


def build_compound_annotations(
    input_df: pd.DataFrame,
    compiled_patterns: dict[str, Chem.Mol],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    parsed_rows: list[dict[str, Any]] = []

    unique_pairs = (
        input_df[["inchikey", "smiles"]]
        .drop_duplicates()
        .reset_index(drop=True)
    )

    total_pairs = len(unique_pairs)
    for pair_index, (_, row) in enumerate(unique_pairs.iterrows(), start=1):
        parsed = parse_structure(row["smiles"])
        parsed_rows.append(
            {
                "inchikey": normalize_inchikey(row["inchikey"]),
                **asdict(parsed),
            }
        )
        if pair_index % 500 == 0 or pair_index == total_pairs:
            print(
                f"[RDKit] Parsed {pair_index:,}/{total_pairs:,} unique "
                "InChIKey-SMILES pairs."
            )

    parsed_df = pd.DataFrame(parsed_rows)
    compound_records: list[dict[str, Any]] = []
    conflict_records: list[dict[str, Any]] = []

    input_counts = input_df.groupby("inchikey", dropna=False).size().to_dict()

    for inchikey, group in parsed_df.groupby("inchikey", sort=False, dropna=False):
        key = normalize_inchikey(inchikey)
        representative, conflict, structures = representative_structure(group, key)
        valid_group = group.loc[group["parse_status"] == "valid"]

        base_record: dict[str, Any] = {
            "inchikey": key,
            "annotation_status": (
                "valid_conflicting_structures"
                if conflict
                else "valid"
                if not valid_group.empty
                else "no_valid_structure"
            ),
            "n_input_rows": int(input_counts.get(key, 0)),
            "n_input_smiles": int(group["input_smiles"].nunique(dropna=True)),
            "n_valid_smiles": int((group["parse_status"] == "valid").sum()),
            "n_unique_standardized_smiles": int(
                valid_group["standardized_smiles"].nunique(dropna=True)
            ),
            "structure_conflict": bool(conflict),
            "representative_input_smiles": representative["input_smiles"],
            "canonical_smiles": representative.get("input_canonical_smiles"),
            "standardized_smiles": representative.get("standardized_smiles"),
            "rdkit_inchikey": representative.get("rdkit_inchikey"),
            "input_fragment_count": representative.get("input_fragment_count"),
            "stereo_fallback_used": bool(
                representative.get("stereo_fallback_used", False)
            ),
            "annotation_error_message": None,
        }
        rdkit_key = clean_text(base_record["rdkit_inchikey"]).upper()
        base_record["inchikey_exact_match"] = bool(rdkit_key and rdkit_key == key)
        base_record["inchikey_connectivity_match"] = bool(
            rdkit_key and key and rdkit_key[:14] == key[:14]
        )

        mol = representative.get("mol")
        if isinstance(mol, Chem.Mol):
            try:
                base_record.update(annotate_molecule(mol, compiled_patterns))
            except Exception as exc:
                base_record["annotation_status"] = "annotation_error"
                base_record["annotation_error_message"] = str(exc)

        compound_records.append(base_record)

        if conflict:
            conflict_records.append(
                {
                    "inchikey": key,
                    "n_unique_standardized_smiles": len(structures),
                    "selected_standardized_smiles": representative[
                        "standardized_smiles"
                    ],
                    "all_standardized_smiles": " | ".join(structures),
                    "selection_rule": (
                        "exact InChIKey match, then connectivity match, then "
                        "most frequent standardized SMILES"
                    ),
                }
            )

    compounds = pd.DataFrame(compound_records)
    conflict_columns = [
        "inchikey",
        "n_unique_standardized_smiles",
        "selected_standardized_smiles",
        "all_standardized_smiles",
        "selection_rule",
    ]
    conflicts = pd.DataFrame(conflict_records, columns=conflict_columns)

    invalid_pairs = parsed_df.loc[parsed_df["parse_status"] != "valid"].copy()
    no_valid_keys = set(
        compounds.loc[
            compounds["annotation_status"] == "no_valid_structure", "inchikey"
        ]
    )
    invalid_pairs["compound_has_no_valid_structure"] = invalid_pairs[
        "inchikey"
    ].isin(no_valid_keys)

    annotation_errors = compounds.loc[
        compounds["annotation_status"] == "annotation_error",
        ["inchikey", "representative_input_smiles", "annotation_error_message"],
    ].copy()
    if not annotation_errors.empty:
        annotation_errors = annotation_errors.rename(
            columns={
                "representative_input_smiles": "input_smiles",
                "annotation_error_message": "error_message",
            }
        )
        annotation_errors["parse_status"] = "annotation_error"
        annotation_errors["compound_has_no_valid_structure"] = False
        invalid_pairs = pd.concat(
            [invalid_pairs, annotation_errors], ignore_index=True, sort=False
        )

    invalid_columns = [
        "inchikey",
        "input_smiles",
        "parse_status",
        "error_message",
        "input_canonical_smiles",
        "standardized_smiles",
        "rdkit_inchikey",
        "input_fragment_count",
        "stereo_fallback_used",
        "compound_has_no_valid_structure",
    ]
    for column in invalid_columns:
        if column not in invalid_pairs.columns:
            invalid_pairs[column] = pd.NA
    invalid_pairs = invalid_pairs[invalid_columns].drop_duplicates()

    return compounds, invalid_pairs, conflicts


def compatibility_targets(input_csv: Path, project_root: Path) -> list[Path]:
    targets = [project_root / DEFAULT_OCCURRENCE_OUTPUT]
    targets.append(input_csv.with_name(DEFAULT_OCCURRENCE_OUTPUT))
    if input_csv.parent.name == "PartII_ALL":
        targets.append(input_csv.parent.parent / DEFAULT_OCCURRENCE_OUTPUT)

    unique: list[Path] = []
    seen: set[Path] = set()
    for target in targets:
        resolved = target.resolve(strict=False)
        if resolved not in seen:
            seen.add(resolved)
            unique.append(target)
    return unique


def write_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def make_metadata(
    input_csv: Path,
    output_dir: Path,
    input_df: pd.DataFrame,
    compounds: pd.DataFrame,
    occurrences: pd.DataFrame,
    invalid: pd.DataFrame,
    conflicts: pd.DataFrame,
    structural_qc: pd.DataFrame,
) -> dict[str, Any]:
    valid_compounds = int(
        compounds["annotation_status"].astype(str).str.startswith("valid").sum()
    )
    return {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_file": str(input_csv.resolve()),
        "input_sha256": sha256_file(input_csv),
        "output_directory": str(output_dir.resolve()),
        "software": {
            "python": platform.python_version(),
            "pandas": pd.__version__,
            "rdkit": rdBase.rdkitVersion,
        },
        "standardization": {
            "cleanup": True,
            "fragment_parent": True,
            "tautomer_canonicalization": False,
            "charge_neutralization": False,
        },
        "motif_smarts": SMARTS_PATTERNS,
        "structural_qc": {
            "primary_minimum_carbons": 14,
            "primary_requires_oxygen": True,
            "primary_excludes_substantial_multicomponent_inputs": True,
            "nitrogen_sulfur_and_halogens": "flagged for review, not automatically excluded",
            "interpretation": "QC flags support manual curation and sensitivity analyses; they do not redefine NPClassifier classes",
        },
        "probable_sugar_rule": {
            "ring_sizes": [5, 6],
            "allowed_ring_atoms": ["C", "O"],
            "non_aromatic": True,
            "minimum_total_noncarbonyl_oxygens": 3,
            "minimum_exocyclic_oxygens": 1,
            "interpretation": "probable sugar-like ring, not definitive glycosylation",
        },
        "counts": {
            "input_rows": int(len(input_df)),
            "input_unique_inchikeys": int(input_df["inchikey"].nunique()),
            "compound_rows": int(len(compounds)),
            "valid_compounds": valid_compounds,
            "compounds_without_valid_structure": int(len(compounds) - valid_compounds),
            "invalid_structure_records": int(len(invalid)),
            "conflicting_inchikeys": int(len(conflicts)),
            "inchikey_exact_matches": int(compounds["inchikey_exact_match"].fillna(False).sum()),
            "inchikey_connectivity_matches": int(compounds["inchikey_connectivity_match"].fillna(False).sum()),
            "stereo_fallback_compounds": int(
                compounds["stereo_fallback_used"].fillna(False).sum()
            ),
            "primary_structure_eligible": int(
                compounds["primary_structure_eligible"].fillna(False).sum()
            ),
            "structural_qc_exclude_primary": int(
                (compounds["structural_qc_status"] == "exclude_primary").sum()
            ),
            "structural_qc_review": int(
                (compounds["structural_qc_status"] == "review").sum()
            ),
            "zero_oxygen_compounds": int(
                compounds["qc_zero_oxygen"].fillna(False).sum()
            ),
            "carbon_below_14_compounds": int(
                compounds["qc_carbon_count_below_14"].fillna(False).sum()
            ),
            "nitrogen_containing_compounds": int(
                compounds["qc_contains_nitrogen"].fillna(False).sum()
            ),
            "sulfur_containing_compounds": int(
                compounds["qc_contains_sulfur"].fillna(False).sum()
            ),
            "halogen_containing_compounds": int(
                compounds["qc_contains_halogen"].fillna(False).sum()
            ),
            "substantial_multicomponent_inputs": int(
                compounds["substantial_multicomponent_input"].fillna(False).sum()
            ),
            "structural_qc_rows": int(len(structural_qc)),
            "occurrence_output_rows": int(len(occurrences)),
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Annotate flavonoid structures with RDKit."
    )
    parser.add_argument(
        "--input",
        help="Path to the Part II *__flavonoids_for_rdkit.csv file.",
    )
    parser.add_argument(
        "--output-dir",
        help="Directory for the complete RDKit output set.",
    )
    parser.add_argument(
        "--project-root",
        default=str(Path.cwd()),
        help="Project root used for automatic input detection and compatibility output.",
    )
    parser.add_argument(
        "--no-compatibility-copies",
        action="store_true",
        help="Do not copy the occurrence annotation file to legacy lookup locations.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    project_root = Path(args.project_root).expanduser().resolve(strict=False)
    input_csv = resolve_input_csv(args.input, project_root)
    output_dir = (
        Path(args.output_dir).expanduser().resolve(strict=False)
        if args.output_dir
        else default_output_dir(input_csv)
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    print("------------------------------------------------------------")
    print("RDKit flavonoid structural annotation")
    print(f"Input: {input_csv}")
    print(f"Output directory: {output_dir}")
    print(f"RDKit version: {rdBase.rdkitVersion}")
    print("------------------------------------------------------------")

    input_df = pd.read_csv(input_csv, dtype=str, keep_default_na=True)
    required = {"inchikey", "smiles"}
    missing = sorted(required.difference(input_df.columns))
    if missing:
        raise KeyError(
            "Input CSV is missing required column(s): " + ", ".join(missing)
        )

    input_df = input_df.copy()
    input_df["inchikey"] = input_df["inchikey"].map(normalize_inchikey)
    input_df["smiles"] = input_df["smiles"].map(clean_text)
    input_df["source_row_id"] = range(1, len(input_df) + 1)

    empty_keys = input_df["inchikey"].eq("")
    if empty_keys.any():
        # Preserve missing-key rows without collapsing unrelated records.
        input_df.loc[empty_keys, "inchikey"] = input_df.loc[
            empty_keys, "source_row_id"
        ].map(lambda value: f"MISSING_INCHIKEY_ROW_{value}")

    compiled_patterns = compile_smarts_patterns(SMARTS_PATTERNS)
    compounds, invalid, conflicts = build_compound_annotations(
        input_df, compiled_patterns
    )
    compounds, structural_qc = add_structural_qc(compounds, input_df)

    occurrences = input_df.merge(
        compounds,
        how="left",
        on="inchikey",
        validate="many_to_one",
    )

    compound_path = output_dir / DEFAULT_COMPOUND_OUTPUT
    occurrence_path = output_dir / DEFAULT_OCCURRENCE_OUTPUT
    invalid_path = output_dir / DEFAULT_INVALID_OUTPUT
    conflict_path = output_dir / DEFAULT_CONFLICT_OUTPUT
    qc_path = output_dir / DEFAULT_QC_OUTPUT
    metadata_path = output_dir / DEFAULT_METADATA_OUTPUT

    write_csv(compounds, compound_path)
    write_csv(occurrences, occurrence_path)
    write_csv(invalid, invalid_path)
    write_csv(conflicts, conflict_path)
    write_csv(structural_qc, qc_path)

    metadata = make_metadata(
        input_csv=input_csv,
        output_dir=output_dir,
        input_df=input_df,
        compounds=compounds,
        occurrences=occurrences,
        invalid=invalid,
        conflicts=conflicts,
        structural_qc=structural_qc,
    )
    metadata_path.write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    if not args.no_compatibility_copies:
        for target in compatibility_targets(input_csv, project_root):
            target.parent.mkdir(parents=True, exist_ok=True)
            if target.resolve(strict=False) != occurrence_path.resolve(strict=False):
                shutil.copy2(occurrence_path, target)
                print(f"[COMPATIBILITY COPY] {target}")

    counts = metadata["counts"]
    print("------------------------------------------------------------")
    print("[RDKit SUMMARY]")
    print(f"Input occurrence rows: {counts['input_rows']}")
    print(f"Unique InChIKeys: {counts['input_unique_inchikeys']}")
    print(f"Valid compound annotations: {counts['valid_compounds']}")
    print(
        "Compounds without valid structures: "
        f"{counts['compounds_without_valid_structure']}"
    )
    print(f"InChIKeys with structural conflicts: {counts['conflicting_inchikeys']}")
    print(
        "Primary structure-eligible compounds: "
        f"{counts['primary_structure_eligible']}"
    )
    print(
        "Structural QC review/exclusion rows: "
        f"{counts['structural_qc_rows']}"
    )
    print(f"Occurrence rows exported: {counts['occurrence_output_rows']}")
    print(f"Compound table: {compound_path}")
    print(f"Occurrence table: {occurrence_path}")
    print(f"Structural QC table: {qc_path}")
    print(f"Audit metadata: {metadata_path}")
    print("------------------------------------------------------------")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
