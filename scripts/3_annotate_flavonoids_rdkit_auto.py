"""RDKit annotation workflow for flavonoid structural features.

This version auto-detects the CSV exported by Part II of the R pipeline,
so the user does not need to rename files or edit the R scripts.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, Optional, Union
import shutil
import sys

import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


DEFAULT_INPUT_NAME = "lotus_flavonoids_for_rdkit.csv"
DEFAULT_OUTPUT_NAME = "lotus_flavonoids_rdkit_annotations.csv"

SMARTS_PATTERNS = {
    "has_phenolic_OH": "[$([OX2H][cX3]:[cX3])]",
    "has_methoxy_aryl": "[c:1][OX2][CH3]",
    "has_prenyl_like": "C=C(C)C",
    "has_conj_carbonyl": "[CX3](=O)[CX3]=[CX3]",
}


def smiles_to_mol(smiles: str) -> Optional[Chem.Mol]:
    """Convert a SMILES string to an RDKit molecule object."""
    try:
        return Chem.MolFromSmiles(smiles)
    except Exception:
        return None


def compile_smarts_patterns(patterns: Dict[str, str]) -> Dict[str, Chem.Mol]:
    """Compile SMARTS queries and fail early if any pattern is invalid."""
    compiled: Dict[str, Chem.Mol] = {}
    for name, smarts in patterns.items():
        query = Chem.MolFromSmarts(smarts)
        if query is None:
            raise ValueError(f"Invalid SMARTS pattern for '{name}': {smarts}")
        compiled[name] = query
    return compiled


def has_probable_sugar(mol: Chem.Mol, min_ring_oxygens: int = 3) -> bool:
    """Detect non-aromatic oxygen-rich 5- or 6-membered rings consistent with glycosides."""
    for ring in mol.GetRingInfo().AtomRings():
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

        oxygen_ids = set()

        for atom in ring_atoms:
            if atom.GetSymbol() == "O":
                is_carbonyl_oxygen = any(
                    bond.GetBondType() == Chem.rdchem.BondType.DOUBLE
                    and bond.GetOtherAtom(atom).GetSymbol() == "C"
                    for bond in atom.GetBonds()
                )
                if not is_carbonyl_oxygen:
                    oxygen_ids.add(atom.GetIdx())

            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() != "O":
                    continue

                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                    continue

                is_carbonyl_oxygen = any(
                    other_bond.GetBondType() == Chem.rdchem.BondType.DOUBLE
                    and other_bond.GetOtherAtom(neighbor).GetSymbol() == "C"
                    for other_bond in neighbor.GetBonds()
                )
                if not is_carbonyl_oxygen:
                    oxygen_ids.add(neighbor.GetIdx())

        if len(oxygen_ids) >= min_ring_oxygens:
            return True

    return False


def annotate_molecule(
    mol: Chem.Mol, compiled_patterns: Dict[str, Chem.Mol]
) -> Dict[str, Union[int, bool]]:
    """Compute structural annotations and global RDKit descriptors for one molecule."""
    annotations = {
        name: mol.HasSubstructMatch(pattern)
        for name, pattern in compiled_patterns.items()
    }
    annotations["has_probable_sugar"] = has_probable_sugar(mol)
    annotations["num_HBD"] = rdMolDescriptors.CalcNumHBD(mol)
    annotations["num_HBA"] = rdMolDescriptors.CalcNumHBA(mol)
    annotations["num_rings"] = rdMolDescriptors.CalcNumRings(mol)
    return annotations


def candidate_input_paths() -> Iterable[Path]:
    """Yield likely input CSV paths in priority order."""
    cwd = Path.cwd()

    explicit = [
        cwd / DEFAULT_INPUT_NAME,
        Path(__file__).resolve().parent / DEFAULT_INPUT_NAME,
    ]
    for path in explicit:
        yield path

    for pattern in (
        "outputs/**/PartII_ALL/*__flavonoids_for_rdkit.csv",
        "**/PartII_ALL/*__flavonoids_for_rdkit.csv",
        "**/*__flavonoids_for_rdkit.csv",
    ):
        matches = sorted(
            cwd.glob(pattern),
            key=lambda p: p.stat().st_mtime,
            reverse=True,
        )
        for path in matches:
            yield path


def resolve_input_csv() -> Path:
    """Find the best available input CSV automatically."""
    seen = set()
    for path in candidate_input_paths():
        try:
            resolved = path.resolve()
        except FileNotFoundError:
            resolved = path.absolute()
        if resolved in seen:
            continue
        seen.add(resolved)
        if path.exists() and path.is_file():
            return path

    raise FileNotFoundError(
        "No RDKit input CSV was found. Expected either '\n"
        f"- {DEFAULT_INPUT_NAME} in the current folder, or\n"
        "- a file matching outputs/**/PartII_ALL/*__flavonoids_for_rdkit.csv\n"
        "Run Part II first."
    )


def infer_output_targets(input_csv: Path) -> list[Path]:
    """Build output paths that are compatible with Part III search logic."""
    cwd = Path.cwd()
    targets: list[Path] = [cwd / DEFAULT_OUTPUT_NAME]

    # Save beside the detected input for traceability.
    targets.append(input_csv.with_name(DEFAULT_OUTPUT_NAME))

    # If the input is inside OUT_DIR/PartII_ALL, also save in OUT_DIR root,
    # because Part III explicitly searches there.
    if input_csv.parent.name == "PartII_ALL":
        targets.append(input_csv.parent.parent / DEFAULT_OUTPUT_NAME)

    unique_targets: list[Path] = []
    seen = set()
    for path in targets:
        resolved = path.resolve(strict=False)
        if resolved in seen:
            continue
        seen.add(resolved)
        unique_targets.append(path)
    return unique_targets


def save_outputs(df: pd.DataFrame, targets: list[Path]) -> None:
    """Write the annotated CSV to all relevant target paths."""
    primary = targets[0]
    primary.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(primary, index=False)

    for extra in targets[1:]:
        extra.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(primary, extra)


def main() -> None:
    input_csv = resolve_input_csv()
    output_targets = infer_output_targets(input_csv)

    print(f"Input detected: {input_csv}")
    print("Outputs to be written:")
    for target in output_targets:
        print(f" - {target}")

    df = pd.read_csv(input_csv)
    if "smiles" not in df.columns:
        raise KeyError("Input CSV does not contain the required 'smiles' column.")

    df["mol"] = df["smiles"].apply(smiles_to_mol)
    df = df.loc[df["mol"].notna()].copy()

    compiled_patterns = compile_smarts_patterns(SMARTS_PATTERNS)
    annotations = pd.DataFrame.from_records(
        [annotate_molecule(mol, compiled_patterns) for mol in df["mol"]]
    )

    annotated_df = pd.concat([df.drop(columns=["mol"]), annotations], axis=1)
    save_outputs(annotated_df, output_targets)

    print(f"RDKit flavonoid annotations saved successfully. Rows exported: {len(annotated_df)}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
