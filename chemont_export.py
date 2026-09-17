"""Helpers for carrying ChemOnt structural-taxonomy rows into workbook exports."""

from __future__ import annotations

from typing import Iterable

import pandas as pd


CHEMONT_EXPORT_COLUMNS: tuple[str, ...] = (
    "InChIKey",
    "Standardized SMILES",
    "ChemOnt_Kingdom",
    "ChemOnt_Superclass",
    "ChemOnt_Class",
    "ChemOnt_Subclass",
    "ChemOnt_Direct_Parent",
    "ChemOnt_Molecular_Framework",
    "ChemOnt_Substituents",
    "ChemOnt_Classification_Version",
    "ChemOnt_Source",
    "ChemOnt_Retrieval_Status",
)

CHEMONT_DISPLAY_LABELS: dict[str, str] = {
    "InChIKey": "InChIKey",
    "Standardized SMILES": "Standardized SMILES",
    "ChemOnt_Kingdom": "ChemOnt Kingdom",
    "ChemOnt_Superclass": "ChemOnt Superclass",
    "ChemOnt_Class": "ChemOnt Class",
    "ChemOnt_Subclass": "ChemOnt Subclass",
    "ChemOnt_Direct_Parent": "ChemOnt Direct Parent",
    "ChemOnt_Molecular_Framework": "ChemOnt Molecular Framework",
    "ChemOnt_Substituents": "ChemOnt Substituents",
    "ChemOnt_Classification_Version": "ChemOnt Version",
    "ChemOnt_Source": "ChemOnt Source",
    "ChemOnt_Retrieval_Status": "ChemOnt Retrieval Status",
}


def build_chemont_bundle(
    inchikey: str = "",
    standardized_smiles: str = "",
    source: str = "Unavailable",
    retrieval_status: str = "Unavailable",
    **values,
) -> dict:
    """Build a normalized ChemOnt bundle for API payloads and workbook rows."""
    bundle = {
        "InChIKey": inchikey or "",
        "Standardized SMILES": standardized_smiles or "",
        "ChemOnt_Kingdom": "",
        "ChemOnt_Superclass": "",
        "ChemOnt_Class": "",
        "ChemOnt_Subclass": "",
        "ChemOnt_Direct_Parent": "",
        "ChemOnt_Molecular_Framework": "",
        "ChemOnt_Substituents": "",
        "ChemOnt_Classification_Version": "",
        "ChemOnt_Source": source,
        "ChemOnt_Retrieval_Status": retrieval_status,
    }
    for key in CHEMONT_EXPORT_COLUMNS:
        if key in values and values[key] is not None:
            bundle[key] = values[key]
    return bundle


def attach_chemont_classifications(df: pd.DataFrame, columns: Iterable[str] | None = None) -> pd.DataFrame:
    """Store the present ChemOnt columns on ``df.attrs`` for workbook export."""
    result = df.copy()
    chemont_columns = tuple(columns or CHEMONT_EXPORT_COLUMNS)
    present = [column for column in chemont_columns if column in result.columns]
    if present:
        result.attrs["chemont_classifications"] = result[present].copy()
    return result


def chemont_fields_present(df: pd.DataFrame, columns: Iterable[str] | None = None) -> list[str]:
    chemont_columns = tuple(columns or CHEMONT_EXPORT_COLUMNS)
    return [column for column in chemont_columns if column in df.columns]


def chemont_display_frame(df: pd.DataFrame) -> pd.DataFrame:
    """Return a copy with visible ChemOnt headers standardized for reports."""
    result = df.copy()
    rename_map = {column: label for column, label in CHEMONT_DISPLAY_LABELS.items() if column in result.columns}
    return result.rename(columns=rename_map)
