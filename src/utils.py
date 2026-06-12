from __future__ import annotations

from pathlib import Path
from typing import Literal, Optional

import anndata as ad
import numpy as np
import pandas as pd

DATA_ROOT = Path(__file__).resolve().parent.parent / "datasets"


def _read_expression_matrix(path: Path) -> pd.DataFrame:
    """Read a whitespace-delimited expression matrix (rows=genes, cols=cells)."""
    return pd.read_csv(path, sep=r"\s+", index_col=0)


def load_smartseq(
    cell_line: Literal["MCF7", "HCC1806"],
    variant: Literal["filtered_normalised", "filtered", "unfiltered"] = "filtered_normalised",
    split: Literal["train", "test"] = "train",
    data_root: Path | None = None,
) -> tuple[pd.DataFrame, pd.Series | None]:
    """Return (expression_df, labels) for a SmartSeq dataset.

    expression_df has shape (n_cells, n_genes) — cells are rows, genes are columns
    (transposed from the on-disk layout) so it matches sklearn / scanpy conventions.
    Labels are 'Hypo'/'Norm' or None for the anonymized test split.
    """
    root = Path(data_root) if data_root else DATA_ROOT
    base = root / "SmartSeq"
    stem = f"{cell_line}_SmartS"

    if variant == "filtered_normalised":
        if split == "train":
            fname = f"{stem}_Filtered_Normalised_3000_Data_train.txt"
        else:
            fname = f"{stem}_Filtered_Normalised_3000_Data_test_anonim.txt"
    elif variant == "filtered":
        if split != "train":
            raise ValueError("filtered variant only has a single (train-like) file")
        fname = f"{stem}_Filtered_Data.txt"
    elif variant == "unfiltered":
        if split != "train":
            raise ValueError("unfiltered variant only has a single (train-like) file")
        fname = f"{stem}_Unfiltered_Data.txt"
    else:
        raise ValueError(f"unknown variant: {variant}")

    expr = _read_expression_matrix(base / fname).T  # cells x genes

    if split == "test":
        return expr, None

    meta = pd.read_csv(base / f"{stem}_MetaData.tsv", sep="\t", index_col=0)
    # Map BAM filename -> Condition (Hypo/Norm). HCC1806 metadata uses 'Normo'
    # while MCF7 uses 'Norm' — collapse to 'Norm' so labels are consistent
    # across cell lines (and align with the DropSeq mapping below).
    raw = meta["Condition"].reindex(expr.index)
    if raw.isna().any():
        missing = raw.index[raw.isna()].tolist()[:3]
        raise ValueError(f"SmartSeq cells missing from metadata, e.g. {missing}")
    labels = raw.replace({"Normo": "Norm"})
    return expr, labels


def load_dropseq(
    cell_line: Literal["MCF7", "HCC1806"],
    split: Literal["train", "test"] = "train",
    data_root: Path | None = None,
) -> tuple[pd.DataFrame, pd.Series | None]:
    """Return (expression_df, labels) for a DropSeq dataset (filtered+normalised)."""
    root = Path(data_root) if data_root else DATA_ROOT
    base = root / "DropSeq"

    if split == "train":
        fname = f"{cell_line}_Filtered_Normalised_3000_Data_train.txt"
    else:
        fname = f"{cell_line}_Filtered_Normalised_3000_Data_test_anonim.txt"

    expr = _read_expression_matrix(base / fname).T  # cells x genes

    if split == "test":
        return expr, None

    # Suffix after the last underscore: 'Normoxia' or 'Hypoxia'
    suffixes = expr.index.to_series().str.rsplit("_", n=1).str[-1]
    mapping = {"Normoxia": "Norm", "Hypoxia": "Hypo"}
    if not set(suffixes.unique()).issubset(mapping):
        raise ValueError(f"unexpected DropSeq label suffixes: {suffixes.unique()}")
    labels = suffixes.map(mapping)
    labels.name = "Condition"
    return expr, labels


def to_anndata(
    expr: pd.DataFrame,
    labels: pd.Series | None = None,
    *,
    dataset: str = "",
) -> ad.AnnData:
    """Wrap an (cells x genes) DataFrame as an AnnData object for scanpy."""
    adata = ad.AnnData(
        X=expr.values.astype(np.float32),
        obs=pd.DataFrame(index=expr.index.astype(str)),
        var=pd.DataFrame(index=expr.columns.astype(str)),
    )
    adata.obs_names_make_unique()
    if labels is not None:
        adata.obs["condition"] = pd.Categorical(labels.values, categories=["Norm", "Hypo"])
    if dataset:
        adata.uns["dataset"] = dataset
    return adata


# Convenience: all four train datasets in one dict
def load_all_train():
    """Load all four (cell_line, technology) train datasets."""
    return {
        "MCF7_SmartSeq": load_smartseq("MCF7", split="train"),
        "HCC1806_SmartSeq": load_smartseq("HCC1806", split="train"),
        "MCF7_DropSeq": load_dropseq("MCF7", split="train"),
        "HCC1806_DropSeq": load_dropseq("HCC1806", split="train"),
    }
