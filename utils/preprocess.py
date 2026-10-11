"""Perturb-seq downsampling helpers (dictionary side only).

Pure / side-effect-light functions used by pipeline/00a_preprocess_perturbseq.py
to shrink a ~7M-cell Perturb-seq AnnData into an object small enough to convert
to a Seurat .rds and feed LearnDictionary(). These operate on the DICTIONARY
object only — the query object is used as-is.
"""
from __future__ import annotations

import numpy as np


def cap_mask(labels, cap: int, seed: int = 1448145):
    """Boolean mask that caps each label group to at most `cap` members.

    `labels` is an array-like of group labels (one per cell). Groups with <= cap
    members are kept whole. Deterministic given `seed`. Works on a backed
    AnnData's .obs column (reads obs only, no matrix load).
    """
    labels = np.asarray(labels).astype(str)
    rng = np.random.default_rng(seed)
    keep = np.zeros(labels.shape[0], dtype=bool)
    for grp in np.unique(labels):
        idx = np.flatnonzero(labels == grp)
        if idx.size > cap:
            idx = rng.choice(idx, size=cap, replace=False)
        keep[idx] = True
    return keep


def cap_cells_per_perturbation(adata, perturbation_col: str, cap: int,
                               seed: int = 1448145):
    """Subsample each perturbation group down to at most `cap` cells.

    Returns a new (view-materialized) AnnData. Groups with <= cap cells are kept
    whole. Deterministic given `seed`.
    """
    if perturbation_col not in adata.obs:
        raise KeyError(f"perturbation_col '{perturbation_col}' not in .obs")
    keep = cap_mask(adata.obs[perturbation_col].to_numpy(), cap, seed)
    return adata[keep].copy()


def drop_small_perturbations(adata, perturbation_col: str, min_cells: int,
                             control_class: str | None = None):
    """Drop perturbation groups with fewer than `min_cells` cells.

    The control class is always retained regardless of size.
    """
    counts = adata.obs[perturbation_col].astype(str).value_counts()
    keep_groups = set(counts[counts >= min_cells].index)
    if control_class is not None:
        keep_groups.add(control_class)
    mask = adata.obs[perturbation_col].astype(str).isin(keep_groups).to_numpy()
    return adata[mask].copy()


def summarize_perturbations(adata, perturbation_col: str, control_class: str) -> dict:
    """Quick report: #perturbations, #cells, control presence. No side effects."""
    vc = adata.obs[perturbation_col].astype(str).value_counts()
    return {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_perturbations": int(vc.size),
        "control_present": control_class in set(vc.index),
        "control_n_cells": int(vc.get(control_class, 0)),
        "median_cells_per_perturbation": float(vc.median()),
    }
