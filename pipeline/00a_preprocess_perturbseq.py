#!/usr/bin/env python
"""00a — Shrink the ~7M-cell mouse Perturb-seq h5ad into a Seurat .rds.

Dictionary side only. Reads the big h5ad backed (obs-only), caps cells per
perturbation, drops tiny perturbation groups, materializes the subset, and
writes it as a Seurat v5 .rds via BRISC so the R step can LearnDictionary().

Run on a node with enough RAM for the *subset* (not the full 7M cells).
BRISC .rds writing needs the `ryp` R bridge (use the brisc conda env).

Usage:
    python pipeline/00a_preprocess_perturbseq.py [--config path/to/config.conf]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from utils.config import load_config, get_path  # noqa: E402
from utils import preprocess as pp  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=None)
    args = ap.parse_args()
    cfg = load_config(args.config)

    import anndata as ad
    from brisc import SingleCell

    h5ad = get_path(cfg, "perturbseq", "h5ad")
    perturbation_col = cfg.get("perturbseq", "perturbation_col")
    control_class = cfg.get("perturbseq", "control_class")
    cap = cfg.getint("perturbseq", "cells_per_perturbation")
    min_cells = cfg.getint("perturbseq", "min_cells")
    reduced_rds = get_path(cfg, "perturbseq", "reduced_rds")
    reduced_rds.parent.mkdir(parents=True, exist_ok=True)

    print(f"[00a] reading (backed) {h5ad}")
    adata = ad.read_h5ad(h5ad, backed="r")
    if perturbation_col not in adata.obs:
        raise KeyError(f"[perturbseq] perturbation_col '{perturbation_col}' not in .obs "
                       f"(have: {list(adata.obs.columns)})")

    keep = pp.cap_mask(adata.obs[perturbation_col].to_numpy(), cap)
    print(f"[00a] capping to <= {cap} cells/perturbation: {adata.n_obs} -> {int(keep.sum())} cells")
    adata = adata[keep].to_memory()
    adata = pp.drop_small_perturbations(adata, perturbation_col, min_cells, control_class)

    report = pp.summarize_perturbations(adata, perturbation_col, control_class)
    print(f"[00a] reduced object: {report}")
    if not report["control_present"]:
        raise ValueError(f"control_class '{control_class}' absent after subsetting — "
                         f"check [perturbseq].control_class")

    # Intermediate h5ad next to the rds, then convert via BRISC.
    tmp_h5ad = reduced_rds.with_suffix(".h5ad")
    adata.write_h5ad(tmp_h5ad)
    print(f"[00a] wrote {tmp_h5ad}; converting to Seurat .rds via BRISC")
    SingleCell(str(tmp_h5ad)).save(str(reduced_rds), assay="RNA", X_key="counts",
                                   overwrite=True)
    print(f"[00a] done -> {reduced_rds}")


if __name__ == "__main__":
    main()
