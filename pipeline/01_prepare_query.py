#!/usr/bin/env python
"""01 — Relabel a query expression object's genes to MOUSE symbols.

The query is used AS-IS (no QC/gene filtering). We only cross-map its genes to
the dictionary species (mouse) via g:Profiler g:Orth (utils/orthology.py,
adapted from PROTON), keep unambiguous 1:1 orthologs, and write the relabeled
object as a Seurat .rds via BRISC for the R fingerprinting step.

The ortholog cache is built here (login node, needs internet) and reused
offline later. Reports mapping rate; nothing is silently dropped.

Usage:
    python pipeline/01_prepare_query.py [--config config.conf] [--cohort rat]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from utils.config import load_config, get_path, query_section  # noqa: E402
from utils import orthology as orth  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=None)
    ap.add_argument("--cohort", default=None, help="query cohort (default: [run].cohort)")
    args = ap.parse_args()
    cfg = load_config(args.config)
    qs = query_section(cfg, args.cohort)

    import anndata as ad
    from brisc import SingleCell

    h5ad = get_path(cfg, qs, "h5ad")
    species = cfg.get(qs, "species")
    symbol_col = cfg.get(qs, "gene_symbol_col")
    target = cfg.get("orthology", "target_organism")
    base_url = cfg.get("orthology", "gprofiler_base_url", fallback="").strip() or None
    cache = get_path(cfg, qs, "ortholog_cache")
    mapped_rds = get_path(cfg, qs, "mapped_rds")
    mapped_rds.parent.mkdir(parents=True, exist_ok=True)
    cache.parent.mkdir(parents=True, exist_ok=True)

    print(f"[01] reading query {h5ad}")
    adata = ad.read_h5ad(h5ad)
    if symbol_col not in adata.var:
        raise KeyError(f"[{qs}] gene_symbol_col '{symbol_col}' not in var "
                       f"(have: {list(adata.var.columns)})")
    src_symbols = adata.var[symbol_col].astype(str).tolist()

    print(f"[01] mapping {species} -> {target} via g:Profiler (cache: {cache})")
    mapping = orth.load_or_build_mapping(src_symbols, cache, source_organism=species,
                                         target_organism=target, base_url=base_url)
    report = orth.summarize_mapping(mapping, src_symbols)
    print(f"[01] mapping: {report['n_mapped']}/{report['n_queried']} "
          f"(rate {report['match_rate']}); ambiguous source={len(report['source_multi'])}, "
          f"target-collisions={len(report['target_multi'])}")
    one2one = orth.best_one_to_one(mapping)
    s2t = dict(zip(one2one["source_symbol"], one2one["target_symbol"]))
    print(f"[01] kept {len(s2t)} unambiguous 1:1 orthologs")

    # Keep only genes with a 1:1 mouse ortholog; relabel var_names to mouse symbol.
    adata.var["mouse_symbol"] = [s2t.get(s, "") for s in src_symbols]
    keep = adata.var["mouse_symbol"].astype(bool).to_numpy()
    adata = adata[:, keep].copy()
    adata.var_names = adata.var["mouse_symbol"].astype(str).values
    # Drop any residual duplicate mouse symbols (keep first occurrence).
    dup = adata.var_names.duplicated()
    if dup.any():
        print(f"[01] dropping {int(dup.sum())} duplicate mouse symbols (keep first)")
        adata = adata[:, ~dup].copy()
    print(f"[01] relabeled query: {adata.n_obs} cells x {adata.n_vars} mouse-symbol genes")

    tmp_h5ad = mapped_rds.with_suffix(".h5ad")
    adata.write_h5ad(tmp_h5ad)
    print(f"[01] wrote {tmp_h5ad}; converting to Seurat .rds via BRISC")
    SingleCell(str(tmp_h5ad)).save(str(mapped_rds), assay="RNA", X_key="counts",
                                   overwrite=True)
    print(f"[01] done -> {mapped_rds}")


if __name__ == "__main__":
    main()
