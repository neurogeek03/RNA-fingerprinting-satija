"""Cross-species gene ortholog mapping via the g:Profiler g:Orth API.

Adapted from PROTON/utils/orthologs.py. Generalized so the TARGET species is
not assumed to be human: here we map a query species (rat `rnorvegicus`, or
human `hsapiens`) onto the Perturb-seq dictionary species (mouse `mmusculus`).
Column names are therefore `target_symbol` / `target_ensg` rather than
`human_*`.

Design (reproducible + HPC-friendly):
  * The API call (`build_mapping`) needs internet -> run it on the LOGIN node.
  * It writes a cache CSV; `load_or_build_mapping` reuses that cache offline on
    compute nodes, so the pipeline never needs network at run time.
  * Pin `base_url` to a g:Profiler archive release for full reproducibility; the
    build date is recorded in the cache regardless.

g:Orth (`gp.orth(organism=<source>, query=[...], target=<target>)`) returns:
  incoming, converted, ortholog_ensg, n_incoming, n_converted, n_result,
  name (= target gene symbol), description, namespaces.
No-ortholog rows come back with 'N/A' in name/ortholog_ensg.
"""
from __future__ import annotations

import datetime as _dt
from pathlib import Path

import pandas as pd

_NA = {"N/A", "None", "nan", ""}


def build_mapping(symbols: list[str], source_organism: str,
                  target_organism: str = "mmusculus",
                  base_url: str | None = None) -> pd.DataFrame:
    """Call g:Profiler g:Orth for `symbols`. Needs internet (login node).

    Returns a tidy frame: source_symbol, target_symbol, target_ensg, built_on,
    source_organism, target_organism (one row per returned ortholog; unmapped
    queries included with blank target_symbol so nothing is silently dropped).
    """
    from gprofiler import GProfiler  # imported here so callers w/o mapping don't need it

    gp = GProfiler(return_dataframe=True, base_url=base_url) if base_url \
        else GProfiler(return_dataframe=True)
    res = gp.orth(organism=source_organism, query=list(dict.fromkeys(symbols)),
                  target=target_organism)

    tidy = (
        res[["incoming", "name", "ortholog_ensg"]]
        .rename(columns={"incoming": "source_symbol",
                         "name": "target_symbol",
                         "ortholog_ensg": "target_ensg"})
        .copy()
    )
    # normalize g:Profiler 'N/A' sentinels to empty (= unmapped)
    for c in ("target_symbol", "target_ensg"):
        tidy[c] = tidy[c].astype(str).where(~tidy[c].astype(str).isin(_NA), "")
    tidy["built_on"] = _dt.date.today().isoformat()
    tidy["source_organism"] = source_organism
    tidy["target_organism"] = target_organism
    return tidy.reset_index(drop=True)


def load_or_build_mapping(symbols: list[str], cache_path: str | Path,
                          source_organism: str, target_organism: str = "mmusculus",
                          base_url: str | None = None,
                          rebuild: bool = False) -> pd.DataFrame:
    """Return the ortholog mapping, building+caching via API only if needed.

    On a compute node the cache must already exist (rebuild would need internet).
    """
    cache_path = Path(cache_path)
    if cache_path.is_file() and not rebuild:
        return pd.read_csv(cache_path, dtype=str).fillna("")
    mapping = build_mapping(symbols, source_organism, target_organism, base_url)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    mapping.to_csv(cache_path, index=False)
    return mapping


def summarize_mapping(mapping: pd.DataFrame, queried: list[str]) -> dict:
    """Match-rate + ambiguity report. No side effects."""
    mapped = mapping[mapping["target_symbol"] != ""]
    mapped_src = set(mapped["source_symbol"])
    queried_set = set(queried)
    # one source -> many target, and one target <- many source
    per_src = mapped.groupby("source_symbol")["target_symbol"].nunique()
    per_tgt = mapped.groupby("target_symbol")["source_symbol"].nunique()
    return {
        "n_queried": len(queried_set),
        "n_mapped": len(mapped_src),
        "match_rate": round(len(mapped_src) / max(len(queried_set), 1), 4),
        "dropped": sorted(queried_set - mapped_src),
        "source_multi": sorted(per_src[per_src > 1].index),   # ambiguous source genes
        "target_multi": sorted(per_tgt[per_tgt > 1].index),   # many-source collisions
    }


def best_one_to_one(mapping: pd.DataFrame) -> pd.DataFrame:
    """Keep only unambiguous 1:1 source<->target pairs (drop many-to-* rows).

    Conservative: a source or target symbol appearing in >1 pair is excluded.
    Callers should log what `summarize_mapping` flagged as ambiguous.
    """
    mapped = mapping[mapping["target_symbol"] != ""].drop_duplicates(
        subset=["source_symbol", "target_symbol"])
    src_ok = mapped.groupby("source_symbol")["target_symbol"].transform("nunique") == 1
    tgt_ok = mapped.groupby("target_symbol")["source_symbol"].transform("nunique") == 1
    return mapped[src_ok & tgt_ok].reset_index(drop=True)
