# CLAUDE.md — RNA Fingerprinting Pipeline

## What this is
Prioritize/nominate **causal genes** behind a disease transcriptional signature by mapping it onto a
**mouse Perturb-seq perturbation dictionary** with the Satija lab
[`rna-fingerprinting`](https://github.com/satijalab/rna-fingerprinting) R package
(`LearnDictionary` → `FingerprintCells` → `SummarizeResults`, SuSiE regression; perturbed genes in the
credible sets are the nominations).

One **fixed** mouse dictionary; **swappable** query cohorts (first = rat Slide-tags OIL vs CORT; later = human).
Each query is **ortholog-mapped to mouse** before fingerprinting. Built on this cluster, then moved to the
cluster holding the Perturb-seq data + more compute.

## Layout (per global rules: utils/ pure+reusable, pipeline/ orchestration, all knobs in config.conf)
- `config.conf` — single source of truth (paths, species, columns, caps, thresholds). Read by both Python
  (`utils/config.py`, configparser) and R (`utils/config.R`, base-R INI parser).
- `utils/` — `orthology.py` (g:Profiler g:Orth, adapted from PROTON; rat/human → mouse),
  `preprocess.py` (scanpy downsampling), `io.R` (Seurat/dictionary/edgeR I/O), `dictionary.R`
  (ChooseFeatures/LearnDictionary/GWPS wrappers), `fingerprint.R` (per-cell-type FingerprintCells loop),
  `results.R` (map nominations back + join DE).
- `pipeline/` — `00a_preprocess_perturbseq.py` → `00b_build_dictionary.R` → `01_prepare_query.py`
  → `02_fingerprint.R` → `03_annotate_results.R`.
- `docker/` — Docker only (`Dockerfile`, R stack).
- Python steps (`00a`, `01`) run in the existing `brisc` conda env
  (`/scratch/mfafouti/miniforge3/envs/brisc`, has brisc+ryp+anndata) plus `gprofiler-official`.

## Key facts (verified)
- **Query (rat):** `/scratch/mfafouti/Mommybrain-PPD/Slide_tags/EdgeR/data/DE_after_mt_filter_108123_...h5ad`
  (symlink → `/scratch/mfafouti/Mommybrain/...`). 108,123 cells × 21,674 genes; **raw int UMI counts in `.X`**;
  var indexed by rat Ensembl with `gene_symbols` + already-present mouse ortholog columns; obs has
  `treatment` (CORT vs OIL, balanced 3v3) and `subclass_name` (40 subclasses). Used **as-is** (no QC/gene filter).
- **Dictionary:** mouse Perturb-seq (~7M cells, **fill `[perturbseq].h5ad`**) + precomputed GWPS.
- **Orthology:** `utils/orthology.py` adapted from `/scratch/mfafouti/Mommybrain-PPD/PROTON/utils/orthologs.py`
  (set `target_organism=mmusculus`). Static cross-check table:
  `/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/genes_mouse_rename/out/rat_to_mouse_filtered.csv`.
- **h5ad → Seurat .rds:** BRISC, `SingleCell(h5ad).save(rds, assay='RNA', X_key='counts')` (Seurat v5; needs
  `ryp`). Existing env: `/scratch/mfafouti/miniforge3/envs/brisc`.

## Running (see docs/ plan for detail)
Login node (internet): `00a`, `01`, and the GWPS install in `00b`. Compute node (Apptainer): `02`, `03`.
```
python pipeline/00a_preprocess_perturbseq.py        # shrink 7M-cell dict -> reduced .rds
Rscript pipeline/00b_build_dictionary.R             # LearnDictionary + GWPS -> cached dicts
python pipeline/01_prepare_query.py --cohort rat    # ortholog-map query -> mouse .rds
Rscript pipeline/02_fingerprint.R config.conf rat   # FingerprintCells per subclass
Rscript pipeline/03_annotate_results.R config.conf rat
```

## Conventions
- Never hard-code paths/tunables — add to `config.conf`. Add a cohort by adding a `[query.<name>]` section and
  setting `[run].cohort`.
- `utils/` never imports from `pipeline/`. Keep cohort-specific maps as arguments, not literals.
- Network steps on the login node only (Trillium compute has no internet).
- R packages via Docker → push → Apptainer; Python via conda/uv (ask before installing).
