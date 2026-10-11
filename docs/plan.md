# RNA Fingerprinting Pipeline — Causal Gene Prioritization

## Context
Prioritize/nominate **causal genes** behind a disease transcriptional signature by mapping it onto a mouse Perturb-seq
**perturbation dictionary**, using the Satija lab [`rna-fingerprinting`](https://github.com/satijalab/rna-fingerprinting)
R package. The method (`LearnDictionary` → `FingerprintCells` → `SummarizeResults`) runs SuSiE regression of a query's
condition-vs-control response onto perturbation fingerprints; perturbed genes in the credible sets are the nominated
causal genes.

**Goal:** one reusable repo — *fixed* mouse dictionary, *swappable* query cohorts (first = rat Slide-tags OIL vs CORT;
later = a human cohort). Dictionary is mouse, so each query is **ortholog-mapped to mouse** before mapping. Built here,
then moved to the cluster that holds the Perturb-seq data + more compute.

## Decisions (confirmed with user)
- **Query input:** full expression object, **used as-is — no added QC/gene filtering**. Rat query = the Slide-tags
  DE-input h5ad (`.../Slide_tags/EdgeR/data/DE_after_mt_filter_108123_mincells_10_in_2_samples_slide_tags.h5ad`,
  symlink → `/scratch/mfafouti/Mommybrain/Slide_tags/EdgeR/...`). **Verified contents:** 108,123 cells × 21,674 genes;
  **raw int UMI counts in `.X`** (CSR) → feeds SCTransform directly; var indexed by rat Ensembl (`ENSRNOG…`) with cols
  `gene_symbols` **and already-mapped mouse orthologs** (`mouse_gene_name`, `mouse_gene_stable_ID`,
  `mouse_orthology_confidence_0_1`); obs has `treatment` (CORT vs OIL, balanced 3v3 by `sample`) and `subclass_name`
  (40 subclasses). Only "filtering" applied for fingerprinting = ortholog restriction (a mapping) + optional
  config-driven min-cells-per-subclass. *Human cohort will also need an expression object.*
- **Dictionary:** build **both** — our own mouse Perturb-seq (`LearnDictionary`) **and** precomputed GWPS
  (`RNAFingerprintingData`) — and compare.
- **Perturb-seq is a ~7M-cell h5ad (dictionary side ONLY):** too big for a dense Seurat/rds. **Preprocess in
  Python/scanpy first** (cap cells per perturbation, select features, drop low-count) → reduced h5ad → `.rds`. This
  subsetting applies to the dictionary object, **not** the query.
- **Orthology:** reuse existing work, do not start from scratch. **Rat cohort:** reuse the mouse-ortholog columns already
  in the query `var`. **Human cohort / cross-check:** PROTON `utils/orthologs.py` (g:Profiler, `target="mmusculus"`).
- **Env:** build Docker **locally → push to Trillium** (user's standard flow); run via Apptainer on HPC.

## Reuse (don't reimplement)
- **Orthology — rat cohort:** the query `var` **already contains** mouse orthologs (`mouse_gene_name`,
  `mouse_gene_stable_ID`, `mouse_orthology_confidence_0_1`). Reuse these directly: keep 1:1 + high-confidence, drop
  unmapped. No API call needed.
- **Orthology — human cohort / cross-check:** `PROTON/utils/orthologs.py` —
  `build_mapping/load_or_build_mapping/best_one_to_one`, g:Profiler g:Orth, cached, species-agnostic; call with
  `target_organism="mmusculus"`. Path `/scratch/mfafouti/Mommybrain-PPD/PROTON/utils/orthologs.py`.
  - Static rat→mouse cross-check table: `/scratch/mfafouti/Mommybrain/Slide_seq/RCTD/genes_mouse_rename/out/rat_to_mouse_filtered.csv`
    (22,781 ENSRNOG↔ENSMUSG pairs; join on rat Ensembl ID).
- **Config pattern:** mirror `PROTON/config/proton.conf` + `PROTON/utils/config.py` (INI loader) for the Python side.
- **Method package:** `LearnDictionary`, `ChooseFeatures(method='de')`, `FingerprintCells`, `SummarizeResults`,
  `ConvertGenes`, `InstallPrecomputedDictionary('gwps')`, `LoadPrecomputedDictionary(gene_key='SYMBOL')`,
  `ProjectFingerprints` (viz). h5ad→rds via **BRISC** (brisc.run).

## Repo structure (`/scratch/mfafouti/RNA-fingerprinting-satija`)
```
config.conf              # INI: paths, species, column names, caps/thresholds, mc.cores, gene_key
CLAUDE.md                # project context (create)
env/Dockerfile           # R + Seurat + RNAFingerprinting(+Data) + BRISC  (built local -> Trillium)
env/environment.yml      # conda/uv: scanpy + anndata + gprofiler-official (preprocess + orthology, login node)
utils/
  config.py / config.R   # parse config.conf on each side
  orthology.py           # thin wrapper importing PROTON utils/orthologs.py (rat/human -> mouse)
  preprocess.py          # scanpy: subset/cap cells per perturbation, feature select
  io.R                   # h5ad->Seurat via BRISC; load/save dictionary (.rds)
  dictionary.R           # wrappers: ChooseFeatures, LearnDictionary; load precomputed GWPS
  fingerprint.R          # wrappers: FingerprintCells + SummarizeResults, looped per cell type
  results.R              # tidy results, map nominated genes back to query species, join DE stats
pipeline/
  00a_preprocess_perturbseq.py  # Python: shrink 7M-cell DICTIONARY h5ad -> reduced h5ad (cap cells/perturbation)
  00b_build_dictionary.R        # reduced h5ad -> rds (BRISC); LearnDictionary; + load GWPS; cache both
  01_prepare_query.py           # relabel query var genes -> mouse symbols (rat: reuse var cols; human: PROTON); keep 1:1
  02_fingerprint.R              # query AS-IS: SCTransform; per subclass x per dictionary FingerprintCells -> SummarizeResults
  03_annotate_results.R         # nominations -> query-species symbols, join DE logFC/FDR, write tables
results/                 # (gitignored) cached dictionaries + per-cohort outputs
```
`utils/` is pure/reusable (species passed as args, mappings via config); `pipeline/` is orchestration only.
All network steps (Docker push, g:Profiler, InstallPrecomputedDictionary) run on the login node.

## Approach
1. **Env:** Dockerfile (R stack) built locally, pushed to Trillium, run as Apptainer; conda/uv env for the Python steps.
2. **Config:** every path/species/column/cap/threshold in `config.conf` — `[perturbseq]` (h5ad, reduced h5ad, dict cache,
   perturbation col, control class, cells_per_perturbation cap, feature method/n), `[gwps]` (gene_key), `[orthology]`
   (source_organism per cohort, target=mmusculus, cache, gprofiler_base_url), `[query.<cohort>]` (h5ad, species,
   condition col, control class, cell-type col), `[run]` (mc.cores, group_level, FDR, out dir).
3. **00a preprocess (dictionary only):** scanpy-load the 7M-cell perturb-seq, cap cells per perturbation + select
   features → small h5ad.
4. **00b dictionary:** reduced h5ad → rds (BRISC); `ChooseFeatures`+`LearnDictionary` → cache `.rds`; also
   `InstallPrecomputedDictionary('gwps')` + load. One-time; reused across cohorts. *(If the 7M set is literally GWPS,
   lean on the precomputed dict and skip rebuild.)*
5. **01 prepare query (relabel only, no QC filtering):** rat = relabel query genes to mouse symbols using the existing
   `var` columns, keep high-confidence 1:1, drop unmapped; human = PROTON g:Profiler `target="mmusculus"`. Report
   mapping rate + feature overlap with each dictionary.
6. **02 fingerprint:** query used as-is → `SCTransform(return.only.var.genes=FALSE)` on raw `.X`; loop subclasses
   (config cell-type col `subclass_name`, condition `treatment`, control `OIL`; optional min-cells-per-subclass);
   per dictionary `FingerprintCells` (group_level TRUE & FALSE) → `SummarizeResults` (credible sets, lBF, signs).
7. **03 annotate:** map nominated mouse genes back to query species; join each subclass's edgeR logFC/FDR; write
   per-subclass + combined causal-gene tables + own-dict vs GWPS comparison.

## Verification
- Smoke test on an interactive allocation (`n`): tiny cells/perturbation cap + feature subset, run step 02 on one
  subclass (e.g. `037_DG_Glut`); confirm `SummarizeResults` returns non-empty `top_credible_set`/`top_lbf`.
- Orthology check: report rat→mouse 1:1 rate; cross-check g:Profiler vs the static reciprocal-best table; fail loudly if
  feature overlap with the dictionary is too low.
- Sanity biology: glucocorticoid/HPA genes present in the rat DE (`Fkbp5`, `Nr3c1`, `Sgk1`) should plausibly surface in
  CORT nominations.
- Full run via Apptainer on the Perturb-seq cluster; outputs under `results/<cohort>/`.

## Open items (not blocking)
- Perturb-seq h5ad path + its perturbation/guide metadata column name; confirm whether it is the GWPS dataset.
- Confirm the human cohort can supply an expression object.
- (Post-plan) Save to memory: Docker built locally → pushed to Trillium.
