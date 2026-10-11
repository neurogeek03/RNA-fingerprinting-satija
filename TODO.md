# TODO — until we can run the pipeline

Status legend: `[ ]` todo · `[~]` in progress · `[x]` done · `[?]` needs a decision

## A. Inputs & environment (blocking)
- [ ] **Stage the mouse Perturb-seq** (~7M cells) onto the run cluster; fill `config.conf [perturbseq].h5ad`.
- [ ] Confirm `[perturbseq].perturbation_col` / `control_class` (defaults `gene` / `non-targeting` — likely wrong for your file).
- [?] **Is the 7M set actually GWPS?** If yes, skip the own-dictionary build (`00a`/`00b`) and use the precomputed one.
- [ ] **Build the Docker image locally** from `docker/Dockerfile` → push → `apptainer pull` on Trillium (you do this).
- [ ] **Python env:** ensure the `brisc` conda env (`/scratch/mfafouti/miniforge3/envs/brisc`) has `scanpy`/`anndata`
      + `gprofiler-official`; `ryp`/R must be visible for BRISC `.rds` writing. (Ask before installing.)
- [ ] Confirm the query h5ad path in `config.conf [query.rat].h5ad` resolves on the run cluster (it's a symlink here).

## B. Code fixes still open (flagged during package verification)
- [?] **GWPS is human, not mouse.** Precomputed GWPS + `ConvertGenes` key off `org.Hs.eg.db`. So "map query→mouse,
      compare to GWPS" won't align — GWPS needs the query mapped to **human**, own dict to **mouse**. Decide:
      (a) drop GWPS and keep only the mouse dictionary, or (b) make target species **per-dictionary**
      (rat→mouse for own, rat→human for GWPS). Option (b) touches `config.conf`, `pipeline/01`, `pipeline/02`, `utils/results.R`.
- [ ] **Remove invented helper** `utils/dictionary.R::dictionary_feature_overlap()` — `FingerprintCells` already handles
      missing features (warns >25%, stops >95%). Follow the vignette instead of reading `dict$fingerprints$Lambdas`.
- [ ] **Align to real `SummarizeResults` output.** `tidy=TRUE` returns `query, cs, lbf, sign` (not
      `top_credible_set/top_lbf`). Fix `utils/results.R` (`split_credible_set` should parse `cs`) and the column
      expectations in `pipeline/03`.
- [?] **Fingerprinting granularity.** `utils/fingerprint.R` loops per subclass (my design choice); the vignette runs on
      the whole object with `condition_meta`/`group_level`. Confirm per-cell-type is what you want.

## C. Decisions to confirm
- [?] Orthology for rat: currently g:Profiler via `utils/orthology.py` (per your "adapt PROTON orthologs.py"). The query
      `var` also already has mouse-ortholog columns — keep g:Profiler, or use those as a cross-check?
- [?] Run both `group_level = TRUE` and `FALSE`, or just one? (`config.conf [run].group_level`).

## D. Run order (once A–C are settled)
Login node (internet): `00a`, `01`, GWPS install in `00b`. Compute node (Apptainer): `02`, `03`.
- [ ] `python pipeline/00a_preprocess_perturbseq.py` → reduced `.rds`
- [ ] `Rscript pipeline/00b_build_dictionary.R` → cached dictionary(ies)
- [ ] `python pipeline/01_prepare_query.py --cohort rat` → mouse-symbol query `.rds`
- [ ] `Rscript pipeline/02_fingerprint.R config.conf rat` → `results/query_rat/fingerprint_results_*.csv`
- [ ] `Rscript pipeline/03_annotate_results.R config.conf rat` → annotated nominations

## E. Verification (smoke test)
- [ ] On an interactive alloc (`n`): tiny `cells_per_perturbation` + one subclass (`037_DG_Glut`); confirm
      `SummarizeResults` returns non-empty credible sets.
- [ ] Report rat→mouse 1:1 mapping rate + dictionary feature overlap.
- [ ] Biology sanity: glucocorticoid/HPA genes (`Fkbp5`, `Nr3c1`, `Sgk1`) plausibly surface in CORT nominations.

## Later
- [ ] Human cohort: fill `[query.human]` once its expression object exists.
- [ ] Optional SLURM `submit-job` wrappers for `02`/`03`.
