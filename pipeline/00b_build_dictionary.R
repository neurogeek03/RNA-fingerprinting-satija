#!/usr/bin/env Rscript
# 00b — Learn the mouse Perturb-seq dictionary(ies).
#
# Loads the reduced Seurat .rds (from 00a), SCTransforms, and runs
# LearnDictionary -> cached .rds. Optionally also installs + loads the
# precomputed GWPS dictionary (login node; InstallPrecomputedDictionary needs
# internet). One-time step; the cached dictionaries are reused across cohorts.
#
# Usage: Rscript pipeline/00b_build_dictionary.R [config.conf]

# --- locate repo + load utils -------------------------------------------------
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(f)) return(dirname(normalizePath(f)))
  getwd()
}
REPO <- dirname(get_script_dir())
source(file.path(REPO, "utils", "config.R"))
source(file.path(REPO, "utils", "io.R"))
source(file.path(REPO, "utils", "dictionary.R"))

args <- commandArgs(trailingOnly = TRUE)
cfg <- load_config(if (length(args) >= 1) args[1] else file.path(REPO, "config.conf"))

mc_cores <- cfg_int(cfg, "run", "mc_cores", required = FALSE, default = 1L)

# --- own dictionary from the reduced Perturb-seq ------------------------------
reduced_rds    <- cfg_get(cfg, "perturbseq", "reduced_rds")
dictionary_rds <- cfg_get(cfg, "perturbseq", "dictionary_rds")
pert_col       <- cfg_get(cfg, "perturbseq", "perturbation_col")
control_class  <- cfg_get(cfg, "perturbseq", "control_class")
feature_method <- cfg_get(cfg, "perturbseq", "feature_method", required = FALSE, default = "de")
min_cells      <- cfg_int(cfg, "perturbseq", "min_cells", required = FALSE, default = 10L)

message("[00b] loading reduced Perturb-seq: ", reduced_rds)
obj <- load_seurat(reduced_rds, assay = "RNA")
message("[00b] SCTransform (all genes)")
obj <- run_sct(obj, assay = "RNA")

message("[00b] LearnDictionary (method=", feature_method, ")")
dict <- learn_dictionary(obj, perturbation_col = pert_col, control_class = control_class,
                         feature_method = feature_method, min_cells = min_cells,
                         mc.cores = mc_cores)
save_dictionary(dict, dictionary_rds)
message("[00b] saved own dictionary -> ", dictionary_rds)

# --- precomputed GWPS dictionary ---------------------------------------------
if (isTRUE(cfg_bool(cfg, "gwps", "use", required = FALSE, default = FALSE))) {
  gwps_name     <- cfg_get(cfg, "gwps", "name", required = FALSE, default = "gwps")
  gwps_key      <- cfg_get(cfg, "gwps", "gene_key", required = FALSE, default = "SYMBOL")
  gwps_rds      <- cfg_get(cfg, "gwps", "dictionary_rds")
  message("[00b] installing + loading precomputed GWPS (needs internet / login node)")
  gwps <- load_gwps_dictionary(name = gwps_name, gene_key = gwps_key, install = TRUE)
  save_dictionary(gwps, gwps_rds)
  message("[00b] saved GWPS dictionary -> ", gwps_rds)
}

message("[00b] done")
