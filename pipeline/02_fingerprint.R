#!/usr/bin/env Rscript
# 02 — Map the (mouse-relabeled) query onto the dictionary(ies), per cell type.
#
# Loads the query .rds from 01, SCTransforms, then for each configured dictionary
# runs FingerprintCells + SummarizeResults over every cell type. Writes a tidy
# results table per dictionary under results/<cohort>/.
#
# Usage: Rscript pipeline/02_fingerprint.R [config.conf] [cohort]

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
source(file.path(REPO, "utils", "fingerprint.R"))

args    <- commandArgs(trailingOnly = TRUE)
cfg     <- load_config(if (length(args) >= 1) args[1] else file.path(REPO, "config.conf"))
cohort  <- if (length(args) >= 2) args[2] else cfg_get(cfg, "run", "cohort")
qs      <- query_section(cfg, cohort)

condition_col <- cfg_get(cfg, qs, "condition_col")
control_class <- cfg_get(cfg, qs, "control_class")
celltype_col  <- cfg_get(cfg, qs, "celltype_col")
mapped_rds    <- cfg_get(cfg, qs, "mapped_rds")
group_level   <- cfg_bool(cfg, "run", "group_level", required = FALSE, default = TRUE)
mc_cores      <- cfg_int(cfg, "run", "mc_cores", required = FALSE, default = 1L)
min_cells     <- cfg_int(cfg, "run", "min_cells_per_celltype", required = FALSE, default = 50L)

out_dir <- file.path(cfg_get(cfg, "paths", "results_dir"), paste0("query_", cohort))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("[02] loading query: ", mapped_rds)
obj <- load_seurat(mapped_rds, assay = "RNA")
message("[02] SCTransform (all genes)")
obj <- run_sct(obj, assay = "RNA")
query_genes <- rownames(GetAssayData(obj, assay = "SCT", layer = "scale.data"))

# Dictionaries to run: own + (optionally) GWPS.
dicts <- list(own = cfg_get(cfg, "perturbseq", "dictionary_rds"))
if (isTRUE(cfg_bool(cfg, "gwps", "use", required = FALSE, default = FALSE)))
  dicts[["gwps"]] <- cfg_get(cfg, "gwps", "dictionary_rds")

for (dname in names(dicts)) {
  message("[02] dictionary: ", dname)
  dict <- load_dictionary(dicts[[dname]])
  ov <- dictionary_feature_overlap(dict, query_genes, min_overlap = 50)
  message(sprintf("[02]   feature overlap: %d of %d dict features in query",
                  ov$n_overlap, ov$n_features))
  res <- fingerprint_all_celltypes(obj, celltype_col = celltype_col,
                                   condition_col = condition_col,
                                   control_class = control_class, dictionary = dict,
                                   dictionary_name = dname, group_level = group_level,
                                   min_cells = min_cells, mc.cores = mc_cores)
  out_csv <- file.path(out_dir, sprintf("fingerprint_results_%s.csv", dname))
  write.csv(res, out_csv, row.names = FALSE)
  message("[02]   wrote ", out_csv, " (", nrow(res), " rows)")
}

message("[02] done")
