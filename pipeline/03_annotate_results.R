#!/usr/bin/env Rscript
# 03 — Annotate nominated causal genes and join cohort DE statistics.
#
# For each dictionary's fingerprint results (from 02): map nominated MOUSE genes
# back to the query species via the inverted ortholog map, join the cohort's
# edgeR logFC/FDR per cell type, and write annotated tables + a cross-dictionary
# comparison of nominations.
#
# Usage: Rscript pipeline/03_annotate_results.R [config.conf] [cohort]

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(f)) return(dirname(normalizePath(f)))
  getwd()
}
REPO <- dirname(get_script_dir())
source(file.path(REPO, "utils", "config.R"))
source(file.path(REPO, "utils", "io.R"))
source(file.path(REPO, "utils", "results.R"))

args   <- commandArgs(trailingOnly = TRUE)
cfg    <- load_config(if (length(args) >= 1) args[1] else file.path(REPO, "config.conf"))
cohort <- if (length(args) >= 2) args[2] else cfg_get(cfg, "run", "cohort")
qs     <- query_section(cfg, cohort)

out_dir        <- file.path(cfg_get(cfg, "paths", "results_dir"), paste0("query_", cohort))
ortholog_cache <- cfg_get(cfg, qs, "ortholog_cache")
edger_dir      <- cfg_get(cfg, qs, "edger_dir", required = FALSE, default = NULL)

message("[03] inverting ortholog map: ", ortholog_cache)
t2s <- invert_ortholog_map(ortholog_cache)

edger_long <- NULL
if (!is.null(edger_dir) && dir.exists(edger_dir)) {
  message("[03] reading edgeR DE tables: ", edger_dir)
  edger_long <- read_edger_dir(edger_dir)
}

dict_names <- c("own")
if (isTRUE(cfg_bool(cfg, "gwps", "use", required = FALSE, default = FALSE)))
  dict_names <- c(dict_names, "gwps")

annotated <- list()
for (dname in dict_names) {
  in_csv <- file.path(out_dir, sprintf("fingerprint_results_%s.csv", dname))
  if (!file.exists(in_csv)) { message("[03] skip missing ", in_csv); next }
  res <- read.csv(in_csv, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"celltype" %in% names(res)) stop("expected a 'celltype' column in ", in_csv)
  res$celltype_file <- gsub(" ", "_", res$celltype)  # match edgeR filename key
  ann <- annotate_nominations(res, t2s, edger_long = edger_long,
                              celltype_join_col = "celltype_file")
  out_csv <- file.path(out_dir, sprintf("nominations_annotated_%s.csv", dname))
  write.csv(ann, out_csv, row.names = FALSE)
  message("[03] wrote ", out_csv, " (", nrow(ann), " rows)")
  ann$dictionary <- dname
  annotated[[dname]] <- ann
}

# Cross-dictionary comparison: which (celltype, query gene) nominations are shared.
if (length(annotated) >= 2) {
  keys <- lapply(annotated, function(a)
    unique(paste(a$celltype, a$nominated_query_gene, sep = "||")))
  shared <- Reduce(intersect, keys)
  cmp <- data.frame(celltype_gene = shared)
  if (nrow(cmp) > 0) {
    split_key <- do.call(rbind, strsplit(cmp$celltype_gene, "||", fixed = TRUE))
    cmp$celltype <- split_key[, 1]; cmp$nominated_query_gene <- split_key[, 2]
    cmp$celltype_gene <- NULL
  }
  out_cmp <- file.path(out_dir, "nominations_shared_own_vs_gwps.csv")
  write.csv(cmp, out_cmp, row.names = FALSE)
  message("[03] wrote ", out_cmp, " (", nrow(cmp), " shared nominations)")
}

message("[03] done")
