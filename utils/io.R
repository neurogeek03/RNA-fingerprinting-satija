# I/O helpers for the RNA fingerprinting pipeline (R side).
#
# Pure/side-effect-light: read Seurat objects (written as .rds by BRISC in the
# Python steps), persist/load learned dictionaries, and read edgeR DE tables.
# Imported by pipeline/ scripts; never imports from pipeline/.

suppressPackageStartupMessages({
  library(Seurat)
})

#' Load a Seurat object from a .rds (produced by BRISC SingleCell$save(...)).
#' Validates that the requested assay and raw counts are present.
load_seurat <- function(rds_path, assay = "RNA") {
  if (!file.exists(rds_path)) stop(sprintf("Seurat rds not found: %s", rds_path))
  obj <- readRDS(rds_path)
  if (!inherits(obj, "Seurat")) stop(sprintf("%s is not a Seurat object", rds_path))
  if (!assay %in% Seurat::Assays(obj))
    stop(sprintf("assay '%s' not in object (have: %s)", assay,
                 paste(Seurat::Assays(obj), collapse = ", ")))
  DefaultAssay(obj) <- assay
  obj
}

#' SCTransform on raw counts, returning ALL genes (required by fingerprinting,
#' which matches query genes against dictionary features by rowname).
run_sct <- function(obj, assay = "RNA", seed.use = 1448145) {
  SCTransform(obj, assay = assay, return.only.var.genes = FALSE,
              verbose = FALSE, seed.use = seed.use)
}

#' Persist / load a learned dictionary object.
save_dictionary <- function(dictionary, rds_path) {
  dir.create(dirname(rds_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(dictionary, rds_path)
  invisible(rds_path)
}

load_dictionary <- function(rds_path) {
  if (!file.exists(rds_path)) stop(sprintf("dictionary not found: %s", rds_path))
  readRDS(rds_path)
}

#' Read all per-cell-type edgeR result TSVs in a directory into one long frame.
#' Gene id is the unlabeled first column (edgeR topTags row index).
#' celltype is parsed from the filename (strips the trailing '_edgeR_results').
read_edger_dir <- function(edger_dir, suffix = "_edgeR_results.tsv") {
  files <- list.files(edger_dir, pattern = paste0(suffix, "$"), full.names = TRUE)
  if (length(files) == 0) stop(sprintf("no '%s' files in %s", suffix, edger_dir))
  frames <- lapply(files, function(f) {
    df <- read.delim(f, check.names = FALSE)
    # first column is the gene id (unnamed in the file -> becomes 'X' or V1)
    names(df)[1] <- "gene"
    df$celltype_file <- sub(suffix, "", basename(f), fixed = TRUE)
    df
  })
  do.call(rbind, frames)
}
