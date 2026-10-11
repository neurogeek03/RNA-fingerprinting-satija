# Fingerprinting wrappers: map a query onto a dictionary, per cell type.
#
# Parameterised so the same code serves any cohort/dictionary. The query object
# is used AS-IS (already SCT-normalized upstream); we only subset by cell type.

suppressPackageStartupMessages({
  library(Seurat)
  library(RNAFingerprinting)
})

#' Run fingerprinting for a single already-subset object (one cell type).
#' Returns the SummarizeResults tidy data.frame (credible sets, lBFs, signs).
fingerprint_one <- function(obj, condition_col, control_class, dictionary,
                            group_level = TRUE, mc.cores = 1) {
  obj <- FingerprintCells(obj, condition_meta = condition_col,
                          control_class = control_class, dictionary = dictionary,
                          group_level = group_level, mc.cores = mc.cores)
  SummarizeResults(obj, group_level = group_level, condition_meta = condition_col,
                   control_class = control_class, tidy = TRUE)
}

#' Loop fingerprinting over every cell type in `celltype_col`.
#'
#' Skips cell types with fewer than `min_cells` cells, or that lack both the
#' control class and >=1 treated cell. Returns one long data.frame with a
#' `celltype` column; the per-cell-type result also carries the dictionary name.
fingerprint_all_celltypes <- function(obj, celltype_col, condition_col,
                                      control_class, dictionary,
                                      dictionary_name = "dictionary",
                                      group_level = TRUE, min_cells = 50,
                                      mc.cores = 1) {
  meta <- obj[[]]
  celltypes <- sort(unique(as.character(meta[[celltype_col]])))
  out <- list()
  for (ct in celltypes) {
    cells <- rownames(meta)[as.character(meta[[celltype_col]]) == ct]
    if (length(cells) < min_cells) {
      message(sprintf("skip %s: %d cells < min %d", ct, length(cells), min_cells))
      next
    }
    sub <- subset(obj, cells = cells)
    conds <- unique(as.character(sub[[]][[condition_col]]))
    if (!(control_class %in% conds) || length(setdiff(conds, control_class)) == 0) {
      message(sprintf("skip %s: needs control + >=1 treated condition (have: %s)",
                      ct, paste(conds, collapse = ",")))
      next
    }
    res <- tryCatch(
      fingerprint_one(sub, condition_col, control_class, dictionary,
                      group_level = group_level, mc.cores = mc.cores),
      error = function(e) { message(sprintf("fail %s: %s", ct, conditionMessage(e))); NULL }
    )
    if (is.null(res) || nrow(res) == 0) next
    res$celltype <- ct
    res$dictionary <- dictionary_name
    out[[ct]] <- res
  }
  if (length(out) == 0) stop("no cell type produced a fingerprint result")
  do.call(rbind, out)
}
