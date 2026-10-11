# Result annotation: map nominated (mouse) causal genes back to the query
# species and join the cohort's edgeR DE statistics.
#
# Dictionary nominations are MOUSE perturbation gene symbols. To interpret them
# in the query species we invert the ortholog map built in pipeline/01.

#' Build a mouse(target) -> query(source) symbol lookup from the ortholog map
#' CSV written by utils/orthology.py (cols: source_symbol, target_symbol, ...).
invert_ortholog_map <- function(ortholog_csv) {
  m <- read.csv(ortholog_csv, stringsAsFactors = FALSE)
  m <- m[nzchar(m$target_symbol) & nzchar(m$source_symbol), c("source_symbol", "target_symbol")]
  m <- m[!duplicated(m$target_symbol), ]
  setNames(m$source_symbol, m$target_symbol)
}

#' Split a credible-set string into individual gene symbols (best-effort).
#' Handles ',', ';', '|' separators and strips leading +/- sign annotations.
split_credible_set <- function(x) {
  if (is.na(x) || !nzchar(x)) return(character(0))
  parts <- unlist(strsplit(x, "[,;|]"))
  parts <- trimws(parts)
  parts <- sub("^[+-]\\s*", "", parts)
  parts[nzchar(parts)]
}

#' Expand a SummarizeResults tidy frame to one row per (celltype, condition,
#' nominated mouse gene) using the top credible set, add the query-species symbol,
#' and left-join the cohort's edgeR logFC/FDR for that gene within that cell type.
#'
#' @param results   fingerprint results (needs `celltype`, `top_credible_set`,
#'                   and a condition column; `top_lbf` carried through if present).
#' @param target_to_source named vector from invert_ortholog_map().
#' @param edger_long edgeR long frame (gene, logFC, FDR, celltype_file) from
#'                   io.R::read_edger_dir(); may be NULL to skip the DE join.
#' @param celltype_join_col column in `results` matching edgeR `celltype_file`.
annotate_nominations <- function(results, target_to_source, edger_long = NULL,
                                 celltype_join_col = "celltype_file") {
  rows <- list()
  for (i in seq_len(nrow(results))) {
    genes <- split_credible_set(results$top_credible_set[i])
    if (length(genes) == 0) next
    base <- results[i, , drop = FALSE]
    for (g in genes) {
      r <- base
      r$nominated_mouse_gene <- g
      r$nominated_query_gene <- unname(target_to_source[g])
      rows[[length(rows) + 1]] <- r
    }
  }
  if (length(rows) == 0) return(results[0, , drop = FALSE])
  long <- do.call(rbind, rows)

  if (!is.null(edger_long) && celltype_join_col %in% names(long)) {
    keep <- c("gene", "logFC", "FDR", "celltype_file")
    de <- edger_long[, intersect(keep, names(edger_long)), drop = FALSE]
    names(de)[names(de) == "gene"] <- "nominated_query_gene"
    names(de)[names(de) == "logFC"] <- "de_logFC"
    names(de)[names(de) == "FDR"]   <- "de_FDR"
    long <- merge(long, de,
                  by.x = c("nominated_query_gene", celltype_join_col),
                  by.y = c("nominated_query_gene", "celltype_file"),
                  all.x = TRUE, sort = FALSE)
  }
  long
}
