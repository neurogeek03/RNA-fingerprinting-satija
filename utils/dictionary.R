# Dictionary-building wrappers around the RNAFingerprinting package.
#
# Thin, parameterised wrappers so pipeline/ scripts don't call the package
# directly with magic strings. All tunables arrive as arguments (from config).

suppressPackageStartupMessages({
  library(Seurat)
  library(RNAFingerprinting)
})

#' Learn a perturbation dictionary from a (subsetted, SCT-normalized) Perturb-seq
#' Seurat object. `feature_method` is passed to ChooseFeatures ('de' or 'expression').
#'
#' Returns the learned dictionary object (cache it with io.R::save_dictionary).
learn_dictionary <- function(obj, perturbation_col, control_class,
                             feature_method = "de", min_cells = 10,
                             mc.cores = 1, seed.use = 1448145) {
  features <- ChooseFeatures(obj, method = feature_method,
                             perturbation_meta = perturbation_col,
                             control_class = control_class)
  # Keep only features that exist in the SCT scale.data (fingerprinting needs them).
  sct <- GetAssayData(obj, assay = "SCT", layer = "scale.data")
  features <- intersect(features, rownames(sct))
  if (length(features) == 0) stop("no usable features after SCT intersection")
  LearnDictionary(obj, perturbation_meta = perturbation_col,
                  control_class = control_class, features = features,
                  min_cells = min_cells, mc.cores = mc.cores, seed.use = seed.use)
}

#' Load the precomputed genome-wide Perturb-seq dictionary (RNAFingerprintingData).
#' `InstallPrecomputedDictionary` needs internet -> run on the login node first.
load_gwps_dictionary <- function(name = "gwps", gene_key = "SYMBOL",
                                 install = FALSE) {
  if (install) {
    suppressPackageStartupMessages(library(RNAFingerprintingData))
    InstallPrecomputedDictionary(name = name)
  }
  LoadPrecomputedDictionary(name = name, gene_key = gene_key)
}

#' Report overlap between a dictionary's features and a query's genes (rownames).
#' Fails loudly (if `min_overlap` set) when cross-species mapping left too few.
dictionary_feature_overlap <- function(dictionary, query_genes, min_overlap = NULL) {
  feats <- rownames(dictionary)
  if (is.null(feats)) feats <- dictionary$features  # fallback if list-shaped
  ov <- intersect(feats, query_genes)
  if (!is.null(min_overlap) && length(ov) < min_overlap)
    stop(sprintf("only %d dictionary features overlap the query (min %d) — check orthology",
                 length(ov), min_overlap))
  list(n_features = length(feats), n_query = length(query_genes), n_overlap = length(ov))
}
