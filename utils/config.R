# Config loading for the RNA fingerprinting pipeline (R side).
#
# Reads the same config.conf as utils/config.py so Python and R share one config.
# Dependency-free base-R INI parser: sections -> named list of key=value strings.
# Scripts in pipeline/ source this instead of hard-coding paths or tunables.

# Locate config.conf one level up from this file (<repo>/utils/config.R -> <repo>/config.conf).
.config_default_path <- function() {
  # When sourced, sys.frame(1)$ofile or the sourced path isn't always available;
  # fall back to a RNAFP_CONFIG env var, then the repo-root guess.
  env <- Sys.getenv("RNAFP_CONFIG", unset = "")
  if (nzchar(env)) return(normalizePath(env, mustWork = TRUE))
  this <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
  if (!is.na(this)) return(file.path(dirname(dirname(this)), "config.conf"))
  "config.conf"
}

#' Load config.conf into a nested named list: cfg[[section]][[key]] = value (string).
load_config <- function(path = NULL) {
  if (is.null(path)) path <- .config_default_path()
  if (!file.exists(path)) stop(sprintf("config not found: %s", path))
  lines <- readLines(path, warn = FALSE)
  cfg <- list()
  section <- NULL
  for (ln in lines) {
    s <- trimws(ln)
    if (s == "" || startsWith(s, "#") || startsWith(s, ";")) next
    if (startsWith(s, "[") && endsWith(s, "]")) {
      section <- substr(s, 2, nchar(s) - 1)
      if (is.null(cfg[[section]])) cfg[[section]] <- list()
      next
    }
    eq <- regexpr("=", s, fixed = TRUE)
    if (eq < 0 || is.null(section)) next
    key <- trimws(substr(s, 1, eq - 1))
    val <- trimws(substr(s, eq + 1, nchar(s)))
    cfg[[section]][[key]] <- val
  }
  cfg
}

#' Fetch a value; error if required and blank/missing.
cfg_get <- function(cfg, section, key, required = TRUE, default = NULL) {
  val <- tryCatch(cfg[[section]][[key]], error = function(e) NULL)
  if (is.null(val) || !nzchar(val)) {
    if (required) stop(sprintf("[%s] %s is blank in config.conf — fill it in.", section, key))
    return(default)
  }
  val
}

#' Typed getters.
cfg_int  <- function(cfg, section, key, required = TRUE, default = NULL)
  { v <- cfg_get(cfg, section, key, required, default); if (is.null(v)) NULL else as.integer(v) }
cfg_num  <- function(cfg, section, key, required = TRUE, default = NULL)
  { v <- cfg_get(cfg, section, key, required, default); if (is.null(v)) NULL else as.numeric(v) }
cfg_bool <- function(cfg, section, key, required = TRUE, default = NULL)
  { v <- cfg_get(cfg, section, key, required, default); if (is.null(v)) NULL else tolower(v) %in% c("true","1","yes") }

#' Active query section name, e.g. "query.rat". Defaults to [run].cohort.
query_section <- function(cfg, cohort = NULL) {
  if (is.null(cohort)) cohort <- cfg_get(cfg, "run", "cohort")
  section <- paste0("query.", cohort)
  if (is.null(cfg[[section]])) stop(sprintf("no [%s] section in config.conf", section))
  section
}
