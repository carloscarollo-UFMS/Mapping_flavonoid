# ======================================================================
# LOTUS PIPELINE FOR THE ANALYSIS OF MINOR FLAVONOIDS
# Integrated workflow:
#   Part I   - Extraction, curation, and normalization
#   Part II  - Biological metadata integration and iTOL preparation
#   Part III - Statistical analyses and figure generation
# ======================================================================

# 0. ENVIRONMENT SETUP AND DEPENDENCIES
# ----------------------------------------------------------------------
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
if (!requireNamespace("arrow", quietly = TRUE)) install.packages("arrow")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(here)
library(dplyr)
library(arrow)

# Helper for nullable default values.
`%||%` <- function(a, b) if (is.null(a)) b else a

# Standardize analytical objects as base data.frame instances before
# downstream processing.
ensure_runtime_data_frames <- function() {
  target_objects <- c("lin_enriched", "uni_enriched")

  for (object_name in target_objects) {
    if (exists(object_name, envir = .GlobalEnv, inherits = FALSE)) {
      object_value <- get(object_name, envir = .GlobalEnv)
      assign(
        object_name,
        as.data.frame(object_value, stringsAsFactors = FALSE),
        envir = .GlobalEnv
      )
    }
  }

  invisible(NULL)
}

# Ensure the output directory exists.
if (!dir.exists(here("outputs"))) {
  dir.create(here("outputs"), recursive = TRUE)
}

cat("=====================================================\n")
cat(" LOTUS ANALYTICAL PIPELINE INITIALIZED\n")
cat(" Execution mode: portable\n")
cat(" Project root: ", here(), "\n", sep = "")
cat("=====================================================\n")

# ----------------------------------------------------------------------
# 1. CONFIGURATION
# ----------------------------------------------------------------------

# Analysis scenario: taxonomic subset.
cfg_mode <- "kingdom"
cfg_values <- "Plantae"
cfg_level <- "family"

cfg <- list(
  ## --- TARGET TAXONOMIC SCOPE -------------------------
  taxon_mode = tolower(cfg_mode),
  taxon_values = cfg_values,
  analysis_tax_level = tolower(cfg_level),

  ## --- MODULE EXECUTION SWITCHES ----------------------
  run_module1 = FALSE,   # Part I: MongoDB extraction, curation, and normalization
  run_module2 = FALSE,   # Part II: biological metadata and iTOL preparation
  run_module3 = TRUE,    # Part III: statistical analyses and figure generation

  ## --- DATABASE CONNECTION PARAMETERS ----------------
  # The connection string must be defined in .Renviron as LOTUS_MONGO_URL.
  mongo_url = Sys.getenv("LOTUS_MONGO_URL"),
  db_name = "lotusdb",
  coll_name = "subset_minor_flavonoids_V",
  mongo_opts = '{"allowDiskUse": true, "batchSize": 5000}',

  ## --- ANALYTICAL FILTERS AND DISPLAY PARAMETERS -----
  analysis_top_taxa = 200L,
  analysis_min_compounds_per_taxon = 5L,

  ## --- OUTPUT AND IDENTIFICATION SETTINGS ------------
  out_dir_base = here("outputs"),
  run_tag_date = Sys.Date(),
  prefix_base_tag = NULL,
  custom_tag_suffix = NULL,

  ## --- ADDITIONAL OPTIONS ----------------------------
  use_WFO_normalization = TRUE,
  export_excel = TRUE,
  export_parquet = TRUE,
  verbose = TRUE
)

# Validate database connectivity settings.
if (cfg$mongo_url == "") {
  warning(
    paste(
      "'LOTUS_MONGO_URL' was not detected in .Renviron.",
      "A localhost fallback will be used for the current session.",
      "Define this variable to ensure reproducibility and portability."
    )
  )
  cfg$mongo_url <- "mongodb://127.0.0.1:27017/?socketTimeoutMS=3600000&connectTimeoutMS=300000"
}

# ----------------------------------------------------------------------
# 2. PATH DEFINITION AND RUN IDENTIFICATION
# ----------------------------------------------------------------------

cat(" Database name: ", cfg$db_name, "\n", sep = "")
cat(" Collection name: ", cfg$coll_name, "\n", sep = "")
cat(" Taxonomic level for analysis: ", cfg$analysis_tax_level, "\n", sep = "")
cat("-----------------------------------------------------\n")

# Generate a stable identifier for files produced during the analysis.
safe_tag_fn <- function(mode, values, run_date, suffix = NULL) {
  tag <- paste(mode, paste(values, collapse = "-"), format(run_date, "%Y%m%d"), sep = "_")

  if (!is.null(suffix) && nzchar(suffix)) {
    tag <- paste0(tag, "_", gsub("[^A-Za-z0-9._-]+", "_", suffix))
  }

  gsub("[^A-Za-z0-9._-]+", "_", tag)
}

tag_base_load <- if (!is.null(cfg$prefix_base_tag) && nzchar(cfg$prefix_base_tag)) {
  cfg$prefix_base_tag
} else {
  safe_tag_fn(
    cfg$taxon_mode,
    cfg$taxon_values,
    cfg$run_tag_date,
    cfg$custom_tag_suffix %||% NULL
  )
}

load_dir <- here("outputs", paste0("lotus_", tag_base_load))

# Register global runtime variables shared across pipeline modules.
runtime <- list(
  OUT_DIR = load_dir,
  base_tag = paste0("lotus_", tag_base_load),
  TAXON_MODE = cfg$taxon_mode,
  TAXON_VALUES = cfg$taxon_values
)

# ----------------------------------------------------------------------
# 3. PIPELINE EXECUTION
# ----------------------------------------------------------------------

# --- PART I: EXTRACTION, CURATION, AND NORMALIZATION ---
if (isTRUE(cfg$run_module1)) {
  cat("\n[1/3] Executing Part I: extraction, curation, and normalization...\n")

  script_path <- here("scripts", "Part I_Extraction.R")

  if (!file.exists(script_path)) {
    stop("Part I script not found at: ", script_path)
  }

  assign("OUT_DIR", runtime$OUT_DIR, envir = .GlobalEnv)
  assign("base_tag", runtime$base_tag, envir = .GlobalEnv)
  assign("cfg", cfg, envir = .GlobalEnv)

  source(script_path)
  ensure_runtime_data_frames()
  cat("[OK] Part I completed successfully.\n")
}

# --- AUTOMATED DATA LOADING FOR DOWNSTREAM MODULES ---
# If Part I is skipped and any downstream module is active, the required
# datasets must be retrieved from the serialized outputs.
need_data <- (cfg$run_module2 || cfg$run_module3)
data_in_memory <- (
  exists("lin_enriched", envir = .GlobalEnv, inherits = FALSE) &&
    exists("uni_enriched", envir = .GlobalEnv, inherits = FALSE)
)

if (need_data && data_in_memory) {
  ensure_runtime_data_frames()
}

if (need_data && !data_in_memory) {
  cat("\n[Auto-load] Downstream modules detected. Searching for Parquet datasets...\n")

  pq_lin <- file.path(load_dir, paste0("lotus_", tag_base_load, "_lin_enriched.parquet"))
  pq_uni <- file.path(load_dir, paste0("lotus_", tag_base_load, "_uni_enriched.parquet"))

  if (file.exists(pq_lin) && file.exists(pq_uni)) {
    cat(" Loading dataset: ", basename(pq_lin), "\n", sep = "")
    lin_enriched <<- as.data.frame(
      arrow::read_parquet(pq_lin, as_data_frame = TRUE),
      stringsAsFactors = FALSE
    )
    uni_enriched <<- as.data.frame(
      arrow::read_parquet(pq_uni, as_data_frame = TRUE),
      stringsAsFactors = FALSE
    )
    cat("[OK] Analytical datasets loaded successfully as data.frame objects.\n")
  } else {
    stop(
      paste0(
        "Required analytical datasets were not found in:\n",
        load_dir,
        "\nSet 'run_module1 = TRUE' to generate the input files before executing downstream modules."
      )
    )
  }
}

# --- PART II: BIOLOGICAL METADATA AND iTOL PREPARATION ---
if (isTRUE(cfg$run_module2)) {
  cat("\n[2/3] Executing Part II: biological metadata integration and iTOL preparation...\n")

  script_path <- here("scripts", "Part II - Bio_iTOL_Prep.R")
  if (!file.exists(script_path)) {
    stop("Part II script not found at: ", script_path)
  }

  assign("OUT_DIR", runtime$OUT_DIR, envir = .GlobalEnv)
  assign("base_tag", runtime$base_tag, envir = .GlobalEnv)
  assign("cfg", cfg, envir = .GlobalEnv)

  ensure_runtime_data_frames()
  source(script_path)
  cat("[OK] Part II completed successfully.\n")
}

# --- PART III: STATISTICAL ANALYSES AND FIGURE GENERATION ---
if (isTRUE(cfg$run_module3)) {
  cat("\n[3/3] Executing Part III: statistical analyses and figure generation...\n")

  script_path <- here("scripts", "Part III - Figures_Stats.R")
  if (!file.exists(script_path)) {
    stop("Part III script not found at: ", script_path)
  }

  assign("OUT_DIR", runtime$OUT_DIR, envir = .GlobalEnv)
  assign("base_tag", runtime$base_tag, envir = .GlobalEnv)
  assign("cfg", cfg, envir = .GlobalEnv)

  ensure_runtime_data_frames()
  source(script_path)
  cat("[OK] Part III completed successfully.\n")
}
