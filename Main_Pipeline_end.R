PROJECT_ROOT <- normalizePath(
  Sys.getenv("MAPPING_FLAVONOID_ROOT", unset = "."),
  winslash = "/",
  mustWork = TRUE
)

setwd(PROJECT_ROOT)

SCRIPTS_DIR <- file.path(PROJECT_ROOT, "scripts")
OUTPUTS_DIR <- file.path(PROJECT_ROOT, "outputs")
DATA_DIR <- file.path(PROJECT_ROOT, "data")
INPUTS_DIR <- file.path(PROJECT_ROOT, "inputs")

dir.create(OUTPUTS_DIR, recursive = TRUE, showWarnings = FALSE)

# Warn about the accidental outputs folder previously created inside scripts.
ACCIDENTAL_OUTPUTS_DIR <- file.path(SCRIPTS_DIR, "outputs")
if (dir.exists(ACCIDENTAL_OUTPUTS_DIR)) {
  warning(
    "An old accidental outputs folder exists inside scripts: ",
    ACCIDENTAL_OUTPUTS_DIR,
    "\nThe validated pipeline writes to: ", OUTPUTS_DIR
  )
}

required_packages <- c("arrow", "dplyr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages)) {
  stop(
    "Missing required packages: ", paste(missing_packages, collapse = ", "),
    "\nInstall them before running the pipeline."
  )
}

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
})

`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0L) b else a
}

ensure_runtime_data_frames <- function() {
  for (object_name in c("lin_enriched", "uni_enriched")) {
    if (exists(object_name, envir = .GlobalEnv, inherits = FALSE)) {
      assign(
        object_name,
        as.data.frame(
          get(object_name, envir = .GlobalEnv, inherits = FALSE),
          stringsAsFactors = FALSE
        ),
        envir = .GlobalEnv
      )
    }
  }
  invisible(NULL)
}

resolve_script <- function(candidates, label) {
  paths <- file.path(SCRIPTS_DIR, candidates)
  hit <- paths[file.exists(paths)]
  if (!length(hit)) {
    stop(
      label, " script not found. Tried:\n- ",
      paste(paths, collapse = "\n- ")
    )
  }
  normalizePath(hit[1], winslash = "/", mustWork = TRUE)
}

safe_tag_fn <- function(mode, values, run_date, suffix = NULL) {
  tag <- paste(
    mode,
    paste(values, collapse = "-"),
    format(run_date, "%Y%m%d"),
    sep = "_"
  )
  if (!is.null(suffix) && nzchar(suffix)) {
    tag <- paste0(tag, "_", gsub("[^A-Za-z0-9._-]+", "_", suffix))
  }
  gsub("[^A-Za-z0-9._-]+", "_", tag)
}

cfg_mode <- "kingdom"
cfg_values <- "Plantae"
cfg_level <- "family"

# Select the modules to run.
# For a clean restart, begin with Part I only.

cfg <- list(
  project_root = PROJECT_ROOT,
  taxon_mode = tolower(cfg_mode),
  taxon_values = cfg_values,
  analysis_tax_level = tolower(cfg_level),

  run_module1 = FALSE,    # Part I: extraction, curation, and normalization
  run_module2 = FALSE,   # Part II: biological metadata integration
  run_module3 = TRUE,   # Part III: statistics and figures

  mongo_url = Sys.getenv("LOTUS_MONGO_URL"),
  db_name = "lotusdb",
  coll_name = "subset_minor_flavonoids_V",
  mongo_opts = '{"allowDiskUse": true, "batchSize": 5000}',
  mongo_batch_size = 2000L,

  use_WFO_normalization = TRUE,
  wfo_csv_path = file.path(DATA_DIR, "wfo", "classification.tsv"),

  flav_classes = c(
    "Flavanones", "Isoflavones", "Flavan-3-ols", "Chalcones",
    "Dihydroflavonols", "Pterocarpan", "Proanthocyanins",
    "Isoflavanones", "2-arylbenzofurans", "Flavandiols",
    "Rotenoids", "Aurones", "Flavans", "Coumestan",
    "Anthocyanidins", "Neoflavonoids", "Flavonolignans"
  ),

  scaffold_field = "murko_framework",
  scaffold_definition = "Bemis-Murcko",

  activity_cutoff_nm = 10000,
  chembl_valid_types = c("IC50", "Ki", "EC50", "Kd", "MIC", "AC50", "GI50"),
  chembl_exact_relations = "=",
  chembl_use_unichem_fallback = FALSE,

  analysis_top_taxa = 200L,
  # Validated primary family threshold for Part III.
  analysis_min_compounds_per_taxon = 10L,
  analysis_min_species_per_taxon = 2L,
  # Sensitivity analyses around the validated primary threshold.
  analysis_compound_thresholds = c(5L, 10L, 20L),
  run_sensitivity_analysis = TRUE,

  chembl_min_records_per_taxon = 20L,
  chembl_min_compounds_per_taxon = 5L,
  chembl_min_targets_per_taxon = 3L,

  out_dir_base = OUTPUTS_DIR,
  run_tag_date = as.Date("2026-06-22"),
  prefix_base_tag = NULL,
  custom_tag_suffix = NULL,

  export_excel = TRUE,
  export_parquet = TRUE,
  auto_install_missing_packages = FALSE,
  auto_install_pkgs = FALSE,
  verbose = TRUE
)

if (!nzchar(cfg$mongo_url)) {
  warning(
    "LOTUS_MONGO_URL was not detected. Using the local MongoDB fallback."
  )
  cfg$mongo_url <- paste0(
    "mongodb://127.0.0.1:27017/",
    "?socketTimeoutMS=3600000&connectTimeoutMS=300000"
  )
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

load_dir <- file.path(OUTPUTS_DIR, paste0("lotus_", tag_base_load))
base_tag <- paste0("lotus_", tag_base_load)

cfg$path_lotus <- file.path(load_dir, paste0(base_tag, ".xlsx"))
cfg$OUT_DIR <- load_dir
cfg$base_tag <- base_tag

runtime <- list(
  OUT_DIR = load_dir,
  base_tag = base_tag,
  TAXON_MODE = cfg$taxon_mode,
  TAXON_VALUES = cfg$taxon_values
)

cat("------------------------------------------------------------\n")
cat("Project root: ", PROJECT_ROOT, "\n", sep = "")
cat("Scripts:      ", SCRIPTS_DIR, "\n", sep = "")
cat("Output:       ", load_dir, "\n", sep = "")
cat("------------------------------------------------------------\n")
cat("Expected standardized module names:\n")
cat("  Part I:   Part_I_Extraction_v2_1.R\n")
cat("  Part II:  Part_II_Bio_iTOL_Prep_v2_1.R\n")
cat("  Part III: Part_III_Figures_Stats_v2_1.R\n")
cat("------------------------------------------------------------\n")

assign("cfg", cfg, envir = .GlobalEnv)
assign("runtime", runtime, envir = .GlobalEnv)
assign("OUT_DIR", runtime$OUT_DIR, envir = .GlobalEnv)
assign("base_tag", runtime$base_tag, envir = .GlobalEnv)

cat("Modules:      Part I=", cfg$run_module1,
    " | Part II=", cfg$run_module2,
    " | Part III=", cfg$run_module3, "\n", sep = "")

if (isTRUE(cfg$run_module1)) {
  required_part1 <- c(
    "here", "mongolite", "jsonlite", "dplyr", "progress",
    "stringr", "writexl", "readr", "stringi", "arrow"
  )
  missing_part1 <- required_part1[
    !vapply(required_part1, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
  ]
  if (length(missing_part1)) {
    stop(
      "Missing packages required by Part I: ", paste(missing_part1, collapse = ", "),
      "\nInstall them before rerunning."
    )
  }
  if (isTRUE(cfg$use_WFO_normalization) && !file.exists(cfg$wfo_csv_path)) {
    stop("WFO classification file not found at: ", cfg$wfo_csv_path)
  }
  cat("\n[1/3] Executing Part I: extraction, curation, and normalization...\n")
  script_path <- resolve_script(
    c("Part_I_Extraction_v2_1.R", "Part I_Extraction.R"),
    "Part I"
  )
  source(script_path, local = .GlobalEnv, chdir = FALSE)
  ensure_runtime_data_frames()
  cat("[OK] Part I completed successfully.\n")
}

need_data <- isTRUE(cfg$run_module2) || isTRUE(cfg$run_module3)
data_in_memory <-
  exists("lin_enriched", envir = .GlobalEnv, inherits = FALSE) &&
  exists("uni_enriched", envir = .GlobalEnv, inherits = FALSE)

if (need_data && data_in_memory) {
  ensure_runtime_data_frames()
}

if (need_data && !data_in_memory) {
  cat("\n[Auto-load] Searching for validated Part I Parquet datasets...\n")

  pq_lin <- file.path(load_dir, paste0(base_tag, "_lin_enriched.parquet"))
  pq_uni <- file.path(load_dir, paste0(base_tag, "_uni_enriched.parquet"))

  if (!file.exists(pq_lin) || !file.exists(pq_uni)) {
    stop(
      "Required analytical datasets were not found in:\n", load_dir,
      "\nRun Part I first."
    )
  }

  lin_enriched <<- as.data.frame(
    arrow::read_parquet(pq_lin, as_data_frame = TRUE),
    stringsAsFactors = FALSE
  )
  uni_enriched <<- as.data.frame(
    arrow::read_parquet(pq_uni, as_data_frame = TRUE),
    stringsAsFactors = FALSE
  )
  cat("[OK] Part I datasets loaded.\n")
}

if (isTRUE(cfg$run_module2)) {
  cat("\n[2/3] Executing Part II: biological metadata integration and RDKit preparation...\n")
  script_path <- resolve_script(
    c("Part_II_Bio_iTOL_Prep_v2_1.R", "Part II - Bio_iTOL_Prep.R"),
    "Part II"
  )
  ensure_runtime_data_frames()
  source(script_path, local = .GlobalEnv, chdir = FALSE)
  cat("[OK] Part II completed successfully.\n")
  cat("[NEXT] Run the revised RDKit Python script before Part III.\n")
}

if (isTRUE(cfg$run_module3)) {
  # Enforce the validated manuscript settings before Part III.
  cfg$analysis_min_compounds_per_taxon <- 10L
  cfg$analysis_compound_thresholds <- c(5L, 10L, 20L)

  stopifnot(
    identical(cfg$analysis_min_compounds_per_taxon, 10L),
    identical(cfg$analysis_compound_thresholds, c(5L, 10L, 20L))
  )

  assign("cfg", cfg, envir = .GlobalEnv)

  cat(
    "\n[3/3] Executing Part III: statistical analyses and figure generation...\n",
    "Primary family threshold: n >= ",
    cfg$analysis_min_compounds_per_taxon,
    "\nSensitivity thresholds: ",
    paste(cfg$analysis_compound_thresholds, collapse = ", "),
    "\n",
    sep = ""
  )
  script_path <- resolve_script(
    c("Part_III_Figures_Stats_v2_1.R", "Part III - Figures_Stats.R", "Part_III_Figures_Stats_v3_1.R"),
    "Part III"
  )
  ensure_runtime_data_frames()
  source(script_path, local = .GlobalEnv, chdir = FALSE)
  cat("[OK] Part III completed successfully.\n")
}
