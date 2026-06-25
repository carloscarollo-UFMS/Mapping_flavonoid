# Part_III v2.13 FINAL QC: corrected clade-domain quartiles and coverage-filtered main network
`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0L) return(b)
  if (is.character(a) && length(a) == 1L && (is.na(a) || !nzchar(a))) return(b)
  a
}

as_base_df <- function(x) as.data.frame(x, stringsAsFactors = FALSE)
trim_chr <- function(x) {
  x <- trimws(as.character(x))
  x[!nzchar(x)] <- NA_character_
  x
}

normalize_path_safe <- function(x) normalizePath(x, winslash = "/", mustWork = FALSE)

normalize_name_vector <- function(n) {
  n <- iconv(as.character(n), from = "", to = "ASCII//TRANSLIT")
  n <- tolower(gsub("[^A-Za-z0-9]+", "_", n))
  n <- gsub("^_+|_+$", "", n)
  make.unique(n, sep = "_")
}

clean_names_simple <- function(x) {
  names(x) <- normalize_name_vector(names(x))
  x
}

find_first_existing <- function(paths) {
  paths <- unique(as.character(paths))
  paths <- paths[!is.na(paths) & nzchar(paths)]
  hit <- paths[file.exists(paths)]
  if (length(hit)) normalize_path_safe(hit[1]) else NA_character_
}

resolve_file <- function(label, candidates, required = TRUE) {
  hit <- find_first_existing(candidates)
  if (is.na(hit) && isTRUE(required)) {
    stop(
      "Required file not found [", label, "]. Tried:\n- ",
      paste(candidates, collapse = "\n- ")
    )
  }
  hit
}

require_pkgs <- function(pkgs, auto_install = FALSE) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (!length(missing)) return(invisible(TRUE))
  if (isTRUE(auto_install)) {
    install.packages(missing, dependencies = TRUE)
    missing <- missing[!vapply(missing, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  }
  if (length(missing)) {
    stop(
      "Missing required packages: ", paste(missing, collapse = ", "), "\n",
      "Install them with:\ninstall.packages(c(",
      paste(sprintf("\"%s\"", missing), collapse = ", "),
      "), dependencies = TRUE)"
    )
  }
  invisible(TRUE)
}

safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

safe_logical <- function(x) {
  if (is.logical(x)) return(x)
  y <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(y))
  out[y %in% c("true", "t", "1", "yes", "y")] <- TRUE
  out[y %in% c("false", "f", "0", "no", "n")] <- FALSE
  out
}

is_binomial_name <- function(x) {
  x <- trim_chr(x)
  !is.na(x) & grepl("^[A-Z][A-Za-z-]+\\s+[a-z][A-Za-z-]+", x)
}

first_nonempty <- function(x) {
  x <- trim_chr(x)
  x <- x[!is.na(x)]
  if (length(x)) x[1] else NA_character_
}

collapse_unique <- function(x, sep = ";") {
  x <- sort(unique(trim_chr(x)))
  x <- x[!is.na(x)]
  if (length(x)) paste(x, collapse = sep) else NA_character_
}

resolve_column <- function(df, candidates, label, required = TRUE) {
  df_names <- names(df)
  candidates <- normalize_name_vector(candidates)


  exact_hit <- candidates[candidates %in% df_names]
  if (length(exact_hit)) return(exact_hit[1])


  compact_df <- gsub("_", "", df_names, fixed = TRUE)
  compact_candidates <- unique(gsub("_", "", candidates, fixed = TRUE))

  for (candidate_compact in compact_candidates) {
    idx <- which(compact_df == candidate_compact)
    if (length(idx) == 1L) {
      message(
        "[INFO] Resolved column [", label, "] by compact-name matching: ",
        df_names[idx]
      )
      return(df_names[idx])
    }
    if (length(idx) > 1L) {
      stop(
        "Ambiguous compact-name match [", label, "]: ",
        paste(df_names[idx], collapse = ", ")
      )
    }
  }

  if (isTRUE(required)) {
    stop(
      "Column not found [", label, "]. Candidates: ",
      paste(candidates, collapse = ", "),
      "
Available columns: ", paste(df_names, collapse = ", ")
    )
  }
  NA_character_
}

read_csv_character <- function(path, label = basename(path), audit_dir = NULL) {
  dat <- readr::read_csv(
    path,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE,
    progress = FALSE,
    trim_ws = TRUE
  )
  probs <- readr::problems(dat)
  if (nrow(probs) > 0L) {
    message(
      "[WARN] Parsing issues detected in ", label, ": ", nrow(probs),
      " row(s). Values were imported as character; see audit file."
    )
    if (!is.null(audit_dir) && dir.exists(audit_dir)) {
      safe_label <- gsub("[^A-Za-z0-9]+", "_", label)
      readr::write_csv(
        as.data.frame(probs, stringsAsFactors = FALSE),
        file.path(audit_dir, paste0("parsing_problems_", safe_label, ".csv")),
        na = ""
      )
    }
  }
  as_base_df(dat)
}

write_csv_safe <- function(df, path) {
  readr::write_csv(as_base_df(df), path, na = "")
  invisible(path)
}

write_json_safe <- function(x, path) {
  jsonlite::write_json(x, path, pretty = TRUE, auto_unbox = TRUE, na = "null")
  invisible(path)
}

scale_safe <- function(x) {
  x <- safe_numeric(x)
  if (sum(is.finite(x)) < 2L || stats::sd(x, na.rm = TRUE) == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}

wilson_ci <- function(k, n, z = 1.96) {
  if (!is.finite(k) || !is.finite(n) || n <= 0) return(c(NA_real_, NA_real_))
  p <- k / n
  den <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / den
  half <- z * sqrt((p * (1 - p) + z^2 / (4 * n)) / n) / den
  c(max(0, center - half), min(1, center + half))
}

weighted_mean_safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  stats::weighted.mean(x[ok], w[ok])
}

mean_finite_safe <- function(x) {
  x <- safe_numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

median_finite_safe <- function(x) {
  x <- safe_numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  stats::median(x)
}

max_finite_safe <- function(x) {
  x <- safe_numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  max(x)
}

pseudo_log10 <- scales::pseudo_log_trans(base = 10)


cfg <- if (exists("cfg", inherits = TRUE)) get("cfg", inherits = TRUE) else list()
runtime <- if (exists("runtime", inherits = TRUE)) get("runtime", inherits = TRUE) else list()

message(">>> PART III v2.13 FINAL QC: n >= 10, corrected quartiles, and Fig. 5 main network n >= 5")

OUT_DIR <- runtime$OUT_DIR %||% cfg$OUT_DIR %||% NA_character_
if (is.na(OUT_DIR) || !nzchar(OUT_DIR) || !dir.exists(OUT_DIR)) {
  stop("OUT_DIR was not found. Run this script through Main_Pipeline_end_v2_1.R or define cfg$OUT_DIR.")
}
OUT_DIR <- normalize_path_safe(OUT_DIR)
base_tag <- runtime$base_tag %||% cfg$base_tag %||% basename(OUT_DIR)
project_root <- normalize_path_safe(cfg$project_root %||% getwd())

cfg$auto_install_pkgs <- isTRUE(cfg$auto_install_pkgs %||% FALSE)
cfg$analysis_min_compounds_per_taxon <- as.integer(cfg$analysis_min_compounds_per_taxon %||% 10L)
if (length(cfg$analysis_min_compounds_per_taxon) != 1L || is.na(cfg$analysis_min_compounds_per_taxon)) {
  stop("cfg$analysis_min_compounds_per_taxon must be one integer value.")
}
if (cfg$analysis_min_compounds_per_taxon != 10L) {
  stop(
    "The validated manuscript analysis requires ",
    "cfg$analysis_min_compounds_per_taxon = 10L; received ",
    cfg$analysis_min_compounds_per_taxon, "."
  )
}
cfg$analysis_compound_thresholds <- as.integer(cfg$analysis_compound_thresholds %||% c(5L, 10L, 20L))
cfg$analysis_compound_thresholds <- sort(unique(cfg$analysis_compound_thresholds[cfg$analysis_compound_thresholds > 0]))
if (!10L %in% cfg$analysis_compound_thresholds) {
  cfg$analysis_compound_thresholds <- sort(unique(c(cfg$analysis_compound_thresholds, 10L)))
}
cfg$permutations <- as.integer(cfg$permutations %||% 999L)
cfg$seed <- as.integer(cfg$seed %||% 20260622L)
cfg$functional_min_compounds_domain <- as.integer(cfg$functional_min_compounds_domain %||% 10L)
cfg$functional_min_records_domain <- as.integer(cfg$functional_min_records_domain %||% 20L)
cfg$motif_model_min_compounds <- as.integer(cfg$motif_model_min_compounds %||% 60L)
cfg$fsp3_gam_min_compounds <- as.integer(cfg$fsp3_gam_min_compounds %||% 50L)
cfg$functional_primary_type_aggregation <- as.character(
  cfg$functional_primary_type_aggregation %||% "equal_mean"
)
cfg$functional_binary_primary_rule <- as.character(
  cfg$functional_binary_primary_rule %||% "concordant_only"
)
cfg$functional_run_record_weight_sensitivity <- isTRUE(
  cfg$functional_run_record_weight_sensitivity %||% TRUE
)
if (!cfg$functional_primary_type_aggregation %in% c("equal_mean", "equal_median")) {
  stop(
    "cfg$functional_primary_type_aggregation must be 'equal_mean' or 'equal_median'."
  )
}
if (!identical(cfg$functional_binary_primary_rule, "concordant_only")) {
  stop(
    "This validated script requires cfg$functional_binary_primary_rule = 'concordant_only'."
  )
}
set.seed(cfg$seed)

pkgs <- c(
  "arrow", "readr", "dplyr", "tidyr", "stringr", "purrr", "tibble",
  "ggplot2", "ggrepel", "patchwork", "scales", "vegan", "mgcv",
  "rstatix", "multcompView", "openxlsx", "broom", "sandwich", "lmtest",
  "jsonlite", "ComplexHeatmap", "circlize", "ggridges", "igraph",
  "tidygraph", "ggraph", "gridExtra"
)
require_pkgs(pkgs, auto_install = cfg$auto_install_pkgs)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})


PARTIII_DIR <- file.path(OUT_DIR, "PartIII_ALL")
TABLE_DIR <- file.path(PARTIII_DIR, "tables")
FIG_DIR <- file.path(PARTIII_DIR, "figures")
ITOL_DIR <- file.path(PARTIII_DIR, "itol")
AUDIT_DIR <- file.path(PARTIII_DIR, "audit")
MODEL_DIR <- file.path(PARTIII_DIR, "models")
for (d in c(PARTIII_DIR, TABLE_DIR, FIG_DIR, ITOL_DIR, AUDIT_DIR, MODEL_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

message("----------------------------------------------------------------")
message(">>> PART III v4.1.2: final tree inputs, n >= 10 primary threshold, and functional-bias controls")
message("----------------------------------------------------------------")

part2_dir <- file.path(OUT_DIR, "PartII_ALL")
rdkit_dir <- file.path(OUT_DIR, "RDKit_ALL")

# The validated family tree may be stored outside OUT_DIR. Priority is:
# 1) cfg$tree_output_dir, 2) project_root/phylo_outputs_FINAL,
# 3) the legacy OUT_DIR/Tree_APG_iTOL directory.
tree_dir_candidates <- unique(c(
  cfg$tree_output_dir %||% NA_character_,
  file.path(project_root, "phylo_outputs_FINAL"),
  file.path(OUT_DIR, "Tree_APG_iTOL")
))
tree_dir_candidates <- tree_dir_candidates[
  !is.na(tree_dir_candidates) & nzchar(tree_dir_candidates) & dir.exists(tree_dir_candidates)
]
if (!length(tree_dir_candidates)) {
  stop(
    "No tree output directory was found. Set cfg$tree_output_dir explicitly."
  )
}
tree_dir <- normalize_path_safe(tree_dir_candidates[[1]])
message("[INPUT] selected tree output directory: ", tree_dir)

path_lin <- resolve_file(
  "Part I lin_enriched parquet",
  c(
    cfg$path_lin_enriched %||% NA_character_,
    file.path(OUT_DIR, paste0(base_tag, "_lin_enriched.parquet")),
    Sys.glob(file.path(OUT_DIR, "*_lin_enriched.parquet"))
  ),
  required = !exists("lin_enriched", inherits = TRUE)
)

path_rdkit_compounds <- resolve_file(
  "RDKit compound table",
  c(
    cfg$path_rdkit_compounds %||% NA_character_,
    file.path(rdkit_dir, "lotus_flavonoids_rdkit_compounds.csv"),
    file.path(OUT_DIR, "lotus_flavonoids_rdkit_compounds.csv")
  )
)
path_rdkit_qc <- resolve_file(
  "RDKit structural QC",
  c(
    cfg$path_rdkit_qc %||% NA_character_,
    file.path(rdkit_dir, "lotus_flavonoids_rdkit_structural_qc.csv"),
    file.path(OUT_DIR, "lotus_flavonoids_rdkit_structural_qc.csv")
  )
)
path_taxonomy_index <- resolve_file(
  "compound taxonomy index",
  c(
    cfg$path_taxonomy_index %||% NA_character_,
    file.path(part2_dir, paste0(base_tag, "_compound_taxonomy_index.csv.gz")),
    Sys.glob(file.path(part2_dir, "*_compound_taxonomy_index.csv.gz"))
  )
)
path_bio_summary <- resolve_file(
  "compound bioactivity summary",
  c(
    cfg$path_bio_summary %||% NA_character_,
    file.path(part2_dir, paste0(base_tag, "_compound_bioactivity_summary.csv.gz")),
    Sys.glob(file.path(part2_dir, "*_compound_bioactivity_summary.csv.gz"))
  )
)
path_domain_primary <- resolve_file(
  "compound-domain primary table",
  c(
    cfg$path_domain_primary %||% NA_character_,
    file.path(part2_dir, paste0(base_tag, "_compound_domain_primary.csv.gz")),
    Sys.glob(file.path(part2_dir, "*_compound_domain_primary.csv.gz"))
  )
)
path_domain_extended <- resolve_file(
  "compound-domain extended table",
  c(
    cfg$path_domain_extended %||% NA_character_,
    file.path(part2_dir, paste0(base_tag, "_compound_domain_extended.csv.gz")),
    Sys.glob(file.path(part2_dir, "*_compound_domain_extended.csv.gz"))
  )
)
path_domain_binary <- resolve_file(
  "compound-domain binary table",
  c(
    cfg$path_domain_binary %||% NA_character_,
    file.path(part2_dir, paste0(base_tag, "_compound_domain_binary.csv.gz")),
    Sys.glob(file.path(part2_dir, "*_compound_domain_binary.csv.gz"))
  )
)
path_target_primary <- resolve_file(
  "compound-target primary table",
  c(
    cfg$path_target_primary %||% NA_character_,
    file.path(part2_dir, paste0(base_tag, "_compound_target_primary.csv.gz")),
    Sys.glob(file.path(part2_dir, "*_compound_target_primary.csv.gz"))
  )
)
path_clades <- resolve_file(
  "APG clade assignments",
  c(
    cfg$path_apg_clades %||% NA_character_,
    file.path(tree_dir, "APG_clade_assignments.csv")
  )
)
path_tree_metadata <- resolve_file(
  "tree metadata",
  c(
    cfg$path_tree_metadata %||% NA_character_,
    file.path(tree_dir, "tree_apg_itol_run_metadata.json")
  ),
  required = FALSE
)
path_nonvascular <- resolve_file(
  "nonvascular exclusion audit",
  c(
    cfg$path_nonvascular %||% NA_character_,
    file.path(tree_dir, "phylogeny_excluded_nonvascular.csv")
  ),
  required = FALSE
)

input_paths <- list(
  lin_enriched = path_lin,
  rdkit_compounds = path_rdkit_compounds,
  rdkit_qc = path_rdkit_qc,
  taxonomy_index = path_taxonomy_index,
  bioactivity_summary = path_bio_summary,
  domain_primary = path_domain_primary,
  domain_extended = path_domain_extended,
  domain_binary = path_domain_binary,
  target_primary = path_target_primary,
  apg_clades = path_clades,
  tree_metadata = path_tree_metadata,
  nonvascular_audit = path_nonvascular
)

for (nm in names(input_paths)) {
  if (!is.na(input_paths[[nm]] %||% NA_character_)) message("[INPUT] ", nm, ": ", input_paths[[nm]])
}


message("----------------------------------------------------------------")
message(">>> [1/10] Loading and harmonizing occurrence, structure, and clade data")
message("----------------------------------------------------------------")

if (exists("lin_enriched", inherits = TRUE)) {
  occurrence <- as_base_df(get("lin_enriched", inherits = TRUE))
  occurrence_source <- "in_memory:lin_enriched"
} else {
  occurrence <- as_base_df(arrow::read_parquet(path_lin, as_data_frame = TRUE))
  occurrence_source <- path_lin
}
occurrence <- clean_names_simple(occurrence)

col_inchikey <- resolve_column(occurrence, c("inchikey", "inchi_key"), "occurrence InChIKey")
col_family <- resolve_column(occurrence, c("family", "family_name"), "occurrence family")
col_genus <- resolve_column(occurrence, c("genus", "genus_name"), "occurrence genus")
col_species <- resolve_column(occurrence, c("species", "accepted_name"), "occurrence species")
col_ref <- resolve_column(occurrence, c("ref_id", "reference_id", "reference", "doi"), "occurrence reference", required = FALSE)
col_class <- resolve_column(
  occurrence,
  c("chemical_taxonomy_npclassifier_class", "chemicaltaxonomy_npclassifier_class", "chemicaltaxonomynpclassifierclass", "class_np", "npclassifier_class"),
  "NPClassifier class"
)

occurrence_std <- data.frame(
  inchikey = trim_chr(occurrence[[col_inchikey]]),
  family_original = trim_chr(occurrence[[col_family]]),
  genus = trim_chr(occurrence[[col_genus]]),
  species = trim_chr(occurrence[[col_species]]),
  ref_id = if (!is.na(col_ref)) trim_chr(occurrence[[col_ref]]) else NA_character_,
  class_np = trim_chr(occurrence[[col_class]]),
  stringsAsFactors = FALSE
)
occurrence_std$family <- occurrence_std$family_original


correction_mask <- !is.na(occurrence_std$family) & !is.na(occurrence_std$genus) &
  occurrence_std$family == "Capparaceae" & occurrence_std$genus == "Cleome"
occurrence_std$family[correction_mask] <- "Cleomaceae"

correction_audit <- occurrence_std[correction_mask, , drop = FALSE] |>
  dplyr::group_by(family_original, genus, family) |>
  dplyr::summarise(
    n_occurrence_rows = dplyr::n(),
    n_compounds = dplyr::n_distinct(inchikey),
    n_references = dplyr::n_distinct(ref_id[!is.na(ref_id)]),
    reason = "Cleome is treated in Cleomaceae by the APG reference and phylogenetic backbone",
    .groups = "drop"
  )
write_csv_safe(correction_audit, file.path(AUDIT_DIR, "taxonomic_family_corrections.csv"))

occurrence_std <- occurrence_std |>
  dplyr::filter(
    !is.na(inchikey), !is.na(family), !is.na(class_np),
    nzchar(inchikey), nzchar(family), nzchar(class_np)
  ) |>
  dplyr::distinct(inchikey, family, genus, species, ref_id, class_np, .keep_all = TRUE)

rdkit <- read_csv_character(
  path_rdkit_compounds,
  label = "rdkit_compounds",
  audit_dir = AUDIT_DIR
) |>
  clean_names_simple()
rdkit$inchikey <- trim_chr(rdkit$inchikey)


rdkit_context_columns <- intersect(
  c(
    "class_np", "family", "family_original", "genus", "species", "ref_id",
    "fine_clade", "analysis_clade", "clade_color"
  ),
  names(rdkit)
)
if (length(rdkit_context_columns) > 0L) {
  message(
    "[INFO] Removing contextual columns from RDKit compound table before joins: ",
    paste(rdkit_context_columns, collapse = ", ")
  )
  rdkit <- rdkit |>
    dplyr::select(-dplyr::all_of(rdkit_context_columns))
}

qc <- read_csv_character(
  path_rdkit_qc,
  label = "rdkit_structural_qc",
  audit_dir = AUDIT_DIR
) |>
  clean_names_simple()
qc$inchikey <- trim_chr(qc$inchikey)


if (!"structural_qc_status" %in% names(rdkit)) {
  rdkit <- rdkit |>
    dplyr::left_join(
      qc |>
        dplyr::select(
          inchikey,
          structural_qc_status,
          structural_qc_reasons,
          primary_structure_eligible
        ),
      by = "inchikey"
    )
}
rdkit$structural_qc_status <- dplyr::coalesce(trim_chr(rdkit$structural_qc_status), "pass")
rdkit$primary_structure_eligible <- safe_logical(rdkit$primary_structure_eligible)
rdkit$primary_structure_eligible[is.na(rdkit$primary_structure_eligible)] <- TRUE

logical_rdkit <- c(
  "has_phenolic_oh", "has_methoxy_aryl", "has_prenyl_like",
  "has_probable_sugar", "has_probable_sugar_ring",
  "has_alpha_beta_unsaturated_carbonyl", "is_simple_selected_motifs"
)
for (nm in intersect(logical_rdkit, names(rdkit))) rdkit[[nm]] <- safe_logical(rdkit[[nm]])

numeric_rdkit <- c(
  "mw", "exact_mw", "logp", "tpsa", "fsp3", "num_hbd", "num_hba",
  "num_rings", "num_aromatic_rings", "num_rotatable_bonds", "heavy_atom_number",
  "selected_motif_count"
)
for (nm in intersect(numeric_rdkit, names(rdkit))) rdkit[[nm]] <- safe_numeric(rdkit[[nm]])

required_rdkit <- c(
  "inchikey", "bemis_murcko_scaffold_smiles", "mw", "logp", "tpsa", "fsp3",
  "has_methoxy_aryl", "has_prenyl_like", "has_probable_sugar"
)
missing_rdkit <- setdiff(required_rdkit, names(rdkit))
if (length(missing_rdkit)) stop("RDKit compound table is missing: ", paste(missing_rdkit, collapse = ", "))

compound_class <- occurrence_std |>
  dplyr::count(inchikey, class_np, name = "n_occurrences") |>
  dplyr::arrange(inchikey, dplyr::desc(n_occurrences), class_np) |>
  dplyr::group_by(inchikey) |>
  dplyr::slice(1) |>
  dplyr::ungroup() |>
  dplyr::select(inchikey, class_np)

class_conflicts <- occurrence_std |>
  dplyr::group_by(inchikey) |>
  dplyr::summarise(
    n_classes = dplyr::n_distinct(class_np),
    classes = collapse_unique(class_np),
    .groups = "drop"
  ) |>
  dplyr::filter(n_classes > 1)
write_csv_safe(class_conflicts, file.path(AUDIT_DIR, "compound_class_conflicts.csv"))

compound_master <- rdkit |>
  dplyr::left_join(compound_class, by = "inchikey") |>
  dplyr::mutate(
    structural_set_full = TRUE,
    structural_set_primary = primary_structure_eligible,
    structural_set_strict = structural_qc_status == "pass",
    motif_methoxy = as.integer(dplyr::coalesce(has_methoxy_aryl, FALSE)),
    motif_prenyl = as.integer(dplyr::coalesce(has_prenyl_like, FALSE)),
    motif_sugar = as.integer(dplyr::coalesce(has_probable_sugar, FALSE)),
    motif_combination = dplyr::case_when(
      motif_methoxy + motif_prenyl + motif_sugar == 0 ~ "Simple",
      TRUE ~ paste0(
        ifelse(motif_methoxy == 1, "Methoxy+", ""),
        ifelse(motif_prenyl == 1, "Prenyl+", ""),
        ifelse(motif_sugar == 1, "Sugar+", "")
      )
    ),
    motif_combination = sub("\\+$", "", motif_combination)
  )

clades <- read_csv_character(
  path_clades,
  label = "apg_clade_assignments",
  audit_dir = AUDIT_DIR
) |>
  clean_names_simple()
clade_family_col <- resolve_column(clades, c("id", "family"), "clade family")
clade_final_col <- resolve_column(clades, c("final_clade", "finalclade", "macro_group"), "final clade")
clade_analysis_col <- resolve_column(
  clades,
  c("analysis_clade", "analysisclade"),
  "analysis/visual clade",
  required = FALSE
)
clade_color_col <- resolve_column(clades, c("color", "colour"), "clade color", required = FALSE)

clades_std <- data.frame(
  family = trim_chr(clades[[clade_family_col]]),
  fine_clade = trim_chr(clades[[clade_final_col]]),
  visual_clade = if (!is.na(clade_analysis_col)) trim_chr(clades[[clade_analysis_col]]) else NA_character_,
  clade_color = if (!is.na(clade_color_col)) trim_chr(clades[[clade_color_col]]) else NA_character_,
  stringsAsFactors = FALSE
) |>
  dplyr::distinct(family, .keep_all = TRUE)


legacy_visual_clade <- function(x) {
  dplyr::case_when(
    x == "Basal Angiosperms" ~ "ANA",
    x %in% c("Caryophyllales", "Santalales") ~ "Caryophy_Santalales",
    TRUE ~ x
  )
}
clades_std$visual_clade <- dplyr::coalesce(
  clades_std$visual_clade,
  legacy_visual_clade(clades_std$fine_clade)
)

collapse_analysis_clade <- function(x) {
  dplyr::case_when(
    x %in% c("Lycophytes", "Ferns") ~ "Non-seed vascular plants",
    x == "Gymnosperms" ~ "Gymnosperms",
    x %in% c("ANA", "Basal Angiosperms", "Chloranthales", "Magnoliids") ~ "Early angiosperms & magnoliids",
    x %in% c("Monocots", "Monocots (Commelinids)") ~ "Monocots",
    x %in% c("Basal Eudicots", "Other Eudicots", "Saxifragales") ~ "Early eudicots & Saxifragales",
    x %in% c("Rosids (Basal)", "Rosids (Fabids)", "Rosids (Malvids)") ~ "Rosids",
    x %in% c("Caryophy_Santalales", "Caryophyllales", "Santalales") ~ "Caryophyllales & Santalales",
    x %in% c("Asterids (Basal)", "Asterids (Lamiids)", "Asterids (Campanulids)") ~ "Asterids",
    TRUE ~ x
  )
}
clades_std$analysis_clade <- collapse_analysis_clade(clades_std$visual_clade)

analysis_clade_levels <- c(
  "Non-seed vascular plants", "Gymnosperms",
  "Early angiosperms & magnoliids", "Monocots",
  "Early eudicots & Saxifragales", "Rosids",
  "Caryophyllales & Santalales", "Asterids"
)
clades_std$analysis_clade <- factor(clades_std$analysis_clade, levels = analysis_clade_levels)

analysis_palette <- c(
  "Non-seed vascular plants" = "#505050",
  "Gymnosperms" = "#000000",
  "Early angiosperms & magnoliids" = "#D9A441",
  "Monocots" = "#7570B3",
  "Early eudicots & Saxifragales" = "#1B9E77",
  "Rosids" = "#66A61E",
  "Caryophyllales & Santalales" = "#E31A1C",
  "Asterids" = "#377EB8"
)

visual_clade_levels <- c(
  "ANA", "Chloranthales", "Magnoliids", "Monocots", "Monocots (Commelinids)",
  "Basal Eudicots", "Saxifragales", "Rosids (Basal)", "Rosids (Fabids)",
  "Rosids (Malvids)", "Caryophy_Santalales", "Asterids (Basal)",
  "Asterids (Lamiids)", "Asterids (Campanulids)", "Gymnosperms",
  "Ferns", "Lycophytes"
)
if (length(visual_clade_levels) != 17L || !"Chloranthales" %in% visual_clade_levels) {
  stop("Internal configuration error: visual_clade_levels must contain 17 categories including Chloranthales.")
}
message("[CHECK] Visual macroclades: ", length(visual_clade_levels), "; Chloranthales explicit: TRUE")

visual_palette_default <- c(
  "ANA" = "#8E0152", "Chloranthales" = "#B35806", "Magnoliids" = "#FFD92F",
  "Monocots" = "#7570B3", "Monocots (Commelinids)" = "#542788",
  "Basal Eudicots" = "#FC8D62", "Saxifragales" = "#1B9E77",
  "Rosids (Basal)" = "#8FBC8F", "Rosids (Fabids)" = "#66A61E",
  "Rosids (Malvids)" = "#E7298A", "Caryophy_Santalales" = "#E31A1C",
  "Asterids (Basal)" = "#FF7F00", "Asterids (Lamiids)" = "#377EB8",
  "Asterids (Campanulids)" = "#984EA3", "Gymnosperms" = "#000000",
  "Ferns" = "#505050", "Lycophytes" = "#7F7F7F"
)
visual_color_from_tree <- clades_std |>
  dplyr::filter(!is.na(visual_clade), !is.na(clade_color)) |>
  dplyr::count(visual_clade, clade_color, sort = TRUE) |>
  dplyr::group_by(visual_clade) |>
  dplyr::slice(1) |>
  dplyr::ungroup()
visual_palette <- visual_palette_default
if (nrow(visual_color_from_tree)) {
  visual_palette[visual_color_from_tree$visual_clade] <- visual_color_from_tree$clade_color
}
clades_std$visual_clade <- factor(clades_std$visual_clade, levels = visual_clade_levels)

# Fail early if the validated tree classification was not loaded.
chloranth_guard <- clades_std |>
  dplyr::filter(family == "Chloranthaceae")
if (
  nrow(chloranth_guard) != 1L ||
    as.character(chloranth_guard$visual_clade[[1]]) != "Chloranthales"
) {
  stop(
    "The loaded APG assignment is not the validated final version: ",
    "Chloranthaceae must have AnalysisClade/visual_clade = 'Chloranthales'. ",
    "Loaded file: ", path_clades
  )
}
if (any(is.na(clades_std$visual_clade)) || any(is.na(clades_std$analysis_clade))) {
  bad_families <- clades_std$family[
    is.na(clades_std$visual_clade) | is.na(clades_std$analysis_clade)
  ]
  stop(
    "Families without a valid visual or broad analysis clade: ",
    paste(sort(unique(bad_families)), collapse = ", ")
  )
}

compound_family_full <- occurrence_std |>
  dplyr::distinct(inchikey, family) |>
  dplyr::inner_join(clades_std, by = "family")

compound_clade_map <- compound_family_full |>
  dplyr::distinct(inchikey, analysis_clade) |>
  dplyr::group_by(inchikey) |>
  dplyr::mutate(n_analysis_clades = dplyr::n_distinct(analysis_clade)) |>
  dplyr::ungroup()


primary_keys <- compound_master$inchikey[compound_master$structural_set_primary]
strict_keys <- compound_master$inchikey[compound_master$structural_set_strict]
full_keys <- compound_master$inchikey

occurrence_primary <- occurrence_std |>
  dplyr::filter(inchikey %in% primary_keys) |>
  dplyr::inner_join(clades_std, by = "family")

compound_family_primary <- occurrence_primary |>
  dplyr::distinct(inchikey, family, fine_clade, analysis_clade)

input_audit <- data.frame(
  metric = c(
    "occurrence_rows_after_standardization", "unique_compounds_full",
    "unique_compounds_primary", "unique_compounds_strict",
    "vascular_families_primary", "taxonomic_correction_rows",
    "class_conflicts"
  ),
  value = c(
    nrow(occurrence_std), length(unique(full_keys)), length(unique(primary_keys)),
    length(unique(strict_keys)), dplyr::n_distinct(compound_family_primary$family),
    sum(correction_mask), nrow(class_conflicts)
  )
)
write_csv_safe(input_audit, file.path(AUDIT_DIR, "part3_input_audit.csv"))


message("----------------------------------------------------------------")
message(">>> [2/10] Building family-level structural and sampling profiles")
message("----------------------------------------------------------------")

compound_structure_master <- compound_master |>
  dplyr::select(
    -dplyr::any_of(
      c(
        "family", "family_original", "genus", "species", "ref_id",
        "fine_clade", "analysis_clade", "clade_color"
      )
    )
  )

family_compound_struct <- compound_family_primary |>
  dplyr::select(inchikey, family, fine_clade, analysis_clade) |>
  dplyr::inner_join(compound_structure_master, by = "inchikey")

family_sampling <- occurrence_primary |>
  dplyr::group_by(family, fine_clade, analysis_clade) |>
  dplyr::summarise(
    n_occurrence_rows = dplyr::n(),
    n_compounds = dplyr::n_distinct(inchikey),
    n_species = dplyr::n_distinct(species[is_binomial_name(species)]),
    n_genera = dplyr::n_distinct(genus[!is.na(genus)]),
    n_references = dplyr::n_distinct(ref_id[!is.na(ref_id)]),
    n_classes = dplyr::n_distinct(class_np),
    .groups = "drop"
  )

scaffold_family_frequency <- family_compound_struct |>
  dplyr::filter(!is.na(bemis_murcko_scaffold_smiles), nzchar(bemis_murcko_scaffold_smiles)) |>
  dplyr::distinct(family, bemis_murcko_scaffold_smiles) |>
  dplyr::count(bemis_murcko_scaffold_smiles, name = "n_families")

family_scaffold_counts <- family_compound_struct |>
  dplyr::filter(!is.na(bemis_murcko_scaffold_smiles), nzchar(bemis_murcko_scaffold_smiles)) |>
  dplyr::distinct(family, inchikey, bemis_murcko_scaffold_smiles) |>
  dplyr::count(family, bemis_murcko_scaffold_smiles, name = "n_compounds_scaffold")

rarefy_family <- function(df, sample_size) {
  counts <- safe_numeric(df$n_compounds_scaffold)
  counts <- counts[is.finite(counts) & counts > 0]
  total <- sum(counts)
  sample_size <- as.integer(sample_size)

  if (!length(counts) || !is.finite(total) || total < sample_size || sample_size < 1L) {
    return(NA_real_)
  }


  log_denominator <- lchoose(total, sample_size)
  probability_absent <- vapply(
    counts,
    function(n_i) {
      remaining <- total - n_i
      if (remaining < sample_size) return(0)
      exp(lchoose(remaining, sample_size) - log_denominator)
    },
    FUN.VALUE = numeric(1)
  )

  sum(1 - probability_absent)
}

rarefaction <- family_scaffold_counts |>
  dplyr::group_by(family) |>
  tidyr::nest() |>
  dplyr::mutate(
    rarefied_scaffolds_5 = purrr::map_dbl(data, rarefy_family, sample_size = 5L),
    rarefied_scaffolds_10 = purrr::map_dbl(data, rarefy_family, sample_size = 10L),
    rarefied_scaffolds_20 = purrr::map_dbl(data, rarefy_family, sample_size = 20L)
  ) |>
  dplyr::select(-data)

family_metrics <- family_compound_struct |>
  dplyr::group_by(family, fine_clade, analysis_clade) |>
  dplyr::summarise(
    n_compounds_structural = dplyr::n_distinct(inchikey),
    n_scaffolds = dplyr::n_distinct(bemis_murcko_scaffold_smiles[!is.na(bemis_murcko_scaffold_smiles)]),
    median_mw = stats::median(mw, na.rm = TRUE),
    median_logp = stats::median(logp, na.rm = TRUE),
    median_tpsa = stats::median(tpsa, na.rm = TRUE),
    median_fsp3 = stats::median(fsp3, na.rm = TRUE),
    prevalence_phenolic = mean(dplyr::coalesce(has_phenolic_oh, FALSE), na.rm = TRUE),
    prevalence_methoxy = mean(dplyr::coalesce(has_methoxy_aryl, FALSE), na.rm = TRUE),
    prevalence_prenyl = mean(dplyr::coalesce(has_prenyl_like, FALSE), na.rm = TRUE),
    prevalence_sugar = mean(dplyr::coalesce(has_probable_sugar, FALSE), na.rm = TRUE),
    prevalence_enone = mean(dplyr::coalesce(has_alpha_beta_unsaturated_carbonyl, FALSE), na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::left_join(family_sampling, by = c("family", "fine_clade", "analysis_clade")) |>
  dplyr::left_join(rarefaction, by = "family")

family_exclusive <- family_scaffold_counts |>
  dplyr::left_join(scaffold_family_frequency, by = "bemis_murcko_scaffold_smiles") |>
  dplyr::group_by(family) |>
  dplyr::summarise(
    n_exclusive_scaffolds = dplyr::n_distinct(bemis_murcko_scaffold_smiles[n_families == 1]),
    .groups = "drop"
  )
family_metrics <- family_metrics |>
  dplyr::left_join(family_exclusive, by = "family") |>
  dplyr::mutate(
    n_exclusive_scaffolds = dplyr::coalesce(n_exclusive_scaffolds, 0L),
    novelty_ratio = n_scaffolds / pmax(n_compounds_structural, 1),
    exclusive_scaffold_share = n_exclusive_scaffolds / pmax(n_scaffolds, 1),
    analysis_clade = factor(analysis_clade, levels = analysis_clade_levels)
  )


novelty_model_data <- family_metrics |>
  dplyr::filter(n_compounds_structural > 0, n_scaffolds > 0)
novelty_k <- min(6L, max(3L, floor(nrow(novelty_model_data) / 30L)))
novelty_model <- mgcv::gam(
  log1p(n_scaffolds) ~ s(log1p(n_compounds_structural), k = novelty_k),
  data = novelty_model_data,
  method = "REML"
)
novelty_model_data$expected_log_scaffolds <- as.numeric(stats::predict(novelty_model, type = "response"))
novelty_model_data$expected_scaffolds <- pmax(0, expm1(novelty_model_data$expected_log_scaffolds))
novelty_model_data$effort_adjusted_novelty <-
  log1p(novelty_model_data$n_scaffolds) - novelty_model_data$expected_log_scaffolds
novelty_model_data$effort_adjusted_novelty_z <- scale_safe(novelty_model_data$effort_adjusted_novelty)

family_metrics <- family_metrics |>
  dplyr::left_join(
    novelty_model_data |>
      dplyr::select(
        family, expected_scaffolds,
        effort_adjusted_novelty, effort_adjusted_novelty_z
      ),
    by = "family"
  )

write_csv_safe(family_metrics, file.path(TABLE_DIR, "family_structural_metrics_primary.csv"))
write_csv_safe(scaffold_family_frequency, file.path(TABLE_DIR, "scaffold_family_frequency.csv"))

novelty_model_diag <- data.frame(
  n_families = nrow(novelty_model_data),
  smooth_k = novelty_k,
  adjusted_r_squared = summary(novelty_model)$r.sq,
  deviance_explained = summary(novelty_model)$dev.expl,
  scale_estimate = summary(novelty_model)$scale,
  stringsAsFactors = FALSE
)
write_csv_safe(novelty_model_diag, file.path(MODEL_DIR, "sampling_adjusted_novelty_model.csv"))
saveRDS(novelty_model, file.path(MODEL_DIR, "sampling_adjusted_novelty_gam.rds"))

sampling_cor_vars <- family_metrics |>
  dplyr::select(
    n_compounds_structural, n_scaffolds, n_occurrence_rows,
    n_species, n_references, novelty_ratio, effort_adjusted_novelty
  )
sampling_cor <- stats::cor(
  sampling_cor_vars,
  method = "spearman",
  use = "pairwise.complete.obs"
)
write_csv_safe(
  tibble::rownames_to_column(as.data.frame(sampling_cor), "variable"),
  file.path(TABLE_DIR, "sampling_effort_spearman_matrix.csv")
)


compound_class_primary <- compound_master |>
  dplyr::filter(structural_set_primary, !is.na(class_np)) |>
  dplyr::select(inchikey, class_np)

family_class_long <- compound_family_primary |>
  dplyr::select(family, fine_clade, analysis_clade, inchikey) |>
  dplyr::inner_join(compound_class_primary, by = "inchikey") |>
  dplyr::distinct(family, fine_clade, analysis_clade, inchikey, class_np) |>
  dplyr::count(family, fine_clade, analysis_clade, class_np, name = "n_compounds_class") |>
  dplyr::left_join(
    family_metrics |> dplyr::select(family, n_compounds_structural),
    by = "family"
  ) |>
  dplyr::mutate(class_prevalence = n_compounds_class / n_compounds_structural)
write_csv_safe(family_class_long, file.path(TABLE_DIR, "family_class_profile_primary.csv"))


message("----------------------------------------------------------------")
message(">>> [3/10] Running coverage-filtered univariate comparisons")
message("----------------------------------------------------------------")

run_univariate_metric <- function(df, metric, threshold) {
  dat <- df |>
    dplyr::filter(
      n_compounds_structural >= threshold,
      !is.na(analysis_clade),
      is.finite(.data[[metric]])
    ) |>
    dplyr::group_by(analysis_clade) |>
    dplyr::filter(dplyr::n() >= 3) |>
    dplyr::ungroup() |>
    dplyr::mutate(analysis_clade = droplevels(analysis_clade))

  if (nrow(dat) < 10 || dplyr::n_distinct(dat$analysis_clade) < 2) {
    return(list(
      global = data.frame(
        metric = metric, threshold = threshold, n_families = nrow(dat),
        n_groups = dplyr::n_distinct(dat$analysis_clade),
        statistic = NA_real_, df = NA_real_, p_value = NA_real_, epsilon_squared = NA_real_
      ),
      posthoc = data.frame()
    ))
  }

  f <- stats::as.formula(paste(metric, "~ analysis_clade"))
  kw <- stats::kruskal.test(f, data = dat)
  k <- dplyr::n_distinct(dat$analysis_clade)
  eps2 <- max(0, (unname(kw$statistic) - k + 1) / (nrow(dat) - k))

  post <- tryCatch(
    rstatix::dunn_test(dat, formula = f, p.adjust.method = "BH") |>
      dplyr::mutate(metric = metric, threshold = threshold),
    error = function(e) data.frame()
  )

  list(
    global = data.frame(
      metric = metric, threshold = threshold, n_families = nrow(dat),
      n_groups = k, statistic = unname(kw$statistic), df = unname(kw$parameter),
      p_value = kw$p.value, epsilon_squared = eps2
    ),
    posthoc = post
  )
}

univariate_metrics <- c(
  "n_compounds_structural", "n_scaffolds", "rarefied_scaffolds_5",
  "effort_adjusted_novelty_z", "exclusive_scaffold_share",
  "median_mw", "median_logp", "median_tpsa", "median_fsp3",
  "prevalence_methoxy", "prevalence_prenyl", "prevalence_sugar"
)

univ_results <- purrr::map(
  univariate_metrics,
  ~ run_univariate_metric(family_metrics, .x, cfg$analysis_min_compounds_per_taxon)
)
univ_global <- dplyr::bind_rows(purrr::map(univ_results, "global")) |>
  dplyr::mutate(p_adjust_bh = stats::p.adjust(p_value, method = "BH"))
univ_posthoc <- dplyr::bind_rows(purrr::map(univ_results, "posthoc"))
write_csv_safe(univ_global, file.path(TABLE_DIR, "univariate_clade_global_tests.csv"))
write_csv_safe(univ_posthoc, file.path(TABLE_DIR, "univariate_clade_dunn_tests.csv"))


message("----------------------------------------------------------------")
message(">>> [4/10] Running multivariate scaffold, class, and structural-profile analyses")
message("----------------------------------------------------------------")

eligible_families_primary <- family_metrics |>
  dplyr::filter(n_compounds_structural >= cfg$analysis_min_compounds_per_taxon) |>
  dplyr::group_by(analysis_clade) |>
  dplyr::filter(dplyr::n() >= 3) |>
  dplyr::ungroup() |>
  dplyr::pull(family)

family_meta <- family_metrics |>
  dplyr::filter(family %in% eligible_families_primary) |>
  dplyr::select(family, fine_clade, analysis_clade, n_compounds_structural) |>
  dplyr::arrange(family)

scaffold_incidence <- family_compound_struct |>
  dplyr::filter(family %in% eligible_families_primary) |>
  dplyr::filter(!is.na(bemis_murcko_scaffold_smiles), nzchar(bemis_murcko_scaffold_smiles)) |>
  dplyr::distinct(family, bemis_murcko_scaffold_smiles) |>
  dplyr::mutate(value = 1L) |>
  tidyr::pivot_wider(
    names_from = bemis_murcko_scaffold_smiles,
    values_from = value,
    values_fill = 0L
  )

class_incidence <- family_class_long |>
  dplyr::filter(family %in% eligible_families_primary) |>
  dplyr::mutate(value = as.integer(n_compounds_class > 0)) |>
  dplyr::select(family, class_np, value) |>
  tidyr::pivot_wider(names_from = class_np, values_from = value, values_fill = 0L)

structural_profile_cols <- c(
  "median_mw", "median_logp", "median_tpsa", "median_fsp3",
  "prevalence_phenolic", "prevalence_methoxy", "prevalence_prenyl",
  "prevalence_sugar", "prevalence_enone"
)
structural_profile <- family_metrics |>
  dplyr::filter(family %in% eligible_families_primary) |>
  dplyr::select(family, dplyr::all_of(structural_profile_cols))

prepare_matrix <- function(
  tbl,
  meta,
  binary = FALSE,
  prevalence_min = 0,
  analysis_name = "matrix"
) {
  tbl <- as_base_df(tbl)
  rn <- trim_chr(tbl$family)
  mat <- as.matrix(tbl[, setdiff(names(tbl), "family"), drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- rn

  common <- intersect(meta$family, rownames(mat))
  mat <- mat[common, , drop = FALSE]
  meta2 <- meta[match(common, meta$family), , drop = FALSE]
  meta2$analysis_clade <- droplevels(factor(meta2$analysis_clade))

  audit <- data.frame(
    analysis = character(0),
    family = character(0),
    reason = character(0),
    stringsAsFactors = FALSE
  )

  add_audit <- function(families, reason) {
    families <- trim_chr(families)
    families <- families[!is.na(families)]
    if (!length(families)) return(invisible(NULL))
    audit <<- dplyr::bind_rows(
      audit,
      data.frame(
        analysis = analysis_name,
        family = families,
        reason = reason,
        stringsAsFactors = FALSE
      )
    )
    invisible(NULL)
  }

  if (!nrow(mat) || !ncol(mat)) {
    return(list(matrix = mat, meta = meta2, audit = audit))
  }

  if (binary) {
    mat[!is.finite(mat)] <- 0
    mat <- (mat > 0) * 1

    minimum_frequency <- max(2L, ceiling(prevalence_min * nrow(mat)))
    frequency <- colSums(mat > 0, na.rm = TRUE)
    keep_columns <- frequency >= minimum_frequency & frequency < nrow(mat)
    mat <- mat[, keep_columns, drop = FALSE]

    if (ncol(mat)) {
      empty_rows <- rowSums(mat > 0, na.rm = TRUE) == 0
      add_audit(rownames(mat)[empty_rows], "empty_after_binary_feature_filter")
      mat <- mat[!empty_rows, , drop = FALSE]
      meta2 <- meta2[!empty_rows, , drop = FALSE]
    }


    if (nrow(mat) && ncol(mat)) {
      frequency <- colSums(mat > 0, na.rm = TRUE)
      keep_columns <- frequency >= 2L & frequency < nrow(mat)
      mat <- mat[, keep_columns, drop = FALSE]
    }

    if (ncol(mat)) {
      empty_rows <- rowSums(mat > 0, na.rm = TRUE) == 0
      add_audit(rownames(mat)[empty_rows], "empty_after_binary_refilter")
      mat <- mat[!empty_rows, , drop = FALSE]
      meta2 <- meta2[!empty_rows, , drop = FALSE]
    }
  } else {
    finite_rows <- apply(mat, 1, function(z) all(is.finite(z)))
    add_audit(rownames(mat)[!finite_rows], "nonfinite_structural_profile")
    mat <- mat[finite_rows, , drop = FALSE]
    meta2 <- meta2[finite_rows, , drop = FALSE]

    if (nrow(mat) && ncol(mat)) {
      keep_columns <- vapply(
        seq_len(ncol(mat)),
        function(j) {
          s <- stats::sd(mat[, j], na.rm = TRUE)
          is.finite(s) && s > 0
        },
        FUN.VALUE = logical(1)
      )
      mat <- mat[, keep_columns, drop = FALSE]
    }

    if (nrow(mat) && ncol(mat)) {
      mat <- scale(mat)
      finite_rows <- apply(mat, 1, function(z) all(is.finite(z)))
      add_audit(rownames(mat)[!finite_rows], "nonfinite_after_scaling")
      mat <- mat[finite_rows, , drop = FALSE]
      meta2 <- meta2[finite_rows, , drop = FALSE]
    }
  }

  meta2 <- meta2[match(rownames(mat), meta2$family), , drop = FALSE]
  meta2$analysis_clade <- droplevels(factor(meta2$analysis_clade))

  list(matrix = mat, meta = meta2, audit = audit)
}

run_multivariate <- function(
  name,
  prepared,
  distance_method,
  binary = FALSE,
  permutations = 999L,
  run_pairwise = TRUE
) {
  mat <- prepared$matrix
  meta <- prepared$meta

  if (nrow(mat) < 5L || ncol(mat) < 1L) {
    stop("Insufficient non-empty matrix for ", name, ": ", nrow(mat), " families x ", ncol(mat), " features.")
  }
  if (dplyr::n_distinct(meta$analysis_clade) < 2L) {
    stop("Fewer than two analytical clades remained for ", name, ".")
  }

  if (binary) {
    empty_rows <- rowSums(mat > 0, na.rm = TRUE) == 0
    if (any(empty_rows)) {
      stop("Internal error: empty binary rows remained in ", name, ".")
    }
  }

  d <- if (distance_method == "euclidean") {
    stats::dist(mat, method = "euclidean")
  } else {
    vegan::vegdist(mat, method = distance_method, binary = binary, na.rm = FALSE)
  }

  distance_values <- as.numeric(d)
  if (!length(distance_values) || any(!is.finite(distance_values))) {
    stop("Non-finite distances remained in ", name, " after matrix sanitization.")
  }
  if (all(abs(distance_values) < .Machine$double.eps^0.5)) {
    stop("All pairwise distances are zero for ", name, ".")
  }

  ad <- vegan::adonis2(d ~ analysis_clade, data = meta, permutations = permutations)
  global <- data.frame(
    analysis = name,
    distance = distance_method,
    binary = binary,
    n_families = nrow(mat),
    n_features = ncol(mat),
    pseudo_f = ad$F[1],
    r_squared = ad$R2[1],
    p_value = ad$`Pr(>F)`[1],
    stringsAsFactors = FALSE
  )


  bd <- vegan::betadisper(
    d,
    group = meta$analysis_clade,
    bias.adjust = TRUE,
    add = "lingoes"
  )
  bd_perm <- vegan::permutest(bd, permutations = permutations)
  dispersion <- data.frame(
    analysis = name,
    distance_correction = "lingoes",
    statistic = bd_perm$tab[1, "F"],
    p_value = bd_perm$tab[1, "Pr(>F)"],
    stringsAsFactors = FALSE
  )

  pairwise <- data.frame()
  if (isTRUE(run_pairwise)) {
    clade_levels <- levels(droplevels(meta$analysis_clade))
    pairs <- if (length(clade_levels) >= 2L) {
      utils::combn(clade_levels, 2, simplify = FALSE)
    } else {
      list()
    }

    pairwise <- purrr::map_dfr(pairs, function(pair) {
      idx <- meta$analysis_clade %in% pair
      sub_meta <- droplevels(meta[idx, , drop = FALSE])
      sub_mat <- mat[idx, , drop = FALSE]

      if (binary && ncol(sub_mat)) {
        frequency <- colSums(sub_mat > 0, na.rm = TRUE)
        sub_mat <- sub_mat[, frequency > 0 & frequency < nrow(sub_mat), drop = FALSE]
        if (ncol(sub_mat)) {
          nonempty <- rowSums(sub_mat > 0, na.rm = TRUE) > 0
          sub_mat <- sub_mat[nonempty, , drop = FALSE]
          sub_meta <- sub_meta[nonempty, , drop = FALSE]
        }
      } else if (ncol(sub_mat)) {
        keep_columns <- vapply(
          seq_len(ncol(sub_mat)),
          function(j) {
            s <- stats::sd(sub_mat[, j], na.rm = TRUE)
            is.finite(s) && s > 0
          },
          FUN.VALUE = logical(1)
        )
        sub_mat <- sub_mat[, keep_columns, drop = FALSE]
      }

      group_counts <- table(sub_meta$analysis_clade)
      base_row <- data.frame(
        analysis = name,
        group1 = pair[1],
        group2 = pair[2],
        n1 = sum(sub_meta$analysis_clade == pair[1], na.rm = TRUE),
        n2 = sum(sub_meta$analysis_clade == pair[2], na.rm = TRUE),
        pseudo_f = NA_real_,
        r_squared = NA_real_,
        p_value = NA_real_,
        status = "skipped",
        stringsAsFactors = FALSE
      )

      if (nrow(sub_mat) < 4L || ncol(sub_mat) < 1L || any(group_counts < 2L)) {
        base_row$status <- "insufficient_nonempty_pairwise_matrix"
        return(base_row)
      }

      sub_d <- if (distance_method == "euclidean") {
        stats::dist(sub_mat, method = "euclidean")
      } else {
        vegan::vegdist(sub_mat, method = distance_method, binary = binary, na.rm = FALSE)
      }
      sub_values <- as.numeric(sub_d)
      if (!length(sub_values) || any(!is.finite(sub_values)) || all(abs(sub_values) < .Machine$double.eps^0.5)) {
        base_row$status <- "noninformative_or_nonfinite_distance"
        return(base_row)
      }

      fit <- vegan::adonis2(sub_d ~ analysis_clade, data = sub_meta, permutations = permutations)
      base_row$pseudo_f <- fit$F[1]
      base_row$r_squared <- fit$R2[1]
      base_row$p_value <- fit$`Pr(>F)`[1]
      base_row$status <- "ok"
      base_row
    }) |>
      dplyr::mutate(p_adjust_bh = stats::p.adjust(p_value, method = "BH"))
  }

  pc <- tryCatch(
    stats::cmdscale(d, k = 2, eig = TRUE, add = TRUE),
    error = function(e) {
      message("[WARN] PCoA failed for ", name, ": ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(pc)) {
    pcoa <- meta |>
      dplyr::mutate(
        Axis1 = NA_real_,
        Axis2 = NA_real_,
        axis1_percent = NA_real_,
        axis2_percent = NA_real_
      )
  } else {
    pcoa <- as.data.frame(pc$points)
    names(pcoa)[1:2] <- c("Axis1", "Axis2")
    pcoa$family <- rownames(pc$points)
    pcoa <- pcoa |>
      dplyr::left_join(meta, by = "family")
    positive_eig <- pc$eig[pc$eig > 0]
    variance <- if (length(positive_eig) >= 2L) {
      100 * positive_eig[1:2] / sum(positive_eig)
    } else {
      c(NA_real_, NA_real_)
    }
    pcoa$axis1_percent <- variance[1]
    pcoa$axis2_percent <- variance[2]
  }

  list(
    global = global,
    dispersion = dispersion,
    pairwise = pairwise,
    pcoa = pcoa,
    distance = d,
    audit = prepared$audit %||% data.frame()
  )
}

prepared_scaffold <- prepare_matrix(
  scaffold_incidence,
  family_meta,
  binary = TRUE,
  analysis_name = "Bemis-Murcko scaffold incidence"
)
prepared_class <- prepare_matrix(
  class_incidence,
  family_meta,
  binary = TRUE,
  analysis_name = "NPClassifier class incidence"
)
prepared_struct <- prepare_matrix(
  structural_profile,
  family_meta,
  binary = FALSE,
  analysis_name = "Descriptor and motif profile"
)

multi_scaffold <- run_multivariate("Bemis-Murcko scaffold incidence", prepared_scaffold, "jaccard", TRUE, cfg$permutations)
multi_class <- run_multivariate("NPClassifier class incidence", prepared_class, "jaccard", TRUE, cfg$permutations)
multi_struct <- run_multivariate("Descriptor and motif profile", prepared_struct, "euclidean", FALSE, cfg$permutations)

multi_global <- dplyr::bind_rows(multi_scaffold$global, multi_class$global, multi_struct$global) |>
  dplyr::mutate(p_adjust_bh = stats::p.adjust(p_value, method = "BH"))
multi_disp <- dplyr::bind_rows(multi_scaffold$dispersion, multi_class$dispersion, multi_struct$dispersion) |>
  dplyr::mutate(p_adjust_bh = stats::p.adjust(p_value, method = "BH"))
multi_pair <- dplyr::bind_rows(multi_scaffold$pairwise, multi_class$pairwise, multi_struct$pairwise)
write_csv_safe(multi_global, file.path(TABLE_DIR, "multivariate_permanova_global.csv"))
write_csv_safe(multi_disp, file.path(TABLE_DIR, "multivariate_dispersion_tests.csv"))
write_csv_safe(multi_pair, file.path(TABLE_DIR, "multivariate_pairwise_permanova.csv"))
write_csv_safe(multi_scaffold$pcoa, file.path(TABLE_DIR, "pcoa_scaffolds.csv"))
write_csv_safe(multi_class$pcoa, file.path(TABLE_DIR, "pcoa_classes.csv"))
write_csv_safe(multi_struct$pcoa, file.path(TABLE_DIR, "pcoa_structural_profile.csv"))
multivariate_exclusion_audit <- dplyr::bind_rows(
  prepared_scaffold$audit,
  prepared_class$audit,
  prepared_struct$audit
) |>
  dplyr::distinct()
write_csv_safe(
  multivariate_exclusion_audit,
  file.path(AUDIT_DIR, "multivariate_excluded_families.csv")
)


build_scaffold_sensitivity <- function(keys, set_name, threshold) {
  dat <- compound_family_full |>
    dplyr::filter(inchikey %in% keys) |>
    dplyr::select(inchikey, family, analysis_clade) |>
    dplyr::inner_join(
      compound_master |>
        dplyr::select(inchikey, bemis_murcko_scaffold_smiles),
      by = "inchikey"
    ) |>
    dplyr::filter(!is.na(bemis_murcko_scaffold_smiles), nzchar(bemis_murcko_scaffold_smiles))

  coverage <- dat |>
    dplyr::group_by(family, analysis_clade) |>
    dplyr::summarise(n_compounds = dplyr::n_distinct(inchikey), .groups = "drop") |>
    dplyr::filter(n_compounds >= threshold) |>
    dplyr::group_by(analysis_clade) |>
    dplyr::filter(dplyr::n() >= 3) |>
    dplyr::ungroup()

  if (nrow(coverage) < 10 || dplyr::n_distinct(coverage$analysis_clade) < 2) {
    return(data.frame(
      structural_set = set_name,
      threshold = threshold,
      n_families = nrow(coverage),
      n_scaffolds = NA_real_,
      r_squared = NA_real_,
      p_value = NA_real_,
      dispersion_p = NA_real_,
      status = "insufficient_family_or_clade_coverage",
      stringsAsFactors = FALSE
    ))
  }

  inc <- dat |>
    dplyr::filter(family %in% coverage$family) |>
    dplyr::distinct(family, bemis_murcko_scaffold_smiles) |>
    dplyr::mutate(value = 1L) |>
    tidyr::pivot_wider(names_from = bemis_murcko_scaffold_smiles, values_from = value, values_fill = 0L)
  prep <- prepare_matrix(
    inc,
    coverage |> dplyr::select(family, analysis_clade),
    binary = TRUE,
    analysis_name = paste0("sensitivity_", set_name, "_", threshold)
  )

  res <- tryCatch(
    run_multivariate(
      paste(set_name, threshold),
      prep,
      "jaccard",
      TRUE,
      cfg$permutations,
      run_pairwise = FALSE
    ),
    error = function(e) e
  )

  if (inherits(res, "error")) {
    return(data.frame(
      structural_set = set_name,
      threshold = threshold,
      n_families = nrow(prep$matrix),
      n_scaffolds = ncol(prep$matrix),
      r_squared = NA_real_,
      p_value = NA_real_,
      dispersion_p = NA_real_,
      status = paste0("skipped: ", conditionMessage(res)),
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    structural_set = set_name,
    threshold = threshold,
    n_families = res$global$n_families,
    n_scaffolds = res$global$n_features,
    r_squared = res$global$r_squared,
    p_value = res$global$p_value,
    dispersion_p = res$dispersion$p_value,
    status = "ok",
    stringsAsFactors = FALSE
  )
}

sensitivity_grid <- tidyr::crossing(
  structural_set = c("full", "primary", "strict"),
  threshold = cfg$analysis_compound_thresholds
)
sensitivity_scaffold <- purrr::pmap_dfr(sensitivity_grid, function(structural_set, threshold) {
  keys <- switch(structural_set, full = full_keys, primary = primary_keys, strict = strict_keys)
  build_scaffold_sensitivity(keys, structural_set, threshold)
}) |>
  dplyr::mutate(p_adjust_bh = stats::p.adjust(p_value, method = "BH"))
write_csv_safe(sensitivity_scaffold, file.path(TABLE_DIR, "sensitivity_scaffold_permanova.csv"))


message("----------------------------------------------------------------")
message(">>> [5/10] Building ChEMBL functional and coverage layers")
message("----------------------------------------------------------------")

bio_summary <- read_csv_character(
  path_bio_summary,
  label = "compound_bioactivity_summary",
  audit_dir = AUDIT_DIR
) |>
  clean_names_simple()
domain_primary <- read_csv_character(
  path_domain_primary,
  label = "compound_domain_primary",
  audit_dir = AUDIT_DIR
) |>
  clean_names_simple()
domain_extended <- read_csv_character(
  path_domain_extended,
  label = "compound_domain_extended",
  audit_dir = AUDIT_DIR
) |>
  clean_names_simple()
domain_binary <- read_csv_character(
  path_domain_binary,
  label = "compound_domain_binary",
  audit_dir = AUDIT_DIR
) |>
  clean_names_simple()
target_primary <- read_csv_character(
  path_target_primary,
  label = "compound_target_primary",
  audit_dir = AUDIT_DIR
) |>
  clean_names_simple()

for (obj_name in c("bio_summary", "domain_primary", "domain_extended", "domain_binary", "target_primary")) {
  obj <- get(obj_name)
  obj$inchikey <- trim_chr(obj$inchikey)
  assign(obj_name, obj)
}

for (nm in intersect(c("mapped_to_chembl", "any_definite_active_10um", "any_definite_above_10um"), names(bio_summary))) {
  bio_summary[[nm]] <- safe_logical(bio_summary[[nm]])
}
for (nm in intersect(c("n_primary_records", "n_quality_records", "n_targets", "n_domains", "median_pchembl"), names(bio_summary))) {
  bio_summary[[nm]] <- safe_numeric(bio_summary[[nm]])
}
for (nm in intersect(c("any_definite_active", "any_definite_above"), names(domain_binary))) {
  domain_binary[[nm]] <- safe_logical(domain_binary[[nm]])
}
for (nm in intersect(c("n_quality_records", "n_definite_active", "n_definite_above", "n_indeterminate"), names(domain_binary))) {
  domain_binary[[nm]] <- safe_numeric(domain_binary[[nm]])
}
for (nm in intersect(c("n_targets", "n_target_summaries", "n_records", "n_assays", "median_target_pchembl", "q25_target_pchembl", "q75_target_pchembl"), names(domain_primary))) {
  domain_primary[[nm]] <- safe_numeric(domain_primary[[nm]])
}

# Binary evidence is classified at the compound x functional-domain level.
# The primary analysis uses only concordant evidence and keeps mixed evidence
# separate, preventing a single positive record from overriding contrary data.
domain_binary <- domain_binary |>
  dplyr::mutate(
    n_quality_records = dplyr::coalesce(n_quality_records, 0),
    n_definite_active = dplyr::coalesce(n_definite_active, 0),
    n_definite_above = dplyr::coalesce(n_definite_above, 0),
    n_indeterminate = dplyr::coalesce(n_indeterminate, 0),
    n_classifiable = n_definite_active + n_definite_above,
    active_fraction_classifiable = dplyr::if_else(
      n_classifiable > 0,
      n_definite_active / n_classifiable,
      NA_real_
    ),
    binary_evidence_status = dplyr::case_when(
      n_definite_active > 0 & n_definite_above == 0 ~ "active_only",
      n_definite_active == 0 & n_definite_above > 0 ~ "above_cutoff_only",
      n_definite_active > 0 & n_definite_above > 0 ~ "mixed_evidence",
      TRUE ~ "indeterminate"
    ),
    binary_majority_status = dplyr::case_when(
      n_classifiable == 0 ~ "indeterminate",
      n_definite_active > n_definite_above ~ "active_majority",
      n_definite_above > n_definite_active ~ "above_cutoff_majority",
      TRUE ~ "tie"
    )
  )

binary_evidence_status_summary <- domain_binary |>
  dplyr::count(binary_evidence_status, name = "n_compound_domain_combinations") |>
  dplyr::mutate(
    fraction = n_compound_domain_combinations / sum(n_compound_domain_combinations)
  )
write_csv_safe(
  binary_evidence_status_summary,
  file.path(TABLE_DIR, "binary_evidence_status_summary.csv")
)

binary_evidence_domain_summary <- domain_binary |>
  dplyr::group_by(target_domain_macro, binary_evidence_status) |>
  dplyr::summarise(
    n_compounds = dplyr::n_distinct(inchikey),
    n_compound_domain_combinations = dplyr::n(),
    n_quality_records = sum(n_quality_records, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::group_by(target_domain_macro) |>
  dplyr::mutate(
    fraction_within_domain = n_compound_domain_combinations / sum(n_compound_domain_combinations)
  ) |>
  dplyr::ungroup()
write_csv_safe(
  binary_evidence_domain_summary,
  file.path(TABLE_DIR, "binary_evidence_status_by_domain.csv")
)

compound_binary_status <- domain_binary |>
  dplyr::group_by(inchikey) |>
  dplyr::summarise(
    n_classifiable_domains = sum(n_classifiable > 0, na.rm = TRUE),
    n_concordant_domains = sum(binary_evidence_status %in% c("active_only", "above_cutoff_only"), na.rm = TRUE),
    n_active_only_domains = sum(binary_evidence_status == "active_only", na.rm = TRUE),
    n_above_cutoff_only_domains = sum(binary_evidence_status == "above_cutoff_only", na.rm = TRUE),
    n_mixed_evidence_domains = sum(binary_evidence_status == "mixed_evidence", na.rm = TRUE),
    any_active_only = any(binary_evidence_status == "active_only", na.rm = TRUE),
    any_above_cutoff_only = any(binary_evidence_status == "above_cutoff_only", na.rm = TRUE),
    any_mixed_evidence = any(binary_evidence_status == "mixed_evidence", na.rm = TRUE),
    .groups = "drop"
  )
write_csv_safe(
  compound_binary_status,
  file.path(TABLE_DIR, "compound_binary_evidence_status.csv")
)


primary_enriched <- domain_primary |>
  dplyr::inner_join(
    compound_master |>
      dplyr::filter(structural_set_primary) |>
      dplyr::select(
        inchikey, class_np, mw, logp, tpsa, fsp3,
        motif_methoxy, motif_prenyl, motif_sugar, motif_combination,
        structural_qc_status
      ),
    by = "inchikey"
  ) |>
  dplyr::filter(
    !is.na(median_target_pchembl),
    !is.na(target_domain_macro),
    target_domain_macro != "Miscellaneous"
  ) |>
  dplyr::group_by(target_domain_macro, standard_type) |>
  dplyr::mutate(
    type_group_n = dplyr::n(),
    z_potency_type = ifelse(type_group_n >= 5, scale_safe(median_target_pchembl), NA_real_)
  ) |>
  dplyr::ungroup()

compound_domain_primary <- primary_enriched |>
  dplyr::group_by(
    inchikey, target_domain_macro, class_np, motif_methoxy, motif_prenyl,
    motif_sugar, motif_combination, mw, logp, tpsa, fsp3,
    structural_qc_status
  ) |>
  dplyr::summarise(
    median_pchembl = median_finite_safe(median_target_pchembl),
    z_potency_equal_mean = mean_finite_safe(z_potency_type),
    z_potency_equal_median = median_finite_safe(z_potency_type),
    z_potency_record_weighted = weighted_mean_safe(z_potency_type, pmax(n_records, 1)),
    n_records = sum(n_records, na.rm = TRUE),
    n_targets = sum(n_targets, na.rm = TRUE),
    n_assays = sum(n_assays, na.rm = TRUE),
    n_standard_types = dplyr::n_distinct(standard_type[is.finite(z_potency_type)]),
    standard_types = collapse_unique(standard_type[is.finite(z_potency_type)]),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    z_potency = if (identical(cfg$functional_primary_type_aggregation, "equal_median")) {
      z_potency_equal_median
    } else {
      z_potency_equal_mean
    },
    delta_equal_mean_minus_record_weighted = z_potency_equal_mean - z_potency_record_weighted,
    delta_equal_median_minus_equal_mean = z_potency_equal_median - z_potency_equal_mean
  ) |>
  dplyr::filter(is.finite(z_potency))

write_csv_safe(compound_domain_primary, file.path(TABLE_DIR, "compound_domain_primary_standardized.csv"))
write_csv_safe(
  compound_domain_primary |>
    dplyr::select(
      inchikey, target_domain_macro, n_standard_types, standard_types,
      n_records, n_targets, n_assays,
      z_potency_equal_mean, z_potency_equal_median,
      z_potency_record_weighted,
      delta_equal_mean_minus_record_weighted,
      delta_equal_median_minus_equal_mean
    ),
  file.path(TABLE_DIR, "compound_domain_potency_aggregation_sensitivity.csv")
)

family_bio_coverage <- compound_family_primary |>
  dplyr::select(family, analysis_clade, inchikey) |>
  dplyr::distinct() |>
  dplyr::left_join(
    bio_summary |>
      dplyr::select(
        inchikey, mapped_to_chembl, n_primary_records, n_quality_records,
        n_targets, n_domains, median_pchembl
      ),
    by = "inchikey"
  ) |>
  dplyr::left_join(compound_binary_status, by = "inchikey") |>
  dplyr::group_by(family, analysis_clade) |>
  dplyr::summarise(
    n_compounds = dplyr::n_distinct(inchikey),
    n_mapped_compounds = dplyr::n_distinct(inchikey[dplyr::coalesce(mapped_to_chembl, FALSE)]),
    n_primary_compounds = dplyr::n_distinct(inchikey[dplyr::coalesce(n_primary_records, 0) > 0]),
    n_binary_evidence_compounds = dplyr::n_distinct(inchikey[dplyr::coalesce(n_classifiable_domains, 0) > 0]),
    n_concordant_binary_compounds = dplyr::n_distinct(inchikey[dplyr::coalesce(n_concordant_domains, 0) > 0]),
    n_active_only_compounds_10um = dplyr::n_distinct(inchikey[dplyr::coalesce(any_active_only, FALSE)]),
    n_above_cutoff_only_compounds_10um = dplyr::n_distinct(inchikey[dplyr::coalesce(any_above_cutoff_only, FALSE)]),
    n_mixed_evidence_compounds_10um = dplyr::n_distinct(inchikey[dplyr::coalesce(any_mixed_evidence, FALSE)]),
    n_primary_records = sum(dplyr::coalesce(n_primary_records, 0), na.rm = TRUE),
    n_targets_sum = sum(dplyr::coalesce(n_targets, 0), na.rm = TRUE),
    median_compound_pchembl = median_finite_safe(median_pchembl),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    mapped_fraction = n_mapped_compounds / pmax(n_compounds, 1),
    primary_fraction = n_primary_compounds / pmax(n_compounds, 1),
    active_fraction_all_compounds_10um = n_active_only_compounds_10um / pmax(n_compounds, 1),
    active_rate_concordant_10um = n_active_only_compounds_10um / pmax(n_concordant_binary_compounds, 1),
    mixed_evidence_fraction_binary = n_mixed_evidence_compounds_10um / pmax(n_binary_evidence_compounds, 1),
    # Backward-compatible field, now based on the concordant-only primary rule.
    active_fraction_10um = active_rate_concordant_10um,
    n_active_compounds_10um = n_active_only_compounds_10um
  )
family_bio_coverage$median_compound_pchembl[!is.finite(family_bio_coverage$median_compound_pchembl)] <- NA_real_
write_csv_safe(family_bio_coverage, file.path(TABLE_DIR, "family_bioactivity_coverage.csv"))

clade_bio_coverage <- compound_clade_map |>
  dplyr::filter(inchikey %in% primary_keys) |>
  dplyr::select(inchikey, analysis_clade) |>
  dplyr::distinct() |>
  dplyr::left_join(
    bio_summary |>
      dplyr::select(inchikey, mapped_to_chembl, n_primary_records),
    by = "inchikey"
  ) |>
  dplyr::group_by(analysis_clade) |>
  dplyr::summarise(
    n_compounds = dplyr::n_distinct(inchikey),
    n_mapped_compounds = dplyr::n_distinct(inchikey[dplyr::coalesce(mapped_to_chembl, FALSE)]),
    n_primary_compounds = dplyr::n_distinct(inchikey[dplyr::coalesce(n_primary_records, 0) > 0]),
    weighted_mapped_fraction = n_mapped_compounds / pmax(n_compounds, 1),
    weighted_primary_fraction = n_primary_compounds / pmax(n_compounds, 1),
    .groups = "drop"
  )
write_csv_safe(clade_bio_coverage, file.path(TABLE_DIR, "clade_bioactivity_coverage.csv"))


single_clade_map <- compound_clade_map |>
  dplyr::filter(n_analysis_clades == 1) |>
  dplyr::select(inchikey, analysis_clade)

compound_domain_clade_single <- compound_domain_primary |>
  dplyr::inner_join(single_clade_map, by = "inchikey")

clade_domain_summary <- compound_domain_clade_single |>
  dplyr::group_by(analysis_clade, target_domain_macro) |>
  dplyr::summarise(
    n_compounds = dplyr::n_distinct(inchikey),
    n_records = sum(n_records, na.rm = TRUE),
    n_targets = sum(n_targets, na.rm = TRUE),
    q25_pchembl = if (any(is.finite(median_pchembl))) {
      as.numeric(stats::quantile(
        median_pchembl[is.finite(median_pchembl)],
        probs = 0.25,
        names = FALSE
      ))
    } else NA_real_,
    q75_pchembl = if (any(is.finite(median_pchembl))) {
      as.numeric(stats::quantile(
        median_pchembl[is.finite(median_pchembl)],
        probs = 0.75,
        names = FALSE
      ))
    } else NA_real_,
    median_z_potency = median_finite_safe(z_potency),
    median_pchembl = median_finite_safe(median_pchembl),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    sufficient_coverage = n_compounds >= cfg$functional_min_compounds_domain &
      n_records >= cfg$functional_min_records_domain,
    analysis_clade = factor(analysis_clade, levels = analysis_clade_levels)
  )
write_csv_safe(clade_domain_summary, file.path(TABLE_DIR, "clade_domain_primary_summary_single_clade.csv"))

compound_domain_clade_all <- compound_domain_primary |>
  dplyr::inner_join(
    compound_clade_map,
    by = "inchikey",
    relationship = "many-to-many"
  ) |>
  dplyr::mutate(clade_weight = 1 / n_analysis_clades)
write_csv_safe(compound_domain_clade_all, file.path(TABLE_DIR, "compound_domain_clade_all_sensitivity.csv"))


message("----------------------------------------------------------------")
message(">>> [6/10] Evaluating nonlinear Fsp3-potency relationships")
message("----------------------------------------------------------------")

fit_fsp3_domain <- function(dat, domain_name, use_record_weights = FALSE) {
  dat <- dat |>
    dplyr::filter(
      is.finite(z_potency), is.finite(fsp3), is.finite(logp),
      is.finite(tpsa), is.finite(mw)
    ) |>
    dplyr::mutate(
      z_logp = scale_safe(logp),
      z_tpsa = scale_safe(tpsa),
      z_mw = scale_safe(mw),
      model_weight = if (isTRUE(use_record_weights)) sqrt(pmax(n_records, 1)) else 1
    )

  if (nrow(dat) < cfg$fsp3_gam_min_compounds || dplyr::n_distinct(dat$fsp3) < 10) {
    return(list(summary = data.frame(
      target_domain_macro = domain_name, n_compounds = nrow(dat), edf = NA,
      smooth_p = NA, adjusted_r_squared = NA, deviance_explained = NA,
      aic_gam = NA, aic_linear = NA, delta_aic_linear_minus_gam = NA,
      weighting = if (isTRUE(use_record_weights)) "sqrt_n_records" else "equal_compound"
    ), prediction = data.frame(), model = NULL))
  }

  k <- min(6L, max(4L, dplyr::n_distinct(dat$fsp3) - 1L))
  gam_fit <- mgcv::gam(
    z_potency ~ s(fsp3, k = k) + z_logp + z_tpsa + z_mw,
    data = dat,
    weights = model_weight,
    method = "REML"
  )
  linear_fit <- stats::lm(
    z_potency ~ fsp3 + z_logp + z_tpsa + z_mw,
    data = dat,
    weights = model_weight
  )

  sm <- summary(gam_fit)
  s_table <- sm$s.table
  grid <- data.frame(
    fsp3 = seq(stats::quantile(dat$fsp3, 0.025), stats::quantile(dat$fsp3, 0.975), length.out = 100),
    z_logp = 0, z_tpsa = 0, z_mw = 0
  )
  pred <- stats::predict(gam_fit, newdata = grid, se.fit = TRUE, type = "response")
  grid$fit <- as.numeric(pred$fit)
  grid$se <- as.numeric(pred$se.fit)
  grid$lower <- grid$fit - 1.96 * grid$se
  grid$upper <- grid$fit + 1.96 * grid$se
  grid$target_domain_macro <- domain_name

  list(
    summary = data.frame(
      target_domain_macro = domain_name,
      n_compounds = dplyr::n_distinct(dat$inchikey),
      edf = s_table[1, "edf"],
      smooth_p = s_table[1, "p-value"],
      adjusted_r_squared = sm$r.sq,
      deviance_explained = sm$dev.expl,
      aic_gam = stats::AIC(gam_fit),
      aic_linear = stats::AIC(linear_fit),
      delta_aic_linear_minus_gam = stats::AIC(linear_fit) - stats::AIC(gam_fit),
      weighting = if (isTRUE(use_record_weights)) "sqrt_n_records" else "equal_compound"
    ),
    prediction = grid,
    model = gam_fit
  )
}

fsp3_split <- split(compound_domain_primary, compound_domain_primary$target_domain_macro)
fsp3_results <- purrr::imap(fsp3_split, fit_fsp3_domain, use_record_weights = FALSE)
fsp3_summary <- dplyr::bind_rows(purrr::map(fsp3_results, "summary")) |>
  dplyr::mutate(smooth_p_adjust_bh = stats::p.adjust(smooth_p, method = "BH"))
fsp3_predictions <- dplyr::bind_rows(purrr::map(fsp3_results, "prediction"))
write_csv_safe(fsp3_summary, file.path(TABLE_DIR, "fsp3_domain_gam_summary.csv"))
write_csv_safe(fsp3_predictions, file.path(TABLE_DIR, "fsp3_domain_gam_predictions.csv"))
for (nm in names(fsp3_results)) {
  if (!is.null(fsp3_results[[nm]]$model)) {
    safe_nm <- gsub("[^A-Za-z0-9]+", "_", nm)
    saveRDS(fsp3_results[[nm]]$model, file.path(MODEL_DIR, paste0("fsp3_gam_equal_compound_", safe_nm, ".rds")))
  }
}

fsp3_results_record_weighted <- list()
fsp3_summary_record_weighted <- data.frame()
if (isTRUE(cfg$functional_run_record_weight_sensitivity)) {
  fsp3_results_record_weighted <- purrr::imap(
    fsp3_split, fit_fsp3_domain, use_record_weights = TRUE
  )
  fsp3_summary_record_weighted <- dplyr::bind_rows(
    purrr::map(fsp3_results_record_weighted, "summary")
  ) |>
    dplyr::mutate(smooth_p_adjust_bh = stats::p.adjust(smooth_p, method = "BH"))
  write_csv_safe(
    fsp3_summary_record_weighted,
    file.path(TABLE_DIR, "fsp3_domain_gam_summary_record_weighted_sensitivity.csv")
  )
  for (nm in names(fsp3_results_record_weighted)) {
    if (!is.null(fsp3_results_record_weighted[[nm]]$model)) {
      safe_nm <- gsub("[^A-Za-z0-9]+", "_", nm)
      saveRDS(
        fsp3_results_record_weighted[[nm]]$model,
        file.path(MODEL_DIR, paste0("fsp3_gam_record_weighted_sensitivity_", safe_nm, ".rds"))
      )
    }
  }
}


message("----------------------------------------------------------------")
message(">>> [7/10] Estimating independent motif effects and combinations")
message("----------------------------------------------------------------")

fit_motif_domain <- function(dat, domain_name, use_record_weights = FALSE) {
  dat <- dat |>
    dplyr::filter(
      is.finite(z_potency), is.finite(fsp3), is.finite(logp),
      is.finite(tpsa), is.finite(mw), !is.na(class_np)
    ) |>
    dplyr::mutate(
      z_fsp3 = scale_safe(fsp3), z_logp = scale_safe(logp),
      z_tpsa = scale_safe(tpsa), z_mw = scale_safe(mw),
      model_weight = if (isTRUE(use_record_weights)) sqrt(pmax(n_records, 1)) else 1
    )

  if (nrow(dat) < cfg$motif_model_min_compounds) {
    return(list(coefficients = data.frame(), diagnostics = data.frame(
      target_domain_macro = domain_name, n_compounds = nrow(dat), status = "insufficient_n"
    ), model = NULL))
  }

  motif_terms <- c("motif_methoxy", "motif_prenyl", "motif_sugar")
  motif_terms <- motif_terms[vapply(motif_terms, function(v) {
    sum(dat[[v]] == 1, na.rm = TRUE) >= 10 && sum(dat[[v]] == 0, na.rm = TRUE) >= 10
  }, logical(1))]
  if (!length(motif_terms)) {
    return(list(coefficients = data.frame(), diagnostics = data.frame(
      target_domain_macro = domain_name, n_compounds = nrow(dat), status = "insufficient_motif_contrast"
    ), model = NULL))
  }

  class_counts <- table(dat$class_np)
  dat$class_model <- ifelse(dat$class_np %in% names(class_counts[class_counts >= 10]), dat$class_np, "Other")
  dat$class_model <- factor(dat$class_model)

  covariates <- c("z_fsp3", "z_logp", "z_tpsa", "z_mw")
  covariates <- covariates[vapply(covariates, function(v) stats::sd(dat[[v]], na.rm = TRUE) > 0, logical(1))]
  rhs <- c(motif_terms, covariates)
  if (nlevels(dat$class_model) >= 2 && nlevels(dat$class_model) <= max(8, floor(nrow(dat) / 12))) rhs <- c(rhs, "class_model")
  formula <- stats::as.formula(paste("z_potency ~", paste(rhs, collapse = " + ")))

  fit <- stats::lm(formula, data = dat, weights = model_weight)
  robust <- lmtest::coeftest(fit, vcov. = sandwich::vcovHC(fit, type = "HC3"))
  coef_df <- data.frame(
    term = rownames(robust),
    estimate = robust[, 1],
    std_error_hc3 = robust[, 2],
    statistic = robust[, 3],
    p_value = robust[, 4],
    stringsAsFactors = FALSE
  ) |>
    dplyr::filter(term %in% motif_terms) |>
    dplyr::mutate(
      target_domain_macro = domain_name,
      n_compounds = dplyr::n_distinct(dat$inchikey),
      conf_low = estimate - 1.96 * std_error_hc3,
      conf_high = estimate + 1.96 * std_error_hc3,
      motif = dplyr::recode(
        term,
        motif_methoxy = "Methoxy",
        motif_prenyl = "Prenyl-like",
        motif_sugar = "Probable sugar"
      ),
      weighting = if (isTRUE(use_record_weights)) "sqrt_n_records" else "equal_compound"
    )

  diag <- data.frame(
    target_domain_macro = domain_name,
    n_compounds = dplyr::n_distinct(dat$inchikey),
    n_rows = nrow(dat),
    n_parameters = length(stats::coef(fit)),
    r_squared = summary(fit)$r.squared,
    adjusted_r_squared = summary(fit)$adj.r.squared,
    condition_number = kappa(stats::model.matrix(fit)),
    weighting = if (isTRUE(use_record_weights)) "sqrt_n_records" else "equal_compound",
    status = "fitted",
    stringsAsFactors = FALSE
  )
  list(coefficients = coef_df, diagnostics = diag, model = fit)
}

motif_split <- split(compound_domain_primary, compound_domain_primary$target_domain_macro)
motif_results <- purrr::imap(motif_split, fit_motif_domain, use_record_weights = FALSE)
motif_coefficients <- dplyr::bind_rows(purrr::map(motif_results, "coefficients")) |>
  dplyr::mutate(p_adjust_bh = stats::p.adjust(p_value, method = "BH"))
motif_diagnostics <- dplyr::bind_rows(purrr::map(motif_results, "diagnostics"))
write_csv_safe(motif_coefficients, file.path(TABLE_DIR, "motif_independent_effects_by_domain.csv"))
write_csv_safe(motif_diagnostics, file.path(MODEL_DIR, "motif_model_diagnostics.csv"))
for (nm in names(motif_results)) {
  if (!is.null(motif_results[[nm]]$model)) {
    safe_nm <- gsub("[^A-Za-z0-9]+", "_", nm)
    saveRDS(motif_results[[nm]]$model, file.path(MODEL_DIR, paste0("motif_lm_equal_compound_", safe_nm, ".rds")))
  }
}

motif_results_record_weighted <- list()
motif_coefficients_record_weighted <- data.frame()
motif_diagnostics_record_weighted <- data.frame()
if (isTRUE(cfg$functional_run_record_weight_sensitivity)) {
  motif_results_record_weighted <- purrr::imap(
    motif_split, fit_motif_domain, use_record_weights = TRUE
  )
  motif_coefficients_record_weighted <- dplyr::bind_rows(
    purrr::map(motif_results_record_weighted, "coefficients")
  ) |>
    dplyr::mutate(p_adjust_bh = stats::p.adjust(p_value, method = "BH"))
  motif_diagnostics_record_weighted <- dplyr::bind_rows(
    purrr::map(motif_results_record_weighted, "diagnostics")
  )
  write_csv_safe(
    motif_coefficients_record_weighted,
    file.path(TABLE_DIR, "motif_independent_effects_record_weighted_sensitivity.csv")
  )
  write_csv_safe(
    motif_diagnostics_record_weighted,
    file.path(MODEL_DIR, "motif_model_diagnostics_record_weighted_sensitivity.csv")
  )
  for (nm in names(motif_results_record_weighted)) {
    if (!is.null(motif_results_record_weighted[[nm]]$model)) {
      safe_nm <- gsub("[^A-Za-z0-9]+", "_", nm)
      saveRDS(
        motif_results_record_weighted[[nm]]$model,
        file.path(MODEL_DIR, paste0("motif_lm_record_weighted_sensitivity_", safe_nm, ".rds"))
      )
    }
  }
}

motif_combination_summary <- compound_domain_primary |>
  dplyr::group_by(target_domain_macro, motif_combination) |>
  dplyr::summarise(
    n_compounds = dplyr::n_distinct(inchikey),
    n_records = sum(n_records, na.rm = TRUE),
    q25_pchembl = if (any(is.finite(median_pchembl))) {
      as.numeric(stats::quantile(median_pchembl[is.finite(median_pchembl)], 0.25, names = FALSE))
    } else NA_real_,
    q75_pchembl = if (any(is.finite(median_pchembl))) {
      as.numeric(stats::quantile(median_pchembl[is.finite(median_pchembl)], 0.75, names = FALSE))
    } else NA_real_,
    median_pchembl = median_finite_safe(median_pchembl),
    median_z_potency = median_finite_safe(z_potency),
    .groups = "drop"
  ) |>
  dplyr::group_by(target_domain_macro) |>
  dplyr::mutate(
    simple_median_pchembl = median_pchembl[motif_combination == "Simple"][1],
    simple_median_z = median_z_potency[motif_combination == "Simple"][1],
    delta_pchembl_vs_simple_descriptive = median_pchembl - simple_median_pchembl,
    delta_standardized_potency_vs_simple = median_z_potency - simple_median_z
  ) |>
  dplyr::ungroup()
write_csv_safe(motif_combination_summary, file.path(TABLE_DIR, "motif_combination_domain_summary.csv"))

binary_enriched <- domain_binary |>
  dplyr::inner_join(
    compound_master |>
      dplyr::filter(structural_set_primary) |>
      dplyr::select(inchikey, motif_combination),
    by = "inchikey"
  ) |>
  dplyr::filter(target_domain_macro != "Miscellaneous")

motif_binary_status_counts <- binary_enriched |>
  dplyr::group_by(target_domain_macro, motif_combination, binary_evidence_status) |>
  dplyr::summarise(
    n_compounds = dplyr::n_distinct(inchikey),
    n_quality_records = sum(n_quality_records, na.rm = TRUE),
    .groups = "drop"
  )
write_csv_safe(
  motif_binary_status_counts,
  file.path(TABLE_DIR, "motif_combination_binary_evidence_status_counts.csv")
)

# Primary binary analysis: only concordant compound-domain evidence.
motif_active_rates <- binary_enriched |>
  dplyr::group_by(target_domain_macro, motif_combination) |>
  dplyr::summarise(
    n_concordant_compounds = dplyr::n_distinct(
      inchikey[binary_evidence_status %in% c("active_only", "above_cutoff_only")]
    ),
    n_active = dplyr::n_distinct(inchikey[binary_evidence_status == "active_only"]),
    n_above_cutoff = dplyr::n_distinct(inchikey[binary_evidence_status == "above_cutoff_only"]),
    n_mixed_excluded = dplyr::n_distinct(inchikey[binary_evidence_status == "mixed_evidence"]),
    n_indeterminate = dplyr::n_distinct(inchikey[binary_evidence_status == "indeterminate"]),
    .groups = "drop"
  ) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    n_compounds = n_concordant_compounds,
    active_rate = dplyr::if_else(n_concordant_compounds > 0, n_active / n_concordant_compounds, NA_real_),
    ci = list(wilson_ci(n_active, n_concordant_compounds)),
    ci_low = ci[[1]],
    ci_high = ci[[2]],
    binary_rule = "concordant_only"
  ) |>
  dplyr::ungroup() |>
  dplyr::select(-ci)
write_csv_safe(
  motif_active_rates,
  file.path(TABLE_DIR, "motif_combination_active_rates_10uM.csv")
)
write_csv_safe(
  motif_active_rates,
  file.path(TABLE_DIR, "motif_combination_active_rates_10uM_concordant_only.csv")
)

# Sensitivity 1: majority of classifiable records within each compound-domain.
motif_active_rates_majority <- binary_enriched |>
  dplyr::filter(binary_majority_status %in% c("active_majority", "above_cutoff_majority")) |>
  dplyr::group_by(target_domain_macro, motif_combination) |>
  dplyr::summarise(
    n_compounds = dplyr::n_distinct(inchikey),
    n_active = dplyr::n_distinct(inchikey[binary_majority_status == "active_majority"]),
    n_above_cutoff = dplyr::n_distinct(inchikey[binary_majority_status == "above_cutoff_majority"]),
    .groups = "drop"
  ) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    active_rate = dplyr::if_else(n_compounds > 0, n_active / n_compounds, NA_real_),
    ci = list(wilson_ci(n_active, n_compounds)),
    ci_low = ci[[1]],
    ci_high = ci[[2]],
    binary_rule = "record_majority_excluding_ties"
  ) |>
  dplyr::ungroup() |>
  dplyr::select(-ci)
write_csv_safe(
  motif_active_rates_majority,
  file.path(TABLE_DIR, "motif_combination_active_rates_10uM_majority_sensitivity.csv")
)

# Sensitivity 2: legacy any-active rule retained only for transparent comparison.
motif_active_rates_legacy_any <- binary_enriched |>
  dplyr::filter(n_classifiable > 0) |>
  dplyr::group_by(target_domain_macro, motif_combination) |>
  dplyr::summarise(
    n_compounds = dplyr::n_distinct(inchikey),
    n_active = dplyr::n_distinct(inchikey[n_definite_active > 0]),
    .groups = "drop"
  ) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    active_rate = dplyr::if_else(n_compounds > 0, n_active / n_compounds, NA_real_),
    ci = list(wilson_ci(n_active, n_compounds)),
    ci_low = ci[[1]],
    ci_high = ci[[2]],
    binary_rule = "legacy_any_active_sensitivity"
  ) |>
  dplyr::ungroup() |>
  dplyr::select(-ci)
write_csv_safe(
  motif_active_rates_legacy_any,
  file.path(TABLE_DIR, "motif_combination_active_rates_10uM_legacy_any_active_sensitivity.csv")
)

message("----------------------------------------------------------------")
message(">>> [8/10] Exporting the original manuscript panels with corrected data")
message("[VERSION] Loaded Part_III_Figures_Stats_v2_9_functional_bias_controls")
message("[INFO] Figure 4 mapping: A = compound_visual_map with n >= 10; B/C = single_visual_map")
message("----------------------------------------------------------------")

pdf_device <- if (.Platform$OS.type == "windows") {
  grDevices::pdf
} else if (capabilities("cairo")) {
  grDevices::cairo_pdf
} else {
  grDevices::pdf
}
message("[INFO] PDF device: ", if (identical(pdf_device, grDevices::pdf)) "pdf" else "cairo_pdf")

save_heatmap_pdf <- function(path, heatmap_object, width, height, title = NULL) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  device_open <- FALSE
  tryCatch(
    {
      grDevices::pdf(
        file = path,
        width = width,
        height = height,
        onefile = TRUE,
        useDingbats = FALSE
      )
      device_open <- TRUE

      draw_args <- list(
        object = heatmap_object,
        heatmap_legend_side = "right",
        annotation_legend_side = "right",
        merge_legend = TRUE,
        padding = grid::unit(c(12, 10, 10, 12), "mm")
      )
      if (!is.null(title) && nzchar(title)) {
        draw_args$column_title <- title
        draw_args$column_title_gp <- grid::gpar(fontsize = 15, fontface = "bold")
      }

      do.call(ComplexHeatmap::draw, draw_args)
      grDevices::dev.off()
      device_open <- FALSE
    },
    error = function(e) {
      if (isTRUE(device_open) && grDevices::dev.cur() > 1L) {
        try(grDevices::dev.off(), silent = TRUE)
      }
      stop(
        "Failed to export heatmap PDF [", basename(path), "]: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop("Heatmap PDF was not created correctly: ", path, call. = FALSE)
  }

  invisible(path)
}


legacy_figure_outputs <- c(
  "Fig_2_Heatmap_Combined_Analysis.pdf",
  "Fig2A_Clade_Dendrogram_Family_Counts.pdf",
  "Fig2B_Physicochemical_Properties.pdf",
  "Fig2C_Structural_Features.pdf",
  "Fig2D_Chemical_Classes.pdf",
  "Fig3A_Volcano_Final_Clades.pdf",
  "Fig3B_Bias_Gradient.pdf",
  "Fig3C_Constellation_AllLabels.pdf",
  "Fig3C_Constellation_NoCorrelation.pdf",
  "Fig4A_Phylo_Heatmap.pdf",
  "Fig4B_Chemical_Drivers.pdf",
  "Fig4C_Potency_Ridges.pdf",
  "Figure4_Full_Composite_A4_Fixed.pdf",
  "Fig5_SAR_Network_Final_Colors.pdf",
  "Fig5_SAR_Network_GradientEdges.pdf"
)

family_visual_map <- clades_std |>
  dplyr::select(family, visual_clade, clade_color) |>
  dplyr::distinct(family, .keep_all = TRUE)

family_metrics_visual <- family_metrics |>
  dplyr::select(-dplyr::any_of(c("visual_clade", "clade_color"))) |>
  dplyr::left_join(family_visual_map, by = "family") |>
  dplyr::mutate(visual_clade = factor(visual_clade, levels = visual_clade_levels))

compound_visual_map <- compound_family_full |>
  dplyr::select(inchikey, family) |>
  dplyr::distinct() |>
  dplyr::left_join(family_visual_map, by = "family") |>
  dplyr::filter(!is.na(visual_clade)) |>
  dplyr::distinct(inchikey, visual_clade)

compound_visual_count <- compound_visual_map |>
  dplyr::count(inchikey, name = "n_visual_clades")

single_visual_map <- compound_visual_map |>
  dplyr::left_join(compound_visual_count, by = "inchikey") |>
  dplyr::filter(n_visual_clades == 1L) |>
  dplyr::select(inchikey, visual_clade)

visual_compound_structure <- compound_master |>
  dplyr::filter(structural_set_primary) |>
  dplyr::inner_join(single_visual_map, by = "inchikey") |>
  dplyr::mutate(visual_clade = factor(visual_clade, levels = visual_clade_levels))


clade_family_counts <- family_metrics_visual |>
  dplyr::filter(!is.na(visual_clade)) |>
  dplyr::count(visual_clade, name = "n_families")

clade_props <- visual_compound_structure |>
  dplyr::group_by(visual_clade) |>
  dplyr::summarise(
    LogP = stats::median(logp, na.rm = TRUE),
    MW = stats::median(mw, na.rm = TRUE),
    Rings = stats::median(num_rings, na.rm = TRUE),
    HBA = stats::median(num_hba, na.rm = TRUE),
    TPSA = stats::median(tpsa, na.rm = TRUE),
    HBD = stats::median(num_hbd, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::filter(dplyr::if_all(-visual_clade, is.finite))

clade_struct <- visual_compound_structure |>
  dplyr::group_by(visual_clade) |>
  dplyr::summarise(
    Prenyl = mean(dplyr::coalesce(has_prenyl_like, FALSE)),
    Methoxy = mean(dplyr::coalesce(has_methoxy_aryl, FALSE)),
    Conj_CO = mean(dplyr::coalesce(has_alpha_beta_unsaturated_carbonyl, FALSE)),
    Sugar = mean(dplyr::coalesce(has_probable_sugar, FALSE)),
    Phenol_OH = mean(dplyr::coalesce(has_phenolic_oh, FALSE)),
    .groups = "drop"
  )

clade_classes <- visual_compound_structure |>
  dplyr::filter(!is.na(class_np)) |>
  dplyr::distinct(visual_clade, inchikey, class_np) |>
  dplyr::count(visual_clade, class_np, name = "n") |>
  dplyr::group_by(visual_clade) |>
  dplyr::mutate(prop = n / sum(n)) |>
  dplyr::ungroup() |>
  dplyr::select(visual_clade, class_np, prop) |>
  tidyr::pivot_wider(names_from = class_np, values_from = prop, values_fill = 0)

common_visual <- Reduce(
  intersect,
  list(
    as.character(clade_props$visual_clade),
    as.character(clade_struct$visual_clade),
    as.character(clade_classes$visual_clade)
  )
)
common_visual <- visual_clade_levels[visual_clade_levels %in% common_visual]
if (length(common_visual) < 2L) stop("Insufficient visual clades for Figure 2.")

matrix_from_clade_table <- function(df, row_col = "visual_clade") {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  rownames(df) <- as.character(df[[row_col]])
  df[[row_col]] <- NULL
  as.matrix(df[common_visual, , drop = FALSE])
}

mat_props_raw <- matrix_from_clade_table(clade_props)
mat_struct_raw <- matrix_from_clade_table(clade_struct)
mat_classes_raw <- matrix_from_clade_table(clade_classes)

zscore_columns <- function(mat) {
  mat <- as.matrix(mat)
  out <- vapply(
    seq_len(ncol(mat)),
    function(j) scale_safe(mat[, j]),
    FUN.VALUE = numeric(nrow(mat))
  )
  if (is.vector(out)) {
    out <- matrix(out, ncol = 1)
  }
  dimnames(out) <- dimnames(mat)
  out[!is.finite(out)] <- 0
  out
}
mat_props_z <- zscore_columns(mat_props_raw)
mat_struct_z <- zscore_columns(mat_struct_raw)
mat_classes_z <- zscore_columns(mat_classes_raw)

struct_label_map <- c(
  Conj_CO = "Conjugated carbonyl",
  Prenyl = "Prenylation",
  Methoxy = "Methoxylation",
  Sugar = "Glycosylation",
  Phenol_OH = "Phenolic OH"
)
struct_labels <- unname(struct_label_map[colnames(mat_struct_z)])
struct_labels[is.na(struct_labels)] <- colnames(mat_struct_z)[is.na(struct_labels)]
colnames(mat_struct_z) <- struct_labels
colnames(mat_struct_raw) <- struct_labels

class_label_map <- c(
  "Flavandiols (Leucoanthocyanidins)" = "Flavandiols (leucoanthocyanidins)",
  "Flavan−3−ols" = "Flavan-3-ols",
  "2−arylbenzofurans" = "2-arylbenzofurans"
)
class_labels <- unname(class_label_map[colnames(mat_classes_z)])
class_labels[is.na(class_labels)] <- colnames(mat_classes_z)[is.na(class_labels)]
colnames(mat_classes_z) <- class_labels
colnames(mat_classes_raw) <- class_labels

display_clade_map <- stats::setNames(common_visual, common_visual)
display_clade_map["Caryophy_Santalales"] <- "Caryophyllales + Santalales"
display_clade_names <- unname(display_clade_map[common_visual])
rownames(mat_props_z) <- display_clade_names
rownames(mat_struct_z) <- display_clade_names
rownames(mat_classes_z) <- display_clade_names
rownames(mat_props_raw) <- display_clade_names
rownames(mat_struct_raw) <- display_clade_names
rownames(mat_classes_raw) <- display_clade_names

validate_figure2_matrix <- function(mat, label) {
  if (!is.matrix(mat) || nrow(mat) < 2L || ncol(mat) < 1L) {
    stop(
      "Figure 2 matrix [", label, "] has invalid dimensions: ",
      paste(dim(mat), collapse = " x ")
    )
  }
  if (!all(is.finite(mat))) {
    stop("Figure 2 matrix [", label, "] still contains non-finite values after standardization.")
  }
  invisible(TRUE)
}
validate_figure2_matrix(mat_props_z, "physicochemical properties")
validate_figure2_matrix(mat_struct_z, "structural features")
validate_figure2_matrix(mat_classes_z, "chemical classes")
message(
  "[INFO] Figure 2 matrices: properties ", nrow(mat_props_z), "x", ncol(mat_props_z),
  "; features ", nrow(mat_struct_z), "x", ncol(mat_struct_z),
  "; classes ", nrow(mat_classes_z), "x", ncol(mat_classes_z)
)


cluster_props <- mat_props_z / sqrt(ncol(mat_props_z))
cluster_struct <- mat_struct_z / sqrt(ncol(mat_struct_z))
cluster_classes <- mat_classes_z / sqrt(ncol(mat_classes_z))
colnames(cluster_props) <- paste0("Properties::", colnames(cluster_props))
colnames(cluster_struct) <- paste0("Features::", colnames(cluster_struct))
colnames(cluster_classes) <- paste0("Classes::", colnames(cluster_classes))
mat_integrated_cluster <- cbind(cluster_props, cluster_struct, cluster_classes)

row_hclust <- stats::hclust(
  stats::dist(mat_integrated_cluster, method = "euclidean"),
  method = "complete"
)
row_order <- rownames(mat_integrated_cluster)[row_hclust$order]

fam_counts <- clade_family_counts$n_families[
  match(common_visual, as.character(clade_family_counts$visual_clade))
]
fam_counts[is.na(fam_counts)] <- 0
names(fam_counts) <- display_clade_names

use_colors <- visual_palette[common_visual]
use_colors[is.na(use_colors)] <- "grey80"

symmetric_heatmap_limit <- function(mat, minimum = 2) {
  value <- suppressWarnings(max(abs(mat), na.rm = TRUE))
  if (!is.finite(value)) value <- minimum
  max(minimum, ceiling(value * 2) / 2)
}
struct_limit <- symmetric_heatmap_limit(mat_struct_z)
class_limit <- symmetric_heatmap_limit(mat_classes_z)

col_props <- circlize::colorRamp2(c(-2, 0, 2), c("#2166AC", "#F7F7F7", "#B2182B"))
col_struct <- circlize::colorRamp2(
  c(-struct_limit, 0, struct_limit),
  c("#D95F0E", "#F7F7F7", "#1B9E77")
)
col_classes <- circlize::colorRamp2(
  c(-class_limit, 0, class_limit),
  c("#2C7FB8", "#F7F7F7", "#88419D")
)

row_anot <- ComplexHeatmap::HeatmapAnnotation(
  `Families (n)` = ComplexHeatmap::anno_barplot(
    fam_counts,
    width = grid::unit(2.5, "cm"),
    gp = grid::gpar(fill = "grey40", col = NA),
    axis_param = list(side = "top", gp = grid::gpar(fontsize = 8))
  ),
  Lineage = factor(common_visual, levels = common_visual),
  col = list(Lineage = use_colors),
  show_legend = TRUE,
  annotation_name_side = "top",
  which = "row"
)

ht_props <- ComplexHeatmap::Heatmap(
  mat_props_z,
  name = "Properties\n(z-score)",
  col = col_props,
  cluster_rows = row_hclust,
  cluster_columns = TRUE,
  row_names_side = "left",
  row_dend_side = "left",
  row_dend_width = grid::unit(1.5, "cm"),
  rect_gp = grid::gpar(col = "white", lwd = 0.5),
  border = TRUE,
  row_names_gp = grid::gpar(fontsize = 11, fontface = "bold"),
  column_names_gp = grid::gpar(fontsize = 9),
  column_names_rot = 45,
  column_title = "Physicochemical\nProperties",
  column_title_gp = grid::gpar(fontsize = 10, fontface = "bold")
)

ht_struct <- ComplexHeatmap::Heatmap(
  mat_struct_z,
  name = "Motif prevalence\n(z-score)",
  col = col_struct,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  rect_gp = grid::gpar(col = "white", lwd = 0.5),
  border = TRUE,
  show_row_names = FALSE,
  column_names_gp = grid::gpar(fontsize = 9),
  column_names_rot = 45,
  column_title = "Structural\nFeatures",
  column_title_gp = grid::gpar(fontsize = 10, fontface = "bold")
)

ht_classes <- ComplexHeatmap::Heatmap(
  mat_classes_z,
  name = "Class prevalence\n(z-score)",
  col = col_classes,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  rect_gp = grid::gpar(col = "white", lwd = 0.5),
  border = TRUE,
  show_row_names = FALSE,
  column_names_gp = grid::gpar(fontsize = 8, fontface = "italic"),
  column_names_rot = 45,
  right_annotation = row_anot,
  column_title = "Chemical\nClasses",
  column_title_gp = grid::gpar(fontsize = 10, fontface = "bold")
)

ht_list <- ht_props + ht_struct + ht_classes

figure2_combined_paths <- c(
  file.path(FIG_DIR, "Fig_2_Heatmap_Combined_Analysis.pdf"),
  file.path(FIG_DIR, "Heatmap_Combined_Analysis.pdf")
)
for (figure2_path in figure2_combined_paths) {
  grDevices::pdf(figure2_path, width = 18, height = 9, useDingbats = FALSE)
  ComplexHeatmap::draw(
    ht_list,
    column_title = "Integrated Chemical Analysis: Properties, Features & Classes",
    column_title_gp = grid::gpar(fontsize = 16, fontface = "bold"),
    gap = grid::unit(3, "mm"),
    padding = grid::unit(c(10, 10, 10, 25), "mm"),
    merge_legend = TRUE,
    heatmap_legend_side = "right"
  )
  grDevices::dev.off()
}
message("[OK] Figure 2 exported with equally weighted integrated row clustering and zero-centered scales.")

write_csv_safe(
  data.frame(visual_clade = row_order, row_order = seq_along(row_order)),
  file.path(TABLE_DIR, "Fig2_clade_row_order.csv")
)
write_csv_safe(
  data.frame(visual_clade = rownames(mat_props_z), mat_props_z, check.names = FALSE),
  file.path(TABLE_DIR, "Fig2B_Physicochemical_Properties_matrix.csv")
)
write_csv_safe(
  data.frame(visual_clade = rownames(mat_struct_z), mat_struct_z, check.names = FALSE),
  file.path(TABLE_DIR, "Fig2C_Structural_Features_matrix.csv")
)
write_csv_safe(
  data.frame(visual_clade = rownames(mat_classes_z), mat_classes_z, check.names = FALSE),
  file.path(TABLE_DIR, "Fig2D_Chemical_Classes_matrix.csv")
)
write_csv_safe(
  data.frame(visual_clade = rownames(mat_integrated_cluster), mat_integrated_cluster, check.names = FALSE),
  file.path(TABLE_DIR, "Fig2_integrated_clustering_matrix_equal_block_weight.csv")
)
write_csv_safe(
  data.frame(
    block = c("Physicochemical properties", "Structural features", "Chemical classes"),
    n_variables = c(ncol(mat_props_z), ncol(mat_struct_z), ncol(mat_classes_z)),
    clustering_multiplier = c(
      1 / sqrt(ncol(mat_props_z)),
      1 / sqrt(ncol(mat_struct_z)),
      1 / sqrt(ncol(mat_classes_z))
    )
  ),
  file.path(TABLE_DIR, "Fig2_integrated_clustering_block_weights.csv")
)


fig3_min_compounds <- 10L
fig3_sensitivity_thresholds <- c(5L, 10L, 20L)

fig3_safe_spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 3 || dplyr::n_distinct(x) < 2 || dplyr::n_distinct(y) < 2) {
    return(list(estimate = NA_real_, p.value = NA_real_))
  }
  suppressWarnings(stats::cor.test(x, y, method = "spearman", exact = FALSE))
}

fig3_format_p <- function(p_value) {
  if (!is.finite(p_value)) return("= NA")
  if (p_value < 0.001) return("< 0.001")
  paste0("= ", format.pval(p_value, digits = 2))
}

fig3_data <- family_metrics_visual |>
  dplyr::left_join(
    family_bio_coverage |>
      dplyr::select(family, primary_fraction, active_fraction_10um),
    by = "family"
  ) |>
  dplyr::mutate(
    opportunity_gap = 1 - dplyr::coalesce(primary_fraction, 0),
    eligible_5 = n_compounds_structural >= 5,
    eligible_10 = n_compounds_structural >= 10,
    eligible_20 = n_compounds_structural >= 20,
    eligible_primary = n_compounds_structural >= fig3_min_compounds
  )

fig3_primary <- fig3_data |>
  dplyr::filter(
    eligible_primary,
    is.finite(n_compounds_structural),
    is.finite(n_scaffolds),
    n_scaffolds > 0
  )

if (nrow(fig3_primary) < 10 || dplyr::n_distinct(fig3_primary$n_compounds_structural) < 4) {
  stop("Figure 3 requires at least 10 eligible families and four distinct sampling-effort values.")
}

fig3_richness_k <- max(
  3L,
  min(
    6L,
    dplyr::n_distinct(fig3_primary$n_compounds_structural) - 1L,
    max(3L, floor(nrow(fig3_primary) / 20L))
  )
)

fig3_richness_model <- mgcv::gam(
  log1p(n_scaffolds) ~ s(log1p(n_compounds_structural), k = fig3_richness_k),
  data = fig3_primary,
  method = "REML"
)

fig3_data$fig3_expected_log_scaffolds <- as.numeric(
  stats::predict(fig3_richness_model, newdata = fig3_data, type = "response")
)
fig3_data$fig3_effort_adjusted_scaffold_richness <-
  log1p(fig3_data$n_scaffolds) - fig3_data$fig3_expected_log_scaffolds
fig3_center <- mean(
  fig3_data$fig3_effort_adjusted_scaffold_richness[fig3_data$eligible_primary],
  na.rm = TRUE
)
fig3_scale <- stats::sd(
  fig3_data$fig3_effort_adjusted_scaffold_richness[fig3_data$eligible_primary],
  na.rm = TRUE
)
if (!is.finite(fig3_scale) || fig3_scale == 0) fig3_scale <- 1
fig3_data$fig3_effort_adjusted_z <-
  (fig3_data$fig3_effort_adjusted_scaffold_richness - fig3_center) / fig3_scale

fig3_highlight_reasons <- dplyr::bind_rows(
  fig3_data |>
    dplyr::filter(eligible_primary, is.finite(fig3_effort_adjusted_z)) |>
    dplyr::slice_max(fig3_effort_adjusted_z, n = 8, with_ties = FALSE) |>
    dplyr::transmute(family, highlight_reason = "high adjusted scaffold richness"),
  fig3_data |>
    dplyr::filter(
      eligible_primary,
      fig3_effort_adjusted_z > 0,
      is.finite(exclusive_scaffold_share)
    ) |>
    dplyr::slice_max(exclusive_scaffold_share, n = 6, with_ties = FALSE) |>
    dplyr::transmute(family, highlight_reason = "high exclusive scaffold fraction"),
  fig3_data |>
    dplyr::filter(
      eligible_primary,
      fig3_effort_adjusted_z > 0,
      is.finite(opportunity_gap)
    ) |>
    dplyr::slice_max(opportunity_gap, n = 6, with_ties = FALSE) |>
    dplyr::transmute(family, highlight_reason = "high bioprospecting opportunity"),
  fig3_data |>
    dplyr::filter(eligible_primary) |>
    dplyr::slice_max(n_compounds_structural, n = 6, with_ties = FALSE) |>
    dplyr::transmute(family, highlight_reason = "high sampling effort")
) |>
  dplyr::distinct() |>
  dplyr::group_by(family) |>
  dplyr::summarise(
    highlight_reason = paste(sort(unique(highlight_reason)), collapse = "; "),
    .groups = "drop"
  )

fig3_data <- fig3_data |>
  dplyr::left_join(fig3_highlight_reasons, by = "family") |>
  dplyr::mutate(
    highlighted = !is.na(highlight_reason),
    label = ifelse(highlighted, family, NA_character_)
  )

fig3_primary_plot <- fig3_data |>
  dplyr::filter(
    eligible_primary,
    is.finite(exclusive_scaffold_share),
    is.finite(fig3_effort_adjusted_z),
    is.finite(opportunity_gap)
  )

p3a <- ggplot2::ggplot(
  fig3_primary_plot,
  ggplot2::aes(exclusive_scaffold_share, fig3_effort_adjusted_z)
) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  ggplot2::geom_vline(
    xintercept = stats::median(fig3_primary_plot$exclusive_scaffold_share, na.rm = TRUE),
    linetype = "dotted", color = "grey55"
  ) +
  ggplot2::geom_point(
    ggplot2::aes(size = n_compounds_structural, fill = opportunity_gap),
    shape = 21, color = "grey35", stroke = 0.3, alpha = 0.62
  ) +
  ggplot2::geom_point(
    data = fig3_primary_plot |> dplyr::filter(highlighted),
    ggplot2::aes(size = n_compounds_structural, fill = opportunity_gap),
    shape = 21, color = "black", stroke = 0.95, alpha = 1
  ) +
  ggrepel::geom_text_repel(
    data = fig3_primary_plot |> dplyr::filter(highlighted),
    ggplot2::aes(label = family), size = 2.7, max.overlaps = Inf,
    box.padding = 0.45, point.padding = 0.25, min.segment.length = 0,
    seed = 20260623
  ) +
  ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::scale_fill_gradientn(
    colors = c("#4575B4", "#FFFFBF", "#D73027"),
    labels = scales::percent_format(), name = "Untested rate\n(opportunity gap)"
  ) +
  ggplot2::scale_size_continuous(trans = "log10", range = c(2, 8), name = "Compounds (n)") +
  ggplot2::labs(
    title = "Scaffold Specialization and Bioprospecting Opportunity",
    subtitle = "Bemis-Murcko scaffolds; primary inference restricted to families with >=10 compounds",
    x = "Exclusive scaffold fraction",
    y = "Sampling-adjusted scaffold richness (z)"
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

ggplot2::ggsave(
  file.path(FIG_DIR, "Fig3A_Volcano_Final_Clades.pdf"),
  p3a, width = 9.5, height = 7.2, device = pdf_device
)

fit_fig3_exclusive_model <- function(dat) {
  dat <- dat |>
    dplyr::filter(
      is.finite(n_compounds_structural),
      is.finite(n_scaffolds),
      is.finite(n_exclusive_scaffolds),
      n_scaffolds > 0,
      n_exclusive_scaffolds >= 0,
      n_exclusive_scaffolds <= n_scaffolds
    )
  if (nrow(dat) < 10 || dplyr::n_distinct(dat$n_compounds_structural) < 4) return(NULL)
  k_value <- max(
    3L,
    min(
      6L,
      dplyr::n_distinct(dat$n_compounds_structural) - 1L,
      max(3L, floor(nrow(dat) / 20L))
    )
  )
  mgcv::gam(
    cbind(n_exclusive_scaffolds, n_scaffolds - n_exclusive_scaffolds) ~
      s(log1p(n_compounds_structural), k = k_value),
    family = quasibinomial(),
    data = dat,
    method = "REML"
  )
}

fig3_exclusive_model <- fit_fig3_exclusive_model(fig3_primary_plot)
if (is.null(fig3_exclusive_model)) {
  stop("The primary Figure 3B exclusive-scaffold model could not be fitted.")
}

fig3_curve_grid <- data.frame(
  n_compounds_structural = exp(seq(
    log(min(fig3_primary_plot$n_compounds_structural, na.rm = TRUE)),
    log(max(fig3_primary_plot$n_compounds_structural, na.rm = TRUE)),
    length.out = 250
  ))
)
fig3_curve_prediction <- stats::predict(
  fig3_exclusive_model,
  newdata = fig3_curve_grid,
  type = "link",
  se.fit = TRUE
)
fig3_curve_grid$expected_exclusive_fraction <- stats::plogis(fig3_curve_prediction$fit)
fig3_curve_grid$lower_95 <- stats::plogis(fig3_curve_prediction$fit - 1.96 * fig3_curve_prediction$se.fit)
fig3_curve_grid$upper_95 <- stats::plogis(fig3_curve_prediction$fit + 1.96 * fig3_curve_prediction$se.fit)

rho_exclusive_primary <- fig3_safe_spearman(
  fig3_primary_plot$n_compounds_structural,
  fig3_primary_plot$exclusive_scaffold_share
)

fig3_background_plot <- fig3_data |>
  dplyr::filter(
    !eligible_primary,
    n_compounds_structural > 0,
    is.finite(exclusive_scaffold_share)
  )

p3b <- ggplot2::ggplot() +
  ggplot2::geom_point(
    data = fig3_background_plot,
    ggplot2::aes(n_compounds_structural, exclusive_scaffold_share),
    color = "grey82", size = 1.7, alpha = 0.55
  ) +
  ggplot2::geom_ribbon(
    data = fig3_curve_grid,
    ggplot2::aes(
      x = n_compounds_structural,
      ymin = lower_95,
      ymax = upper_95
    ),
    fill = "grey75", alpha = 0.35
  ) +
  ggplot2::geom_line(
    data = fig3_curve_grid,
    ggplot2::aes(n_compounds_structural, expected_exclusive_fraction),
    color = "black", linewidth = 0.9
  ) +
  ggplot2::geom_point(
    data = fig3_primary_plot,
    ggplot2::aes(
      n_compounds_structural,
      exclusive_scaffold_share,
      size = n_compounds_structural,
      fill = exclusive_scaffold_share
    ),
    shape = 21, color = "grey35", stroke = 0.3, alpha = 0.68
  ) +
  ggplot2::geom_point(
    data = fig3_primary_plot |> dplyr::filter(highlighted),
    ggplot2::aes(
      n_compounds_structural,
      exclusive_scaffold_share,
      size = n_compounds_structural,
      fill = exclusive_scaffold_share
    ),
    shape = 21, color = "black", stroke = 0.95, alpha = 1
  ) +
  ggplot2::geom_vline(
    xintercept = c(5, 20),
    linetype = "dotted", color = "grey70", linewidth = 0.35
  ) +
  ggplot2::geom_vline(
    xintercept = fig3_min_compounds,
    linetype = "dashed", color = "black", linewidth = 0.5
  ) +
  ggrepel::geom_text_repel(
    data = fig3_primary_plot |> dplyr::filter(highlighted),
    ggplot2::aes(
      n_compounds_structural,
      exclusive_scaffold_share,
      label = family
    ),
    size = 2.6, max.overlaps = Inf, box.padding = 0.35,
    point.padding = 0.25, min.segment.length = 0, seed = 20260623
  ) +
  ggplot2::scale_x_continuous(trans = "log10") +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1.05)
  ) +
  ggplot2::scale_fill_gradientn(
    colors = c("#4575B4", "#FFFFBF", "#D73027"),
    labels = scales::percent_format(), name = "Exclusive scaffold\nfraction"
  ) +
  ggplot2::scale_size_continuous(trans = "log10", range = c(2, 7), guide = "none") +
  ggplot2::labs(
    title = "Sampling Effort versus Exclusive Scaffold Fraction",
    subtitle = paste0(
      "Primary analysis: n >=10; Spearman rho = ",
      round(unname(rho_exclusive_primary$estimate), 2),
      ", p ",
      fig3_format_p(rho_exclusive_primary$p.value),
      ". Black curve and 95% CI are fitted for n >=10."
    ),
    x = "Sampling effort (unique compounds; log scale)",
    y = "Exclusive scaffold fraction"
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

ggplot2::ggsave(
  file.path(FIG_DIR, "Fig3B_Bias_Gradient.pdf"),
  p3b, width = 9.5, height = 7.2, device = pdf_device
)

fig3_threshold_sensitivity <- purrr::map_dfr(
  fig3_sensitivity_thresholds,
  function(threshold_value) {
    dat <- fig3_data |>
      dplyr::filter(
        n_compounds_structural >= threshold_value,
        is.finite(exclusive_scaffold_share),
        n_scaffolds > 0
      )
    cor_result <- fig3_safe_spearman(
      dat$n_compounds_structural,
      dat$exclusive_scaffold_share
    )
    model_result <- fit_fig3_exclusive_model(dat)
    model_summary <- if (is.null(model_result)) NULL else summary(model_result)
    data.frame(
      minimum_compounds = threshold_value,
      n_families = nrow(dat),
      spearman_rho = unname(cor_result$estimate),
      spearman_p = cor_result$p.value,
      median_compounds = stats::median(dat$n_compounds_structural, na.rm = TRUE),
      median_exclusive_scaffold_fraction = stats::median(dat$exclusive_scaffold_share, na.rm = TRUE),
      gam_adjusted_r_squared = if (is.null(model_summary)) NA_real_ else model_summary$r.sq,
      gam_deviance_explained = if (is.null(model_summary)) NA_real_ else model_summary$dev.expl,
      stringsAsFactors = FALSE
    )
  }
)

fit_fig3_richness_threshold <- function(threshold_value) {
  dat <- fig3_data |>
    dplyr::filter(
      n_compounds_structural >= threshold_value,
      is.finite(n_scaffolds),
      n_scaffolds > 0
    )
  if (nrow(dat) < 10 || dplyr::n_distinct(dat$n_compounds_structural) < 4) return(data.frame())
  k_value <- max(
    3L,
    min(
      6L,
      dplyr::n_distinct(dat$n_compounds_structural) - 1L,
      max(3L, floor(nrow(dat) / 20L))
    )
  )
  fit <- mgcv::gam(
    log1p(n_scaffolds) ~ s(log1p(n_compounds_structural), k = k_value),
    data = dat,
    method = "REML"
  )
  dat$adjusted_residual <- log1p(dat$n_scaffolds) - as.numeric(stats::predict(fit, type = "response"))
  dat |>
    dplyr::arrange(dplyr::desc(adjusted_residual)) |>
    dplyr::mutate(
      minimum_compounds = threshold_value,
      adjusted_richness_rank = dplyr::row_number()
    ) |>
    dplyr::slice_head(n = 15) |>
    dplyr::select(
      minimum_compounds, adjusted_richness_rank, family,
      n_compounds_structural, n_scaffolds, n_exclusive_scaffolds,
      exclusive_scaffold_share, adjusted_residual
    )
}

fig3_threshold_top_families <- purrr::map_dfr(
  fig3_sensitivity_thresholds,
  fit_fig3_richness_threshold
)

write_csv_safe(
  fig3_threshold_sensitivity,
  file.path(TABLE_DIR, "Fig3_threshold_sensitivity_n5_n10_n20.csv")
)
write_csv_safe(
  fig3_threshold_top_families,
  file.path(TABLE_DIR, "Fig3_top_adjusted_families_threshold_sensitivity.csv")
)
write_csv_safe(
  fig3_data |>
    dplyr::filter(eligible_primary) |>
    dplyr::select(
      family, fine_clade, analysis_clade, n_compounds_structural,
      n_scaffolds, n_exclusive_scaffolds, exclusive_scaffold_share,
      fig3_expected_log_scaffolds, fig3_effort_adjusted_scaffold_richness,
      fig3_effort_adjusted_z, opportunity_gap, highlighted, highlight_reason
    ) |>
    dplyr::arrange(dplyr::desc(highlighted), dplyr::desc(fig3_effort_adjusted_z)),
  file.path(TABLE_DIR, "Fig3A_family_metrics_primary_n10.csv")
)
write_csv_safe(
  fig3_highlight_reasons |>
    dplyr::left_join(
      fig3_data |>
        dplyr::select(
          family, n_compounds_structural, n_scaffolds,
          n_exclusive_scaffolds, exclusive_scaffold_share,
          fig3_effort_adjusted_z, opportunity_gap
        ),
      by = "family"
    ) |>
    dplyr::arrange(dplyr::desc(fig3_effort_adjusted_z)),
  file.path(TABLE_DIR, "Fig3_highlighted_families_n10.csv")
)
saveRDS(fig3_richness_model, file.path(MODEL_DIR, "Fig3A_scaffold_richness_GAM_n10.rds"))
saveRDS(fig3_exclusive_model, file.path(MODEL_DIR, "Fig3B_exclusive_scaffold_GAM_n10.rds"))

compound_overall_potency <- compound_domain_primary |>
  dplyr::group_by(inchikey) |>
  dplyr::summarise(
    median_standardized_potency = stats::median(z_potency, na.rm = TRUE),
    median_pchembl = stats::median(median_pchembl, na.rm = TRUE),
    n_domains = dplyr::n_distinct(target_domain_macro),
    .groups = "drop"
  )

class_constellation <- compound_master |>
  dplyr::filter(structural_set_primary, !is.na(class_np)) |>
  dplyr::inner_join(compound_overall_potency, by = "inchikey") |>
  dplyr::group_by(class_np) |>
  dplyr::summarise(
    n_compounds = dplyr::n_distinct(inchikey),
    median_fsp3 = stats::median(fsp3, na.rm = TRUE),
    median_logp = stats::median(logp, na.rm = TRUE),
    median_standardized_potency = stats::median(median_standardized_potency, na.rm = TRUE),
    median_pchembl = stats::median(median_pchembl, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::filter(n_compounds >= 5, is.finite(median_fsp3), is.finite(median_logp),
                is.finite(median_standardized_potency))

potency_breaks <- stats::quantile(
  class_constellation$median_standardized_potency,
  probs = c(1 / 3, 2 / 3), na.rm = TRUE, names = FALSE
)
class_constellation <- class_constellation |>
  dplyr::mutate(
    potency_group = dplyr::case_when(
      median_standardized_potency <= potency_breaks[1] ~ "Lower reported",
      median_standardized_potency >= potency_breaks[2] ~ "Higher reported",
      TRUE ~ "Intermediate"
    ),
    potency_group = factor(potency_group, levels = c("Lower reported", "Intermediate", "Higher reported"))
  )

rho_class <- suppressWarnings(stats::cor.test(
  class_constellation$median_fsp3,
  class_constellation$median_standardized_potency,
  method = "spearman", exact = FALSE
))

make_constellation <- function(show_correlation = TRUE) {
  subtitle <- if (show_correlation) {
    paste0(
      "Class-level Spearman rho = ", round(unname(rho_class$estimate), 2),
      " (p = ", format.pval(rho_class$p.value, digits = 2),
      "); compound-level nonlinear GAMs are reported separately."
    )
  } else {
    "Higher reported relative potency occurs across both planar and more saturated flavonoid classes."
  }
  ggplot2::ggplot(
    class_constellation,
    ggplot2::aes(median_fsp3, median_logp)
  ) +
    ggplot2::annotate(
      "rect", xmin = 0.20, xmax = 0.40, ymin = -Inf, ymax = Inf,
      fill = "grey90", alpha = 0.45
    ) +
    ggplot2::geom_vline(
      xintercept = stats::median(class_constellation$median_fsp3),
      linetype = "dashed", color = "grey60"
    ) +
    ggplot2::geom_hline(
      yintercept = stats::median(class_constellation$median_logp),
      linetype = "dashed", color = "grey60"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(size = n_compounds, fill = potency_group),
      shape = 21, color = "grey20", stroke = 0.5, alpha = 0.9
    ) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = class_np), size = 3, max.overlaps = Inf,
      box.padding = 0.45, min.segment.length = 0
    ) +
    ggplot2::scale_fill_manual(
      values = c("Lower reported" = "#4575B4", "Intermediate" = "#FDAE61", "Higher reported" = "#D73027")
    ) +
    ggplot2::scale_size_continuous(trans = "sqrt", range = c(3, 10), name = "Compounds (n)") +
    ggplot2::labs(
      title = "Chemotaxonomic Constellation",
      subtitle = subtitle,
      x = "Complexity (median Fsp3)",
      y = "Lipophilicity (median LogP)",
      fill = "Reported relative potency"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}

p3c_all <- make_constellation(TRUE)
p3c_no <- make_constellation(FALSE)
ggplot2::ggsave(file.path(FIG_DIR, "Fig3C_Constellation_AllLabels.pdf"), p3c_all,
                width = 9.5, height = 7.2, device = pdf_device)
ggplot2::ggsave(file.path(FIG_DIR, "Fig3C_Constellation_NoCorrelation.pdf"), p3c_no,
                width = 9.5, height = 7.2, device = pdf_device)


# Single-lineage mapping retained exclusively for Figures 4B and 4C.
# compound_domain_primary contains repeated InChIKeys across target domains,
# whereas single_visual_map has one lineage per retained InChIKey.
compound_domain_visual_4bc <- compound_domain_primary |>
  dplyr::inner_join(
    single_visual_map,
    by = "inchikey",
    relationship = "many-to-one"
  ) |>
  dplyr::mutate(visual_clade = factor(visual_clade, levels = visual_clade_levels))

# Figure 4A atlas: compounds reported in multiple lineages are counted once
# within each lineage x target-domain combination.
compound_domain_visual_4a <- compound_domain_primary |>
  dplyr::inner_join(
    compound_visual_map,
    by = "inchikey",
    relationship = "many-to-many"
  ) |>
  dplyr::distinct(inchikey, visual_clade, target_domain_macro, .keep_all = TRUE) |>
  dplyr::mutate(visual_clade = factor(visual_clade, levels = visual_clade_levels))

visual_domain_summary <- compound_domain_visual_4a |>
  dplyr::group_by(visual_clade, target_domain_macro) |>
  dplyr::summarise(
    n_compounds = dplyr::n_distinct(inchikey),
    n_records = sum(n_records, na.rm = TRUE),
    median_z_potency = stats::median(z_potency, na.rm = TRUE),
    median_pchembl = stats::median(median_pchembl, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    # Figure 4A coverage criterion: at least 10 distinct compounds
    # per visual lineage x bioactivity domain. Record count is retained
    # for auditing, but it is not used as an additional display filter.
    sufficient_coverage = n_compounds >= cfg$functional_min_compounds_domain
  )

domain_levels_4a <- c(
  "Infectious Diseases",
  "Cell Signaling and Survival",
  "Inflammation and Proteolysis",
  "Phenotypic: Cancer and Cells",
  "Gene Regulation",
  "Metabolic Enzymes",
  "Neuroscience Targets",
  "General Assays and Safety",
  "Transporters and Ion Channels"
)

heat_grid_4a <- tidyr::expand_grid(
  visual_clade = visual_clade_levels,
  target_domain_macro = domain_levels_4a
) |>
  dplyr::left_join(
    visual_domain_summary |>
      dplyr::mutate(
        visual_clade = as.character(visual_clade),
        target_domain_macro = as.character(target_domain_macro)
      ),
    by = c("visual_clade", "target_domain_macro")
  ) |>
  dplyr::mutate(
    sufficient_coverage = dplyr::coalesce(sufficient_coverage, FALSE),
    value = dplyr::if_else(sufficient_coverage, median_z_potency, NA_real_)
  )

heat_visual <- heat_grid_4a |>
  dplyr::select(visual_clade, target_domain_macro, value) |>
  tidyr::pivot_wider(
    names_from = target_domain_macro,
    values_from = value,
    names_expand = TRUE
  ) |>
  as.data.frame(stringsAsFactors = FALSE)

rownames(heat_visual) <- heat_visual$visual_clade
heat_visual$visual_clade <- NULL
mat_heat <- as.matrix(heat_visual)
storage.mode(mat_heat) <- "double"
mat_heat[!is.finite(mat_heat)] <- NA_real_
mat_heat <- mat_heat[visual_clade_levels, domain_levels_4a, drop = FALSE]

mat_heat_cluster <- mat_heat
for (j in seq_len(ncol(mat_heat_cluster))) {
  finite_j <- is.finite(mat_heat_cluster[, j])
  fill_j <- if (any(finite_j)) stats::median(mat_heat_cluster[finite_j, j], na.rm = TRUE) else 0
  mat_heat_cluster[!finite_j, j] <- fill_j
}
mat_heat_cluster[!is.finite(mat_heat_cluster)] <- 0

safe_hclust <- function(x, margin = c("rows", "columns")) {
  margin <- match.arg(margin)
  y <- if (margin == "rows") x else t(x)
  if (nrow(y) < 2) return(FALSE)
  d <- stats::dist(y, method = "euclidean")
  if (!length(d) || any(!is.finite(d))) return(FALSE)
  tryCatch(stats::hclust(d, method = "ward.D2"), error = function(e) FALSE)
}

col_cluster_4a <- safe_hclust(mat_heat_cluster, "columns")

lineage_compound_counts_4a <- compound_domain_visual_4a |>
  dplyr::distinct(visual_clade, inchikey) |>
  dplyr::mutate(visual_clade = as.character(visual_clade)) |>
  dplyr::count(visual_clade, name = "n_primary_bioactivity_compounds")

hit_vector <- lineage_compound_counts_4a$n_primary_bioactivity_compounds[
  match(rownames(mat_heat), lineage_compound_counts_4a$visual_clade)
]
hit_vector[is.na(hit_vector)] <- 0

row_hit_annotation <- ComplexHeatmap::rowAnnotation(
  `Compounds (n)` = ComplexHeatmap::anno_barplot(
    hit_vector,
    gp = grid::gpar(fill = "grey45", col = NA),
    axis_param = list(side = "bottom", gp = grid::gpar(fontsize = 7))
  ),
  width = grid::unit(25, "mm")
)

row_labels_4a <- setNames(
  visual_clade_levels,
  visual_clade_levels
)
row_labels_4a["Caryophy_Santalales"] <- "Caryophyllales + Santalales"

ht4a <- ComplexHeatmap::Heatmap(
  mat_heat,
  name = "Standardized\npotency",
  col = circlize::colorRamp2(c(-1.5, 0, 1.5), c("#3B4CC0", "#F7F7F7", "#B40426")),
  na_col = "grey90",
  cluster_rows = FALSE,
  cluster_columns = col_cluster_4a,
  row_order = visual_clade_levels,
  row_labels = unname(row_labels_4a[rownames(mat_heat)]),
  row_names_side = "right",
  row_names_gp = grid::gpar(fontsize = 8.5),
  column_names_rot = 35,
  column_names_gp = grid::gpar(fontsize = 8),
  right_annotation = row_hit_annotation,
  border = TRUE,
  rect_gp = grid::gpar(col = "white", lwd = 0.45),
  column_title = "Functional bioactivity domains (clustered)",
  row_title = "Plant lineages"
)

utils::write.csv(
  heat_grid_4a |>
    dplyr::select(
      visual_clade, target_domain_macro, n_compounds, n_records,
      median_z_potency, sufficient_coverage, value
    ),
  file.path(TABLE_DIR, "Fig4A_lineage_domain_coverage_complete.csv"),
  row.names = FALSE
)
utils::write.csv(
  lineage_compound_counts_4a,
  file.path(TABLE_DIR, "Fig4A_unique_primary_bioactivity_compounds_by_lineage.csv"),
  row.names = FALSE
)

message("[VERSION] Part III v2.9 — functional bias controls and inline Figure 4A export")

fig4a_pdf_path <- file.path(FIG_DIR, "Fig4A_Phylo_Heatmap.pdf")
dir.create(dirname(fig4a_pdf_path), recursive = TRUE, showWarnings = FALSE)

fig4a_device_open <- FALSE
tryCatch(
  {
    grDevices::pdf(
      file = fig4a_pdf_path,
      width = 11.5,
      height = 8.4,
      onefile = TRUE,
      useDingbats = FALSE
    )
    fig4a_device_open <- TRUE

    do.call(
      ComplexHeatmap::draw,
      list(
        object = ht4a,
        heatmap_legend_side = "right",
        annotation_legend_side = "right",
        merge_legend = TRUE,
        padding = grid::unit(c(12, 10, 10, 12), "mm"),
        column_title = "A. Lineage-level Bioactivity Atlas",
        column_title_gp = grid::gpar(fontsize = 15, fontface = "bold")
      )
    )

    grDevices::dev.off()
    fig4a_device_open <- FALSE
  },
  error = function(e) {
    if (isTRUE(fig4a_device_open) && grDevices::dev.cur() > 1L) {
      try(grDevices::dev.off(), silent = TRUE)
    }
    stop(
      "Failed to export Fig4A_Phylo_Heatmap.pdf: ",
      conditionMessage(e),
      call. = FALSE
    )
  }
)

if (!file.exists(fig4a_pdf_path) || is.na(file.info(fig4a_pdf_path)$size) ||
    file.info(fig4a_pdf_path)$size <= 0) {
  stop("Figure 4A PDF was not created correctly: ", fig4a_pdf_path, call. = FALSE)
}
message("[OK] Figure 4A exported: ", fig4a_pdf_path)


class_domain_visual <- compound_domain_visual_4bc |>
  dplyr::filter(!is.na(class_np)) |>
  dplyr::group_by(class_np, target_domain_macro) |>
  dplyr::summarise(
    n_compounds = dplyr::n_distinct(inchikey),
    n_records = sum(n_records, na.rm = TRUE),
    median_z_potency = stats::median(z_potency, na.rm = TRUE),
    .groups = "drop"
  )

dominant_visual <- compound_domain_visual_4bc |>
  dplyr::filter(!is.na(class_np)) |>
  dplyr::group_by(class_np, target_domain_macro, visual_clade) |>
  dplyr::summarise(n_lineage_compounds = dplyr::n_distinct(inchikey), .groups = "drop") |>
  dplyr::arrange(class_np, target_domain_macro, dplyr::desc(n_lineage_compounds), visual_clade) |>
  dplyr::group_by(class_np, target_domain_macro) |>
  dplyr::slice(1) |>
  dplyr::ungroup() |>
  dplyr::select(class_np, target_domain_macro, source_lineage = visual_clade)

class_domain_visual <- class_domain_visual |>
  dplyr::left_join(dominant_visual, by = c("class_np", "target_domain_macro")) |>
  dplyr::filter(n_compounds >= 3)

keep_classes_4b <- class_domain_visual |>
  dplyr::group_by(class_np) |>
  dplyr::summarise(total_compounds = sum(n_compounds), .groups = "drop") |>
  dplyr::slice_max(total_compounds, n = 14, with_ties = FALSE) |>
  dplyr::pull(class_np)

p4b <- ggplot2::ggplot(
  class_domain_visual |> dplyr::filter(class_np %in% keep_classes_4b),
  ggplot2::aes(target_domain_macro, class_np)
) +
  ggplot2::geom_point(
    ggplot2::aes(size = n_compounds, color = source_lineage, alpha = median_z_potency)
  ) +
  ggplot2::scale_color_manual(values = visual_palette, drop = FALSE) +
  ggplot2::scale_size_continuous(range = c(2, 10), name = "Compounds (n)") +
  ggplot2::scale_alpha_continuous(range = c(0.45, 1), name = "Median standardized\npotency") +
  ggplot2::labs(title = "B. Chemical Correlates of Reported Potency", x = NULL, y = NULL, color = "Source lineage") +
  ggplot2::theme_classic() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 32, hjust = 1, size = 8),
    axis.text.y = ggplot2::element_text(size = 8),
    plot.title = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(file.path(FIG_DIR, "Fig4B_Chemical_Drivers.pdf"), p4b,
                width = 10.5, height = 7.2, device = pdf_device)

ridge_data <- compound_domain_visual_4bc |>
  dplyr::group_by(visual_clade) |>
  dplyr::filter(dplyr::n_distinct(inchikey) >= 10) |>
  dplyr::ungroup() |>
  dplyr::mutate(visual_clade = stats::reorder(visual_clade, z_potency, stats::median))

p4c <- ggplot2::ggplot(
  ridge_data,
  ggplot2::aes(z_potency, visual_clade, fill = visual_clade)
) +
  ggridges::geom_density_ridges(
    scale = 2, rel_min_height = 0.01, alpha = 0.82,
    color = "white", linewidth = 0.25
  ) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.55) +
  ggplot2::scale_fill_manual(values = visual_palette, guide = "none") +
  ggridges::theme_ridges() +
  ggplot2::labs(
    title = "C. Lineage-level Potency Distributions",
    subtitle = "Standardized within functional domain x measurement type; line marks zero within context",
    x = "Standardized potency", y = NULL
  ) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

ggplot2::ggsave(file.path(FIG_DIR, "Fig4C_Potency_Ridges.pdf"), p4c,
                width = 10.2, height = 6.2, device = pdf_device)


network_data <- compound_domain_primary |>
  dplyr::mutate(
    motif_count = motif_methoxy + motif_prenyl + motif_sugar,
    display_motif = dplyr::case_when(
      motif_count == 0 ~ "Simple",
      motif_count > 1 ~ "Multi-motif",
      motif_prenyl == 1 ~ "Prenyl-like",
      motif_methoxy == 1 ~ "Aryl methoxy",
      motif_sugar == 1 ~ "Sugar-like",
      TRUE ~ "Simple"
    ),
    node_label = stringr::str_squish(paste(display_motif, class_np))
  )

network_edges_complete <- network_data |>
  dplyr::group_by(node_label, target_domain_macro) |>
  dplyr::summarise(
    median_z_potency = stats::median(z_potency, na.rm = TRUE),
    median_pchembl = stats::median(median_pchembl, na.rm = TRUE),
    n_compounds = dplyr::n_distinct(inchikey),
    n_records = sum(n_records, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::filter(n_compounds >= 3, is.finite(median_z_potency)) |>
  dplyr::mutate(from = node_label, to = target_domain_macro)

# Main-text network: require at least five compounds per chemical-node x domain edge.
# The complete n >= 3 network is retained as a supplementary audit export.
network_edges <- network_edges_complete |>
  dplyr::filter(n_compounds >= 5)

network_node_compounds <- network_data |>
  dplyr::group_by(node_label, inchikey) |>
  dplyr::summarise(
    mean_z_potency = mean_finite_safe(z_potency),
    max_pchembl = max_finite_safe(median_pchembl),
    display_motif = dplyr::first(display_motif),
    class_np = dplyr::first(class_np),
    .groups = "drop"
  ) |>
  dplyr::arrange(node_label, dplyr::desc(mean_z_potency), dplyr::desc(max_pchembl))

representative_compounds <- network_node_compounds |>
  dplyr::group_by(node_label) |>
  dplyr::slice(1) |>
  dplyr::ungroup()

chemical_nodes <- network_edges |>
  dplyr::distinct(name = from) |>
  dplyr::left_join(
    representative_compounds |>
      dplyr::select(name = node_label, best_inchikey = inchikey, display_motif, class_np, lead_z = mean_z_potency, lead_pchembl = max_pchembl),
    by = "name"
  ) |>
  dplyr::mutate(type = "Chemical")

target_nodes_network <- network_edges |>
  dplyr::distinct(name = to) |>
  dplyr::mutate(
    best_inchikey = NA_character_, display_motif = "Target", class_np = NA_character_,
    lead_z = NA_real_, lead_pchembl = NA_real_, type = "Target"
  )

network_nodes <- dplyr::bind_rows(chemical_nodes, target_nodes_network) |>
  dplyr::distinct(name, .keep_all = TRUE)

network_graph <- tidygraph::tbl_graph(
  nodes = network_nodes,
  edges = network_edges |> dplyr::select(from, to, median_z_potency, median_pchembl, n_compounds, n_records),
  directed = FALSE
) |>
  tidygraph::activate(nodes) |>
  dplyr::mutate(degree = tidygraph::centrality_degree())

network_nodes_active <- network_graph |>
  tidygraph::activate(nodes) |>
  tidygraph::as_tibble()
top_network_labels <- network_nodes_active |>
  dplyr::filter(type == "Chemical") |>
  dplyr::slice_max(lead_z, n = 15, with_ties = FALSE) |>
  dplyr::pull(name)

network_graph <- network_graph |>
  tidygraph::activate(nodes) |>
  dplyr::mutate(
    final_label = dplyr::case_when(
      type == "Target" ~ name,
      name %in% top_network_labels ~ paste0(name, "\n", substr(best_inchikey, 1, 9), "..."),
      TRUE ~ NA_character_
    )
  )

set.seed(cfg$seed)
p5 <- ggraph::ggraph(network_graph, layout = "fr") +
  ggraph::geom_edge_link(
    ggplot2::aes(width = n_compounds, alpha = n_records, color = median_z_potency),
    show.legend = c(width = FALSE, alpha = FALSE, color = FALSE)
  ) +


  ggraph::scale_edge_color_gradient2(
    low = "#4575B4", mid = "grey88", high = "#D73027", midpoint = 0,
    guide = "none"
  ) +
  ggraph::scale_edge_width(range = c(0.4, 2.6)) +
  ggraph::scale_edge_alpha(range = c(0.25, 0.85)) +
  ggraph::geom_node_point(
    ggplot2::aes(filter = type == "Target", size = degree),
    shape = 22, fill = "#26334D", color = "white", stroke = 1
  ) +
  ggraph::geom_node_point(
    ggplot2::aes(filter = type == "Chemical", fill = lead_z, shape = display_motif, size = degree),
    color = "grey15", stroke = 0.5, alpha = 0.92
  ) +
  ggraph::geom_node_text(
    ggplot2::aes(label = final_label, filter = type == "Target"),
    repel = TRUE, fontface = "bold", size = 3.5
  ) +
  ggraph::geom_node_text(
    ggplot2::aes(label = final_label, filter = type == "Chemical" & !is.na(final_label)),
    repel = TRUE, fontface = "bold", color = "red4", size = 2.6,
    min.segment.length = 0
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#4575B4", mid = "#FDAE61", high = "#D73027", midpoint = 0,
    name = "Representative\nstandardized potency", na.value = "grey90"
  ) +
  ggplot2::scale_shape_manual(values = c(
    "Simple" = 21, "Prenyl-like" = 24, "Aryl methoxy" = 23,
    "Sugar-like" = 25, "Multi-motif" = 22, "Target" = 22
  )) +
  ggplot2::scale_size_continuous(range = c(3, 11), guide = "none") +
  ggplot2::theme_void() +
  ggplot2::labs(title = "Structure-Bioactivity Association Network") +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 15), legend.position = "right")

ggplot2::ggsave(file.path(FIG_DIR, "Fig5_SAR_Network_Final_Colors.pdf"), p5,
                width = 15, height = 12, device = pdf_device)
ggplot2::ggsave(file.path(FIG_DIR, "Fig5_SAR_Network_GradientEdges.pdf"), p5,
                width = 15, height = 12, device = pdf_device)
message("[OK] Figure 5 network exported; edge colours retain standardized potency and the edge legend is left for Affinity.")

lead_metadata <- representative_compounds |>
  dplyr::filter(node_label %in% top_network_labels) |>
  dplyr::left_join(
    compound_master |>
      dplyr::select(inchikey, dplyr::any_of(c("canonical_smiles", "standardized_smiles", "representative_input_smiles"))),
    by = "inchikey"
  )
write_csv_safe(lead_metadata, file.path(PARTIII_DIR, "Fig5_Top_Leads_Metadata.csv"))
write_csv_safe(network_edges, file.path(PARTIII_DIR, "Fig5_Network_Edges_Corrected.csv"))
write_csv_safe(
  network_edges_complete,
  file.path(PARTIII_DIR, "Fig5_Network_Edges_Complete_n3_Supplement.csv")
)
write_csv_safe(network_nodes_active, file.path(PARTIII_DIR, "Fig5_Network_Nodes_Corrected.csv"))


motif_association_compat <- motif_combination_summary |>
  dplyr::filter(motif_combination %in% c("Simple", "Methoxy", "Prenyl", "Sugar")) |>
  dplyr::mutate(
    Motif_Type = dplyr::recode(
      motif_combination,
      Simple = "Simple", Methoxy = "Aryl methoxy",
      Prenyl = "Prenyl-like", Sugar = "Sugar-like"
    ),
    Functional_Domain = target_domain_macro,
    N_Compounds = n_compounds,
    Median_pChEMBL_Descriptive = median_pchembl,
    Standardized_Median = median_z_potency,
    Delta_Standardized_vs_Simple = delta_standardized_potency_vs_simple,
    Delta_pChEMBL_vs_Simple_Descriptive = delta_pchembl_vs_simple_descriptive
  ) |>
  dplyr::select(
    Functional_Domain, Motif_Type, N_Compounds,
    Median_pChEMBL_Descriptive, Standardized_Median,
    Delta_Standardized_vs_Simple,
    Delta_pChEMBL_vs_Simple_Descriptive
  )
write_csv_safe(
  motif_association_compat,
  file.path(PARTIII_DIR, "Table_Statistics_Motif_Bioactivity_Associations.csv")
)
# Legacy filename retained for pipeline compatibility; content now uses non-causal labels.
write_csv_safe(
  motif_association_compat,
  file.path(PARTIII_DIR, "Table_Statistics_Potency_Gain.csv")
)
write_csv_safe(motif_coefficients, file.path(PARTIII_DIR, "Table_Significance_Tests.csv"))


boxplot_data <- network_data |>
  dplyr::mutate(
    display_motif = factor(
      display_motif,
      levels = c("Simple", "Sugar-like", "Aryl methoxy", "Prenyl-like", "Multi-motif")
    )
  )
p_box <- ggplot2::ggplot(boxplot_data, ggplot2::aes(display_motif, z_potency, fill = display_motif)) +
  ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  ggplot2::geom_jitter(width = 0.18, alpha = 0.15, size = 0.45) +
  ggplot2::facet_wrap(~ target_domain_macro, scales = "free_y") +
  ggplot2::scale_fill_manual(values = c(
    "Simple" = "grey75", "Sugar-like" = "#440154",
    "Aryl methoxy" = "#21908C", "Prenyl-like" = "#FDE725",
    "Multi-motif" = "#D95F02"
  )) +
  ggplot2::labs(
    title = "Descriptive potency distributions by predominant motif category",
    subtitle = "Inferential estimates use independent motif flags and robust multivariable models",
    x = NULL, y = "Standardized potency"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1), legend.position = "none")
ggplot2::ggsave(file.path(PARTIII_DIR, "Fig_Statistical_Boxplot.pdf"), p_box,
                width = 12, height = 8.5, device = pdf_device)


message("----------------------------------------------------------------")
message(">>> [9/10] Writing Figure 1 and Supplementary iTOL datasets")
message("----------------------------------------------------------------")

itol_write_heatmap <- function(df, fields, path, label, color_min = "#FFFFFF", color_max = "#000000", color_nan = "#E6E6E6", color_mid = NULL, user_min = NULL, user_mid = NULL, user_max = NULL) {
  required_columns <- c("family", fields)
  missing_columns <- setdiff(required_columns, names(df))
  if (length(missing_columns) > 0L) stop("Missing iTOL heatmap columns: ", paste(missing_columns, collapse = ", "))
  export_df <- as.data.frame(df[, required_columns, drop = FALSE], stringsAsFactors = FALSE)
  export_df$family <- as.character(export_df$family)
  if (anyNA(export_df$family) || any(export_df$family == "")) stop("Empty family identifiers in iTOL heatmap: ", label)
  if (anyDuplicated(export_df$family)) stop("Duplicated family identifiers in iTOL heatmap: ", label)
  for (field in fields) {
    values <- suppressWarnings(as.numeric(export_df[[field]]))
    export_df[[field]] <- ifelse(
      is.na(values) | !is.finite(values),
      "X",
      as.character(signif(values, 15))
    )
  }
  gradient_options <- c(
    paste("COLOR_MIN", color_min, sep = "\t"),
    paste("COLOR_MAX", color_max, sep = "\t")
  )
  if (!is.null(color_mid)) {
    gradient_options <- c(
      gradient_options,
      "USE_MID_COLOR\t1",
      paste("COLOR_MID", color_mid, sep = "\t")
    )
  }
  if (!is.null(user_min)) gradient_options <- c(gradient_options, paste("USER_MIN_VALUE", user_min, sep = "\t"))
  if (!is.null(user_mid)) gradient_options <- c(gradient_options, paste("USER_MID_VALUE", user_mid, sep = "\t"))
  if (!is.null(user_max)) gradient_options <- c(gradient_options, paste("USER_MAX_VALUE", user_max, sep = "\t"))
  header <- c(
    "DATASET_HEATMAP", "SEPARATOR TAB",
    paste("DATASET_LABEL", label, sep = "\t"),
    "COLOR\t#000000",
    paste(c("FIELD_LABELS", fields), collapse = "\t"),
    gradient_options,
    paste("COLOR_NAN", color_nan, sep = "\t"),
    "DATA"
  )
  body <- vapply(seq_len(nrow(export_df)), function(i) paste(unlist(export_df[i, ], use.names = FALSE), collapse = "\t"), character(1))
  observed_fields <- lengths(strsplit(body, "\t", fixed = TRUE)) - 1L
  if (any(observed_fields != length(fields))) {
    bad <- which(observed_fields != length(fields))
    stop("Invalid iTOL heatmap row width in ", label, ": ", paste(export_df$family[bad], collapse = ", "))
  }
  writeLines(c(header, body), path, useBytes = TRUE)
}

itol_write_binary <- function(df, fields, colors, shapes, path, label) {
  header <- c(
    "DATASET_BINARY", "SEPARATOR TAB",
    paste("DATASET_LABEL", label, sep = "\t"),
    "COLOR\t#000000",
    paste(c("FIELD_SHAPES", shapes), collapse = "\t"),
    paste(c("FIELD_LABELS", fields), collapse = "\t"),
    paste(c("FIELD_COLORS", colors), collapse = "\t"),
    "DATA"
  )
  export_df <- df[, c("family", fields), drop = FALSE]
  export_df[] <- lapply(export_df, function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x
  })
  body <- apply(export_df, 1, function(z) paste(z, collapse = "\t"))
  writeLines(c(header, body), path, useBytes = TRUE)
}


class_wide <- family_class_long |>
  dplyr::select(family, class_np, class_prevalence) |>
  tidyr::pivot_wider(names_from = class_np, values_from = class_prevalence, values_fill = 0) |>
  dplyr::filter(family %in% family_visual_map$family)
class_fields <- setdiff(names(class_wide), "family")
if (!is.null(cfg$flav_classes)) {
  preferred <- cfg$flav_classes[cfg$flav_classes %in% class_fields]
  class_fields <- c(preferred, setdiff(class_fields, preferred))
}
itol_write_heatmap(
  class_wide, class_fields,
  file.path(ITOL_DIR, "Fig1_iTOL_01_Flavonoid_Subclasses_Heatmap.txt"),
  "Flavonoid subclasses", "#FFFFFF", "#000000"
)


motif_fig1 <- family_metrics_visual |>
  dplyr::filter(!is.na(visual_clade)) |>
  dplyr::transmute(
    family,
    Glycosylation = ifelse(prevalence_sugar > 0, 1, -1),
    Prenylation = ifelse(prevalence_prenyl > 0, 1, -1),
    Methoxylation = ifelse(prevalence_methoxy > 0, 1, -1),
    Phenolic_OH = ifelse(prevalence_phenolic > 0, 1, -1)
  )
itol_write_binary(
  motif_fig1,
  c("Glycosylation", "Prenylation", "Methoxylation", "Phenolic_OH"),
  c("#6A3D9A", "#E31A1C", "#FF7F00", "#1B9E77"),
  c(3, 4, 2, 1),
  file.path(ITOL_DIR, "Fig1_iTOL_02_Structural_Motifs_Symbols.txt"),
  "Structural motifs"
)


family_context_evidence <- compound_family_primary |>
  dplyr::select(family, inchikey) |>
  dplyr::distinct() |>
  dplyr::inner_join(
    target_primary |>
      dplyr::mutate(
        n_records = safe_numeric(n_records),
        median_pchembl = safe_numeric(median_pchembl),
        functional_context_id = stringr::str_c(
          dplyr::coalesce(target_domain_macro, "Unclassified"),
          dplyr::coalesce(target_chembl_id, "NoTargetID"),
          dplyr::coalesce(target_pref_name_final, "NoTargetName"),
          dplyr::coalesce(target_type_final, "NoTargetType"),
          sep = " | "
        )
      ) |>
      dplyr::filter(target_domain_macro != "Miscellaneous"),
    by = "inchikey",
    relationship = "many-to-many"
  ) |>
  dplyr::group_by(family) |>
  dplyr::summarise(
    n_functional_contexts = dplyr::n_distinct(functional_context_id),
    reported_high_potency = any(median_pchembl >= 7, na.rm = TRUE),
    .groups = "drop"
  )

functional_fig1 <- family_bio_coverage |>
  dplyr::left_join(family_context_evidence, by = "family") |>
  dplyr::transmute(
    family,
    Concordant_active = ifelse(n_active_only_compounds_10um > 0, 1, -1),
    Multi_context = ifelse(dplyr::coalesce(n_functional_contexts, 0L) >= 2, 1, -1),
    Reported_high_potency = ifelse(dplyr::coalesce(reported_high_potency, FALSE), 1, -1)
  )
itol_write_binary(
  functional_fig1,
  c("Concordant_active", "Multi_context", "Reported_high_potency"),
  c("#2CA25F", "#3182BD", "#CB181D"),
  c(1, 2, 3),
  file.path(ITOL_DIR, "Fig1_iTOL_03_Functional_Annotation.txt"),
  "Functional annotation"
)


s1_fields <- c(
  "Infectious Diseases",
  "Cell Signaling and Survival",
  "Inflammation and Proteolysis",
  "Phenotypic: Cancer and Cells",
  "Gene Regulation",
  "Metabolic Enzymes",
  "Neuroscience Targets",
  "General Assays and Safety",
  "Transporters and Ion Channels"
)
family_domain_s1 <- compound_family_primary |>
  dplyr::select(family, inchikey) |>
  dplyr::distinct() |>
  dplyr::inner_join(
    compound_domain_primary,
    by = "inchikey",
    relationship = "many-to-many"
  ) |>
  dplyr::filter(target_domain_macro %in% s1_fields) |>
  dplyr::group_by(family, target_domain_macro) |>
  dplyr::summarise(
    n_compounds = dplyr::n_distinct(inchikey),
    median_z = ifelse(any(is.finite(z_potency)), stats::median(z_potency[is.finite(z_potency)], na.rm = TRUE), NA_real_),
    .groups = "drop"
  ) |>
  dplyr::mutate(value = ifelse(n_compounds >= 1, median_z, NA_real_)) |>
  dplyr::select(family, target_domain_macro, value) |>
  tidyr::pivot_wider(names_from = target_domain_macro, values_from = value, values_fill = NA_real_)
for (field in s1_fields) {
  if (!field %in% names(family_domain_s1)) family_domain_s1[[field]] <- NA_real_
}
family_domain_s1 <- family_visual_map |>
  dplyr::select(family) |>
  dplyr::distinct() |>
  dplyr::left_join(family_domain_s1, by = "family") |>
  dplyr::select(family, dplyr::all_of(s1_fields))
itol_write_heatmap(
  family_domain_s1, s1_fields,
  file.path(ITOL_DIR, "FigS1_iTOL_Functional_Bioactivity_Domains_Heatmap.txt"),
  "Functional bioactivity domains", "#F7FBFF", "#67000D", "#E6E6E6"
)


top30_categories <- target_primary |>
  dplyr::filter(
    target_domain_macro != "Miscellaneous",
    !is.na(target_category_l3),
    target_category_l3 != ""
  ) |>
  dplyr::group_by(target_category_l3) |>
  dplyr::summarise(n_compounds = dplyr::n_distinct(inchikey), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(n_compounds), target_category_l3) |>
  dplyr::slice_head(n = 30) |>
  dplyr::pull(target_category_l3)

if (length(top30_categories) != 30L) stop("Figure S2 requires exactly 30 target categories; found ", length(top30_categories))

s2_family_universe <- family_visual_map |>
  dplyr::filter(!is.na(visual_clade)) |>
  dplyr::select(family) |>
  dplyr::distinct() |>
  dplyr::arrange(family)

s2_primary_evidence <- compound_family_primary |>
  dplyr::select(family, inchikey) |>
  dplyr::distinct() |>
  dplyr::inner_join(
    target_primary |>
      dplyr::filter(
        target_domain_macro != "Miscellaneous",
        !is.na(target_category_l3),
        target_category_l3 != ""
      ) |>
      dplyr::select(inchikey, target_category_l3),
    by = "inchikey",
    relationship = "many-to-many"
  ) |>
  dplyr::distinct(family, inchikey, target_category_l3)

s2_family_denominator <- s2_primary_evidence |>
  dplyr::group_by(family) |>
  dplyr::summarise(
    n_primary_bioactivity_compounds = dplyr::n_distinct(inchikey),
    .groups = "drop"
  )

s2_category_counts <- s2_primary_evidence |>
  dplyr::filter(target_category_l3 %in% top30_categories) |>
  dplyr::group_by(family, target_category_l3) |>
  dplyr::summarise(n_category_compounds = dplyr::n_distinct(inchikey), .groups = "drop")

family_target_s2_long <- tidyr::expand_grid(
  family = s2_family_universe$family,
  target_category_l3 = top30_categories
) |>
  dplyr::left_join(s2_category_counts, by = c("family", "target_category_l3")) |>
  dplyr::left_join(s2_family_denominator, by = "family") |>
  dplyr::mutate(
    value = dplyr::case_when(
      is.na(n_primary_bioactivity_compounds) | n_primary_bioactivity_compounds == 0L ~ NA_real_,
      TRUE ~ dplyr::coalesce(n_category_compounds, 0L) / n_primary_bioactivity_compounds
    ),
    target_category_l3 = factor(target_category_l3, levels = top30_categories)
  )

family_target_s2 <- family_target_s2_long |>
  dplyr::select(family, target_category_l3, value) |>
  tidyr::pivot_wider(
    names_from = target_category_l3,
    values_from = value,
    names_expand = TRUE
  ) |>
  dplyr::select(family, dplyr::all_of(top30_categories)) |>
  dplyr::arrange(family)

s2_fields <- top30_categories

if (nrow(family_target_s2) != nrow(s2_family_universe)) stop("Figure S2 family count mismatch: ", nrow(family_target_s2), " versus ", nrow(s2_family_universe))
if (ncol(family_target_s2) != 31L) stop("Figure S2 must contain family plus 30 target-category columns")
if (!identical(family_target_s2$family, s2_family_universe$family)) stop("Figure S2 family identifiers do not match the analytical tree universe")
if (!identical(names(family_target_s2)[-1], top30_categories)) stop("Figure S2 target-category columns are not in the expected order")

s2_family_audit <- s2_family_universe |>
  dplyr::left_join(s2_family_denominator, by = "family") |>
  dplyr::mutate(
    n_primary_bioactivity_compounds = dplyr::coalesce(n_primary_bioactivity_compounds, 0L),
    export_status = ifelse(n_primary_bioactivity_compounds > 0L, "evaluated", "no_primary_bioactivity_data")
  )

write_csv_safe(s2_family_audit, file.path(AUDIT_DIR, "FigS2_family_bioactivity_coverage.csv"))
write_csv_safe(family_target_s2, file.path(ITOL_DIR, "FigS2_target_category_prevalence_matrix.csv"))

itol_write_heatmap(
  family_target_s2, s2_fields,
  file.path(ITOL_DIR, "FigS2_iTOL_Top30_Target_Categories_Heatmap.txt"),
  "Target-category prevalence", "#FFFFFF", "#000000", "#F2E6C9"
)


s3_family_universe <- family_visual_map |>
  dplyr::filter(!is.na(visual_clade)) |>
  dplyr::select(family) |>
  dplyr::distinct() |>
  dplyr::arrange(family)

s3_motifs <- s3_family_universe |>
  dplyr::left_join(
    family_metrics_visual |>
      dplyr::select(
        family,
        prevalence_prenyl,
        prevalence_methoxy,
        prevalence_sugar,
        prevalence_phenolic
      ),
    by = "family"
  ) |>
  dplyr::transmute(
    family,
    Prenylation = prevalence_prenyl,
    Methoxylation = prevalence_methoxy,
    Glycosylation = prevalence_sugar,
    `Phenolic OH` = prevalence_phenolic
  )

s3_physicochemical <- s3_family_universe |>
  dplyr::left_join(
    family_metrics_visual |>
      dplyr::select(family, median_mw, median_logp, median_tpsa),
    by = "family"
  ) |>
  dplyr::mutate(
    mw_z = pmax(-2, pmin(2, scale_safe(median_mw))),
    logp_z = pmax(-2, pmin(2, scale_safe(median_logp))),
    tpsa_z = pmax(-2, pmin(2, scale_safe(median_tpsa)))
  ) |>
  dplyr::transmute(
    family,
    `MW (z)` = mw_z,
    `LogP (z)` = logp_z,
    `TPSA (z)` = tpsa_z
  )

if (!identical(s3_motifs$family, s3_family_universe$family)) stop("Figure S3 motif family identifiers do not match the analytical tree universe")
if (!identical(s3_physicochemical$family, s3_family_universe$family)) stop("Figure S3 physicochemical family identifiers do not match the analytical tree universe")

write_csv_safe(s3_motifs, file.path(ITOL_DIR, "FigS3_structural_motif_prevalence_matrix.csv"))
write_csv_safe(s3_physicochemical, file.path(ITOL_DIR, "FigS3_standardized_physicochemical_matrix.csv"))

itol_write_heatmap(
  s3_motifs,
  c("Prenylation", "Methoxylation", "Glycosylation", "Phenolic OH"),
  file.path(ITOL_DIR, "FigS3_iTOL_01_Structural_Motif_Prevalence.txt"),
  "Structural motif prevalence",
  "#FFFFFF",
  "#000000",
  "#F2D9A6",
  user_min = 0,
  user_max = 1
)

itol_write_heatmap(
  s3_physicochemical,
  c("MW (z)", "LogP (z)", "TPSA (z)"),
  file.path(ITOL_DIR, "FigS3_iTOL_02_Standardized_Physicochemical_Properties.txt"),
  "Standardized physicochemical properties",
  "#2166AC",
  "#B2182B",
  "#F2D9A6",
  color_mid = "#FFFFFF",
  user_min = -2,
  user_mid = 0,
  user_max = 2
)

itol_manifest <- data.frame(
  figure = c("Fig1", "Fig1", "Fig1", "FigS1", "FigS2", "FigS3", "FigS3"),
  upload_order = c(1, 2, 3, 1, 1, 1, 2),
  file = c(
    "Fig1_iTOL_01_Flavonoid_Subclasses_Heatmap.txt",
    "Fig1_iTOL_02_Structural_Motifs_Symbols.txt",
    "Fig1_iTOL_03_Functional_Annotation.txt",
    "FigS1_iTOL_Functional_Bioactivity_Domains_Heatmap.txt",
    "FigS2_iTOL_Top30_Target_Categories_Heatmap.txt",
    "FigS3_iTOL_01_Structural_Motif_Prevalence.txt",
    "FigS3_iTOL_02_Standardized_Physicochemical_Properties.txt"
  ),
  description = c(
    "17 NPClassifier subclasses; family prevalence 0-1",
    "Independent presence of glycosylation, prenylation, methoxylation, and phenolic OH",
    "Concordant active evidence, multiple functional contexts, and reported high-potency family flags",
    "Median standardized potency for nine functional bioactivity domains",
    "Relative family coverage of the 30 most frequent target categories",
    "Family prevalence of prenylation, methoxylation, glycosylation, and phenolic OH",
    "Family median MW, LogP, and TPSA standardized across families and clipped to -2 to 2 for display"
  )
)
write_csv_safe(itol_manifest, file.path(ITOL_DIR, "iTOL_UPLOAD_MANIFEST.csv"))

write_csv_safe(
  data.frame(file = legacy_figure_outputs, stringsAsFactors = FALSE),
  file.path(AUDIT_DIR, "legacy_figure_output_manifest.csv")
)


message("----------------------------------------------------------------")
message(">>> [10/10] Writing master workbook and reproducibility metadata")
message("----------------------------------------------------------------")

workbook_path <- file.path(PARTIII_DIR, "PartIII_Results_Master.xlsx")
wb <- openxlsx::createWorkbook()
add_sheet <- function(name, df) {
  name <- substr(gsub("[^A-Za-z0-9_]+", "_", name), 1, 31)
  df <- as_base_df(df)
  if (ncol(df) == 0L) df <- data.frame(note = "No records generated for this table", stringsAsFactors = FALSE)
  openxlsx::addWorksheet(wb, name)
  openxlsx::writeDataTable(wb, name, df, withFilter = TRUE)
  openxlsx::freezePane(wb, name, firstRow = TRUE)
}
add_sheet("family_metrics", family_metrics)
add_sheet("family_class_profile", family_class_long)
add_sheet("univariate_global", univ_global)
add_sheet("univariate_dunn", univ_posthoc)
add_sheet("permanova_global", multi_global)
add_sheet("dispersion", multi_disp)
add_sheet("pairwise_permanova", multi_pair)
add_sheet("sensitivity_permanova", sensitivity_scaffold)
add_sheet("family_bio_coverage", family_bio_coverage)
add_sheet("clade_bio_coverage", clade_bio_coverage)
add_sheet("clade_domain", clade_domain_summary)
add_sheet("fsp3_gam", fsp3_summary)
add_sheet("fsp3_weight_sens", fsp3_summary_record_weighted)
add_sheet("motif_effects", motif_coefficients)
add_sheet("motif_model_diag", motif_diagnostics)
add_sheet("motif_weight_sens", motif_coefficients_record_weighted)
add_sheet("motif_combinations", motif_combination_summary)
add_sheet("binary_status", binary_evidence_status_summary)
add_sheet("binary_by_domain", binary_evidence_domain_summary)
add_sheet("motif_active_rates", motif_active_rates)
add_sheet("motif_active_majority", motif_active_rates_majority)
add_sheet("motif_active_legacy", motif_active_rates_legacy_any)
add_sheet("potency_agg_sens", compound_domain_primary |>
  dplyr::select(
    inchikey, target_domain_macro, n_standard_types, n_records,
    z_potency_equal_mean, z_potency_equal_median, z_potency_record_weighted,
    delta_equal_mean_minus_record_weighted
  ))
add_sheet("input_audit", input_audit)
add_sheet("taxonomic_corrections", correction_audit)
openxlsx::saveWorkbook(wb, workbook_path, overwrite = TRUE)

metadata <- list(
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  base_tag = base_tag,
  out_dir = OUT_DIR,
  part3_dir = PARTIII_DIR,
  occurrence_source = occurrence_source,
  input_paths = input_paths,
  primary_analysis = list(
    taxonomic_scope = "vascular families represented in the validated family tree",
    structural_scope = "primary_structure_eligible compounds",
    minimum_compounds_per_family = cfg$analysis_min_compounds_per_taxon,
    figure3_minimum_compounds_per_family = fig3_min_compounds,
    figure3_sensitivity_thresholds = fig3_sensitivity_thresholds,
    scaffold_definition = "RDKit Bemis-Murcko scaffold",
    novelty_inference = "rarefaction plus residual from GAM scaffold-richness expectation; descriptive exclusive scaffolds",
    functional_layer = "primary pChEMBL standardized within functional-domain x standard-type and aggregated with equal weight across measurement types",
    functional_lineage_rule = "single broad analytical clade per compound",
    functional_binary_rule = "concordant-only active versus above-cutoff evidence; mixed evidence excluded from the primary rate analysis",
    functional_model_weighting = "equal weight per compound; sqrt(n_records) retained only as sensitivity",
    permutations = cfg$permutations,
    figure_design = "legacy manuscript layout preserved; corrected data, equal-compound weighting, and mixed-evidence controls",
    visual_macroclades = visual_clade_levels
  ),
  counts = list(
    occurrence_rows = nrow(occurrence_std),
    compounds_full = length(unique(full_keys)),
    compounds_primary = length(unique(primary_keys)),
    compounds_strict = length(unique(strict_keys)),
    families_primary_structural = dplyr::n_distinct(compound_family_primary$family),
    families_primary_threshold = length(eligible_families_primary),
    primary_compound_domain_rows = nrow(compound_domain_primary),
    primary_compounds_with_function = dplyr::n_distinct(compound_domain_primary$inchikey),
    single_clade_functional_compounds = dplyr::n_distinct(compound_domain_clade_single$inchikey),
    taxonomic_correction_rows = sum(correction_mask)
  ),
  software = list(
    R = R.version.string,
    dplyr = as.character(utils::packageVersion("dplyr")),
    vegan = as.character(utils::packageVersion("vegan")),
    mgcv = as.character(utils::packageVersion("mgcv")),
    ggplot2 = as.character(utils::packageVersion("ggplot2")),
    openxlsx = as.character(utils::packageVersion("openxlsx"))
  )
)
write_json_safe(metadata, file.path(AUDIT_DIR, "part3_run_metadata.json"))

capture.output(sessionInfo(), file = file.path(AUDIT_DIR, "sessionInfo.txt"))
saveRDS(
  list(
    family_metrics = family_metrics,
    family_class_long = family_class_long,
    compound_domain_primary = compound_domain_primary,
    compound_binary_status = compound_binary_status,
    binary_evidence_status_summary = binary_evidence_status_summary,
    clade_domain_summary = clade_domain_summary,
    motif_coefficients = motif_coefficients,
    motif_coefficients_record_weighted = motif_coefficients_record_weighted,
    motif_active_rates = motif_active_rates,
    motif_active_rates_majority = motif_active_rates_majority,
    fsp3_summary = fsp3_summary,
    fsp3_summary_record_weighted = fsp3_summary_record_weighted
  ),
  file.path(PARTIII_DIR, "PartIII_analysis_objects.rds")
)

message("------------------------------------------------------------")
message("[PART III SUMMARY]")
message("Occurrence rows: ", nrow(occurrence_std))
message("Full / primary / strict compounds: ", length(unique(full_keys)), " / ", length(unique(primary_keys)), " / ", length(unique(strict_keys)))
message("Families in primary structural data: ", dplyr::n_distinct(compound_family_primary$family))
message("Families meeting primary threshold: ", length(eligible_families_primary))
message("Primary compound-domain functional rows: ", nrow(compound_domain_primary))
message("Master workbook: ", workbook_path)
message("Legacy-preserving figure panels: ", FIG_DIR)
message("iTOL tracks: ", ITOL_DIR)
message("Metadata: ", file.path(AUDIT_DIR, "part3_run_metadata.json"))
message("------------------------------------------------------------")
