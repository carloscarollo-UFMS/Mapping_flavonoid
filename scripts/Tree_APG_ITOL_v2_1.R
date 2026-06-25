# Tree_APG_ITOL_v2_4: robust primary occurrence input selection
# ============================================================================
# Updated from Tree_APG_ITOL_v2_3_chloranthales_fix.R
#
# Main changes:
#   1. Load only the primary *_lin_enriched.parquet occurrence table; Part II
#      summary workbooks are never selected automatically.
#   2. Prefer the current pipeline run and exact base-tag filename using a
#      deterministic rule; tied candidates require an explicit cfg path.
#   3. Preserve the explicit Chloranthales classification and clade guards.
#   4. Validate the selected occurrence table before tree construction.
# ============================================================================


`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0L) return(b)
  if (is.character(a) && length(a) == 1L && (is.na(a) || !nzchar(a))) return(b)
  a
}

as_base_df <- function(x) as.data.frame(x, stringsAsFactors = FALSE)

normalize_path <- function(x) {
  normalizePath(x, winslash = "/", mustWork = FALSE)
}

trim_chr <- function(x) {
  x <- trimws(as.character(x))
  x[!nzchar(x)] <- NA_character_
  x
}

squish_chr <- function(x) {
  x <- gsub("[[:space:]]+", " ", trim_chr(x))
  x[!nzchar(x)] <- NA_character_
  x
}

normalize_taxon_word <- function(x) {
  x <- squish_chr(x)
  keep <- !is.na(x)
  x[keep] <- paste0(
    toupper(substr(x[keep], 1L, 1L)),
    tolower(substr(x[keep], 2L, nchar(x[keep])))
  )
  x
}

normalize_species <- function(x) {
  x <- squish_chr(gsub("_", " ", as.character(x), fixed = TRUE))
  x
}

first_two_words <- function(x) {
  x <- normalize_species(x)
  vapply(
    strsplit(ifelse(is.na(x), "", x), " ", fixed = TRUE),
    function(parts) {
      parts <- parts[nzchar(parts)]
      if (length(parts) < 2L) return(NA_character_)
      paste(parts[1:2], collapse = " ")
    },
    FUN.VALUE = character(1)
  )
}

is_binomial <- function(x) {
  x2 <- first_two_words(x)
  !is.na(x2) & grepl(
    "^[A-Z][[:alpha:]-]+ [a-z][[:alpha:]-]+$",
    x2,
    perl = TRUE
  ) & !grepl(
    " (sp|spp|cf|aff|indet|unknown)$",
    x2,
    ignore.case = TRUE,
    perl = TRUE
  )
}

safe_filename <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
}

find_project_root <- function(start = getwd()) {
  candidates <- unique(normalize_path(c(
    start,
    file.path(start, ".."),
    file.path(start, "../.."),
    file.path(start, "../../..")
  )))

  for (candidate in candidates) {
    if (
      dir.exists(file.path(candidate, "scripts")) &&
      dir.exists(file.path(candidate, "outputs"))
    ) {
      return(candidate)
    }
  }

  normalize_path(start)
}

require_pkgs <- function(pkgs) {
  missing <- pkgs[!vapply(
    pkgs,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )]

  if (length(missing) > 0L) {
    stop(
      "Missing required R packages: ", paste(missing, collapse = ", "), "\n",
      "Install the CRAN packages with install.packages(), and install ",
      "V.PhyloMaker2 from its official repository when necessary."
    )
  }

  invisible(TRUE)
}

find_column <- function(df, candidates, required = TRUE) {
  nms <- names(df)
  idx <- match(tolower(candidates), tolower(nms))
  idx <- idx[!is.na(idx)]

  if (length(idx) > 0L) return(nms[idx[1L]])

  if (isTRUE(required)) {
    stop(
      "Required column not found. Expected one of: ",
      paste(candidates, collapse = ", ")
    )
  }

  NA_character_
}

resolve_single_file <- function(label, explicit = NULL, candidates = character(0)) {
  explicit <- explicit %||% NA_character_
  explicit <- as.character(explicit)

  if (length(explicit) == 1L && !is.na(explicit) && nzchar(explicit)) {
    if (!file.exists(explicit)) {
      stop(label, " was explicitly configured but does not exist:\n", explicit)
    }
    return(normalize_path(explicit))
  }

  candidates <- unique(as.character(candidates))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  candidates <- candidates[file.exists(candidates)]

  if (length(candidates) == 0L) {
    stop(label, " was not found. Configure its path explicitly in cfg.")
  }

  if (length(candidates) > 1L) {
    info <- file.info(candidates)
    candidates <- candidates[order(info$mtime, decreasing = TRUE)]
    warning(
      "Multiple candidates were found for ", label, ". Using the most recently modified file:\n",
      candidates[1L],
      "\nSet the path explicitly in cfg for strict reproducibility."
    )
  }

  normalize_path(candidates[1L])
}

load_occurrence_table <- function(cfg, project_root, runtime_out_dir, runtime_base_tag) {
  # When the tree is executed inside the full pipeline, the in-memory object is
  # the safest and most reproducible input because it is the exact output of
  # Part I used by the current run.
  if (exists("lin_enriched", envir = .GlobalEnv, inherits = FALSE)) {
    message("[INFO] Using in-memory object: lin_enriched")
    return(list(
      data = as_base_df(get("lin_enriched", envir = .GlobalEnv)),
      source = "in_memory:lin_enriched"
    ))
  }

  explicit <- cfg$tree_occurrence_path %||% cfg$path_lin_enriched %||% NA_character_
  explicit <- as.character(explicit)

  # An explicitly configured path always has priority. Unlike cfg$path_lotus,
  # these two parameters are reserved specifically for the occurrence-level
  # analytical table.
  if (length(explicit) == 1L && !is.na(explicit) && nzchar(explicit)) {
    if (!file.exists(explicit)) {
      stop(
        "The explicitly configured occurrence table does not exist:\n",
        explicit
      )
    }
    path_input <- normalize_path(explicit)
  } else {
    # Discover only the primary lin_enriched Parquet. Summary workbooks from
    # Part II must never be considered occurrence inputs.
    parquet_candidates <- character(0)

    if (!is.null(runtime_out_dir) && !is.null(runtime_base_tag)) {
      parquet_candidates <- c(
        parquet_candidates,
        file.path(
          runtime_out_dir,
          paste0(runtime_base_tag, "_lin_enriched.parquet")
        ),
        file.path(
          runtime_out_dir,
          "PartI_ALL",
          paste0(runtime_base_tag, "_lin_enriched.parquet")
        )
      )
    }

    search_roots <- unique(c(
      if (!is.null(runtime_out_dir)) runtime_out_dir else character(0),
      file.path(project_root, "outputs")
    ))
    search_roots <- search_roots[dir.exists(search_roots)]

    for (root in search_roots) {
      parquet_candidates <- c(
        parquet_candidates,
        list.files(
          root,
          pattern = "_lin_enriched\\.parquet$",
          recursive = TRUE,
          full.names = TRUE
        )
      )
    }

    parquet_candidates <- unique(normalize_path(parquet_candidates))
    parquet_candidates <- parquet_candidates[
      !is.na(parquet_candidates) &
        nzchar(parquet_candidates) &
        file.exists(parquet_candidates)
    ]

    # Exclude audit/derived tables that are not the primary occurrence set.
    excluded_name_pattern <- paste(
      c(
        "all_taxonomic_records",
        "compound_species",
        "taxonomy_excluded",
        "excluded_records",
        "uni_enriched"
      ),
      collapse = "|"
    )
    parquet_candidates <- parquet_candidates[
      !grepl(
        excluded_name_pattern,
        basename(parquet_candidates),
        ignore.case = TRUE,
        perl = TRUE
      )
    ]

    if (length(parquet_candidates) == 0L) {
      # A workbook is accepted only when deliberately configured as cfg$path_lotus
      # and only when it actually contains the lin_enriched sheet. There is no
      # recursive XLSX discovery because BIO_A/B/C summary files are not valid
      # occurrence inputs.
      fallback <- cfg$path_lotus %||% NA_character_
      fallback <- as.character(fallback)

      if (length(fallback) == 1L && !is.na(fallback) && nzchar(fallback) && file.exists(fallback)) {
        fallback_ext <- tolower(tools::file_ext(fallback))
        if (fallback_ext %in% c("xlsx", "xls")) {
          require_pkgs("readxl")
          fallback_sheets <- readxl::excel_sheets(fallback)
          if ("lin_enriched" %in% fallback_sheets) {
            path_input <- normalize_path(fallback)
          } else {
            stop(
              "No primary *_lin_enriched.parquet file was found, and the configured ",
              "cfg$path_lotus workbook does not contain a 'lin_enriched' sheet:\n",
              fallback,
              "\nSet cfg$tree_occurrence_path explicitly to the Part I Parquet file."
            )
          }
        } else {
          stop(
            "No primary *_lin_enriched.parquet file was found. The configured ",
            "cfg$path_lotus is not an Excel workbook with a lin_enriched sheet:\n",
            fallback,
            "\nSet cfg$tree_occurrence_path explicitly."
          )
        }
      } else {
        stop(
          "No primary occurrence-level *_lin_enriched.parquet file was found.\n",
          "Set cfg$tree_occurrence_path explicitly to the Part I analytical Parquet."
        )
      }
    } else {
      # Deterministic scoring: current run > exact base-tag filename > PartI_ALL.
      scores <- rep(0L, length(parquet_candidates))

      if (!is.null(runtime_out_dir)) {
        current_root <- paste0(normalize_path(runtime_out_dir), "/")
        scores <- scores + ifelse(
          startsWith(paste0(parquet_candidates, "/"), current_root),
          100L,
          0L
        )
      }

      if (!is.null(runtime_base_tag)) {
        exact_name <- paste0(runtime_base_tag, "_lin_enriched.parquet")
        scores <- scores + ifelse(
          basename(parquet_candidates) == exact_name,
          50L,
          0L
        )
      }

      scores <- scores + ifelse(
        grepl("[/\\\\]PartI_ALL[/\\\\]", parquet_candidates, perl = TRUE),
        20L,
        0L
      )

      best <- which(scores == max(scores))

      if (length(best) > 1L) {
        # Ties are not resolved by modification time because cloud-sync times are
        # not a reliable provenance signal. Require an explicit path instead.
        stop(
          "Multiple equally suitable primary occurrence Parquets were found:\n- ",
          paste(parquet_candidates[best], collapse = "\n- "),
          "\nSet cfg$tree_occurrence_path explicitly to the file from the current run."
        )
      }

      path_input <- parquet_candidates[best]
    }
  }

  ext <- tolower(tools::file_ext(path_input))

  if (ext == "parquet") {
    require_pkgs("arrow")
    dat <- arrow::read_parquet(path_input, as_data_frame = TRUE)
  } else if (ext %in% c("xlsx", "xls")) {
    require_pkgs("readxl")
    sheets <- readxl::excel_sheets(path_input)
    if (!"lin_enriched" %in% sheets) {
      stop("Sheet 'lin_enriched' was not found in: ", path_input)
    }
    dat <- readxl::read_excel(
      path_input,
      sheet = "lin_enriched",
      .name_repair = "minimal"
    )
  } else if (ext %in% c("csv", "gz")) {
    require_pkgs("readr")
    dat <- readr::read_csv(path_input, show_col_types = FALSE)
  } else {
    stop("Unsupported occurrence input format: ", ext)
  }

  # Fail early if a wrong file happens to satisfy the filename rule.
  required_semantic_columns <- c("family", "genus")
  missing_semantic_columns <- required_semantic_columns[
    !tolower(required_semantic_columns) %in% tolower(names(dat))
  ]
  if (length(missing_semantic_columns) > 0L) {
    stop(
      "The selected occurrence input lacks required columns: ",
      paste(missing_semantic_columns, collapse = ", "),
      "\nSelected file: ", path_input
    )
  }

  message("[INFO] Selected occurrence input: ", path_input)
  list(data = as_base_df(dat), source = path_input)
}

parse_logical_flag <- function(x) {
  if (is.logical(x)) return(x)
  y <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(y))
  out[y %in% c("true", "t", "1", "yes", "y")] <- TRUE
  out[y %in% c("false", "f", "0", "no", "n")] <- FALSE
  out
}

load_rdkit_primary_universe <- function(cfg, project_root, runtime_out_dir) {
  require_pkgs("readr")

  compounds_explicit <- cfg$path_rdkit_compounds %||%
    cfg$rdkit_compounds_path %||%
    cfg$phylo_rdkit_compounds_path %||%
    NA_character_
  qc_explicit <- cfg$path_rdkit_qc %||%
    cfg$rdkit_qc_path %||%
    cfg$path_rdkit_structural_qc %||%
    cfg$phylo_rdkit_qc_path %||%
    NA_character_

  compounds_candidates <- c(
    compounds_explicit,
    if (!is.null(runtime_out_dir)) {
      file.path(runtime_out_dir, "RDKit_ALL", "lotus_flavonoids_rdkit_compounds.csv")
    } else {
      character(0)
    },
    file.path(project_root, "lotus_flavonoids_rdkit_compounds.csv")
  )
  qc_candidates <- c(
    qc_explicit,
    if (!is.null(runtime_out_dir)) {
      file.path(runtime_out_dir, "RDKit_ALL", "lotus_flavonoids_rdkit_structural_qc.csv")
    } else {
      character(0)
    },
    file.path(project_root, "lotus_flavonoids_rdkit_structural_qc.csv")
  )

  if (dir.exists(file.path(project_root, "outputs"))) {
    compounds_candidates <- c(
      compounds_candidates,
      list.files(
        file.path(project_root, "outputs"),
        pattern = "^lotus_flavonoids_rdkit_compounds\\.csv$",
        recursive = TRUE,
        full.names = TRUE
      )
    )
    qc_candidates <- c(
      qc_candidates,
      list.files(
        file.path(project_root, "outputs"),
        pattern = "^lotus_flavonoids_rdkit_structural_qc\\.csv$",
        recursive = TRUE,
        full.names = TRUE
      )
    )
  }

  compounds_path <- resolve_single_file(
    "RDKit compound table",
    explicit = compounds_explicit,
    candidates = compounds_candidates
  )
  qc_path <- resolve_single_file(
    "RDKit structural-QC table",
    explicit = qc_explicit,
    candidates = qc_candidates
  )

  compounds <- as_base_df(readr::read_csv(
    compounds_path,
    show_col_types = FALSE,
    progress = FALSE
  ))
  qc <- as_base_df(readr::read_csv(
    qc_path,
    show_col_types = FALSE,
    progress = FALSE
  ))

  compound_key_col <- find_column(
    compounds,
    c("inchikey", "inchi_key", "standard_inchi_key")
  )
  annotation_col <- find_column(
    compounds,
    c("annotation_status"),
    required = FALSE
  )
  qc_key_col <- find_column(
    qc,
    c("inchikey", "inchi_key", "standard_inchi_key")
  )
  eligibility_col <- find_column(
    qc,
    c("primary_structure_eligible"),
    required = FALSE
  )
  status_col <- find_column(
    qc,
    c("structural_qc_status"),
    required = FALSE
  )
  reasons_col <- find_column(
    qc,
    c("structural_qc_reasons"),
    required = FALSE
  )

  compound_keys <- trim_chr(compounds[[compound_key_col]])
  if (!is.na(annotation_col)) {
    annotation_status <- tolower(trim_chr(compounds[[annotation_col]]))
    compound_keys <- compound_keys[
      is.na(annotation_status) | annotation_status == "valid"
    ]
  }
  compound_keys <- sort(unique(compound_keys[!is.na(compound_keys)]))

  qc_keys <- trim_chr(qc[[qc_key_col]])
  eligible_flag <- if (!is.na(eligibility_col)) {
    parse_logical_flag(qc[[eligibility_col]])
  } else {
    rep(NA, nrow(qc))
  }
  qc_status <- if (!is.na(status_col)) {
    tolower(trim_chr(qc[[status_col]]))
  } else {
    rep(NA_character_, nrow(qc))
  }
  qc_reasons <- if (!is.na(reasons_col)) {
    trim_chr(qc[[reasons_col]])
  } else {
    rep(NA_character_, nrow(qc))
  }

  excluded_mask <- (!is.na(eligible_flag) & !eligible_flag) |
    (!is.na(qc_status) & qc_status == "exclude_primary")
  excluded_details <- data.frame(
    inchikey = qc_keys[excluded_mask],
    structural_qc_status = qc_status[excluded_mask],
    structural_qc_reasons = qc_reasons[excluded_mask],
    stringsAsFactors = FALSE
  )
  excluded_details <- excluded_details[
    !is.na(excluded_details$inchikey),
    ,
    drop = FALSE
  ]
  excluded_details <- excluded_details[
    !duplicated(excluded_details$inchikey),
    ,
    drop = FALSE
  ]

  excluded_keys <- sort(unique(excluded_details$inchikey))
  eligible_keys <- setdiff(compound_keys, excluded_keys)

  if (length(compound_keys) == 0L) {
    stop("The RDKit compound table contains no valid InChIKeys.")
  }
  if (length(eligible_keys) == 0L) {
    stop("No primary-structure-eligible compounds remain after RDKit QC.")
  }

  list(
    eligible_keys = eligible_keys,
    all_rdkit_keys = compound_keys,
    excluded_details = excluded_details,
    compounds_source = compounds_path,
    qc_source = qc_path
  )
}

collapse_values <- function(x) {
  x <- sort(unique(trim_chr(x)))
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA_character_)
  paste(x, collapse = ";")
}

filter_primary_structure_occurrences <- function(occ, primary_universe) {
  require_pkgs("dplyr")

  all_keys <- primary_universe$all_rdkit_keys
  eligible_keys <- primary_universe$eligible_keys
  excluded_details <- primary_universe$excluded_details

  dat <- as_base_df(occ)
  dat$primary_structure_eligible <-
    !is.na(dat$inchikey) & dat$inchikey %in% eligible_keys
  dat$primary_structure_filter_reason <- NA_character_
  dat$primary_structure_filter_reason[is.na(dat$inchikey)] <- "missing_inchikey"
  dat$primary_structure_filter_reason[
    !is.na(dat$inchikey) & !dat$inchikey %in% all_keys
  ] <- "not_found_in_rdkit_compound_output"

  if (nrow(excluded_details) > 0L) {
    reason_map <- excluded_details$structural_qc_reasons
    names(reason_map) <- excluded_details$inchikey
    status_map <- excluded_details$structural_qc_status
    names(status_map) <- excluded_details$inchikey
    rejected <- !is.na(dat$inchikey) & dat$inchikey %in% names(reason_map)
    rejection_reason <- unname(reason_map[dat$inchikey[rejected]])
    rejection_status <- unname(status_map[dat$inchikey[rejected]])
    rejection_reason[is.na(rejection_reason) | !nzchar(rejection_reason)] <-
      rejection_status[is.na(rejection_reason) | !nzchar(rejection_reason)]
    rejection_reason[is.na(rejection_reason) | !nzchar(rejection_reason)] <-
      "exclude_primary"
    dat$primary_structure_filter_reason[rejected] <- rejection_reason
  }
  dat$primary_structure_filter_reason[dat$primary_structure_eligible] <- "eligible"

  family_total <- dat |>
    dplyr::group_by(.data$family) |>
    dplyr::summarise(
      n_occurrence_rows_total = dplyr::n(),
      n_compounds_total = dplyr::n_distinct(
        .data$inchikey[!is.na(.data$inchikey)]
      ),
      n_references_total = dplyr::n_distinct(
        .data$ref_id[!is.na(.data$ref_id)]
      ),
      .groups = "drop"
    )

  family_primary <- dat |>
    dplyr::filter(.data$primary_structure_eligible) |>
    dplyr::group_by(.data$family) |>
    dplyr::summarise(
      n_occurrence_rows_primary = dplyr::n(),
      n_primary_compounds = dplyr::n_distinct(.data$inchikey),
      n_primary_references = dplyr::n_distinct(
        .data$ref_id[!is.na(.data$ref_id)]
      ),
      .groups = "drop"
    )

  family_excluded <- dat |>
    dplyr::filter(!.data$primary_structure_eligible) |>
    dplyr::group_by(.data$family) |>
    dplyr::summarise(
      n_occurrence_rows_excluded = dplyr::n(),
      n_excluded_compounds = dplyr::n_distinct(
        .data$inchikey[!is.na(.data$inchikey)]
      ),
      excluded_inchikeys = collapse_values(.data$inchikey),
      exclusion_reasons = collapse_values(
        .data$primary_structure_filter_reason
      ),
      .groups = "drop"
    )

  family_audit <- family_total |>
    dplyr::left_join(family_primary, by = "family") |>
    dplyr::left_join(family_excluded, by = "family")

  numeric_columns <- c(
    "n_occurrence_rows_primary", "n_primary_compounds",
    "n_primary_references", "n_occurrence_rows_excluded",
    "n_excluded_compounds"
  )
  for (column in numeric_columns) {
    family_audit[[column]][is.na(family_audit[[column]])] <- 0L
  }
  family_audit$retained_in_primary_tree <- family_audit$n_primary_compounds > 0L
  family_audit$tree_exclusion_reason <- ifelse(
    family_audit$retained_in_primary_tree,
    NA_character_,
    "no_primary_structure_eligible_compounds"
  )
  family_audit <- family_audit |>
    dplyr::arrange(.data$family)

  excluded_occurrences <- dat |>
    dplyr::filter(!.data$primary_structure_eligible) |>
    dplyr::arrange(.data$family, .data$inchikey)

  filtered <- dat |>
    dplyr::filter(.data$primary_structure_eligible) |>
    dplyr::select(
      -dplyr::all_of(c(
        "primary_structure_eligible",
        "primary_structure_filter_reason"
      ))
    )

  list(
    data = as_base_df(filtered),
    family_audit = as_base_df(family_audit),
    families_without_primary = as_base_df(
      family_audit[!family_audit$retained_in_primary_tree, , drop = FALSE]
    ),
    excluded_occurrences = as_base_df(excluded_occurrences)
  )
}

standardize_occurrence_columns <- function(df) {
  family_col <- find_column(df, c("family", "family_name"))
  genus_col <- find_column(df, c("genus", "genus_name"))
  species_col <- find_column(df, c("species", "species_name"), required = FALSE)
  inchikey_col <- find_column(df, c("inchikey", "inchi_key", "standard_inchi_key"), required = FALSE)
  ref_col <- find_column(df, c("ref_id", "reference", "reference_id", "doi"), required = FALSE)

  out <- data.frame(
    family = normalize_taxon_word(df[[family_col]]),
    genus = normalize_taxon_word(df[[genus_col]]),
    species = if (!is.na(species_col)) normalize_species(df[[species_col]]) else NA_character_,
    inchikey = if (!is.na(inchikey_col)) trim_chr(df[[inchikey_col]]) else NA_character_,
    ref_id = if (!is.na(ref_col)) trim_chr(df[[ref_col]]) else NA_character_,
    stringsAsFactors = FALSE
  )

  out <- out[
    !is.na(out$family) & !is.na(out$genus) &
      nzchar(out$family) & nzchar(out$genus),
    ,
    drop = FALSE
  ]

  out$species_binomial <- first_two_words(out$species)
  out$valid_binomial <- is_binomial(out$species_binomial)
  out$species_genus <- ifelse(
    out$valid_binomial,
    sub(" .*", "", out$species_binomial),
    NA_character_
  )
  out$valid_binomial <- out$valid_binomial & out$species_genus == out$genus

  as_base_df(out)
}

apply_curated_family_corrections <- function(
  occ,
  corrections = character(0),
  correction_reasons = character(0)
) {
  require_pkgs("dplyr")

  corrections <- corrections[
    !is.na(names(corrections)) & nzchar(names(corrections)) &
      !is.na(corrections) & nzchar(corrections)
  ]

  if (length(corrections) == 0L) {
    return(list(
      data = as_base_df(occ),
      audit = data.frame(
        original_family = character(0),
        genus = character(0),
        corrected_family = character(0),
        reason = character(0),
        n_occurrence_rows = integer(0),
        n_compounds = integer(0),
        n_references = integer(0),
        stringsAsFactors = FALSE
      )
    ))
  }

  pair_key <- paste(occ$family, occ$genus, sep = "|")
  matched <- pair_key %in% names(corrections)

  if (!any(matched)) {
    return(list(
      data = as_base_df(occ),
      audit = data.frame(
        original_family = character(0),
        genus = character(0),
        corrected_family = character(0),
        reason = character(0),
        n_occurrence_rows = integer(0),
        n_compounds = integer(0),
        n_references = integer(0),
        stringsAsFactors = FALSE
      )
    ))
  }

  audit_rows <- occ[matched, , drop = FALSE]
  audit_rows$original_family <- audit_rows$family
  audit_rows$corrected_family <- unname(corrections[pair_key[matched]])
  audit_rows$reason <- unname(correction_reasons[pair_key[matched]])
  audit_rows$reason[is.na(audit_rows$reason) | !nzchar(audit_rows$reason)] <-
    "Curated genus-family correction for APG/backbone consistency"

  audit <- audit_rows |>
    dplyr::group_by(
      .data$original_family,
      .data$genus,
      .data$corrected_family,
      .data$reason
    ) |>
    dplyr::summarise(
      n_occurrence_rows = dplyr::n(),
      n_compounds = dplyr::n_distinct(.data$inchikey[!is.na(.data$inchikey)]),
      n_references = dplyr::n_distinct(.data$ref_id[!is.na(.data$ref_id)]),
      .groups = "drop"
    )

  occ$family[matched] <- unname(corrections[pair_key[matched]])

  list(
    data = as_base_df(occ),
    audit = as_base_df(audit)
  )
}

select_family_representatives <- function(occ, forbidden_pairs = character(0)) {
  require_pkgs("dplyr")

  forbidden_pairs <- forbidden_pairs[!is.na(names(forbidden_pairs)) & nzchar(names(forbidden_pairs))]
  pair_key <- paste(occ$family, occ$genus, sep = "|")
  occ$representative_allowed <- !pair_key %in% names(forbidden_pairs)
  occ$representative_exclusion_reason <- unname(forbidden_pairs[pair_key])

  species_candidates <- occ |>
    dplyr::filter(.data$valid_binomial) |>
    dplyr::group_by(.data$family, .data$genus, .data$species_binomial) |>
    dplyr::summarise(
      n_occurrence_rows = dplyr::n(),
      n_compounds = dplyr::n_distinct(.data$inchikey[!is.na(.data$inchikey)]),
      n_references = dplyr::n_distinct(.data$ref_id[!is.na(.data$ref_id)]),
      representative_allowed = all(.data$representative_allowed),
      representative_exclusion_reason = dplyr::first(
        .data$representative_exclusion_reason
      ),
      .groups = "drop"
    ) |>
    dplyr::arrange(
      .data$family,
      dplyr::desc(.data$representative_allowed),
      dplyr::desc(.data$n_compounds),
      dplyr::desc(.data$n_occurrence_rows),
      dplyr::desc(.data$n_references),
      .data$species_binomial
    ) |>
    dplyr::group_by(.data$family) |>
    dplyr::slice_head(n = 1L) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      representative_source = "observed_binomial_species",
      tree_species = gsub(" ", "_", .data$species_binomial, fixed = TRUE),
      representative_taxonomy_status = dplyr::if_else(
        .data$representative_allowed,
        "compatible",
        "forced_incompatible_no_alternative"
      )
    )

  represented_families <- unique(species_candidates$family)

  genus_candidates <- occ |>
    dplyr::filter(!.data$family %in% represented_families) |>
    dplyr::group_by(.data$family, .data$genus) |>
    dplyr::summarise(
      n_occurrence_rows = dplyr::n(),
      n_compounds = dplyr::n_distinct(.data$inchikey[!is.na(.data$inchikey)]),
      n_references = dplyr::n_distinct(.data$ref_id[!is.na(.data$ref_id)]),
      representative_allowed = all(.data$representative_allowed),
      representative_exclusion_reason = dplyr::first(
        .data$representative_exclusion_reason
      ),
      .groups = "drop"
    ) |>
    dplyr::arrange(
      .data$family,
      dplyr::desc(.data$representative_allowed),
      dplyr::desc(.data$n_compounds),
      dplyr::desc(.data$n_occurrence_rows),
      dplyr::desc(.data$n_references),
      .data$genus
    ) |>
    dplyr::group_by(.data$family) |>
    dplyr::slice_head(n = 1L) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      species_binomial = NA_character_,
      representative_source = "genus_level_dummy",
      representative_taxonomy_status = dplyr::if_else(
        .data$representative_allowed,
        "compatible",
        "forced_incompatible_no_alternative"
      )
    )

  if (nrow(genus_candidates) > 0L) {
    genus_candidates$tree_species <- paste0(
      genus_candidates$genus,
      "_sp",
      sprintf("%04d", seq_len(nrow(genus_candidates)))
    )
  } else {
    genus_candidates$tree_species <- character(0)
  }

  out <- dplyr::bind_rows(species_candidates, genus_candidates) |>
    dplyr::arrange(.data$family)

  duplicated_species <- duplicated(out$tree_species) | duplicated(out$tree_species, fromLast = TRUE)
  if (any(duplicated_species)) {
    idx <- which(duplicated_species)
    out$tree_species[idx] <- paste0(
      out$genus[idx],
      "_spdup",
      sprintf("%04d", seq_along(idx))
    )
    out$representative_source[idx] <- "genus_level_dummy_due_to_duplicate_tip"
  }

  if (anyDuplicated(out$family)) stop("Representative selection produced duplicate families.")
  if (anyDuplicated(out$tree_species)) stop("Representative selection produced duplicate tree tips.")

  incompatible <- out$representative_taxonomy_status != "compatible"
  if (any(incompatible)) {
    bad <- paste0(
      out$family[incompatible], " -> ", out$genus[incompatible],
      collapse = "; "
    )
    stop(
      "No APG/backbone-compatible representative was available for: ", bad, ". ",
      "Provide an explicit compatible representative or revise the family assignment."
    )
  }

  as_base_df(out)
}

load_apg_table <- function(path_apg) {
  require_pkgs("readr")
  apg <- as_base_df(readr::read_csv(path_apg, show_col_types = FALSE))

  family_col <- find_column(apg, c("Family", "family"))
  names(apg)[match(family_col, names(apg))] <- "Family"

  for (nm in c("Order", "Node.1", "Node.2", "Node.3", "Node.4")) {
    if (!nm %in% names(apg)) apg[[nm]] <- NA_character_
  }

  apg$Family <- normalize_taxon_word(apg$Family)
  for (nm in c("Order", "Node.1", "Node.2", "Node.3", "Node.4")) {
    apg[[nm]] <- squish_chr(apg[[nm]])
  }

  duplicate_families <- apg$Family[duplicated(apg$Family)]
  if (length(duplicate_families) > 0L) {
    warning(
      "Duplicated families were found in the APG reference table; the first record will be used: ",
      paste(unique(duplicate_families), collapse = ", ")
    )
    apg <- apg[!duplicated(apg$Family), , drop = FALSE]
  }

  apg
}

contains_term <- function(x, pattern) {
  !is.na(x) & grepl(pattern, x, ignore.case = TRUE, perl = TRUE)
}

classify_major_clade <- function(family, order_name, nodes_text, manual_clade_overrides) {
  family <- normalize_taxon_word(family)
  order_name <- squish_chr(order_name)
  nodes_text <- squish_chr(nodes_text)

  if (!is.na(family) && family %in% names(manual_clade_overrides)) {
    return(unname(manual_clade_overrides[[family]]))
  }

  combined <- paste(
    c(order_name, nodes_text),
    collapse = " | "
  )
  combined <- tolower(combined)
  ord <- if (is.na(order_name) || !nzchar(order_name)) "" else tolower(order_name)

 
  if (grepl("bryophy|moss|liverwort|hornwort|marchanti|anthocerot", combined)) {
    return("Bryophytes")
  }

  if (grepl("lycophy|lycopod|selaginell|isoet", combined)) {
    return("Lycophytes")
  }

  if (grepl("fern|monilophy|polypod|equiset|ophiogloss|maratt|psilot", combined)) {
    return("Ferns")
  }

  if (grepl("gymnosperm|conifer|cycad|ginkgo|gnetophy", combined)) {
    return("Gymnosperms")
  }

 
  if (ord == "saxifragales") return("Saxifragales")
  if (ord == "santalales") return("Santalales")
  if (ord == "caryophyllales") return("Caryophyllales")
  if (ord == "vitales") return("Rosids (Basal)")

  if (ord %in% c("amborellales", "nymphaeales", "austrobaileyales")) {
    return("Basal Angiosperms")
  }

  # Chloranthales is retained as an explicit early-diverging mesangiosperm
  # lineage rather than being folded into the ANA grade.
  if (ord == "chloranthales") {
    return("Chloranthales")
  }

  if (ord %in% c("canellales", "laurales", "magnoliales", "piperales")) {
    return("Magnoliids")
  }

  if (grepl("campanulids", combined)) return("Asterids (Campanulids)")
  if (grepl("lamiids", combined)) return("Asterids (Lamiids)")
  if (grepl("fabids", combined)) return("Rosids (Fabids)")
  if (grepl("malvids", combined)) return("Rosids (Malvids)")
  if (grepl("commelinids", combined)) return("Monocots (Commelinids)")

  if (ord %in% c("cornales", "ericales")) return("Asterids (Basal)")

  if (grepl("asterids", combined)) return("Asterids (Basal)")
  if (grepl("rosids", combined)) return("Rosids (Basal)")
  if (grepl("magnoliids", combined)) return("Magnoliids")
  if (grepl("monocots", combined)) return("Monocots")

  if (ord %in% c(
    "acorales", "alismatales", "petrosaviales", "dioscoreales",
    "pandanales", "liliales", "asparagales", "arecales",
    "commelinales", "poales", "zingiberales"
  )) {
    if (ord %in% c("arecales", "commelinales", "poales", "zingiberales")) {
      return("Monocots (Commelinids)")
    }
    return("Monocots")
  }

  if (ord %in% c("ranunculales", "proteales", "trochodendrales", "buxales")) {
    return("Basal Eudicots")
  }

  if (grepl("eudicots", combined)) return("Other Eudicots")
  if (grepl("angiosperm", combined)) return("Other Angiosperms")

  "Other"
}

assign_apg_clades <- function(representatives, apg, manual_clade_overrides, nonvascular_families) {
  idx <- match(representatives$family, apg$Family)

  out <- cbind(
    representatives,
    apg[idx, c("Order", "Node.1", "Node.2", "Node.3", "Node.4"), drop = FALSE]
  )
  out <- as_base_df(out)

  out$apg_match <- !is.na(idx)
  out$nodes_text <- apply(
    out[, c("Node.1", "Node.2", "Node.3", "Node.4"), drop = FALSE],
    1L,
    function(z) paste(z[!is.na(z) & nzchar(z)], collapse = " | ")
  )

  out$FinalClade <- vapply(
    seq_len(nrow(out)),
    function(i) {
      classify_major_clade(
        family = out$family[i],
        order_name = out$Order[i],
        nodes_text = out$nodes_text[i],
        manual_clade_overrides = manual_clade_overrides
      )
    },
    FUN.VALUE = character(1)
  )

  manual_nonvascular <- out$family %in% nonvascular_families
  out$FinalClade[manual_nonvascular] <- "Bryophytes"

  out$is_vascular <- !out$FinalClade %in% c("Bryophytes")
  out$clade_assignment_source <- ifelse(
    out$family %in% names(manual_clade_overrides),
    "manual_family_override",
    ifelse(
      manual_nonvascular,
      "manual_nonvascular_override",
      ifelse(out$apg_match, "APG_reference", "unmatched")
    )
  )

  as_base_df(out)
}

apply_backbone_taxonomy <- function(
  df,
  family_overrides,
  genus_family_overrides,
  available_backbone_families
) {
  out <- df
  out$original_family <- out$family
  out$backbone_family <- out$family
  out$original_family_in_backbone <- out$family %in% available_backbone_families
  out$backbone_override_reason <- NA_character_


  family_hit <- !out$original_family_in_backbone & out$family %in% names(family_overrides)
  out$backbone_family[family_hit] <- unname(family_overrides[out$family[family_hit]])
  out$backbone_override_reason[family_hit] <- "family_override_absent_from_backbone"

  genus_hit <- !out$original_family_in_backbone & out$genus %in% names(genus_family_overrides)
  out$backbone_family[genus_hit] <- unname(genus_family_overrides[out$genus[genus_hit]])
  out$backbone_override_reason[genus_hit] <- ifelse(
    is.na(out$backbone_override_reason[genus_hit]),
    "genus_override_absent_from_backbone",
    paste(
      out$backbone_override_reason[genus_hit],
      "genus_override_absent_from_backbone",
      sep = ";"
    )
  )

  out$backbone_family_in_backbone <- out$backbone_family %in% available_backbone_families
  as_base_df(out)
}

extract_scenario_tree <- function(result, scenario) {
  scenario_number <- sub("^S", "", toupper(scenario))
  key <- paste0("scenario.", scenario_number)
  tree <- result[[key]]
  if (is.null(tree) || !inherits(tree, "phylo")) {
    stop("V.PhyloMaker2 did not return a valid ", key, " phylogeny.")
  }
  tree
}

restore_family_tip_labels <- function(tree, mapping) {
  normalized_tree_tips <- gsub(" ", "_", tree$tip.label, fixed = TRUE)
  normalized_input_tips <- gsub(" ", "_", mapping$tree_species, fixed = TRUE)
  idx <- match(normalized_tree_tips, normalized_input_tips)

  if (anyNA(idx)) {
    missing_tips <- tree$tip.label[is.na(idx)]
    stop(
      "Some phylogeny tips could not be mapped back to analytical families:\n- ",
      paste(missing_tips, collapse = "\n- ")
    )
  }

  restored <- mapping$original_family[idx]
  if (anyNA(restored) || any(!nzchar(restored))) {
    stop("Family restoration produced missing tip labels.")
  }
  if (anyDuplicated(restored)) {
    stop(
      "Family restoration produced duplicate tips: ",
      paste(unique(restored[duplicated(restored)]), collapse = ", ")
    )
  }

  tree$tip.label <- restored
  tree
}


collapse_analysis_clade <- function(final_clade, family = NA_character_) {
  final_clade <- as.character(final_clade)
  family <- as.character(family)
  out <- final_clade

  # ANA is restricted to Amborellales, Nymphaeales, and Austrobaileyales.
  # Chloranthales remains explicit and is only combined with other early
  # angiosperm groups later in broader analytical summaries.
  out[final_clade == "Basal Angiosperms"] <- "ANA"
  out[final_clade %in% c("Basal Eudicots", "Other Eudicots")] <- "Basal Eudicots"
  out[final_clade %in% c("Caryophyllales", "Santalales")] <- "Caryophy_Santalales"
  out
}

get_descendant_tips <- function(tree, node) {
  if (length(node) != 1L || is.na(node)) return(character(0))
  if (node <= ape::Ntip(tree)) return(tree$tip.label[node])
  ape::extract.clade(tree, node = node)$tip.label
}

build_local_neighborhood_audit <- function(tree, assignments) {
  assignments <- as_base_df(assignments)
  rownames(assignments) <- assignments$family
  n_tip <- ape::Ntip(tree)

  out <- lapply(seq_len(n_tip), function(i) {
    tip <- tree$tip.label[i]
    incoming <- which(tree$edge[, 2] == i)
    parent <- if (length(incoming) == 1L) tree$edge[incoming, 1] else NA_integer_
    sister_children <- if (!is.na(parent)) {
      setdiff(tree$edge[tree$edge[, 1] == parent, 2], i)
    } else {
      integer(0)
    }
    sister_tips <- unique(unlist(lapply(sister_children, function(x) {
      get_descendant_tips(tree, x)
    })))

    sister_orders <- unique(assignments[sister_tips, "Order", drop = TRUE])
    sister_orders <- sister_orders[!is.na(sister_orders) & nzchar(sister_orders)]
    sister_clades <- unique(assignments[sister_tips, "FinalClade", drop = TRUE])
    sister_clades <- sister_clades[!is.na(sister_clades) & nzchar(sister_clades)]
    expected_order <- assignments[tip, "Order", drop = TRUE]
    expected_clade <- assignments[tip, "FinalClade", drop = TRUE]

    data.frame(
      family = tip,
      expected_order = expected_order,
      expected_clade = expected_clade,
      n_sister_tips = length(sister_tips),
      sister_tips = paste(sister_tips, collapse = ";"),
      sister_orders = paste(sister_orders, collapse = ";"),
      sister_clades = paste(sister_clades, collapse = ";"),
      same_order_sister = !is.na(expected_order) && expected_order %in% sister_orders,
      same_clade_sister = !is.na(expected_clade) && expected_clade %in% sister_clades,
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(out)
}

group_mrca_purity <- function(tree, tips) {
  tips <- intersect(unique(as.character(tips)), tree$tip.label)
  if (length(tips) < 2L) {
    return(list(node = NA_integer_, n_group = length(tips), n_descendants = NA_integer_, purity = NA_real_))
  }
  node <- ape::getMRCA(tree, tips)
  descendants <- get_descendant_tips(tree, node)
  list(
    node = node,
    n_group = length(tips),
    n_descendants = length(descendants),
    purity = length(tips) / length(descendants)
  )
}

audit_order_monophyly <- function(tree, assignments) {
  assignments <- as_base_df(assignments)
  assignments <- assignments[
    assignments$family %in% tree$tip.label &
      !is.na(assignments$Order) & nzchar(assignments$Order),
    , drop = FALSE
  ]

  orders <- sort(unique(assignments$Order))
  rows <- lapply(orders, function(ord) {
    tips <- assignments$family[assignments$Order == ord]
    base <- group_mrca_purity(tree, tips)

    candidate <- NA_character_
    best_without <- NA_real_
    improvement <- NA_real_
    if (length(tips) >= 3L && is.finite(base$purity) && base$purity < 1) {
      leave_one_out <- vapply(tips, function(tip) {
        x <- group_mrca_purity(tree, setdiff(tips, tip))$purity
        ifelse(is.finite(x), x, NA_real_)
      }, FUN.VALUE = numeric(1))
      if (any(is.finite(leave_one_out))) {
        idx <- which.max(leave_one_out)
        candidate <- tips[idx]
        best_without <- leave_one_out[idx]
        improvement <- best_without - base$purity
      }
    }

    data.frame(
      order = ord,
      n_families = length(tips),
      mrca_descendant_tips = base$n_descendants,
      order_purity = base$purity,
      likely_outlier = candidate,
      purity_without_likely_outlier = best_without,
      purity_improvement = improvement,
      review_flag = is.finite(base$purity) && base$purity < 0.8,
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

apply_topology_corrections <- function(tree, assignments, corrections = list()) {
  if (length(corrections) == 0L) {
    return(list(tree = tree, audit = data.frame()))
  }

  corrected_tree <- tree
  corrected_tree$edge.length <- NULL
  audit_rows <- list()

  for (family in names(corrections)) {
    spec <- corrections[[family]]
    if (!family %in% corrected_tree$tip.label) {
      stop("Topology correction requested for missing tree tip: ", family)
    }

    target_order <- as.character(spec$target_order %||% NA_character_)
    target_families <- as.character(spec$target_families %||% character(0))
    reason <- as.character(spec$reason %||% "Curated topology correction")

    if (!is.na(target_order) && nzchar(target_order)) {
      target_families <- assignments$family[
        assignments$Order == target_order & assignments$family != family
      ]
    }
    target_families <- intersect(unique(target_families), corrected_tree$tip.label)
    target_families <- setdiff(target_families, family)

    if (length(target_families) < 2L) {
      stop(
        "Topology correction for ", family,
        " requires at least two target families; found ", length(target_families), "."
      )
    }

    before_audit <- build_local_neighborhood_audit(corrected_tree, assignments)
    before_row <- before_audit[before_audit$family == family, , drop = FALSE]

    corrected_tree <- ape::drop.tip(corrected_tree, family)
    target_node <- ape::getMRCA(corrected_tree, target_families)
    if (is.null(target_node) || is.na(target_node)) {
      stop("Could not determine the target MRCA for topology correction: ", family)
    }


    donor_tree <- ape::read.tree(text = paste0("(", family, ");"))
    donor_tree$edge.length <- NULL

    corrected_tree <- ape::bind.tree(
      x = corrected_tree,
      y = donor_tree,
      where = target_node,
      position = 0
    )
    corrected_tree <- ape::collapse.singles(corrected_tree)
    corrected_tree <- ape::multi2di(corrected_tree, random = FALSE)
    corrected_tree <- ape::ladderize(corrected_tree, right = FALSE)

    after_audit <- build_local_neighborhood_audit(corrected_tree, assignments)
    after_row <- after_audit[after_audit$family == family, , drop = FALSE]

    audit_rows[[length(audit_rows) + 1L]] <- data.frame(
      family = family,
      target_order = target_order,
      target_families = paste(sort(target_families), collapse = ";"),
      reason = reason,
      sister_orders_before = before_row$sister_orders %||% NA_character_,
      sister_clades_before = before_row$sister_clades %||% NA_character_,
      sister_orders_after = after_row$sister_orders %||% NA_character_,
      sister_clades_after = after_row$sister_clades %||% NA_character_,
      stringsAsFactors = FALSE
    )
  }

  list(tree = corrected_tree, audit = dplyr::bind_rows(audit_rows))
}

kv_tab <- function(key, values) {
  paste(c(key, values), collapse = "\t")
}

write_tree_colors <- function(df, path_out) {
  body_branch <- paste(df$ID, "branch", df$Color, "normal", "2", sep = "\t")
  body_label <- paste(df$ID, "label", df$Color, "bold", "1.5", sep = "\t")
  writeLines(
    c("TREE_COLORS", "SEPARATOR TAB", "DATA", body_branch, body_label),
    path_out
  )
}

write_color_strip <- function(
  df, palette, path_out,
  dataset_label = "APG macroclades",
  legend_title = "Major clades",
  strip_width = 100
) {
  used_labels <- names(palette)[names(palette) %in% unique(df$FinalClade)]
  used_colors <- unname(palette[used_labels])

  header <- c(
    "DATASET_COLORSTRIP",
    "SEPARATOR TAB",
    kv_tab("DATASET_LABEL", dataset_label),
    kv_tab("COLOR", "#000000"),
    kv_tab("STRIP_WIDTH", as.character(strip_width)),
    kv_tab("LEGEND_TITLE", legend_title),
    kv_tab("LEGEND_SHAPES", rep("1", length(used_labels))),
    kv_tab("LEGEND_COLORS", used_colors),
    kv_tab("LEGEND_LABELS", used_labels),
    "DATA"
  )

  body <- paste(df$ID, df$Color, df$FinalClade, sep = "\t")
  writeLines(c(header, body), path_out)
}

write_metadata <- function(metadata, path_out) {
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(metadata, path_out, pretty = TRUE, auto_unbox = TRUE, na = "null")
  } else {
    dput(metadata, file = sub("\\.json$", ".R", path_out))
    warning("Package 'jsonlite' is unavailable; metadata were written as an R object instead.")
  }
}



cfg <- if (exists("cfg", envir = .GlobalEnv, inherits = FALSE)) {
  get("cfg", envir = .GlobalEnv)
} else {
  list()
}

require_pkgs(c("ape", "V.PhyloMaker2", "dplyr", "readr"))

project_root <- normalize_path(cfg$project_root %||% find_project_root())
runtime_out_dir <- if (exists("OUT_DIR", envir = .GlobalEnv, inherits = FALSE)) {
  normalize_path(get("OUT_DIR", envir = .GlobalEnv))
} else {
  NULL
}
runtime_base_tag <- if (exists("base_tag", envir = .GlobalEnv, inherits = FALSE)) {
  as.character(get("base_tag", envir = .GlobalEnv))
} else {
  NULL
}

# Prefer files from the current pipeline run whenever the user has not
# supplied explicit paths in cfg. This prevents an older, more recently
# synchronized file elsewhere under outputs/ from being selected implicitly.
if (!is.null(runtime_out_dir)) {
  default_rdkit_compounds <- file.path(
    runtime_out_dir,
    "RDKit_ALL",
    "lotus_flavonoids_rdkit_compounds.csv"
  )
  default_rdkit_qc <- file.path(
    runtime_out_dir,
    "RDKit_ALL",
    "lotus_flavonoids_rdkit_structural_qc.csv"
  )

  if (
    is.null(cfg$path_rdkit_compounds) &&
      is.null(cfg$rdkit_compounds_path) &&
      is.null(cfg$phylo_rdkit_compounds_path) &&
      file.exists(default_rdkit_compounds)
  ) {
    cfg$path_rdkit_compounds <- normalize_path(default_rdkit_compounds)
  }

  if (
    is.null(cfg$path_rdkit_qc) &&
      is.null(cfg$rdkit_qc_path) &&
      is.null(cfg$path_rdkit_structural_qc) &&
      is.null(cfg$phylo_rdkit_qc_path) &&
      file.exists(default_rdkit_qc)
  ) {
    cfg$path_rdkit_qc <- normalize_path(default_rdkit_qc)
  }

  if (!is.null(runtime_base_tag)) {
    default_occurrence <- file.path(
      runtime_out_dir,
      paste0(runtime_base_tag, "_lin_enriched.parquet")
    )

    if (
      is.null(cfg$tree_occurrence_path) &&
        is.null(cfg$path_lin_enriched) &&
        file.exists(default_occurrence)
    ) {
      cfg$tree_occurrence_path <- normalize_path(default_occurrence)
    }
  }
}

# The current analytical workflow is a vascular-family visualization tree.
# Non-vascular records remain available in the audit outputs.
if (is.null(cfg$phylo_scenario)) cfg$phylo_scenario <- "S3"
if (is.null(cfg$phylo_exclude_nonvascular)) cfg$phylo_exclude_nonvascular <- TRUE

scenario <- toupper(cfg$phylo_scenario %||% "S3")
if (!scenario %in% c("S1", "S2", "S3")) {
  stop("cfg$phylo_scenario must be one of: S1, S2, S3.")
}

exclude_nonvascular <- isTRUE(cfg$phylo_exclude_nonvascular %||% TRUE)

path_apg <- resolve_single_file(
  "APG reference table",
  explicit = cfg$path_apg %||% NA_character_,
  candidates = c(
    file.path(project_root, "scripts", "APGIV_family_order_clades_WorldFlora.csv"),
    file.path(project_root, "APGIV_family_order_clades_WorldFlora.csv")
  )
)

out_dir <- normalize_path(
  cfg$phylo_out_dir %||%
    if (!is.null(runtime_out_dir)) {
      file.path(runtime_out_dir, "Tree_APG_iTOL")
    } else {
      file.path(project_root, "phylo_outputs_FINAL")
    }
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

manual_group_clade_overrides <- c(
  # Lycophytes
  "Lycopodiaceae" = "Lycophytes",
  "Selaginellaceae" = "Lycophytes",
  "Isoetaceae" = "Lycophytes",
  # Ferns and fern allies
  "Ophioglossaceae" = "Ferns",
  "Helminthostachyaceae" = "Ferns",
  "Botrychiaceae" = "Ferns",
  "Psilotaceae" = "Ferns",
  "Equisetaceae" = "Ferns",
  "Marattiaceae" = "Ferns",
  "Osmundaceae" = "Ferns",
  "Hymenophyllaceae" = "Ferns",
  "Gleicheniaceae" = "Ferns",
  "Schizaeaceae" = "Ferns",
  "Lygodiaceae" = "Ferns",
  "Pteridaceae" = "Ferns",
  "Dennstaedtiaceae" = "Ferns",
  "Thelypteridaceae" = "Ferns",
  "Blechnaceae" = "Ferns",
  "Dryopteridaceae" = "Ferns",
  "Polypodiaceae" = "Ferns",
  "Aspleniaceae" = "Ferns",
  "Cyatheaceae" = "Ferns",
  "Dicksoniaceae" = "Ferns",
  "Salviniaceae" = "Ferns",
  "Marsileaceae" = "Ferns",
  # Gymnosperms
  "Pinaceae" = "Gymnosperms",
  "Araucariaceae" = "Gymnosperms",
  "Podocarpaceae" = "Gymnosperms",
  "Sciadopityaceae" = "Gymnosperms",
  "Taxaceae" = "Gymnosperms",
  "Cephalotaxaceae" = "Gymnosperms",
  "Cupressaceae" = "Gymnosperms",
  "Cycadaceae" = "Gymnosperms",
  "Zamiaceae" = "Gymnosperms",
  "Ginkgoaceae" = "Gymnosperms",
  "Gnetaceae" = "Gymnosperms",
  "Ephedraceae" = "Gymnosperms",
  "Welwitschiaceae" = "Gymnosperms",
  # Non-vascular families retained only for auditing
  "Polytrichaceae" = "Bryophytes",
  "Aulacomniaceae" = "Bryophytes",
  "Hypnaceae" = "Bryophytes",
  "Mniaceae" = "Bryophytes",
  "Dicranaceae" = "Bryophytes",
  "Bryaceae" = "Bryophytes",
  "Leucobryaceae" = "Bryophytes",
  "Marchantiaceae" = "Bryophytes",
  "Sphagnaceae" = "Bryophytes",
  "Acrobolbaceae" = "Bryophytes",
  "Corsiniaceae" = "Bryophytes",
  "Lembophyllaceae" = "Bryophytes"
)

manual_clade_overrides <- c(
  manual_group_clade_overrides,
  "Namaceae" = "Asterids (Lamiids)",
  "Heliotropiaceae" = "Asterids (Lamiids)",
  "Cordiaceae" = "Asterids (Lamiids)",
  "Ehretiaceae" = "Asterids (Lamiids)",
  "Viburnaceae" = "Asterids (Campanulids)",
  "Vitaceae" = "Rosids (Basal)",
  "Cynomoriaceae" = "Saxifragales"
)
manual_clade_overrides <- c(
  manual_clade_overrides,
  cfg$phylo_manual_clade_overrides %||% character(0)
)
manual_clade_overrides <- manual_clade_overrides[!duplicated(names(manual_clade_overrides), fromLast = TRUE)]

nonvascular_families <- unique(c(
  "Polytrichaceae", "Aulacomniaceae", "Hypnaceae", "Mniaceae",
  "Dicranaceae", "Bryaceae", "Leucobryaceae", "Marchantiaceae",
  "Sphagnaceae", "Acrobolbaceae", "Corsiniaceae", "Lembophyllaceae",
  cfg$phylo_nonvascular_families %||% character(0)
))


family_backbone_overrides <- c(
  "Helminthostachyaceae" = "Ophioglossaceae",
  "Botrychiaceae" = "Ophioglossaceae",
  "Namaceae" = "Boraginaceae",
  "Heliotropiaceae" = "Boraginaceae",
  "Cordiaceae" = "Boraginaceae",
  "Ehretiaceae" = "Boraginaceae",
  "Cephalotaxaceae" = "Taxaceae",
  "Viburnaceae" = "Adoxaceae"
)
family_backbone_overrides <- c(
  family_backbone_overrides,
  cfg$phylo_family_backbone_overrides %||% character(0)
)
family_backbone_overrides <- family_backbone_overrides[
  !duplicated(names(family_backbone_overrides), fromLast = TRUE)
]

genus_family_overrides <- c(
  "Cephalotaxus" = "Taxaceae",
  "Sambucus" = "Adoxaceae",
  "Viburnum" = "Adoxaceae",
  "Cleome" = "Cleomaceae"
)
genus_family_overrides <- c(
  genus_family_overrides,
  cfg$phylo_genus_family_overrides %||% character(0)
)
genus_family_overrides <- genus_family_overrides[
  !duplicated(names(genus_family_overrides), fromLast = TRUE)
]


curated_family_corrections <- c(
  "Capparaceae|Cleome" = "Cleomaceae"
)
curated_family_corrections <- c(
  curated_family_corrections,
  cfg$phylo_curated_family_corrections %||% character(0)
)
curated_family_corrections <- curated_family_corrections[
  !duplicated(names(curated_family_corrections), fromLast = TRUE)
]

curated_family_correction_reasons <- c(
  "Capparaceae|Cleome" = paste(
    "Cleome droserifolia is treated in Cleomaceae by the APG reference",
    "and the V.PhyloMaker2 backbone; the source family label Capparaceae",
    "is retained only in the correction audit"
  )
)
curated_family_correction_reasons <- c(
  curated_family_correction_reasons,
  cfg$phylo_curated_family_correction_reasons %||% character(0)
)
curated_family_correction_reasons <- curated_family_correction_reasons[
  !duplicated(names(curated_family_correction_reasons), fromLast = TRUE)
]

# Genus-family pairs that must not be used as family representatives after
# curated corrections. Empty by default; users may supply additional pairs.
representative_forbidden_pairs <- character(0)
representative_forbidden_pairs <- c(
  representative_forbidden_pairs,
  cfg$phylo_representative_forbidden_pairs %||% character(0)
)
representative_forbidden_pairs <- representative_forbidden_pairs[
  !duplicated(names(representative_forbidden_pairs), fromLast = TRUE)
]

clade_palette <- c(
  "Bryophytes" = "#8C8C8C",
  "Lycophytes" = "#6B6B6B",
  "Ferns" = "#4D4D4D",
  "Gymnosperms" = "#000000",
  "Basal Angiosperms" = "#E5C494",
  "Chloranthales" = "#B35806",
  "Magnoliids" = "#FFD92F",
  "Monocots" = "#7570B3",
  "Monocots (Commelinids)" = "#54278F",
  "Basal Eudicots" = "#FC8D62",
  "Other Eudicots" = "#BDBDBD",
  "Saxifragales" = "#1B9E77",
  "Rosids (Basal)" = "#8FBC8F",
  "Rosids (Fabids)" = "#66A61E",
  "Rosids (Malvids)" = "#E7298A",
  "Caryophyllales" = "#E31A1C",
  "Santalales" = "#FF7F00",
  "Asterids (Basal)" = "#FDAE6B",
  "Asterids (Lamiids)" = "#377EB8",
  "Asterids (Campanulids)" = "#984EA3",
  "Other Angiosperms" = "#969696",
  "Other" = "#252525"
)



topology_corrections <- list(
  "Cynomoriaceae" = list(
    target_order = "Saxifragales",
    reason = paste(
      "Cynomoriaceae was placed within Rosales by the GBOTB representative",
      "although the APG reference used in this workflow assigns it to Saxifragales;",
      "the family tip is regrafted as sister to the sampled Saxifragales clade"
    )
  )
)
if (length(cfg$phylo_topology_corrections %||% list()) > 0L) {
  topology_corrections <- utils::modifyList(
    topology_corrections,
    cfg$phylo_topology_corrections
  )
}


analysis_clade_palette <- c(
  "ANA" = "#7F007F",
  "Chloranthales" = "#B35806",
  "Magnoliids" = "#FFD92F",
  "Monocots" = "#7570B3",
  "Monocots (Commelinids)" = "#54278F",
  "Basal Eudicots" = "#FC8D62",
  "Saxifragales" = "#1B9E77",
  "Rosids (Basal)" = "#8FBC8F",
  "Rosids (Fabids)" = "#66A61E",
  "Rosids (Malvids)" = "#E7298A",
  "Caryophy_Santalales" = "#E31A1C",
  "Asterids (Basal)" = "#FF7F00",
  "Asterids (Lamiids)" = "#377EB8",
  "Asterids (Campanulids)" = "#984EA3",
  "Gymnosperms" = "#000000",
  "Ferns" = "#4D4D4D",
  "Lycophytes" = "#6B6B6B"
)

message("----------------------------------------------------------------")
message(">>> CONFIGURATION: current-run analytical inputs")
message("----------------------------------------------------------------")
message(
  "[CONFIG] Occurrence input: ",
  cfg$tree_occurrence_path %||% cfg$path_lin_enriched %||%
    if (exists("lin_enriched", envir = .GlobalEnv, inherits = FALSE)) {
      "in_memory:lin_enriched"
    } else {
      "automatic discovery"
    }
)
message(
  "[CONFIG] RDKit compounds: ",
  cfg$path_rdkit_compounds %||% cfg$rdkit_compounds_path %||%
    cfg$phylo_rdkit_compounds_path %||% "automatic discovery"
)
message(
  "[CONFIG] RDKit QC: ",
  cfg$path_rdkit_qc %||% cfg$rdkit_qc_path %||%
    cfg$path_rdkit_structural_qc %||%
    cfg$phylo_rdkit_qc_path %||% "automatic discovery"
)
message("[CONFIG] Phylogeny scenario: ", scenario)
message("[CONFIG] Exclude non-vascular families from tree: ", exclude_nonvascular)

message("----------------------------------------------------------------")
message(">>> [1/6] Loading and standardizing analytical occurrence data")
message("----------------------------------------------------------------")

occ_loaded <- load_occurrence_table(
  cfg = cfg,
  project_root = project_root,
  runtime_out_dir = runtime_out_dir,
  runtime_base_tag = runtime_base_tag
)
occ <- standardize_occurrence_columns(occ_loaded$data)

if (nrow(occ) == 0L) stop("No valid family/genus occurrence records were found.")

family_correction_result <- apply_curated_family_corrections(
  occ = occ,
  corrections = curated_family_corrections,
  correction_reasons = curated_family_correction_reasons
)
occ <- family_correction_result$data
family_correction_audit <- family_correction_result$audit
occ_before_primary_filter <- occ

primary_universe <- load_rdkit_primary_universe(
  cfg = cfg,
  project_root = project_root,
  runtime_out_dir = runtime_out_dir
)
primary_filter_result <- filter_primary_structure_occurrences(
  occ = occ,
  primary_universe = primary_universe
)
occ <- primary_filter_result$data
family_primary_structure_audit <- primary_filter_result$family_audit
families_without_primary <- primary_filter_result$families_without_primary
primary_structure_excluded_occurrences <-
  primary_filter_result$excluded_occurrences

if (nrow(occ) == 0L) {
  stop("No occurrence records remain after primary structural-QC filtering.")
}

message("[INFO] Input source: ", occ_loaded$source)
message("[INFO] RDKit compounds: ", primary_universe$compounds_source)
message("[INFO] RDKit structural QC: ", primary_universe$qc_source)
message("[INFO] Valid occurrence rows before structural QC: ", nrow(occ_before_primary_filter))
message("[INFO] Primary-eligible occurrence rows: ", nrow(occ))
message(
  "[INFO] Curated family corrections applied: ",
  sum(family_correction_audit$n_occurrence_rows)
)
message(
  "[INFO] Families removed for absence of primary-eligible compounds: ",
  nrow(families_without_primary)
)
if (nrow(families_without_primary) > 0L) {
  message(
    "[INFO] Removed families: ",
    paste(families_without_primary$family, collapse = ", ")
  )
}
message("[INFO] Analytical families after structural QC: ", length(unique(occ$family)))

message("----------------------------------------------------------------")
message(">>> [2/6] Selecting one deterministic representative per family")
message("----------------------------------------------------------------")

occ_pair_key <- paste(occ$family, occ$genus, sep = "|")
representative_forbidden_audit <- occ[
  occ_pair_key %in% names(representative_forbidden_pairs),
  ,
  drop = FALSE
] |>
  dplyr::group_by(.data$family, .data$genus) |>
  dplyr::summarise(
    n_occurrence_rows = dplyr::n(),
    n_compounds = dplyr::n_distinct(.data$inchikey[!is.na(.data$inchikey)]),
    n_references = dplyr::n_distinct(.data$ref_id[!is.na(.data$ref_id)]),
    .groups = "drop"
  )
if (nrow(representative_forbidden_audit) > 0L) {
  audit_key <- paste(
    representative_forbidden_audit$family,
    representative_forbidden_audit$genus,
    sep = "|"
  )
  representative_forbidden_audit$reason <- unname(
    representative_forbidden_pairs[audit_key]
  )
}

representatives <- select_family_representatives(
  occ,
  forbidden_pairs = representative_forbidden_pairs
)
message("[INFO] Representatives selected: ", nrow(representatives))
message(
  "[INFO] Excluded incompatible representative pairs encountered: ",
  nrow(representative_forbidden_audit)
)
message(
  "[INFO] Observed binomial representatives: ",
  sum(representatives$representative_source == "observed_binomial_species")
)
message(
  "[INFO] Genus-level dummy representatives: ",
  sum(grepl("^genus_level_dummy", representatives$representative_source))
)

message("----------------------------------------------------------------")
message(">>> [3/6] Assigning APG macroclades and vascular status")
message("----------------------------------------------------------------")

apg <- load_apg_table(path_apg)
classified <- assign_apg_clades(
  representatives = representatives,
  apg = apg,
  manual_clade_overrides = manual_clade_overrides,
  nonvascular_families = nonvascular_families
)

unmatched <- classified[!classified$apg_match, , drop = FALSE]
unresolved <- classified[
  !classified$apg_match & classified$clade_assignment_source == "unmatched",
  ,
  drop = FALSE
]
nonvascular <- classified[!classified$is_vascular, , drop = FALSE]

if (nrow(unresolved) > 0L) {
  stop(
    "APG-unresolved analytical families remain: ",
    paste(sort(unique(unresolved$family)), collapse = ", "),
    ". Add an explicit APG/clade assignment before constructing the tree."
  )
}

# APG-matched families can still fall into generic catch-all labels when the
# reference table lacks sufficient node information. These groups are not ANA
# and must be assigned explicitly before they are used in analysis-clade plots.
analysis_unclassified <- classified[
  classified$is_vascular &
    classified$FinalClade %in% c("Other Angiosperms", "Other"),
  ,
  drop = FALSE
]
if (nrow(analysis_unclassified) > 0L) {
  stop(
    "Vascular families lack a specific analytical clade assignment: ",
    paste(
      paste0(
        sort(unique(analysis_unclassified$family)),
        " [",
        analysis_unclassified$FinalClade[
          match(
            sort(unique(analysis_unclassified$family)),
            analysis_unclassified$family
          )
        ],
        "]"
      ),
      collapse = ", "
    ),
    ". Add entries to cfg$phylo_manual_clade_overrides or to ",
    "manual_clade_overrides before constructing the tree."
  )
}

if (exclude_nonvascular) {
  tree_families <- classified[classified$is_vascular, , drop = FALSE]
} else {
  stop(
    "V.PhyloMaker2 is used here for vascular plants. ",
    "Set cfg$phylo_exclude_nonvascular = TRUE for this workflow."
  )
}

if (nrow(tree_families) < 2L) stop("Fewer than two vascular families are available for tree construction.")

message("[INFO] Vascular families retained for tree: ", nrow(tree_families))
message("[INFO] Non-vascular families retained only in audit: ", nrow(nonvascular))
message("[INFO] APG-unmatched families (including explicit overrides): ", nrow(unmatched))
message("[INFO] APG-unresolved families: ", nrow(unresolved))

message("----------------------------------------------------------------")
message(">>> [4/6] Building the vascular family-level visualization tree")
message("----------------------------------------------------------------")

utils::data("GBOTB.extended.TPL", package = "V.PhyloMaker2", envir = environment())
utils::data("nodes.info.1.TPL", package = "V.PhyloMaker2", envir = environment())

if (!exists("GBOTB.extended.TPL", inherits = FALSE)) {
  stop("Dataset GBOTB.extended.TPL was not loaded from V.PhyloMaker2.")
}
if (!exists("nodes.info.1.TPL", inherits = FALSE)) {
  stop("Dataset nodes.info.1.TPL was not loaded from V.PhyloMaker2.")
}

available_backbone_families <- unique(
  normalize_taxon_word(nodes.info.1.TPL$family[nodes.info.1.TPL$level == "F"])
)
available_backbone_families <- available_backbone_families[
  !is.na(available_backbone_families) & nzchar(available_backbone_families)
]

tree_families <- apply_backbone_taxonomy(
  tree_families,
  family_overrides = family_backbone_overrides,
  genus_family_overrides = genus_family_overrides,
  available_backbone_families = available_backbone_families
)

sp_input <- data.frame(
  species = tree_families$tree_species,
  genus = tree_families$genus,
  family = tree_families$backbone_family,
  stringsAsFactors = FALSE
)

if (anyDuplicated(sp_input$species)) stop("Phylogeny input contains duplicate species labels.")
if (anyNA(sp_input$family) || anyNA(sp_input$genus) || anyNA(sp_input$species)) {
  stop("Phylogeny input contains missing species, genus, or family values.")
}

phylo_result <- V.PhyloMaker2::phylo.maker(
  sp.list = sp_input,
  tree = GBOTB.extended.TPL,
  nodes = nodes.info.1.TPL,
  scenarios = scenario
)

maker_species_list <- as_base_df(phylo_result$species.list %||% data.frame())
if (nrow(maker_species_list) > 0L && "species" %in% names(maker_species_list)) {
  maker_species_list$species_normalized <- gsub(
    " ", "_", maker_species_list$species, fixed = TRUE
  )
  idx_status <- match(
    gsub(" ", "_", tree_families$tree_species, fixed = TRUE),
    maker_species_list$species_normalized
  )
  tree_families$phylo_maker_status <- if (
    "status" %in% names(maker_species_list)
  ) {
    maker_species_list$status[idx_status]
  } else {
    NA_character_
  }
} else {
  tree_families$phylo_maker_status <- NA_character_
}

failed_to_bind <- !is.na(tree_families$phylo_maker_status) &
  grepl("fail", tree_families$phylo_maker_status, ignore.case = TRUE)
if (any(failed_to_bind)) {
  stop(
    "V.PhyloMaker2 failed to bind the following analytical families:
- ",
    paste(tree_families$original_family[failed_to_bind], collapse = "
- ")
  )
}

tree_species <- extract_scenario_tree(phylo_result, scenario)
tree_family_raw <- restore_family_tip_labels(tree_species, tree_families)
tree_family_raw <- ape::ladderize(tree_family_raw, right = FALSE)

# The raw S3 tree is retained for transparency. Topology corrections are applied
# to a branch-length-free copy because this family tree is a visualization tree,
# not a calibrated family-age phylogeny.
assignment_for_tree <- data.frame(
  family = classified$family,
  Order = classified$Order,
  FinalClade = classified$FinalClade,
  stringsAsFactors = FALSE
)

order_audit_before <- audit_order_monophyly(tree_family_raw, assignment_for_tree)
local_audit_before <- build_local_neighborhood_audit(tree_family_raw, assignment_for_tree)

topology_result <- apply_topology_corrections(
  tree = tree_family_raw,
  assignments = assignment_for_tree,
  corrections = topology_corrections
)
tree_family <- topology_result$tree
topology_correction_audit <- topology_result$audit
tree_family <- ape::ladderize(tree_family, right = FALSE)

order_audit_after <- audit_order_monophyly(tree_family, assignment_for_tree)
local_audit_after <- build_local_neighborhood_audit(tree_family, assignment_for_tree)

# Enforce the intended result of each curated order-level correction.
for (family in names(topology_corrections)) {
  target_order <- as.character(topology_corrections[[family]]$target_order %||% NA_character_)
  if (!is.na(target_order) && nzchar(target_order)) {
    row <- local_audit_after[local_audit_after$family == family, , drop = FALSE]
    if (nrow(row) != 1L || !isTRUE(row$same_order_sister)) {
      stop(
        "Topology correction validation failed for ", family,
        ": no immediate sister lineage from order ", target_order, "."
      )
    }
  }
}

# Posidoniaceae can appear visually adjacent to magnoliids after circular
# rotation, but its actual topology must remain with Araceae in Alismatales.
posidonia_check <- local_audit_after[
  local_audit_after$family == "Posidoniaceae",
  , drop = FALSE
]
if (nrow(posidonia_check) == 1L && !isTRUE(posidonia_check$same_order_sister)) {
  stop("Posidoniaceae is not placed with another Alismatales family after tree construction.")
}

expected_families <- sort(tree_families$original_family)
observed_families <- sort(tree_family$tip.label)
if (!identical(expected_families, observed_families)) {
  stop("Exported tree tips do not exactly match the retained vascular families.")
}

message("[OK] Tree tips: ", ape::Ntip(tree_family))
message("[INFO] Internal nodes: ", tree_family$Nnode)
message("[INFO] Curated topology corrections: ", nrow(topology_correction_audit))

message("----------------------------------------------------------------")
message(">>> [5/6] Writing tree, iTOL datasets, and audit tables")
message("----------------------------------------------------------------")

scenario_tag <- safe_filename(tolower(scenario))
path_newick <- file.path(out_dir, paste0("family_tree_", scenario_tag, ".nwk"))
path_newick_legacy <- file.path(out_dir, "newick_end.nwk")
path_newick_raw <- file.path(out_dir, paste0("family_tree_", scenario_tag, "_raw_vphylomaker2.nwk"))
path_newick_grafen <- file.path(out_dir, paste0("family_tree_", scenario_tag, "_grafen.nwk"))

# Default file for iTOL: corrected topology without branch lengths. This avoids
# misleading visual emphasis from backbone time lengths and permits a compact
# circular cladogram. Raw and Grafen-length versions are exported separately.
tree_family_cladogram <- tree_family
tree_family_cladogram$edge.length <- NULL
tree_family_grafen <- ape::compute.brlen(tree_family_cladogram, method = "Grafen", power = 1)
ape::write.tree(tree_family_cladogram, file = path_newick)
ape::write.tree(tree_family_cladogram, file = path_newick_legacy)
ape::write.tree(tree_family_raw, file = path_newick_raw)
ape::write.tree(tree_family_grafen, file = path_newick_grafen)

export_df <- data.frame(ID = tree_family$tip.label, stringsAsFactors = FALSE)
idx_export <- match(export_df$ID, classified$family)
export_df$FinalClade <- classified$FinalClade[idx_export]
export_df$AnalysisClade <- collapse_analysis_clade(
  export_df$FinalClade,
  family = export_df$ID
)
export_df$Order <- classified$Order[idx_export]
export_df$apg_match <- classified$apg_match[idx_export]
export_df$DetailedColor <- unname(clade_palette[export_df$FinalClade])
export_df$DetailedColor[is.na(export_df$DetailedColor)] <- unname(clade_palette[["Other"]])
export_df$Color <- unname(analysis_clade_palette[export_df$AnalysisClade])
if (anyNA(export_df$Color)) {
  missing_groups <- unique(export_df$AnalysisClade[is.na(export_df$Color)])
  stop("Missing analysis-clade colors for: ", paste(missing_groups, collapse = ", "))
}

if (anyNA(export_df$FinalClade)) stop("Missing clade assignments after tree export.")

clade_summary <- export_df |>
  dplyr::count(.data$AnalysisClade, name = "n_families", sort = TRUE)
detailed_clade_summary <- export_df |>
  dplyr::count(.data$FinalClade, name = "n_families", sort = TRUE)

backbone_audit <- tree_families[
  !is.na(tree_families$backbone_override_reason),
  c(
    "original_family", "genus", "tree_species", "backbone_family",
    "original_family_in_backbone", "backbone_family_in_backbone",
    "backbone_override_reason", "FinalClade"
  ),
  drop = FALSE
]

tip_audit <- tree_families[
  ,
  c(
    "original_family", "genus", "species_binomial", "tree_species",
    "representative_source", "representative_taxonomy_status",
    "representative_exclusion_reason", "n_occurrence_rows", "n_compounds",
    "n_references", "backbone_family", "original_family_in_backbone",
    "backbone_family_in_backbone", "backbone_override_reason",
    "phylo_maker_status", "Order", "Node.1", "Node.2", "Node.3", "Node.4", "FinalClade",
    "apg_match", "clade_assignment_source"
  ),
  drop = FALSE
]
tip_audit$AnalysisClade <- collapse_analysis_clade(
  tip_audit$FinalClade,
  family = tip_audit$original_family
)

readr::write_csv(
  family_primary_structure_audit,
  file.path(out_dir, "family_primary_structure_coverage.csv"),
  na = ""
)
readr::write_csv(
  families_without_primary,
  file.path(out_dir, "families_in_tree_without_primary_compounds.csv"),
  na = ""
)
readr::write_csv(
  primary_structure_excluded_occurrences,
  file.path(out_dir, "primary_structure_excluded_occurrences.csv"),
  na = ""
)
readr::write_csv(representatives, file.path(out_dir, "family_representatives.csv"), na = "")
readr::write_csv(
  family_correction_audit,
  file.path(out_dir, "phylogeny_taxonomic_family_corrections.csv"),
  na = ""
)
readr::write_csv(
  representative_forbidden_audit,
  file.path(out_dir, "phylogeny_representative_exclusions.csv"),
  na = ""
)
classified$AnalysisClade <- collapse_analysis_clade(
  classified$FinalClade,
  family = classified$family
)
readr::write_csv(classified, file.path(out_dir, "APG_clade_assignments_all_families.csv"), na = "")
readr::write_csv(export_df, file.path(out_dir, "APG_clade_assignments.csv"), na = "")
readr::write_csv(clade_summary, file.path(out_dir, "APG_clade_summary.csv"), na = "")
readr::write_csv(detailed_clade_summary, file.path(out_dir, "APG_detailed_clade_summary.csv"), na = "")
readr::write_csv(unmatched, file.path(out_dir, "APG_unmatched_families.csv"), na = "")
readr::write_csv(unresolved, file.path(out_dir, "APG_unresolved_families.csv"), na = "")
readr::write_csv(
  analysis_unclassified,
  file.path(out_dir, "analysis_clade_unclassified_families.csv"),
  na = ""
)
readr::write_csv(nonvascular, file.path(out_dir, "phylogeny_excluded_nonvascular.csv"), na = "")
readr::write_csv(backbone_audit, file.path(out_dir, "phylogeny_backbone_overrides.csv"), na = "")
readr::write_csv(tip_audit, file.path(out_dir, "phylogeny_tip_audit.csv"), na = "")
readr::write_csv(sp_input, file.path(out_dir, "phylogeny_input_vphylomaker2.csv"), na = "")
readr::write_csv(topology_correction_audit, file.path(out_dir, "phylogeny_topology_corrections.csv"), na = "")
readr::write_csv(order_audit_before, file.path(out_dir, "phylogeny_order_monophyly_before.csv"), na = "")
readr::write_csv(order_audit_after, file.path(out_dir, "phylogeny_order_monophyly_after.csv"), na = "")
readr::write_csv(local_audit_before, file.path(out_dir, "phylogeny_local_neighborhood_before.csv"), na = "")
readr::write_csv(local_audit_after, file.path(out_dir, "phylogeny_local_neighborhood_after.csv"), na = "")

# Default iTOL files use the manuscript-level analysis macroclades. Detailed APG
# files are also exported for taxonomic auditing.
analysis_export <- data.frame(
  ID = export_df$ID,
  FinalClade = export_df$AnalysisClade,
  Color = export_df$Color,
  stringsAsFactors = FALSE
)
write_tree_colors(analysis_export, file.path(out_dir, "iTOL_tree_colors.txt"))
write_color_strip(
  analysis_export, analysis_clade_palette,
  file.path(out_dir, "iTOL_clade_strip.txt"),
  dataset_label = "Analysis macroclades",
  legend_title = "Analysis macroclades",
  strip_width = 70
)

detailed_export <- data.frame(
  ID = export_df$ID,
  FinalClade = export_df$FinalClade,
  Color = export_df$DetailedColor,
  stringsAsFactors = FALSE
)
write_tree_colors(detailed_export, file.path(out_dir, "iTOL_tree_colors_detailed_APG.txt"))
write_color_strip(
  detailed_export, clade_palette,
  file.path(out_dir, "iTOL_clade_strip_detailed_APG.txt"),
  dataset_label = "Detailed APG clades",
  legend_title = "Detailed APG clades",
  strip_width = 70
)

message("----------------------------------------------------------------")
message(">>> [6/6] Writing reproducibility metadata")
message("----------------------------------------------------------------")

metadata <- list(
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  input_source = occ_loaded$source,
  rdkit_compounds_source = primary_universe$compounds_source,
  rdkit_structural_qc_source = primary_universe$qc_source,
  apg_reference = path_apg,
  output_directory = out_dir,
  phylogeny = list(
    package = "V.PhyloMaker2",
    package_version = as.character(utils::packageVersion("V.PhyloMaker2")),
    backbone = "GBOTB.extended.TPL",
    nodes = "nodes.info.1.TPL",
    scenario = scenario,
    unit = "one deterministic representative per analytical family",
    intended_use = "family-level visualization and annotation; not family crown-age comparative modelling"
  ),
  counts = list(
    occurrence_rows_before_primary_filter = nrow(occ_before_primary_filter),
    occurrence_rows_primary_eligible = nrow(occ),
    rdkit_compounds = length(primary_universe$all_rdkit_keys),
    rdkit_primary_structure_eligible = length(primary_universe$eligible_keys),
    rdkit_excluded_primary = nrow(primary_universe$excluded_details),
    families_without_primary_compounds = nrow(families_without_primary),
    curated_family_correction_pairs = nrow(family_correction_audit),
    curated_family_corrected_rows = sum(family_correction_audit$n_occurrence_rows),
    analytical_families = nrow(classified),
    vascular_tree_families = nrow(tree_families),
    nonvascular_excluded_families = nrow(nonvascular),
    apg_unmatched_families = nrow(unmatched),
    apg_unresolved_families = nrow(unresolved),
    analysis_clade_unclassified_families = nrow(analysis_unclassified),
    observed_binomial_representatives = sum(
      representatives$representative_source == "observed_binomial_species"
    ),
    genus_dummy_representatives = sum(
      grepl("^genus_level_dummy", representatives$representative_source)
    ),
    excluded_incompatible_representative_pairs = nrow(
      representative_forbidden_audit
    ),
    backbone_overrides = nrow(backbone_audit),
    phylo_maker_bound_tips = sum(tree_families$phylo_maker_status == "bind", na.rm = TRUE),
    phylo_maker_pruned_tips = sum(tree_families$phylo_maker_status == "prune", na.rm = TRUE),
    tree_tips = ape::Ntip(tree_family),
    tree_internal_nodes = tree_family$Nnode,
    curated_topology_corrections = nrow(topology_correction_audit),
    order_review_flags_before = sum(order_audit_before$review_flag, na.rm = TRUE),
    order_review_flags_after = sum(order_audit_after$review_flag, na.rm = TRUE),
    analysis_macroclades = length(unique(export_df$AnalysisClade))
  ),
  rules = list(
    family_tree_scope = "families with at least one primary_structure_eligible RDKit compound",
    unlisted_structural_qc_compounds = "eligible when present as valid compounds in the complete RDKit compound table",
    nonvascular_excluded = exclude_nonvascular,
    representative_ranking = c(
      "descending distinct compounds",
      "descending occurrence rows",
      "descending references",
      "alphabetical taxon name"
    ),
    curated_family_corrections = as.list(curated_family_corrections),
    curated_family_correction_reasons = as.list(curated_family_correction_reasons),
    family_backbone_overrides = as.list(family_backbone_overrides),
    genus_family_overrides = as.list(genus_family_overrides),
    representative_forbidden_pairs = as.list(representative_forbidden_pairs),
    manual_clade_overrides = as.list(manual_clade_overrides),
    topology_corrections = topology_corrections,
    default_tree_export = "corrected cladogram without branch lengths",
    detailed_and_analysis_clades_exported = TRUE
  ),
  software = list(
    R = R.version.string,
    ape = as.character(utils::packageVersion("ape")),
    dplyr = as.character(utils::packageVersion("dplyr")),
    readr = as.character(utils::packageVersion("readr"))
  )
)

write_metadata(metadata, file.path(out_dir, "tree_apg_itol_run_metadata.json"))

message("------------------------------------------------------------")
message("[TREE/APG SUMMARY]")
message(
  "Curated family corrections: ",
  nrow(family_correction_audit),
  " pair(s), ",
  sum(family_correction_audit$n_occurrence_rows),
  " occurrence row(s)"
)
message("Analytical families: ", nrow(classified))
message("Vascular families in tree: ", nrow(tree_families))
message("Non-vascular families excluded from tree: ", nrow(nonvascular))
message("APG-unmatched families (including explicit overrides): ", nrow(unmatched))
message("APG-unresolved families: ", nrow(unresolved))
message(
  "Excluded incompatible representative pairs: ",
  nrow(representative_forbidden_audit)
)
message("Backbone substitutions: ", nrow(backbone_audit))
message("Curated topology corrections: ", nrow(topology_correction_audit))
message("Order-level review flags before correction: ", sum(order_audit_before$review_flag, na.rm = TRUE))
message("Order-level review flags after correction: ", sum(order_audit_after$review_flag, na.rm = TRUE))
message("Analysis macroclades: ", length(unique(export_df$AnalysisClade)))
message("Tree (corrected cladogram): ", path_newick)
message(
  "Families excluded for no primary compounds: ",
  file.path(out_dir, "families_in_tree_without_primary_compounds.csv")
)
message("Output directory: ", out_dir)
message("------------------------------------------------------------")
