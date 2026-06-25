
`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0) return(b)
  if (is.character(a) && length(a) == 1 && !nzchar(a)) return(b)
  a
}

as_base_df <- function(x) as.data.frame(x, stringsAsFactors = FALSE)

require_pkgs <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing) > 0) {
    stop(
      "Missing required packages: ", paste(missing, collapse = ", "), "\n",
      "Install them with:\n",
      "install.packages(c(", paste(sprintf('\"%s\"', missing), collapse = ", "), "), dependencies = TRUE)\n"
    )
  }
  invisible(TRUE)
}

required_pkgs <- c(
  "dplyr", "tidyr", "readr", "stringr", "readxl", "writexl",
  "httr", "jsonlite", "progress", "mongolite"
)
require_pkgs(required_pkgs)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(readxl)
  library(writexl)
  library(httr)
  library(jsonlite)
  library(progress)
  library(mongolite)
})

if (!exists("cfg")) stop("Part II requires object 'cfg' created by Main_Pipeline_end.R.")
if (!exists("runtime")) stop("Part II requires object 'runtime' created by Main_Pipeline_end.R.")

OUT_DIR <- runtime$OUT_DIR %||% if (exists("OUT_DIR")) OUT_DIR else NA_character_
base_tag <- runtime$base_tag %||% if (exists("base_tag")) base_tag else NA_character_

if (is.na(OUT_DIR) || !nzchar(OUT_DIR) || !dir.exists(OUT_DIR)) {
  stop("OUT_DIR does not exist: ", OUT_DIR)
}
if (is.na(base_tag) || !nzchar(base_tag)) stop("base_tag is missing or empty.")

out2_dir <- file.path(OUT_DIR, "PartII_ALL")
cache_dir <- file.path(out2_dir, "cache")
dir.create(out2_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

out2_path <- function(x) file.path(out2_dir, x)
cache_path <- function(x) file.path(cache_dir, x)


flav_classes <- cfg$flav_classes %||% c(
  "Flavanones", "Isoflavones", "Flavan-3-ols", "Chalcones",
  "Dihydroflavonols", "Pterocarpan", "Proanthocyanins",
  "Isoflavanones", "2-arylbenzofurans", "Flavandiols", "Rotenoids",
  "Aurones", "Flavans", "Coumestan", "Anthocyanidins",
  "Neoflavonoids", "Flavonolignans"
)

ACTIVITY_CUTOFF_NM <- as.numeric(cfg$activity_cutoff_nm %||% 10000)
VALID_TYPES <- unique(as.character(cfg$chembl_valid_types %||% c(
  "IC50", "Ki", "EC50", "Kd", "MIC", "AC50", "GI50"
)))


PCHEMBL_TYPES <- unique(as.character(cfg$chembl_pchembl_types %||% c(
  "IC50", "XC50", "EC50", "AC50", "Ki", "Kd", "Potency", "ED50"
)))

CHEMBL_BASE_URL <- cfg$chembl_base_url %||% "https://www.ebi.ac.uk/chembl/api/data"
UNICHEM_LEGACY_URL <- cfg$unichem_legacy_url %||% "https://www.ebi.ac.uk/unichem/rest/inchikey"
CHEMBL_BATCH_SIZE <- as.integer(cfg$chembl_batch_size %||% 30L)
CHEMBL_PAGE_LIMIT <- as.integer(cfg$chembl_page_limit %||% 1000L)
CHEMBL_PAUSE_SECONDS <- as.numeric(cfg$chembl_pause_seconds %||% 0.15)
CHEMBL_REFRESH_MAPPING <- isTRUE(cfg$chembl_refresh_mapping %||% FALSE)
CHEMBL_REFRESH_ACTIVITIES <- isTRUE(cfg$chembl_refresh_activities %||% FALSE)
CHEMBL_REFRESH_METADATA <- isTRUE(cfg$chembl_refresh_metadata %||% FALSE)

CHEMBL_USE_UNICHEM_FALLBACK <- isTRUE(cfg$chembl_use_unichem_fallback %||% FALSE)
RUN_GLOBAL_CONTEXT <- isTRUE(cfg$run_global_context %||% TRUE)
GLOBAL_COLLECTION <- cfg$global_context_collection %||% "lotusUniqueNaturalProduct"
ALLOW_VALIDITY_COMMENTS <- as.character(
  cfg$chembl_allowed_validity_comments %||% c(NA_character_, "", "Manually validated")
)

path_tax_excel <- cfg$path_lotus %||% file.path(OUT_DIR, paste0(base_tag, ".xlsx"))


ACTIVITY_API_FIELDS <- c(
  "activity_id", "molecule_chembl_id", "standard_type", "standard_relation",
  "standard_value", "standard_units", "pchembl_value", "data_validity_comment",
  "potential_duplicate", "assay_chembl_id", "assay_description", "assay_type",
  "bao_label", "target_chembl_id", "target_pref_name", "target_organism",
  "target_type", "document_chembl_id", "src_id", "activity_comment"
)
ASSAY_API_FIELDS <- c(
  "assay_chembl_id", "assay_type", "confidence_score", "description",
  "relationship_type", "assay_organism", "bao_format", "cell_chembl_id",
  "tissue_chembl_id", "src_id"
)
DOCUMENT_API_FIELDS <- c(
  "document_chembl_id", "doi", "pubmed_id", "year", "title", "journal",
  "doc_type", "src_id"
)
TARGET_API_FIELDS <- c(
  "target_chembl_id", "pref_name", "target_type", "organism", "tax_id"
)


clean_chr <- function(x) {
  x <- stringr::str_squish(as.character(x))
  x[is.na(x) | x == "" | toupper(x) == "NA"] <- NA_character_
  x
}

safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

collapse_unique <- function(x, sep = ";") {
  x <- clean_chr(x)
  x <- sort(unique(x[!is.na(x)]))
  if (!length(x)) return(NA_character_)
  paste(x, collapse = sep)
}

first_non_missing <- function(x) {
  x <- x[!is.na(x) & nzchar(as.character(x))]
  if (!length(x)) return(NA)
  x[[1]]
}

split_batches <- function(x, size) {
  x <- unique(x[!is.na(x) & nzchar(as.character(x))])
  if (!length(x)) return(list())
  split(x, ceiling(seq_along(x) / size))
}

empty_df <- function(...) data.frame(..., stringsAsFactors = FALSE)

scalarize_for_storage <- function(df) {
  df <- as_base_df(df)
  if (!ncol(df)) return(df)

  for (nm in names(df)) {
    if (is.data.frame(df[[nm]]) || is.list(df[[nm]])) {
      df[[nm]] <- vapply(
        df[[nm]],
        function(value) {
          if (is.null(value) || length(value) == 0) return(NA_character_)
          if (length(value) == 1 && is.atomic(value)) return(as.character(value))
          jsonlite::toJSON(value, auto_unbox = TRUE, null = "null", na = "null")
        },
        FUN.VALUE = character(1)
      )
    }
  }
  df
}


normalize_api_frame <- function(df, fields = NULL) {
  df <- scalarize_for_storage(df)
  n <- nrow(df)

  if (!is.null(fields) && length(fields)) {
    for (nm in setdiff(fields, names(df))) {
      df[[nm]] <- rep(NA_character_, n)
    }
    df <- df[, fields, drop = FALSE]
  }

  for (nm in names(df)) {
    df[[nm]] <- as.character(df[[nm]])
  }

  as_base_df(df)
}

write_table_cache <- function(df, stem) {
  df <- scalarize_for_storage(df)
  csv_path <- cache_path(paste0(stem, ".csv.gz"))
  readr::write_csv(df, csv_path, na = "")

  if (requireNamespace("arrow", quietly = TRUE)) {
    parquet_path <- cache_path(paste0(stem, ".parquet"))
    tryCatch(
      arrow::write_parquet(df, parquet_path),
      error = function(e) {
        warning(
          "Parquet cache was skipped for '", stem, "': ", e$message,
          ". The CSV cache remains available."
        )
        if (file.exists(parquet_path)) unlink(parquet_path)
      }
    )
  }
  invisible(csv_path)
}

read_table_cache <- function(stem) {
  parquet_path <- cache_path(paste0(stem, ".parquet"))
  csv_path <- cache_path(paste0(stem, ".csv.gz"))

  if (file.exists(parquet_path) && requireNamespace("arrow", quietly = TRUE)) {
    parquet_result <- tryCatch(
      as_base_df(arrow::read_parquet(parquet_path, as_data_frame = TRUE)),
      error = function(e) {
        warning("Could not read Parquet cache '", basename(parquet_path),
                "'; falling back to CSV: ", e$message)
        NULL
      }
    )
    if (!is.null(parquet_result)) return(parquet_result)
  }
  if (file.exists(csv_path)) {
    return(as_base_df(readr::read_csv(csv_path, show_col_types = FALSE, progress = FALSE)))
  }
  NULL
}

write_primary_table <- function(df, filename_stem) {
  df <- scalarize_for_storage(df)
  csv_path <- out2_path(paste0(filename_stem, ".csv.gz"))
  readr::write_csv(df, csv_path, na = "")
  if (requireNamespace("arrow", quietly = TRUE)) {
    parquet_path <- out2_path(paste0(filename_stem, ".parquet"))
    tryCatch(
      arrow::write_parquet(df, parquet_path),
      error = function(e) {
        warning(
          "Parquet output was skipped for '", filename_stem, "': ", e$message,
          ". The CSV output remains available."
        )
        if (file.exists(parquet_path)) unlink(parquet_path)
      }
    )
  }
  invisible(csv_path)
}

api_get_json <- function(url, query = list(), timeout_seconds = 90) {
  response <- tryCatch(
    httr::RETRY(
      verb = "GET",
      url = url,
      query = query,
      times = 5,
      pause_base = 1,
      pause_cap = 30,
      terminate_on = c(400, 401, 403, 404),
      httr::timeout(timeout_seconds),
      httr::user_agent("ComparativeFlavonoidPipeline/2.0")
    ),
    error = function(e) e
  )

  if (inherits(response, "error")) {
    return(list(ok = FALSE, status = NA_integer_, error = response$message, json = NULL))
  }
  status <- httr::status_code(response)
  if (status < 200 || status >= 300) {
    return(list(
      ok = FALSE,
      status = status,
      error = httr::content(response, as = "text", encoding = "UTF-8"),
      json = NULL
    ))
  }

  text <- httr::content(response, as = "text", encoding = "UTF-8")
  parsed <- tryCatch(
    jsonlite::fromJSON(text, flatten = TRUE, simplifyVector = TRUE),
    error = function(e) e
  )
  if (inherits(parsed, "error")) {
    return(list(ok = FALSE, status = status, error = parsed$message, json = NULL))
  }
  list(ok = TRUE, status = status, error = NULL, json = parsed)
}

chembl_fetch_pages <- function(resource, object_key, query = list(), limit = CHEMBL_PAGE_LIMIT, only_fields = NULL) {
  endpoint <- paste0(CHEMBL_BASE_URL, "/", resource, ".json")
  offset <- 0L
  pages <- list()
  page_no <- 0L
  total_count <- NA_integer_

  repeat {
    page_no <- page_no + 1L
    q <- c(query, list(limit = limit, offset = offset))
    if (!is.null(only_fields) && length(only_fields)) {
      q$only <- paste(unique(only_fields), collapse = ",")
    }
    result <- api_get_json(endpoint, q)
    if (!isTRUE(result$ok)) {
      stop(
        "ChEMBL request failed [", resource, ", offset=", offset, "]: ",
        result$error %||% paste0("HTTP ", result$status)
      )
    }

    json <- result$json
    objects <- json[[object_key]]
    page_df <- if (is.null(objects) || length(objects) == 0) {
      normalize_api_frame(data.frame(stringsAsFactors = FALSE), only_fields)
    } else {
      normalize_api_frame(as_base_df(objects), only_fields)
    }
    pages[[page_no]] <- page_df

    if (!is.null(json$page_meta$total_count)) {
      total_count <- as.integer(json$page_meta$total_count)
    }

    n_page <- nrow(page_df)
    offset <- offset + n_page
    next_url <- json$page_meta[["next"]]
    has_next <- !is.null(next_url) && length(next_url) == 1L && !is.na(next_url) && nzchar(next_url)

    if (n_page == 0 || !has_next || (!is.na(total_count) && offset >= total_count)) break
    Sys.sleep(CHEMBL_PAUSE_SECONDS)
  }

  normalize_api_frame(dplyr::bind_rows(pages), only_fields)
}


load_lotus_table <- function(object_name, sheet_name) {
  if (exists(object_name, inherits = TRUE)) {
    return(as_base_df(get(object_name, inherits = TRUE)))
  }
  if (!file.exists(path_tax_excel)) {
    stop("LOTUS workbook not found and object '", object_name, "' is absent: ", path_tax_excel)
  }
  sheets <- readxl::excel_sheets(path_tax_excel)
  if (!sheet_name %in% sheets) stop("Sheet '", sheet_name, "' was not found in: ", path_tax_excel)
  as_base_df(readxl::read_excel(path_tax_excel, sheet = sheet_name, .name_repair = "minimal"))
}

lin_source <- load_lotus_table("lin_enriched", "lin_enriched")
uni_source <- load_lotus_table("uni_enriched", "uni_enriched")

required_lin <- c(
  "inchikey", "smiles", "chemicalTaxonomyNPclassifierClass",
  "family", "genus", "species"
)
missing_lin <- setdiff(required_lin, names(lin_source))
if (length(missing_lin)) {
  stop("lin_enriched is missing required columns: ", paste(missing_lin, collapse = ", "))
}


message("----------------------------------------------------------------")
message(">>> STEP A: flavonoid occurrence and compound tables for RDKit")
message("----------------------------------------------------------------")

lin_prepared <- as_base_df(
  lin_source %>%
    dplyr::mutate(
      inchikey = toupper(clean_chr(inchikey)),
      smiles = clean_chr(smiles),
      class_np = clean_chr(chemicalTaxonomyNPclassifierClass),
      class_np = dplyr::recode(
        class_np,
        "Flavandiols (Leucoanthocyanidins)" = "Flavandiols"
      ),
      family = clean_chr(family),
      genus = clean_chr(genus),
      species = clean_chr(species)
    )
)

flav_occurrences <- as_base_df(
  lin_prepared %>%
    dplyr::filter(!is.na(inchikey), !is.na(smiles), !is.na(class_np)) %>%
    dplyr::filter(class_np %in% flav_classes) %>%
    dplyr::select(
      dplyr::any_of(c(
        "lotus_id", "inchikey", "smiles", "class_np", "family", "genus",
        "species", "ref_id", "source", "accepted_rank", "tax_status",
        "original_name", "matched_name", "accepted_name"
      ))
    ) %>%
    dplyr::distinct()
)

if (nrow(flav_occurrences) == 0) {
  stop("The 17-class flavonoid filter returned zero occurrence rows.")
}

rdkit_input_path <- out2_path(paste0(base_tag, "__flavonoids_for_rdkit.csv"))
readr::write_csv(flav_occurrences, rdkit_input_path, na = "")

compound_taxonomy_index <- as_base_df(
  flav_occurrences %>%
    dplyr::group_by(inchikey) %>%
    dplyr::summarise(
      
      n_occurrence_rows = dplyr::n(),
      n_smiles = dplyr::n_distinct(smiles[!is.na(smiles) & nzchar(smiles)]),
      n_classes = dplyr::n_distinct(class_np[!is.na(class_np) & nzchar(class_np)]),
      n_families = dplyr::n_distinct(family[!is.na(family) & nzchar(family)]),
      n_genera = dplyr::n_distinct(genus[!is.na(genus) & nzchar(genus)]),
      n_species = dplyr::n_distinct(
        species[
          !is.na(species) & nzchar(species) &
            stringr::str_detect(species, "\\s")
        ]
      ),
      n_references = if ("ref_id" %in% names(flav_occurrences)) {
        dplyr::n_distinct(ref_id[!is.na(ref_id) & nzchar(ref_id)])
      } else {
        NA_integer_
      },
      smiles = first_non_missing(smiles),
      class_np = collapse_unique(class_np),
      family = collapse_unique(family),
      genus = collapse_unique(genus),
      species = collapse_unique(species),
      ref_ids = if ("ref_id" %in% names(flav_occurrences)) collapse_unique(ref_id, sep = " | ") else NA_character_,
      .groups = "drop"
    )
)

compound_family_map <- as_base_df(
  flav_occurrences %>%
    dplyr::filter(!is.na(family)) %>%
    dplyr::distinct(inchikey, family)
)
compound_class_map <- as_base_df(
  flav_occurrences %>%
    dplyr::filter(!is.na(class_np)) %>%
    dplyr::distinct(inchikey, class_np)
)
compound_species_map <- as_base_df(
  flav_occurrences %>%
    dplyr::filter(!is.na(species), stringr::str_detect(species, "\\s")) %>%
    dplyr::distinct(inchikey, family, genus, species)
)

write_primary_table(compound_taxonomy_index, paste0(base_tag, "_compound_taxonomy_index"))
write_primary_table(compound_family_map, paste0(base_tag, "_compound_family_map"))
write_primary_table(compound_class_map, paste0(base_tag, "_compound_class_map"))
write_primary_table(compound_species_map, paste0(base_tag, "_compound_species_map"))

class_filter_audit <- as_base_df(
  lin_prepared %>%
    dplyr::count(class_np, name = "n_rows", sort = TRUE) %>%
    dplyr::mutate(included_in_minor_flavonoid_set = class_np %in% flav_classes)
)
readr::write_csv(class_filter_audit, out2_path("class_filter_audit.csv"), na = "")

message("[OK] RDKit occurrence input: ", rdkit_input_path)
message("[INFO] Occurrence rows: ", nrow(flav_occurrences))
message("[INFO] Unique InChIKeys: ", nrow(compound_taxonomy_index))


message("----------------------------------------------------------------")
message(">>> STEP B1: InChIKey to ChEMBL mapping")
message("----------------------------------------------------------------")

inchikeys_to_map <- sort(unique(compound_taxonomy_index$inchikey))
map_cache_stem <- "chembl_inchikey_map"
chembl_map_cache <- if (!CHEMBL_REFRESH_MAPPING) read_table_cache(map_cache_stem) else NULL

if (is.null(chembl_map_cache)) {
  chembl_map_cache <- empty_df(
    inchikey = character(),
    chembl_id = character(),
    mapping_method = character(),
    mapping_status = character()
  )
}

chembl_map_cache <- as_base_df(
  chembl_map_cache %>%
    dplyr::mutate(
      inchikey = toupper(clean_chr(inchikey)),
      chembl_id = clean_chr(chembl_id)
    ) %>%
    dplyr::filter(!is.na(inchikey)) %>%
    dplyr::distinct(inchikey, chembl_id, .keep_all = TRUE)
)

map_via_chembl_molecule_endpoint <- function(keys) {
  batches <- split_batches(keys, CHEMBL_BATCH_SIZE)
  out <- list()

  if (!length(batches)) return(empty_df(inchikey = character(), chembl_id = character()))

  for (i in seq_along(batches)) {
    keys_batch <- batches[[i]]
    fetched <- tryCatch(
      chembl_fetch_pages(
        resource = "molecule",
        object_key = "molecules",
        query = list(
          molecule_structures__standard_inchi_key__in = paste(keys_batch, collapse = ",")
        )
      ),
      error = function(e) NULL
    )

    if (!is.null(fetched) && nrow(fetched)) {
      key_candidates <- c(
        "molecule_structures.standard_inchi_key",
        "molecule_structures_standard_inchi_key",
        "standard_inchi_key"
      )
      key_col <- key_candidates[key_candidates %in% names(fetched)]
      if (length(key_col) > 0 && "molecule_chembl_id" %in% names(fetched)) {
        key_col <- key_col[[1]]
        out[[length(out) + 1L]] <- as_base_df(
          fetched %>%
            dplyr::transmute(
              inchikey = toupper(clean_chr(.data[[key_col]])),
              chembl_id = clean_chr(molecule_chembl_id),
              mapping_method = "ChEMBL molecule endpoint",
              mapping_status = "mapped"
            ) %>%
            dplyr::filter(!is.na(inchikey), !is.na(chembl_id))
        )
      }
    }
    Sys.sleep(CHEMBL_PAUSE_SECONDS)
  }
  as_base_df(dplyr::bind_rows(out))
}

map_one_unichem_legacy <- function(inchikey) {
  result <- api_get_json(paste0(UNICHEM_LEGACY_URL, "/", inchikey), timeout_seconds = 30)
  if (!isTRUE(result$ok) || is.null(result$json)) return(NULL)

  json <- result$json
  if (!is.data.frame(json)) {
    json <- tryCatch(as_base_df(json), error = function(e) NULL)
  }
  if (is.null(json) || !nrow(json)) return(NULL)

  source_col <- c("src_id", "source_id")[c("src_id", "source_id") %in% names(json)]
  id_col <- c("src_compound_id", "source_compound_id")[c("src_compound_id", "source_compound_id") %in% names(json)]
  if (!length(source_col) || !length(id_col)) return(NULL)
  source_col <- source_col[[1]]
  id_col <- id_col[[1]]

  hit <- json[as.character(json[[source_col]]) == "1", , drop = FALSE]
  if (!nrow(hit)) return(NULL)
  empty_df(
    inchikey = inchikey,
    chembl_id = as.character(hit[[id_col]][1]),
    mapping_method = "UniChem legacy",
    mapping_status = "mapped"
  )
}

mapped_keys <- unique(chembl_map_cache$inchikey)
missing_keys <- setdiff(inchikeys_to_map, mapped_keys)

if (length(missing_keys)) {
  message("[INFO] Mapping ", length(missing_keys), " previously uncached InChIKeys.")

  direct_map <- map_via_chembl_molecule_endpoint(missing_keys)
  direct_mapped_keys <- unique(direct_map$inchikey)
  missing_after_direct <- setdiff(missing_keys, direct_mapped_keys)

  if (nrow(direct_map)) {
    chembl_map_cache <- as_base_df(dplyr::bind_rows(chembl_map_cache, direct_map))
  }

  fallback_df <- empty_df(
    inchikey = character(),
    chembl_id = character(),
    mapping_method = character(),
    mapping_status = character()
  )

  if (length(missing_after_direct) && isTRUE(CHEMBL_USE_UNICHEM_FALLBACK)) {
    message(
      "[INFO] Optional UniChem fallback enabled for ",
      length(missing_after_direct),
      " InChIKeys."
    )
    pb <- utils::txtProgressBar(min = 0, max = length(missing_after_direct), style = 3)
    fallback_rows <- list()

    for (i in seq_along(missing_after_direct)) {
      key <- missing_after_direct[[i]]
      hit <- map_one_unichem_legacy(key)
      if (!is.null(hit)) fallback_rows[[length(fallback_rows) + 1L]] <- hit

      if (i %% 100L == 0L || i == length(missing_after_direct)) {
        partial <- as_base_df(dplyr::bind_rows(fallback_rows))
        checkpoint <- as_base_df(dplyr::bind_rows(chembl_map_cache, partial))
        write_table_cache(checkpoint, map_cache_stem)
      }

      utils::setTxtProgressBar(pb, i)
      Sys.sleep(CHEMBL_PAUSE_SECONDS)
    }
    close(pb)

    fallback_df <- as_base_df(dplyr::bind_rows(fallback_rows))
    if (nrow(fallback_df)) {
      chembl_map_cache <- as_base_df(dplyr::bind_rows(chembl_map_cache, fallback_df))
    }
  } else if (length(missing_after_direct)) {
    message(
      "[INFO] UniChem fallback disabled. Recording ",
      length(missing_after_direct),
      " compounds as not mapped by exact ChEMBL InChIKey lookup."
    )
  }

  mapped_after_all_methods <- unique(c(direct_mapped_keys, fallback_df$inchikey))
  still_unmapped <- setdiff(missing_keys, mapped_after_all_methods)

  if (length(still_unmapped)) {
    unmapped_rows <- empty_df(
      inchikey = still_unmapped,
      chembl_id = NA_character_,
      mapping_method = "ChEMBL exact Standard InChIKey",
      mapping_status = if (isTRUE(CHEMBL_USE_UNICHEM_FALLBACK)) {
        "not_mapped_after_unichem"
      } else {
        "not_mapped_exact_inchikey"
      }
    )
    chembl_map_cache <- as_base_df(dplyr::bind_rows(chembl_map_cache, unmapped_rows))
  }

  chembl_map_cache <- as_base_df(
    chembl_map_cache %>%
      dplyr::arrange(inchikey, dplyr::desc(!is.na(chembl_id))) %>%
      dplyr::distinct(inchikey, chembl_id, .keep_all = TRUE)
  )
  write_table_cache(chembl_map_cache, map_cache_stem)
}

chembl_map <- as_base_df(
  chembl_map_cache %>%
    dplyr::filter(inchikey %in% inchikeys_to_map, !is.na(chembl_id)) %>%
    dplyr::distinct(inchikey, chembl_id, .keep_all = TRUE)
)

unmapped <- setdiff(inchikeys_to_map, unique(chembl_map$inchikey))
chembl_mapping_full <- as_base_df(
  dplyr::bind_rows(
    chembl_map,
    empty_df(
      inchikey = unmapped,
      chembl_id = NA_character_,
      mapping_method = NA_character_,
      mapping_status = "not_mapped"
    )
  ) %>%
    dplyr::arrange(inchikey, chembl_id)
)

write_table_cache(chembl_mapping_full, map_cache_stem)
readr::write_csv(chembl_mapping_full, out2_path("chembl_inchikey_map.csv"), na = "")

mapping_summary <- empty_df(
  metric = c("input_inchikeys", "mapped_inchikeys", "unmapped_inchikeys", "chembl_ids"),
  n = c(
    length(inchikeys_to_map),
    dplyr::n_distinct(chembl_map$inchikey),
    length(unmapped),
    dplyr::n_distinct(chembl_map$chembl_id)
  )
)
readr::write_csv(mapping_summary, out2_path("chembl_mapping_summary.csv"), na = "")
message("[INFO] Mapped InChIKeys: ", dplyr::n_distinct(chembl_map$inchikey), " / ", length(inchikeys_to_map))


message("----------------------------------------------------------------")
message(">>> STEP B2: ChEMBL activity acquisition")
message("----------------------------------------------------------------")

activity_cache_stem <- "chembl_activities_raw"
raw_activities_cache <- if (!CHEMBL_REFRESH_ACTIVITIES) read_table_cache(activity_cache_stem) else NULL
if (is.null(raw_activities_cache)) raw_activities_cache <- data.frame(stringsAsFactors = FALSE)
raw_activities_cache <- normalize_api_frame(raw_activities_cache, ACTIVITY_API_FIELDS)

chembl_ids_needed <- sort(unique(chembl_map$chembl_id))
chembl_ids_cached <- if (nrow(raw_activities_cache) && "molecule_chembl_id" %in% names(raw_activities_cache)) {
  unique(clean_chr(raw_activities_cache$molecule_chembl_id))
} else {
  character(0)
}
chembl_ids_missing <- setdiff(chembl_ids_needed, chembl_ids_cached)

fetch_activity_batch <- function(ids) {
  chembl_fetch_pages(
    resource = "activity",
    object_key = "activities",
    query = list(
      molecule_chembl_id__in = paste(ids, collapse = ","),
      standard_type__in = paste(VALID_TYPES, collapse = ",")
    ),
    only_fields = ACTIVITY_API_FIELDS
  )
}

if (length(chembl_ids_missing)) {
  message("[INFO] Retrieving activities for ", length(chembl_ids_missing), " uncached ChEMBL IDs.")
  batches <- split_batches(chembl_ids_missing, CHEMBL_BATCH_SIZE)
  new_activity_rows <- list()
  pb <- utils::txtProgressBar(min = 0, max = length(batches), style = 3)

  for (i in seq_along(batches)) {
    batch_df <- tryCatch(
      fetch_activity_batch(batches[[i]]),
      error = function(e) {
        warning("Activity batch ", i, " failed: ", e$message)
        data.frame(stringsAsFactors = FALSE)
      }
    )
    batch_df <- normalize_api_frame(batch_df, ACTIVITY_API_FIELDS)
    if (nrow(batch_df)) new_activity_rows[[length(new_activity_rows) + 1L]] <- batch_df

    if (i %% 5L == 0L || i == length(batches)) {
      new_activity_df <- if (length(new_activity_rows)) {
        normalize_api_frame(dplyr::bind_rows(new_activity_rows), ACTIVITY_API_FIELDS)
      } else {
        normalize_api_frame(data.frame(stringsAsFactors = FALSE), ACTIVITY_API_FIELDS)
      }
      checkpoint <- normalize_api_frame(
        dplyr::bind_rows(raw_activities_cache, new_activity_df),
        ACTIVITY_API_FIELDS
      )
      write_table_cache(checkpoint, activity_cache_stem)
    }
    utils::setTxtProgressBar(pb, i)
    Sys.sleep(CHEMBL_PAUSE_SECONDS)
  }
  close(pb)

  new_activity_df <- if (length(new_activity_rows)) {
    normalize_api_frame(dplyr::bind_rows(new_activity_rows), ACTIVITY_API_FIELDS)
  } else {
    normalize_api_frame(data.frame(stringsAsFactors = FALSE), ACTIVITY_API_FIELDS)
  }
  raw_activities_cache <- normalize_api_frame(
    dplyr::bind_rows(raw_activities_cache, new_activity_df),
    ACTIVITY_API_FIELDS
  )
}

raw_activities <- normalize_api_frame(raw_activities_cache, ACTIVITY_API_FIELDS)
if (nrow(raw_activities)) {
  raw_activities <- as_base_df(
    raw_activities %>%
      dplyr::filter(molecule_chembl_id %in% chembl_ids_needed) %>%
      dplyr::distinct()
  )
}
write_table_cache(raw_activities, activity_cache_stem)
write_primary_table(raw_activities, paste0(base_tag, "_chembl_activities_raw"))
message("[INFO] Raw activity records retained: ", nrow(raw_activities))


message("----------------------------------------------------------------")
message(">>> STEP B3: ChEMBL metadata acquisition")
message("----------------------------------------------------------------")

fetch_metadata_by_ids <- function(resource, object_key, id_field, ids, cache_stem, refresh = FALSE, only_fields = NULL) {
  cache <- if (!refresh) read_table_cache(cache_stem) else NULL
  if (is.null(cache)) cache <- data.frame(stringsAsFactors = FALSE)
  cache <- normalize_api_frame(cache, only_fields)

  cached_ids <- if (nrow(cache) && id_field %in% names(cache)) unique(clean_chr(cache[[id_field]])) else character(0)
  missing <- setdiff(ids, cached_ids)

  if (length(missing)) {
    batches <- split_batches(missing, CHEMBL_BATCH_SIZE)
    new_rows <- list()
    for (i in seq_along(batches)) {
      fetched <- tryCatch(
        chembl_fetch_pages(
          resource = resource,
          object_key = object_key,
          query = stats::setNames(list(paste(batches[[i]], collapse = ",")), paste0(id_field, "__in")),
          only_fields = only_fields
        ),
        error = function(e) {
          warning(resource, " metadata batch ", i, " failed: ", e$message)
          data.frame(stringsAsFactors = FALSE)
        }
      )
      fetched <- normalize_api_frame(fetched, only_fields)
      if (nrow(fetched)) new_rows[[length(new_rows) + 1L]] <- fetched
      Sys.sleep(CHEMBL_PAUSE_SECONDS)
    }
    new_metadata <- if (length(new_rows)) {
      normalize_api_frame(dplyr::bind_rows(new_rows), only_fields)
    } else {
      normalize_api_frame(data.frame(stringsAsFactors = FALSE), only_fields)
    }
    cache <- normalize_api_frame(dplyr::bind_rows(cache, new_metadata), only_fields)
  }

  if (nrow(cache) && id_field %in% names(cache)) {
    cache <- as_base_df(cache %>% dplyr::filter(.data[[id_field]] %in% ids) %>% dplyr::distinct())
  }
  write_table_cache(cache, cache_stem)
  cache
}

assay_ids <- if (nrow(raw_activities) && "assay_chembl_id" %in% names(raw_activities)) {
  sort(unique(clean_chr(raw_activities$assay_chembl_id)))
} else character(0)
document_ids <- if (nrow(raw_activities) && "document_chembl_id" %in% names(raw_activities)) {
  sort(unique(clean_chr(raw_activities$document_chembl_id)))
} else character(0)
target_ids <- if (nrow(raw_activities) && "target_chembl_id" %in% names(raw_activities)) {
  sort(unique(clean_chr(raw_activities$target_chembl_id)))
} else character(0)

assay_metadata <- fetch_metadata_by_ids(
  "assay", "assays", "assay_chembl_id", assay_ids,
  "chembl_assays", CHEMBL_REFRESH_METADATA, ASSAY_API_FIELDS
)
document_metadata <- fetch_metadata_by_ids(
  "document", "documents", "document_chembl_id", document_ids,
  "chembl_documents", CHEMBL_REFRESH_METADATA, DOCUMENT_API_FIELDS
)
target_metadata <- fetch_metadata_by_ids(
  "target", "targets", "target_chembl_id", target_ids,
  "chembl_targets", CHEMBL_REFRESH_METADATA, TARGET_API_FIELDS
)

chembl_status_result <- api_get_json(paste0(CHEMBL_BASE_URL, "/status.json"), timeout_seconds = 30)
chembl_status <- if (isTRUE(chembl_status_result$ok)) chembl_status_result$json else list(
  status_error = chembl_status_result$error %||% "Unavailable"
)
jsonlite::write_json(
  chembl_status,
  out2_path("chembl_api_status.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  na = "null"
)


message("----------------------------------------------------------------")
message(">>> STEP C: ChEMBL curation and potency normalization")
message("----------------------------------------------------------------")

ensure_col <- function(df, name, default = NA) {
  if (!name %in% names(df)) df[[name]] <- rep(default, nrow(df))
  df
}

required_activity_cols <- c(
  "activity_id", "molecule_chembl_id", "standard_type", "standard_relation",
  "standard_value", "standard_units", "pchembl_value", "data_validity_comment",
  "potential_duplicate", "assay_chembl_id", "assay_description", "assay_type",
  "bao_label", "target_chembl_id", "target_pref_name", "target_organism",
  "target_type", "document_chembl_id", "src_id", "activity_comment"
)
for (nm in required_activity_cols) raw_activities <- ensure_col(raw_activities, nm, NA)

normalize_relation <- function(x) {
  x <- clean_chr(x)
  x <- stringr::str_replace_all(x, c("≤" = "<=", "≥" = ">="))
  x
}

normalize_unit <- function(x) {
  x <- clean_chr(x)
  x <- stringr::str_replace_all(x, "μ", "u")
  x <- stringr::str_replace_all(x, "µ", "u")
  x <- stringr::str_replace_all(x, "\\s+", "")
  tolower(x)
}

convert_to_nm <- function(value, unit) {
  value <- safe_numeric(value)
  unit <- normalize_unit(unit)
  multiplier <- dplyr::case_when(
    unit == "fm" ~ 1e-6,
    unit == "pm" ~ 1e-3,
    unit == "nm" ~ 1,
    unit == "um" ~ 1e3,
    unit == "mm" ~ 1e6,
    unit == "m" ~ 1e9,
    TRUE ~ NA_real_
  )
  value * multiplier
}

parse_duplicate_flag <- function(x) {
  x0 <- tolower(clean_chr(x))
  x0 %in% c("1", "true", "t", "yes", "y")
}

validity_allowed <- function(x) {
  x0 <- clean_chr(x)
  is.na(x0) | x0 %in% ALLOW_VALIDITY_COMMENTS
}

classify_cutoff <- function(relation, value_nm, cutoff_nm) {
  relation <- normalize_relation(relation)
  dplyr::case_when(
    is.na(value_nm) | value_nm <= 0 ~ "not_classifiable",
    relation == "=" & value_nm <= cutoff_nm ~ "active_at_or_below_cutoff",
    relation == "=" & value_nm > cutoff_nm ~ "above_cutoff",
    relation %in% c("<", "<=") & value_nm <= cutoff_nm ~ "active_at_or_below_cutoff",
    relation == ">" & value_nm >= cutoff_nm ~ "above_cutoff",
    relation == ">=" & value_nm > cutoff_nm ~ "above_cutoff",
    relation %in% c("<", "<=", ">", ">=") ~ "indeterminate_censored",
    TRUE ~ "not_classifiable"
  )
}


assay_keep <- intersect(
  c(
    "assay_chembl_id", "assay_type", "confidence_score", "description",
    "relationship_type", "assay_organism", "bao_format", "cell_chembl_id",
    "tissue_chembl_id", "src_id"
  ),
  names(assay_metadata)
)
assay_join <- if (length(assay_keep)) {
  assay_metadata[, assay_keep, drop = FALSE]
} else empty_df(assay_chembl_id = character())
names(assay_join)[names(assay_join) != "assay_chembl_id"] <- paste0(
  "assay_meta_", names(assay_join)[names(assay_join) != "assay_chembl_id"]
)

doc_keep <- intersect(
  c("document_chembl_id", "doi", "pubmed_id", "year", "title", "journal", "doc_type", "src_id"),
  names(document_metadata)
)
document_join <- if (length(doc_keep)) {
  document_metadata[, doc_keep, drop = FALSE]
} else empty_df(document_chembl_id = character())
names(document_join)[names(document_join) != "document_chembl_id"] <- paste0(
  "document_", names(document_join)[names(document_join) != "document_chembl_id"]
)

target_keep <- intersect(
  c("target_chembl_id", "pref_name", "target_type", "organism", "tax_id"),
  names(target_metadata)
)
target_join <- if (length(target_keep)) {
  target_metadata[, target_keep, drop = FALSE]
} else empty_df(target_chembl_id = character())
names(target_join)[names(target_join) != "target_chembl_id"] <- paste0(
  "target_meta_", names(target_join)[names(target_join) != "target_chembl_id"]
)


for (nm in c("assay_meta_assay_type", "assay_meta_confidence_score",
             "assay_meta_description", "assay_meta_relationship_type",
             "assay_meta_assay_organism", "assay_meta_bao_format",
             "assay_meta_cell_chembl_id", "assay_meta_tissue_chembl_id",
             "assay_meta_src_id")) {
  assay_join <- ensure_col(assay_join, nm, NA)
}
for (nm in c("document_doi", "document_pubmed_id", "document_year",
             "document_title", "document_journal", "document_doc_type",
             "document_src_id")) {
  document_join <- ensure_col(document_join, nm, NA)
}
for (nm in c("target_meta_pref_name", "target_meta_target_type",
             "target_meta_organism", "target_meta_tax_id")) {
  target_join <- ensure_col(target_join, nm, NA)
}

chembl_curated <- as_base_df(
  raw_activities %>%
    dplyr::mutate(
      molecule_chembl_id = clean_chr(molecule_chembl_id),
      standard_type = clean_chr(standard_type),
      standard_relation = normalize_relation(standard_relation),
      standard_units = clean_chr(standard_units),
      standard_value_numeric = safe_numeric(standard_value),
      value_nm = convert_to_nm(standard_value_numeric, standard_units),
      pchembl_value_numeric = safe_numeric(pchembl_value),
      p_activity_exact = dplyr::if_else(
        standard_relation == "=" & !is.na(value_nm) & value_nm > 0,
        -log10(value_nm * 1e-9),
        NA_real_
      ),
      p_activity_primary = dplyr::if_else(
        standard_type %in% PCHEMBL_TYPES,
        pchembl_value_numeric,
        NA_real_
      ),
      is_potential_duplicate = parse_duplicate_flag(potential_duplicate),
      validity_is_allowed = validity_allowed(data_validity_comment),
      unit_is_convertible = !is.na(value_nm),
      relation_is_exact = standard_relation == "=",
      type_is_selected = standard_type %in% VALID_TYPES,
      activity_class_10uM = classify_cutoff(standard_relation, value_nm, ACTIVITY_CUTOFF_NM),
      definite_active_10uM = activity_class_10uM == "active_at_or_below_cutoff",
      definite_above_10uM = activity_class_10uM == "above_cutoff"
    ) %>%
    dplyr::left_join(chembl_map %>% dplyr::select(inchikey, chembl_id), by = c("molecule_chembl_id" = "chembl_id")) %>%
    dplyr::left_join(assay_join, by = "assay_chembl_id") %>%
    dplyr::left_join(document_join, by = "document_chembl_id") %>%
    dplyr::left_join(target_join, by = "target_chembl_id") %>%
    dplyr::mutate(
      target_pref_name_final = dplyr::coalesce(
        clean_chr(target_pref_name),
        clean_chr(target_meta_pref_name)
      ),
      target_type_final = dplyr::coalesce(
        clean_chr(target_type),
        clean_chr(target_meta_target_type)
      ),
      target_organism_final = dplyr::coalesce(
        clean_chr(target_organism),
        clean_chr(target_meta_organism)
      ),
      assay_type_final = dplyr::coalesce(
        clean_chr(assay_type),
        clean_chr(assay_meta_assay_type)
      ),
      assay_description_final = dplyr::coalesce(
        clean_chr(assay_description),
        clean_chr(assay_meta_description)
      ),
      reference_id = dplyr::case_when(
        !is.na(document_doi) ~ document_doi,
        !is.na(document_pubmed_id) ~ paste0("PMID:", document_pubmed_id),
        TRUE ~ document_chembl_id
      ),
      quality_record = (
        type_is_selected &
          unit_is_convertible &
          !is.na(standard_value_numeric) & standard_value_numeric > 0 &
          !is_potential_duplicate & validity_is_allowed
      ),
      eligible_extended_quantitative = quality_record & relation_is_exact & !is.na(p_activity_exact),
      eligible_primary_quantitative = (
        quality_record & relation_is_exact &
          standard_type %in% PCHEMBL_TYPES & !is.na(p_activity_primary)
      ),
      eligible_binary_10uM = quality_record & activity_class_10uM %in% c(
        "active_at_or_below_cutoff", "above_cutoff"
      ),
      exclusion_reason_primary = dplyr::case_when(
        !type_is_selected ~ "standard_type_not_selected",
        is.na(standard_value_numeric) | standard_value_numeric <= 0 ~ "invalid_or_missing_value",
        !unit_is_convertible ~ "unsupported_unit",
        is_potential_duplicate ~ "potential_duplicate",
        !validity_is_allowed ~ "data_validity_flag",
        !relation_is_exact ~ "non_exact_relation",
        !standard_type %in% PCHEMBL_TYPES ~ "not_pchembl_comparable_type",
        is.na(p_activity_primary) ~ "missing_pchembl_value",
        TRUE ~ "eligible"
      )
    )
)


chembl_curated <- as_base_df(chembl_curated %>% dplyr::distinct())
write_primary_table(chembl_curated, paste0(base_tag, "_chembl_activities_curated"))


message("----------------------------------------------------------------")
message(">>> STEP D: functional target-domain categorization")
message("----------------------------------------------------------------")

rx <- function(pattern) stringr::regex(pattern, ignore_case = TRUE)

categorize_target_domain <- function(df) {
  df <- as_base_df(df)

  target_text <- paste(
    clean_chr(df$target_pref_name_final),
    clean_chr(df$target_type_final),
    clean_chr(df$target_organism_final)
  )
  assay_text <- clean_chr(df$assay_description_final)

  molecular_target_types <- c(
    "SINGLE PROTEIN", "PROTEIN COMPLEX", "PROTEIN FAMILY",
    "PROTEIN COMPLEX GROUP", "PROTEIN-PROTEIN INTERACTION",
    "CHIMERIC PROTEIN", "SELECTIVITY GROUP"
  )
  is_molecular_target <- toupper(clean_chr(df$target_type_final)) %in% molecular_target_types


  classification_text <- ifelse(
    is_molecular_target,
    target_text,
    paste(target_text, assay_text)
  )

  detail <- dplyr::case_when(
    stringr::str_detect(classification_text, rx("A549|NCI-H|Calu|HOP|lung cancer")) ~ "Cells: Lung Cancer",
    stringr::str_detect(classification_text, rx("MCF-?7|MDA-MB|T47D|SK-BR|breast cancer")) ~ "Cells: Breast Cancer",
    stringr::str_detect(classification_text, rx("HCT|HT-29|SW480|SW620|Caco|colorectal|colon cancer")) ~ "Cells: Colorectal Cancer",
    stringr::str_detect(classification_text, rx("K562|HL-60|Jurkat|CEM T-lymphocyte|Molt 4|MOLT-?4|leukemia|lymphoma|myeloma")) ~ "Cells: Hematologic Cancer",
    stringr::str_detect(classification_text, rx("HeLa|SiHa|cervical cancer")) ~ "Cells: Cervical Cancer",
    stringr::str_detect(classification_text, rx("PC-3|DU-145|LNCaP|prostate cancer")) ~ "Cells: Prostate Cancer",
    stringr::str_detect(classification_text, rx("melanoma|SK-MEL|A375|B16")) ~ "Cells: Melanoma",
    stringr::str_detect(classification_text, rx("HepG2|Hep3B|liver cancer|hepatoma")) ~ "Cells: Liver Cancer",
    stringr::str_detect(classification_text, rx("OVCAR|SK-OV|A2780|ovarian cancer")) ~ "Cells: Ovarian Cancer",
    stringr::str_detect(classification_text, rx("PANC-1|MIA PaCa|pancreatic cancer|gastric cancer|AGS")) ~ "Cells: GI/Pancreatic Cancer",
    stringr::str_detect(classification_text, rx("cell line|CELL-LINE|cytotoxicity|antiproliferative|tumou?r cell|cancer cell")) ~ "Cells: Phenotypic/Other",

    
    stringr::str_detect(classification_text, rx("virus|HIV|influenza|dengue|SARS|hepatitis|chikungunya|coronavirus")) ~ "Pathogen: Virus",
    stringr::str_detect(classification_text, rx(
      "bacter|Staphylococcus|Escherichia|Mycobacterium|Pseudomonas|Streptococcus|Bacillus|Enterococcus|Klebsiella|Salmonella|Enterobacter|Proteus|Helicobacter|Clostridium|Listeria|Neisseria|Shigella|Vibrio|Acinetobacter|telomere resolvase|\\bResT\\b"
    )) ~ "Pathogen: Bacteria",
    stringr::str_detect(classification_text, rx(
      "fung|Candida|Aspergillus|Cryptococcus|Trichophyton|Microsporum|Pyricularia|Cladosporium|Alternaria|Botrytis|Fusarium|Saccharomyces"
    )) ~ "Pathogen: Fungus",
    stringr::str_detect(classification_text, rx("parasite|Plasmodium|Trypanosoma|Leishmania|Entamoeba|Taenia|Cruzipain")) ~ "Pathogen: Parasite",

    
    stringr::str_detect(classification_text, rx("cyclooxygenase|COX-|prostaglandin|thromboxane")) ~ "Enzymes: Inflammatory COX/PG",
    stringr::str_detect(classification_text, rx("lipoxygenase|ALOX|LOX-")) ~ "Enzymes: Lipoxygenase",
    stringr::str_detect(classification_text, rx("nitric oxide synthase|\\bNOS\\b")) ~ "Enzymes: Nitric Oxide Synthase",
    stringr::str_detect(classification_text, rx("complement pathway|complement system|NLRP3|NACHT, LRR and PYD|inflammasome")) ~ "Proteins: Inflammasome/Complement",
    stringr::str_detect(classification_text, rx("neutrophil elastase|elastase|stromelysin|matrilysin|metalloproteinase|ADAMTS|matrix metallo|\\bMMP[- ]|collagenase|gelatinase")) ~ "Enzymes: Matrix Metalloproteases",
    stringr::str_detect(classification_text, rx("protease|peptidase|cathepsin|calpain|thrombin|trypsin|caspase|renin|proteasome|papain|kallikrein|chymase|angiotensin-converting enzyme|tissue factor")) ~ "Enzymes: Proteases",
    stringr::str_detect(classification_text, rx("macrophage migration inhibitory factor|high mobility group protein B1|lysozyme C|antiinflamm|ear edema|neutrophil|LPS/IFN|superoxide anion")) ~ "Proteins: Inflammatory Response",

    
    stringr::str_detect(classification_text, rx("cholinesterase|AChE|butyrylcholinesterase")) ~ "Enzymes: Cholinesterase",
    stringr::str_detect(classification_text, rx("monoamine oxidase|\\bMAO[- ]")) ~ "Enzymes: Monoamine Oxidase",
    stringr::str_detect(classification_text, rx("amyloid|BACE|tau|secretase|alpha-synuclein|Alzheimer|transthyretin")) ~ "Proteins: Neurodegeneration",
    stringr::str_detect(classification_text, rx("serotonin|5-HT|hydroxytryptamine receptor|dopamine|adrenergic|histamine receptor|muscarinic|nicotinic|acetylcholine receptor|GABA|opioid|cannabinoid|adenosine( A2)? receptor|glutamate receptor|taste receptor|melatonin receptor|vasopressin V1a receptor|corticotropin-releasing factor")) ~ "Receptors: Neurotransmission",

    
    stringr::str_detect(classification_text, rx("kinase|CDK|Aurora|MAPK|PI3K|Akt|mTOR|GSK|PKC|EGFR|VEGF|FLT3|Src|Abl|JAK")) ~ "Enzymes: Kinases",
    stringr::str_detect(classification_text, rx("phosphatase|synaptojanin")) ~ "Enzymes: Phosphatases",
    stringr::str_detect(classification_text, rx("NF-kappa|NF-kB|tumou?r necrosis factor|\\bTNF\\b|hedgehog|Wnt|beta-catenin|signal transducer and activator of transcription|\\bSTAT[0-9]")) ~ "Proteins: NF-kB/TNF Signaling",
    stringr::str_detect(classification_text, rx("interleukin|cytokine|sclerostin")) ~ "Proteins: Cytokines",
    stringr::str_detect(classification_text, rx("apoptosis|Bcl-2|BCL2|BID|BAD|Bcl-xL|Mcl-1|p53|MDM2|survivin|BAX|programmed cell death protein 1|PD-?1|PD-?L1")) ~ "Proteins: Apoptosis",
    stringr::str_detect(classification_text, rx("HSP90|HSP60|HSP10|heat shock|chaperone")) ~ "Proteins: Heat Shock",
    stringr::str_detect(classification_text, rx("growth factor|angiopoietin|EGF|HGF|IGF|hypoxia-inducible|\\bHIF[- ]")) ~ "Proteins: Growth Factors",
    stringr::str_detect(classification_text, rx("ubiquitin|autophagy|LC3|microtubule|tubulin|cytoskeleton|actin|calmodulin")) ~ "Proteins: Proteostasis/Cytoskeleton",

    
    stringr::str_detect(classification_text, rx("HDAC|sirtuin|histone|methyltransferase|demethylase|epigenetic|polycomb|\\bEED\\b|\\bFTO\\b")) ~ "Epigenetics: Chromatin",
    stringr::str_detect(classification_text, rx("estrogen receptor|androgen receptor|PPAR|peroxisome proliferator-activated receptor|glucocorticoid receptor|progesterone receptor|retinoic|vitamin D receptor|RXR|nuclear receptor|oxysterols receptor|steroid hormone receptor ERR")) ~ "Receptors: Nuclear",
    stringr::str_detect(classification_text, rx("DNA|RNA polymerase|helicase|telomerase|transcription factor|zinc finger|X-box-binding|translation initiation factor|polyadenylate-binding|RNA-binding protein|ribosomal subunit|elongation factor|RuvB-like|CREB-binding|MUS81|poly \\[ADP-ribose\\] polymerase|ELAV-like")) ~ "Proteins: DNA/RNA Regulation",

    
    stringr::str_detect(classification_text, rx("glucosidase|glucuronidase|amylase|maltase|aldose reductase|diabetes")) ~ "Enzymes: Carbohydrate Metabolism",
    stringr::str_detect(classification_text, rx("lipase|fatty acid synthase|HMG-CoA|squalene|lipid metabolism|epoxide hydrolase|sterol O-acyltransferase|diacylglycerol O-acyltransferase|perilipin|ABHD5|hepatic steatosis|serum LDL")) ~ "Enzymes: Lipid Metabolism",
    stringr::str_detect(classification_text, rx("aromatase|hydroxysteroid|steroid metabolism|steroid 17-alpha-hydroxylase")) ~ "Enzymes: Steroid Metabolism",
    stringr::str_detect(classification_text, rx("cytochrome P450|\\bCYP[0-9]")) ~ "Enzymes: CYP450",
    stringr::str_detect(classification_text, rx("UGT|UDP-glucuronosyltransferase|glucuronosyltransferase|sulfotransferase|glutathione transferase|phase II")) ~ "Enzymes: Phase II Metabolism",
    stringr::str_detect(classification_text, rx("carbonic anhydrase")) ~ "Enzymes: Carbonic Anhydrase",
    stringr::str_detect(classification_text, rx("phosphodiesterase|\\bPDE[0-9]")) ~ "Enzymes: Phosphodiesterase",
    stringr::str_detect(classification_text, rx("tyrosinase")) ~ "Enzymes: Tyrosinase",
    stringr::str_detect(classification_text, rx("sialidase|neuraminidase")) ~ "Enzymes: Sialidase",
    stringr::str_detect(classification_text, rx("oxidoreductase|dehydrogenase|reductase|oxidase|oxygenase|dioxygenase|decarboxylase|isomerase|esterase|lactoylglutathione lyase|phosphoglycerate mutase|alpha-mannosidase|glycogen phosphorylase|ATP-diphosphohydrolase|thiosulfate sulfurtransferase|liver microsome|drug metabolism")) ~ "Enzymes: General Metabolism",

    
    stringr::str_detect(classification_text, rx(
      "ion channel|channel|transporter|translocase|pump|ATPase|SERCA|CFTR|cystic fibrosis transmembrane conductance regulator|SGLT|solute carrier|\\bSLC[0-9]|sodium channel|potassium channel|calcium channel|P-glycoprotein|ABC transporter|\\bABCB[0-9]|\\bABCG[0-9]|multidrug resistance-associated protein|scavenger receptor class B|sodium/hydrogen exchanger|hERG"
    )) ~ "Transporters and Ion Channels",

    
    stringr::str_detect(classification_text, rx("radical scavenging|DPPH|ABTS|FRAP|antioxidant activity")) ~ "General Assays and Safety",
    stringr::str_detect(classification_text, rx("ADMET|toxicity|safety|non-protein|unchecked|no relevant target|physicochemical|solubility|permeability|vasorelaxant|aortic ring|aorta|phytotoxic|root growth|chloroplast|ATP synthesis|sea urchin|antimitotic|lipid peroxidation|TBARS")) ~ "General Assays and Safety",
    TRUE ~ "Other Targets"
  )

  macro <- dplyr::case_when(
    stringr::str_starts(detail, "Cells:") ~ "Phenotypic: Cancer and Cells",
    stringr::str_starts(detail, "Pathogen:") ~ "Infectious Diseases",
    stringr::str_detect(detail, rx("Kinases|Phosphatases|Apoptosis|Heat Shock|Growth Factors|NF-kB|Cytokines|Proteostasis|Cytoskeleton")) ~ "Cell Signaling and Survival",
    stringr::str_detect(detail, rx("Inflammatory|Inflammation|Lipoxygenase|Nitric Oxide|Proteases|Metalloproteases|Inflammasome|Complement")) ~ "Inflammation and Proteolysis",
    stringr::str_detect(detail, rx("Cholinesterase|Monoamine|Neurodegeneration|Neurotransmission")) ~ "Neuroscience Targets",
    stringr::str_detect(detail, rx("Chromatin|Nuclear|DNA/RNA")) ~ "Gene Regulation",
    stringr::str_detect(detail, rx("Carbohydrate|Lipid|Steroid|CYP450|Phase II|Carbonic|Phosphodiesterase|Tyrosinase|Sialidase|General Metabolism")) ~ "Metabolic Enzymes",
    detail == "Transporters and Ion Channels" ~ "Transporters and Ion Channels",
    detail == "General Assays and Safety" ~ "General Assays and Safety",
    TRUE ~ "Miscellaneous"
  )

  df$target_category_l3 <- detail
  df$target_domain_macro <- macro
  df$target_category_rule <- detail
  df$target_classification_basis <- ifelse(
    is_molecular_target,
    "target_identity",
    "target_and_assay_context"
  )
  df
}

chembl_curated <- categorize_target_domain(chembl_curated)
write_primary_table(chembl_curated, paste0(base_tag, "_chembl_activities_curated"))


message("----------------------------------------------------------------")
message(">>> STEP E: compound-target and compound-domain summaries")
message("----------------------------------------------------------------")

quantile_safe <- function(x, p) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(stats::quantile(x, p, na.rm = TRUE, names = FALSE, type = 7))
}

primary_records <- as_base_df(
  chembl_curated %>%
    dplyr::filter(eligible_primary_quantitative, !is.na(inchikey), !is.na(target_chembl_id))
)
extended_records <- as_base_df(
  chembl_curated %>%
    dplyr::filter(eligible_extended_quantitative, !is.na(inchikey), !is.na(target_chembl_id))
)
binary_records <- as_base_df(
  chembl_curated %>%
    dplyr::filter(quality_record, !is.na(inchikey))
)

compound_target_primary <- as_base_df(
  primary_records %>%
    dplyr::group_by(
      inchikey, target_chembl_id, target_pref_name_final, target_type_final,
      target_organism_final, target_category_l3, target_domain_macro,
      standard_type
    ) %>%
    dplyr::summarise(
      n_records = dplyr::n(),
      n_assays = dplyr::n_distinct(assay_chembl_id, na.rm = TRUE),
      n_documents = dplyr::n_distinct(document_chembl_id, na.rm = TRUE),
      median_pchembl = stats::median(p_activity_primary, na.rm = TRUE),
      q25_pchembl = quantile_safe(p_activity_primary, 0.25),
      q75_pchembl = quantile_safe(p_activity_primary, 0.75),
      min_pchembl = min(p_activity_primary, na.rm = TRUE),
      max_pchembl = max(p_activity_primary, na.rm = TRUE),
      references = collapse_unique(reference_id, sep = " | "),
      .groups = "drop"
    )
)

compound_target_extended <- as_base_df(
  extended_records %>%
    dplyr::group_by(
      inchikey, target_chembl_id, target_pref_name_final, target_type_final,
      target_organism_final, target_category_l3, target_domain_macro,
      standard_type
    ) %>%
    dplyr::summarise(
      n_records = dplyr::n(),
      n_assays = dplyr::n_distinct(assay_chembl_id, na.rm = TRUE),
      n_documents = dplyr::n_distinct(document_chembl_id, na.rm = TRUE),
      median_p_activity = stats::median(p_activity_exact, na.rm = TRUE),
      q25_p_activity = quantile_safe(p_activity_exact, 0.25),
      q75_p_activity = quantile_safe(p_activity_exact, 0.75),
      min_p_activity = min(p_activity_exact, na.rm = TRUE),
      max_p_activity = max(p_activity_exact, na.rm = TRUE),
      references = collapse_unique(reference_id, sep = " | "),
      .groups = "drop"
    )
)


compound_domain_primary <- as_base_df(
  compound_target_primary %>%
    dplyr::group_by(inchikey, target_domain_macro, standard_type) %>%
    dplyr::summarise(
      n_targets = dplyr::n_distinct(target_chembl_id),
      n_target_summaries = dplyr::n(),
      n_records = sum(n_records, na.rm = TRUE),
      n_assays = sum(n_assays, na.rm = TRUE),
      median_target_pchembl = stats::median(median_pchembl, na.rm = TRUE),
      q25_target_pchembl = quantile_safe(median_pchembl, 0.25),
      q75_target_pchembl = quantile_safe(median_pchembl, 0.75),
      .groups = "drop"
    )
)

compound_domain_extended <- as_base_df(
  compound_target_extended %>%
    dplyr::group_by(inchikey, target_domain_macro, standard_type) %>%
    dplyr::summarise(
      n_targets = dplyr::n_distinct(target_chembl_id),
      n_target_summaries = dplyr::n(),
      n_records = sum(n_records, na.rm = TRUE),
      n_assays = sum(n_assays, na.rm = TRUE),
      median_target_p_activity = stats::median(median_p_activity, na.rm = TRUE),
      q25_target_p_activity = quantile_safe(median_p_activity, 0.25),
      q75_target_p_activity = quantile_safe(median_p_activity, 0.75),
      .groups = "drop"
    )
)

compound_binary_summary <- as_base_df(
  binary_records %>%
    dplyr::group_by(inchikey, target_domain_macro) %>%
    dplyr::summarise(
      n_quality_records = dplyr::n(),
      n_definite_active = sum(definite_active_10uM, na.rm = TRUE),
      n_definite_above = sum(definite_above_10uM, na.rm = TRUE),
      n_indeterminate = sum(activity_class_10uM == "indeterminate_censored", na.rm = TRUE),
      any_definite_active = any(definite_active_10uM, na.rm = TRUE),
      any_definite_above = any(definite_above_10uM, na.rm = TRUE),
      .groups = "drop"
    )
)

compound_bioactivity_summary <- as_base_df(
  compound_taxonomy_index %>%
    dplyr::select(inchikey, class_np, family, genus, species, n_occurrence_rows, n_families, n_species) %>%
    dplyr::left_join(
      chembl_map %>%
        dplyr::group_by(inchikey) %>%
        dplyr::summarise(
          chembl_ids = collapse_unique(chembl_id),
          n_chembl_ids = dplyr::n_distinct(chembl_id),
          .groups = "drop"
        ),
      by = "inchikey"
    ) %>%
    dplyr::left_join(
      chembl_curated %>%
        dplyr::group_by(inchikey) %>%
        dplyr::summarise(
          n_records_raw = dplyr::n(),
          n_quality_records = sum(quality_record, na.rm = TRUE),
          n_primary_records = sum(eligible_primary_quantitative, na.rm = TRUE),
          n_extended_records = sum(eligible_extended_quantitative, na.rm = TRUE),
          n_binary_records = sum(eligible_binary_10uM, na.rm = TRUE),
          n_targets = dplyr::n_distinct(target_chembl_id[quality_record], na.rm = TRUE),
          n_domains = dplyr::n_distinct(target_domain_macro[quality_record], na.rm = TRUE),
          any_definite_active_10uM = any(definite_active_10uM & quality_record, na.rm = TRUE),
          any_definite_above_10uM = any(definite_above_10uM & quality_record, na.rm = TRUE),
          median_pchembl = if (any(eligible_primary_quantitative)) stats::median(p_activity_primary[eligible_primary_quantitative], na.rm = TRUE) else NA_real_,
          .groups = "drop"
        ),
      by = "inchikey"
    ) %>%
    dplyr::mutate(
      mapped_to_chembl = !is.na(n_chembl_ids),
      n_chembl_ids = tidyr::replace_na(n_chembl_ids, 0L),
      n_records_raw = tidyr::replace_na(n_records_raw, 0L),
      n_quality_records = tidyr::replace_na(n_quality_records, 0L),
      n_primary_records = tidyr::replace_na(n_primary_records, 0L),
      n_extended_records = tidyr::replace_na(n_extended_records, 0L),
      n_binary_records = tidyr::replace_na(n_binary_records, 0L),
      n_targets = tidyr::replace_na(n_targets, 0L),
      n_domains = tidyr::replace_na(n_domains, 0L),
      any_definite_active_10uM = tidyr::replace_na(any_definite_active_10uM, FALSE),
      any_definite_above_10uM = tidyr::replace_na(any_definite_above_10uM, FALSE)
    )
)

write_primary_table(compound_target_primary, paste0(base_tag, "_compound_target_primary"))
write_primary_table(compound_target_extended, paste0(base_tag, "_compound_target_extended"))
write_primary_table(compound_domain_primary, paste0(base_tag, "_compound_domain_primary"))
write_primary_table(compound_domain_extended, paste0(base_tag, "_compound_domain_extended"))
write_primary_table(compound_binary_summary, paste0(base_tag, "_compound_domain_binary"))
write_primary_table(compound_bioactivity_summary, paste0(base_tag, "_compound_bioactivity_summary"))


message("----------------------------------------------------------------")
message(">>> STEP F: optional global LOTUS family context")
message("----------------------------------------------------------------")

extract_families_recursive <- function(x) {
  if (is.null(x)) return(character(0))
  flat <- tryCatch(unlist(x, recursive = TRUE, use.names = TRUE), error = function(e) NULL)
  if (is.null(flat) || !length(flat)) return(character(0))
  nm <- names(flat)
  if (is.null(nm)) return(character(0))
  idx <- grepl("(^|\\.)family$", nm, ignore.case = TRUE)
  fam <- clean_chr(flat[idx])
  sort(unique(fam[!is.na(fam)]))
}

global_context <- empty_df(
  inchikey = inchikeys_to_map,
  global_family_count = NA_integer_,
  global_family_list = NA_character_,
  global_context_status = "not_run"
)

if (RUN_GLOBAL_CONTEXT) {
  mongo_url <- cfg$mongo_url %||% "mongodb://127.0.0.1:27017"
  mongo_db_name <- cfg$db_name %||% "lotusdb"

  global_context <- tryCatch({
    lotus_global <- mongolite::mongo(
      collection = GLOBAL_COLLECTION,
      db = mongo_db_name,
      url = mongo_url
    )

    chunks <- split_batches(inchikeys_to_map, as.integer(cfg$mongo_batch_size %||% 2000L))
    rows <- list()
    pb <- progress::progress_bar$new(
      format = "   Global context [:bar] :percent | Batch :current/:total | ETA: :eta",
      total = length(chunks), width = 70
    )

    for (i in seq_along(chunks)) {
      pb$tick()
      q <- jsonlite::toJSON(
        list(inchikey = list(`$in` = unname(chunks[[i]]))),
        auto_unbox = TRUE
      )
      batch <- as_base_df(
        lotus_global$find(
          query = q,
          fields = '{"inchikey":1,"taxonomyReferenceObjects":1,"_id":0}'
        )
      )
      if (!nrow(batch)) next

      for (j in seq_len(nrow(batch))) {
        tax_obj <- if (is.list(batch$taxonomyReferenceObjects)) {
          batch$taxonomyReferenceObjects[[j]]
        } else {
          batch$taxonomyReferenceObjects[j, , drop = FALSE]
        }
        fam <- extract_families_recursive(tax_obj)
        if (!length(fam)) fam <- NA_character_
        rows[[length(rows) + 1L]] <- empty_df(
          inchikey = rep(clean_chr(batch$inchikey[j]), length(fam)),
          family = fam
        )
      }
    }

    long <- as_base_df(dplyr::bind_rows(rows))
    if (!nrow(long)) {
      empty_df(
        inchikey = inchikeys_to_map,
        global_family_count = 0L,
        global_family_list = NA_character_,
        global_context_status = "not_found"
      )
    } else {
      summary <- as_base_df(
        long %>%
          dplyr::filter(!is.na(inchikey), !is.na(family)) %>%
          dplyr::distinct(inchikey, family) %>%
          dplyr::group_by(inchikey) %>%
          dplyr::summarise(
            global_family_count = dplyr::n_distinct(family),
            global_family_list = collapse_unique(family),
            .groups = "drop"
          )
      )
      empty_df(inchikey = inchikeys_to_map) %>%
        dplyr::left_join(summary, by = "inchikey") %>%
        dplyr::mutate(
          global_family_count = tidyr::replace_na(global_family_count, 0L),
          global_context_status = dplyr::if_else(global_family_count > 0, "found", "not_found")
        ) %>%
        as_base_df()
    }
  }, error = function(e) {
    warning("Global LOTUS context was not completed: ", e$message)
    empty_df(
      inchikey = inchikeys_to_map,
      global_family_count = NA_integer_,
      global_family_list = NA_character_,
      global_context_status = paste0("error: ", e$message)
    )
  })
}

write_primary_table(global_context, paste0(base_tag, "_BIO_C_Global_Context"))
writexl::write_xlsx(
  list(Global_Context = as_base_df(global_context)),
  out2_path(paste0(base_tag, "_BIO_C_Global_Context.xlsx"))
)


message("----------------------------------------------------------------")
message(">>> STEP G: audit tables and compatibility workbook")
message("----------------------------------------------------------------")

filter_flow <- empty_df(
  stage = c(
    "01_raw_api_records",
    "02_selected_standard_type",
    "03_numeric_positive_value",
    "04_convertible_concentration_unit",
    "05_not_potential_duplicate",
    "06_validity_allowed",
    "07_quality_records",
    "08_extended_exact_quantitative",
    "09_primary_pchembl_quantitative",
    "10_binary_classifiable_at_10uM"
  ),
  n_records = c(
    nrow(chembl_curated),
    sum(chembl_curated$type_is_selected, na.rm = TRUE),
    sum(!is.na(chembl_curated$standard_value_numeric) & chembl_curated$standard_value_numeric > 0, na.rm = TRUE),
    sum(chembl_curated$unit_is_convertible, na.rm = TRUE),
    sum(!chembl_curated$is_potential_duplicate, na.rm = TRUE),
    sum(chembl_curated$validity_is_allowed, na.rm = TRUE),
    sum(chembl_curated$quality_record, na.rm = TRUE),
    sum(chembl_curated$eligible_extended_quantitative, na.rm = TRUE),
    sum(chembl_curated$eligible_primary_quantitative, na.rm = TRUE),
    sum(chembl_curated$eligible_binary_10uM, na.rm = TRUE)
  )
)

unit_summary <- as_base_df(
  chembl_curated %>%
    dplyr::count(standard_units, unit_is_convertible, name = "n_records", sort = TRUE)
)
relation_summary <- as_base_df(
  chembl_curated %>%
    dplyr::count(standard_relation, activity_class_10uM, name = "n_records", sort = TRUE)
)
type_summary <- as_base_df(
  chembl_curated %>%
    dplyr::count(
      standard_type,
      eligible_primary_quantitative,
      eligible_extended_quantitative,
      name = "n_records",
      sort = TRUE
    )
)
validity_summary <- as_base_df(
  chembl_curated %>%
    dplyr::count(data_validity_comment, is_potential_duplicate, name = "n_records", sort = TRUE)
)
target_category_summary <- as_base_df(
  chembl_curated %>%
    dplyr::count(target_domain_macro, target_category_l3, name = "n_records", sort = TRUE)
)
target_category_primary_summary <- as_base_df(
  chembl_curated %>%
    dplyr::filter(eligible_primary_quantitative) %>%
    dplyr::count(target_domain_macro, target_category_l3, name = "n_records", sort = TRUE)
)
unclassified_targets <- as_base_df(
  chembl_curated %>%
    dplyr::filter(target_domain_macro == "Miscellaneous") %>%
    dplyr::count(
      target_chembl_id, target_pref_name_final, target_type_final,
      assay_description_final, target_classification_basis,
      name = "n_records", sort = TRUE
    )
)
unclassified_targets_primary <- as_base_df(
  chembl_curated %>%
    dplyr::filter(
      eligible_primary_quantitative,
      target_domain_macro == "Miscellaneous"
    ) %>%
    dplyr::count(
      target_chembl_id, target_pref_name_final, target_type_final,
      assay_description_final, target_classification_basis,
      name = "n_records", sort = TRUE
    )
)

readr::write_csv(filter_flow, out2_path("chembl_filter_flow.csv"), na = "")
readr::write_csv(unit_summary, out2_path("chembl_unit_summary.csv"), na = "")
readr::write_csv(relation_summary, out2_path("chembl_relation_summary.csv"), na = "")
readr::write_csv(type_summary, out2_path("chembl_standard_type_summary.csv"), na = "")
readr::write_csv(validity_summary, out2_path("chembl_validity_summary.csv"), na = "")
readr::write_csv(target_category_summary, out2_path("chembl_target_category_summary.csv"), na = "")
readr::write_csv(
  target_category_primary_summary,
  out2_path("chembl_target_category_primary_summary.csv"),
  na = ""
)
readr::write_csv(unclassified_targets, out2_path("chembl_unclassified_targets.csv"), na = "")
readr::write_csv(
  unclassified_targets_primary,
  out2_path("chembl_unclassified_targets_primary.csv"),
  na = ""
)


master_data <- as_base_df(
  primary_records %>%
    dplyr::left_join(
      compound_taxonomy_index %>%
        dplyr::select(inchikey, class_np, family, genus, species, n_occurrence_rows, n_families, n_species),
      by = "inchikey"
    ) %>%
    dplyr::left_join(global_context, by = "inchikey") %>%
    dplyr::rename(
      Target_Category_L3 = target_category_l3,
      Target_Category_Macro = target_domain_macro,
      Plot_Class = class_np,
      pActivity = p_activity_primary
    )
)

legacy_workbook <- out2_path(paste0(base_tag, "_Lotus_Final_Database_v9_Fixed.xlsx"))
writexl::write_xlsx(
  list(
    MASTER_DATA = master_data,
    DOMAIN_PRIMARY = compound_domain_primary,
    DOMAIN_EXTENDED = compound_domain_extended,
    TARGET_PRIMARY = compound_target_primary,
    TARGET_EXTENDED = compound_target_extended,
    COMPOUND_SUMMARY = compound_bioactivity_summary,
    BINARY_SUMMARY = compound_binary_summary,
    MAPPING = chembl_mapping_full,
    FILTER_FLOW = filter_flow,
    UNIT_SUMMARY = unit_summary,
    TYPE_SUMMARY = type_summary,
    TARGET_DOMAINS = target_category_summary
  ),
  path = legacy_workbook
)


bio_summary_workbook <- out2_path(paste0(base_tag, "_BIO_A_Summary.xlsx"))
writexl::write_xlsx(
  list(
    Compound_Summary = compound_bioactivity_summary,
    Domain_Primary = compound_domain_primary,
    Domain_Extended = compound_domain_extended,
    Target_Primary = compound_target_primary,
    Mapping = chembl_mapping_full,
    Filter_Flow = filter_flow,
    Units = unit_summary,
    Relations = relation_summary,
    Types = type_summary,
    Validity = validity_summary,
    Target_Domains = target_category_summary
  ),
  path = bio_summary_workbook
)

run_metadata <- list(
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  base_tag = base_tag,
  out_dir = normalizePath(OUT_DIR, winslash = "/", mustWork = FALSE),
  activity_cutoff_nm = ACTIVITY_CUTOFF_NM,
  selected_standard_types = VALID_TYPES,
  pchembl_comparable_types = PCHEMBL_TYPES,
  allowed_validity_comments = ALLOW_VALIDITY_COMMENTS,
  target_compounds = length(inchikeys_to_map),
  mapped_compounds = dplyr::n_distinct(chembl_map$inchikey),
  raw_activity_records = nrow(raw_activities),
  quality_activity_records = sum(chembl_curated$quality_record, na.rm = TRUE),
  primary_quantitative_records = nrow(primary_records),
  extended_quantitative_records = nrow(extended_records),
  global_context_collection = GLOBAL_COLLECTION,
  global_context_requested = RUN_GLOBAL_CONTEXT,
  chembl_api_status = chembl_status,
  software = list(
    R = R.version.string,
    dplyr = as.character(utils::packageVersion("dplyr")),
    httr = as.character(utils::packageVersion("httr")),
    jsonlite = as.character(utils::packageVersion("jsonlite"))
  )
)
jsonlite::write_json(
  run_metadata,
  out2_path("part2_run_metadata.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  na = "null"
)


assign("flav_for_rdkit", flav_occurrences, envir = .GlobalEnv)
assign("compound_taxonomy_index", compound_taxonomy_index, envir = .GlobalEnv)
assign("compound_family_map", compound_family_map, envir = .GlobalEnv)
assign("chembl_mapping_full", chembl_mapping_full, envir = .GlobalEnv)
assign("chembl_raw", raw_activities, envir = .GlobalEnv)
assign("chembl_curated", chembl_curated, envir = .GlobalEnv)
assign("chembl_master", master_data, envir = .GlobalEnv)
assign("chembl_target_primary", compound_target_primary, envir = .GlobalEnv)
assign("chembl_target_extended", compound_target_extended, envir = .GlobalEnv)
assign("chembl_domain_primary", compound_domain_primary, envir = .GlobalEnv)
assign("chembl_domain_extended", compound_domain_extended, envir = .GlobalEnv)
assign("chembl_compound_summary", compound_bioactivity_summary, envir = .GlobalEnv)

message("----------------------------------------------------------------")
message(">>> PART II COMPLETED SUCCESSFULLY")
message("RDKit input: ", rdkit_input_path)
message("Mapped compounds: ", dplyr::n_distinct(chembl_map$inchikey), " / ", length(inchikeys_to_map))
message("Raw ChEMBL records: ", nrow(raw_activities))
message("Primary pChEMBL records: ", nrow(primary_records))
message("Extended exact records: ", nrow(extended_records))
message("Compatibility workbook: ", legacy_workbook)
message("OUT_DIR: ", normalizePath(OUT_DIR, winslash = "/", mustWork = FALSE))
message("----------------------------------------------------------------")
