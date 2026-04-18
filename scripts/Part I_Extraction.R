# Part I: LOTUS extraction, taxonomic normalization, and dataset export.

suppressPackageStartupMessages({
  library(here)
  library(mongolite)
  library(jsonlite)
  library(dplyr)
  library(progress)
  library(stringr)
  library(writexl)
  library(readr)
  library(stringi)
})

options(OutDec = ".", scipen = 999)

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}

as_base_df <- function(x) as.data.frame(x, stringsAsFactors = FALSE)

fix_ref_id <- function(x) {
  x <- as.character(x)
  x <- stringr::str_squish(x)
  x[!nzchar(x)] <- NA_character_
  gsub("\\$x\\$x\\$", ".", x, perl = TRUE)
}

norm_ascii <- function(x) stringi::stri_trans_general(x, "Latin-ASCII")
tidy_space <- function(x) trimws(gsub("\\s+", " ", x))

title_case_1 <- function(x) {
  x <- as.character(x)
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

canon_genus <- function(x) {
  x <- tidy_space(norm_ascii(x))
  x <- sub("\\s+.*$", "", x)
  x <- gsub("[^A-Za-z-]", "", x)
  ifelse(nzchar(x), title_case_1(x), NA_character_)
}

is_binomial <- function(x) {
  x <- tidy_space(as.character(x))
  grepl("^\\S+\\s+\\S+", x)
}

fix_glued_species <- function(genus, species) {
  genus <- as.character(genus)
  species <- as.character(species)
  needs <- !is.na(genus) & !is.na(species) &
    grepl(paste0("^", genus, "[A-Za-z]"), species)

  species[needs] <- sub(
    paste0("^(", genus[needs], ")([A-Za-z])"),
    "\\1 \\2",
    species[needs]
  )
  species
}

normalize_taxon_mode <- function(x) {
  x0 <- tolower(trimws(as.character(x %||% "")))
  map <- c(
    "family" = "family", "families" = "family", "fam" = "family",
    "genus" = "genus", "gen" = "genus",
    "species" = "species", "specie" = "species", "sp" = "species",
    "kingdom" = "kingdom", "kingdoms" = "kingdom"
  )
  x1 <- map[[x0]]
  if (is.null(x1)) {
    stop(
      sprintf("Invalid TAXON_MODE: '%s'. Use: family | genus | species | kingdom.", x0),
      call. = FALSE
    )
  }
  x1
}

TAXON_MODE <- normalize_taxon_mode(cfg$taxon_mode %||% "family")
TAXON_VALUES <- cfg$taxon_values %||% character(0)
TAXON_VALUES <- unique(na.omit(trimws(as.character(TAXON_VALUES))))
if (length(TAXON_VALUES) == 0L) {
  stop("TAXON_VALUES is empty after normalization.", call. = FALSE)
}

MONGO_URL <- Sys.getenv(
  "LOTUS_MONGO_URL",
  unset = (
    cfg$mongo_url %||%
      "mongodb://127.0.0.1:27017/?socketTimeoutMS=3600000&connectTimeoutMS=300000&serverSelectionTimeoutMS=300000"
  )
)

DB_NAME <- cfg$db_name %||% "lotus"
COLL_NAME <- cfg$coll_name %||% "lotusUniqueNaturalProduct"
opts <- cfg$mongo_opts %||% '{"allowDiskUse": true, "batchSize": 5000}'

PAGE <- cfg$page_size_lines %||% 50000L
chunk_size <- cfg$chunk_size_inchikey %||% 1000L

safe_tag <- function(mode, values, run_date, suffix = NULL) {
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

tag_base <- if (!is.null(cfg$prefix_base_tag) && nzchar(cfg$prefix_base_tag)) {
  cfg$prefix_base_tag
} else {
  safe_tag(
    TAXON_MODE,
    TAXON_VALUES,
    cfg$run_tag_date %||% Sys.Date(),
    cfg$custom_tag_suffix %||% NULL
  )
}

out_dir_base <- cfg$out_dir_base %||% getwd()
OUT_DIR_RUN <- file.path(out_dir_base, paste0("lotus_", tag_base))
dir.create(OUT_DIR_RUN, showWarnings = FALSE, recursive = TRUE)

if (exists("runtime") && !is.null(runtime$OUT_DIR) && nzchar(runtime$OUT_DIR)) {
  OUT_DIR <- runtime$OUT_DIR
} else {
  OUT_DIR <- OUT_DIR_RUN
}
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

PART1_DIR <- file.path(OUT_DIR, "PartI_ALL")
dir.create(PART1_DIR, showWarnings = FALSE, recursive = TRUE)

safe_file <- function(base, ext) file.path(OUT_DIR, paste0(base, ext))
base_tag <- paste0("lotus_", tag_base)

PROPS_CORE_FIELDS <- cfg$props_core_fields %||% c(
  "lotus_id", "wikidata_id", "inchikey", "smiles", "iupac_name",
  "molecular_formula", "molecular_weight",
  "xlogp", "alogp", "amralogp", "manholdlogp",
  "topoPSA", "tpsaEfficiency", "fsp3",
  "hBondAcceptorCount", "hBondDonorCount", "LipinskiRuleOf5Failures",
  "contains_sugar", "contains_ring_sugars", "contains_linear_sugars",
  "number_of_carbons", "number_of_oxygens", "number_of_nitrogens",
  "total_atom_number", "heavy_atom_number", "max_number_of_rings", "min_number_of_rings",
  "murko_framework", "ertlFunctionalFragmentsPseudoSmiles",
  "chemicalTaxonomyNPclassifierSuperclass", "chemicalTaxonomyNPclassifierClass",
  "chemicalTaxonomyClassyfireSuperclass", "chemicalTaxonomyClassyfireClass",
  "traditional_name", "allWikidataIds"
)

collapse <- function(x) {
  x <- unique(x)
  x <- x[!is.na(x) & nzchar(as.character(x))]
  if (!length(x)) return(NA_character_)
  paste(x, collapse = ";")
}

safe_first <- function(x) {
  x <- x[!is.na(x) & nzchar(as.character(x))]
  if (!length(x)) return(NA_character_)
  x[1]
}

lotus <- mongo(collection = COLL_NAME, db = DB_NAME, url = MONGO_URL)
if (isTRUE(cfg$verbose %||% TRUE)) {
  cat("[INFO] Connected to MongoDB. Total documents:", lotus$count(), "\n")
  cat(
    sprintf(
      "[INFO] Part I scope: TAXON_MODE=%s | TAXON_VALUES=%s\n",
      TAXON_MODE,
      paste(TAXON_VALUES, collapse = ", ")
    )
  )
}

regex_escape <- function(x) {
  x <- as.character(x)
  x <- stringr::str_squish(x)
  gsub("([\\.^$|()?*+\\[\\]{}-])", "\\\\\\1", x, perl = TRUE)
}

build_taxon_match <- function(prefix, mode, values) {
  vals <- unique(na.omit(trimws(values)))
  vals <- vals[nzchar(vals)]
  if (!length(vals)) stop("Empty or invalid TAXON_VALUES.")

  or_list <- list()
  for (v0 in vals) {
    v <- regex_escape(v0)
    if (mode == "genus") {
      or_list <- c(
        or_list,
        list(setNames(list(list("$regex" = paste0("^", v, "$"), "$options" = "i")), paste0(prefix, "genus"))),
        list(setNames(list(list("$regex" = paste0("^", v, "\\s"), "$options" = "i")), paste0(prefix, "species")))
      )
    } else if (mode == "family") {
      or_list <- c(
        or_list,
        list(setNames(list(list("$regex" = paste0("^", v, "$"), "$options" = "i")), paste0(prefix, "family")))
      )
    } else if (mode == "species") {
      or_list <- c(
        or_list,
        list(setNames(list(list("$regex" = paste0("^", v), "$options" = "i")), paste0(prefix, "species"))),
        list(setNames(list(list("$regex" = paste0("^", v), "$options" = "i")), paste0(prefix, "organism_value"))),
        list(setNames(list(list("$regex" = paste0("^", v), "$options" = "i")), paste0(prefix, "cleaned_organism_id")))
      )
    } else if (mode == "kingdom") {
      or_list <- c(
        or_list,
        list(setNames(list(list("$regex" = paste0("^", v, "$"), "$options" = "i")), paste0(prefix, "kingdom")))
      )
    } else {
      stop("Invalid TAXON_MODE.")
    }
  }
  list("$or" = or_list)
}

pipe_count <- list(
  list("$project" = list(tx1 = list("$objectToArray" = "$taxonomyReferenceObjects"))),
  list("$unwind" = "$tx1"),
  list("$project" = list(tx2 = list("$objectToArray" = "$tx1.v"))),
  list("$unwind" = "$tx2"),
  list("$unwind" = "$tx2.v"),
  list("$match" = build_taxon_match("tx2.v.", TAXON_MODE, TAXON_VALUES)),
  list("$count" = "n")
)

cnt <- lotus$aggregate(jsonlite::toJSON(pipe_count, auto_unbox = TRUE), options = opts)
cnt <- as_base_df(cnt)
total_lines <- if (nrow(cnt)) cnt$n[1] else 0L
if (isTRUE(cfg$verbose %||% TRUE)) {
  cat("[INFO] Matched rows (compound x reference x source x species):", total_lines, "\n")
}
if (total_lines == 0L) {
  stop("No records were retrieved for the selected taxonomic criterion.")
}

n_batches <- ceiling(total_lines / PAGE)
pb <- progress::progress_bar$new(
  format = "Batch :current/:total [:bar] :percent | rows=:rows | :elapsed (ETA :eta)",
  total = n_batches,
  clear = FALSE,
  width = 80
)

processed <- 0L
pages <- vector("list", n_batches)

for (i in seq_len(n_batches)) {
  skip_rows <- (i - 1L) * PAGE
  pipe_page <- list(
    list("$project" = list(
      `_id` = 0,
      lotus_id = 1,
      smiles = 1,
      inchikey = 1,
      iupac_name = 1,
      molecular_formula = 1,
      tx1 = list("$objectToArray" = "$taxonomyReferenceObjects")
    )),
    list("$unwind" = "$tx1"),
    list("$project" = list(
      lotus_id = 1,
      smiles = 1,
      inchikey = 1,
      iupac_name = 1,
      molecular_formula = 1,
      ref_id = "$tx1.k",
      srcMap = "$tx1.v"
    )),
    list("$project" = list(
      lotus_id = 1,
      smiles = 1,
      inchikey = 1,
      iupac_name = 1,
      molecular_formula = 1,
      ref_id = 1,
      src = list("$objectToArray" = "$srcMap")
    )),
    list("$unwind" = "$src"),
    list("$unwind" = "$src.v"),
    list("$match" = build_taxon_match("src.v.", TAXON_MODE, TAXON_VALUES)),
    list("$project" = list(
      lotus_id = 1,
      smiles = 1,
      inchikey = 1,
      iupac_name = 1,
      molecular_formula = 1,
      ref_id = 1,
      source = "$src.k",
      family = "$src.v.family",
      genus = "$src.v.genus",
      species = "$src.v.species"
    )),
    list("$sort" = list("lotus_id" = 1L, "ref_id" = 1L, "source" = 1L, "species" = 1L)),
    list("$skip" = skip_rows),
    list("$limit" = PAGE)
  )

  page <- lotus$aggregate(jsonlite::toJSON(pipe_page, auto_unbox = TRUE), options = opts)
  page <- as_base_df(page)

  if (nrow(page)) {
    pages[[i]] <- page
    processed <- processed + nrow(page)
  } else {
    pages[[i]] <- NULL
  }

  pb$tick(tokens = list(rows = format(processed, big.mark = ",", decimal.mark = ".", scientific = FALSE)))
}

lin <- as_base_df(dplyr::bind_rows(pages))
lin <- as_base_df(dplyr::distinct(
  lin,
  lotus_id,
  ref_id,
  source,
  family,
  genus,
  species,
  inchikey,
  smiles,
  iupac_name,
  molecular_formula,
  .keep_all = TRUE
))

if (isTRUE(cfg$verbose %||% TRUE)) {
  cat("[OK] Rows loaded into memory:", nrow(lin), "\n")
}

lin <- as_base_df(
  lin %>%
    dplyr::mutate(ref_id = fix_ref_id(ref_id))
)

if (identical(TAXON_MODE, "genus")) {
  target_gen <- unique(tolower(trimws(TAXON_VALUES)))
  target_gen <- target_gen[nzchar(target_gen)]
  esc <- function(x) gsub("([\\.^$|()?*+\\[\\]{}-])", "\\\\\\1", x, perl = TRUE)
  pat_gen <- paste0("^(", paste(esc(target_gen), collapse = "|"), ")$")
  pat_pref <- paste0("^(", paste(esc(target_gen), collapse = "|"), ")\\s")

  lin <- as_base_df(
    lin %>%
      dplyr::mutate(
        .g0 = tolower(trimws(genus)),
        .s0 = tolower(trimws(species)),
        .s_ok = !is.na(.s0) & grepl(pat_pref, .s0, perl = TRUE),
        .g_ok = !is.na(.g0) & grepl(pat_gen, .g0, perl = TRUE)
      )
  )

  keep_mask <- with(lin, (.s_ok) | (.g_ok & (is.na(.s0) | .s_ok)))
  drop_df <- as_base_df(lin[!keep_mask, c("lotus_id", "ref_id", "source", "family", "genus", "species"), drop = FALSE])

  if (nrow(drop_df)) {
    utils::write.table(
      drop_df,
      file = file.path(OUT_DIR, paste0(base_tag, "_genus_mode_inconsistent_pairs.tsv")),
      sep = "\t",
      quote = TRUE,
      row.names = FALSE,
      col.names = TRUE
    )

    message(
      sprintf(
        "[FILTER] Genus-mode consistency filtering removed %d records. Output saved to OUT_DIR.",
        nrow(drop_df)
      )
    )
  }

  lin <- as_base_df(lin[keep_mask, , drop = FALSE])
  lin$.g0 <- NULL
  lin$.s0 <- NULL
  lin$.s_ok <- NULL
  lin$.g_ok <- NULL
}

fields_props <- sprintf(
  '{"_id":0,%s}',
  paste(sprintf('"%s":1', PROPS_CORE_FIELDS), collapse = ",")
)

scalarize_field <- function(doc, nm) {
  val <- doc[[nm]]
  if (is.null(val)) return(NA_character_)
  if (is.atomic(val) && length(val) == 1) return(as.character(val))
  v <- tryCatch(unlist(val, use.names = FALSE), error = function(e) val)
  v <- v[!is.na(v)]
  v <- v[is.atomic(v)]
  if (!length(v)) return(NA_character_)
  paste(unique(as.character(v)), collapse = ";")
}

inchis <- unique(stats::na.omit(lin$inchikey))

if ((!length(inchis) || all(is.na(inchis))) && isTRUE(cfg$stop_if_no_inchikey %||% TRUE)) {
  stop("No valid InChIKeys were detected in 'lin'. Property enrichment cannot proceed.")
}

props_list <- vector("list", ceiling(length(inchis) / chunk_size))
for (j in seq_along(props_list)) {
  idx <- ((j - 1L) * chunk_size + 1):min(j * chunk_size, length(inchis))
  q <- jsonlite::toJSON(list(inchikey = list(`$in` = unname(inchis[idx]))), auto_unbox = TRUE)
  it <- lotus$iterate(query = q, fields = fields_props)
  rows <- list()
  k <- 0L

  repeat {
    doc <- it$one()
    if (is.null(doc)) break
    k <- k + 1L
    if (!is.list(doc)) doc <- as.list(doc)
    rows[[k]] <- setNames(
      lapply(PROPS_CORE_FIELDS, function(nm) scalarize_field(doc, nm)),
      PROPS_CORE_FIELDS
    )
  }

  props_list[[j]] <- if (k > 0L) as_base_df(dplyr::bind_rows(rows)) else NULL
}

props_core <- unique(as_base_df(dplyr::bind_rows(props_list)))
if (isTRUE(cfg$verbose %||% TRUE)) {
  cat("[OK] Chemical properties retrieved for unique InChIKeys:", nrow(props_core), "\n")
}

unify_dupes <- function(df, bases = c("lotus_id", "smiles", "iupac_name", "molecular_formula", "molecular_weight")) {
  for (b in bases) {
    x <- paste0(b, ".x")
    y <- paste0(b, ".y")
    has_x <- x %in% names(df)
    has_y <- y %in% names(df)
    if (has_x && has_y) {
      df[[b]] <- dplyr::coalesce(df[[x]], df[[y]])
      df[[x]] <- NULL
      df[[y]] <- NULL
    } else if (has_x && !has_y) {
      df[[b]] <- df[[x]]
      df[[x]] <- NULL
    } else if (!has_x && has_y) {
      df[[b]] <- df[[y]]
      df[[y]] <- NULL
    }
  }
  df
}

lin_enriched <- as_base_df(
  lin %>%
    dplyr::left_join(props_core, by = "inchikey") %>%
    unify_dupes()
)

NUMERIC_PROPS <- cfg$numeric_props %||% c(
  "molecular_weight", "xlogp", "alogp", "amralogp", "manholdlogp",
  "topoPSA", "tpsaEfficiency", "fsp3",
  "hBondAcceptorCount", "hBondDonorCount",
  "number_of_carbons", "number_of_oxygens", "number_of_nitrogens",
  "total_atom_number", "heavy_atom_number",
  "max_number_of_rings", "min_number_of_rings",
  "LipinskiRuleOf5Failures"
)

LOGICAL_PROPS <- cfg$logical_props %||% c(
  "contains_ring_sugars", "contains_linear_sugars", "contains_sugar"
)

lin_enriched <- as_base_df(
  lin_enriched %>%
    dplyr::mutate(dplyr::across(all_of(NUMERIC_PROPS), ~ suppressWarnings(as.numeric(.)))) %>%
    dplyr::mutate(dplyr::across(all_of(LOGICAL_PROPS), ~ {
      if (is.logical(.)) . else tolower(as.character(.)) %in% c("true", "t", "1")
    }))
)

USE_WFO <- isTRUE(cfg$use_WFO_normalization %||% TRUE)

if (USE_WFO) {
  cat("[WFO] Starting taxonomic normalization.\n")

  canon_name <- function(x) {
    x <- tolower(stringr::str_squish(as.character(x)))
    x <- gsub("\\s+\\(.*?\\)", "", x)
    x <- gsub("\\b(ex|sensu|auct\\.|non)\\b.*$", "", x)
    x <- gsub("\\b(subsp\\.|ssp\\.|var\\.|subvar\\.|f\\.|forma|cv\\.|group)\\b", "", x)
    x <- stringr::str_squish(x)
    toks <- strsplit(x, "\\s+")
    vapply(
      toks,
      function(tt) {
        if (length(tt) >= 3) {
          paste(tt[1:3], collapse = " ")
        } else if (length(tt) >= 2) {
          paste(tt[1:2], collapse = " ")
        } else {
          tt[1] %||% NA_character_
        }
      },
      FUN.VALUE = character(1)
    )
  }

  mk_join_key_vec <- function(genus, species) {
    cand <- ifelse(!is.na(species) & nzchar(species), species, paste(genus, species))
    tolower(stringr::str_squish(cand))
  }

  first_non_na_chr <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x)) x[1] else NA_character_
  }

  split_scientific2 <- function(nm) {
    nm <- stringr::str_squish(as.character(nm))
    parts <- strsplit(nm, "\\s+")
    data.frame(
      corrected_genus = vapply(parts, function(p) if (length(p) >= 1) p[1] else NA_character_, ""),
      corrected_specific_epithet = vapply(parts, function(p) if (length(p) >= 2) p[2] else NA_character_, ""),
      corrected_infraspecific = vapply(parts, function(p) {
        if (length(p) >= 3) paste(p[3:length(p)], collapse = " ") else NA_character_
      }, ""),
      stringsAsFactors = FALSE
    )
  }

  WFO_CSV_PATH <- cfg$wfo_csv_path %||% here::here("data", "wfo", "classification.tsv")
  TODAY_STR <- format(Sys.Date())
  col_taxonID <- cfg$wfo_cols$taxonID %||% "taxonID"
  col_sciName <- cfg$wfo_cols$scientificName %||% "scientificName"
  col_status <- cfg$wfo_cols$taxonomicStatus %||% "taxonomicStatus"
  col_accID <- cfg$wfo_cols$acceptedNameUsageID %||% "acceptedNameUsageID"
  col_family <- cfg$wfo_cols$family %||% "family"
  col_genus <- cfg$wfo_cols$genus %||% "genus"

  wfo_raw <- as_base_df(readr::read_tsv(
    file = WFO_CSV_PATH,
    col_types = readr::cols(.default = readr::col_character()),
    progress = TRUE,
    locale = readr::locale(encoding = "UTF-8"),
    na = c("", "NA", "NULL")
  ))
  if (ncol(wfo_raw) == 1L) stop("The WFO file appears to be invalid. Please verify the delimiter.")

  wfo <- as_base_df(
    wfo_raw %>%
      dplyr::transmute(
        taxonID = .data[[col_taxonID]],
        name = .data[[col_sciName]],
        status = tolower(stringr::str_squish(.data[[col_status]])),
        accID = .data[[col_accID]],
        family = .data[[col_family]],
        genus = .data[[col_genus]]
      )
  )

  accepted <- as_base_df(
    wfo %>%
      dplyr::filter(status == "accepted") %>%
      dplyr::transmute(accepted_id = taxonID, accepted_name = name)
  )

  synonyms <- as_base_df(
    wfo %>%
      dplyr::filter(status != "accepted", !is.na(accID), nzchar(accID))
  )

  syn_map <- as_base_df(
    synonyms %>%
      dplyr::left_join(accepted, by = c("accID" = "accepted_id")) %>%
      dplyr::transmute(
        synonym_name = name,
        synonym_id = taxonID,
        accepted_name = accepted_name,
        accepted_id = accID
      )
  )

  acc_key <- as_base_df(
    accepted %>%
      dplyr::transmute(key = canon_name(accepted_name), accepted_name, accepted_id) %>%
      dplyr::filter(!is.na(key) & nzchar(key)) %>%
      dplyr::distinct(key, .keep_all = TRUE)
  )

  syn_key <- as_base_df(
    syn_map %>%
      dplyr::transmute(key = canon_name(synonym_name), accepted_name, accepted_id) %>%
      dplyr::filter(!is.na(key) & nzchar(key)) %>%
      dplyr::distinct(key, .keep_all = TRUE)
  )

  dict <- as_base_df(
    dplyr::bind_rows(acc_key, syn_key) %>%
      dplyr::distinct(key, .keep_all = TRUE)
  )

  if (isTRUE(cfg$verbose %||% TRUE)) {
    cat("[WFO] Dictionary loaded:", nrow(dict), "unique keys.\n")
  }

  lin_join <- mk_join_key_vec(lin$genus, lin$species)
  resolved <- as_base_df(
    data.frame(
      original_name = lin_join,
      key = canon_name(lin_join),
      stringsAsFactors = FALSE
    ) %>%
      dplyr::left_join(dict, by = "key") %>%
      dplyr::mutate(
        tax_provider = ifelse(!is.na(accepted_name), "WFO_offline", NA_character_),
        tax_status = dplyr::case_when(
          !is.na(accepted_name) & canon_name(original_name) != canon_name(accepted_name) ~ "synonym",
          !is.na(accepted_name) ~ "accepted",
          TRUE ~ NA_character_
        ),
        tax_checked_at = TODAY_STR
      )
  )

  crosswalk_wfo <- as_base_df(
    resolved %>%
      dplyr::transmute(original_name, accepted_name, accepted_id, tax_status, tax_provider, tax_checked_at) %>%
      dplyr::group_by(original_name) %>%
      dplyr::summarise(
        accepted_name = first_non_na_chr(accepted_name),
        accepted_id = first_non_na_chr(accepted_id),
        tax_status = first_non_na_chr(tax_status),
        tax_provider = first_non_na_chr(tax_provider),
        tax_checked_at = first_non_na_chr(tax_checked_at),
        .groups = "drop"
      )
  )

  lin_wfo <- as_base_df(
    lin_enriched %>%
      dplyr::mutate(.original_name = mk_join_key_vec(genus, species)) %>%
      dplyr::left_join(crosswalk_wfo, by = c(".original_name" = "original_name"))
  )

  accepted_meta <- as_base_df(
    wfo %>%
      dplyr::filter(status == "accepted") %>%
      dplyr::transmute(
        accepted_id = taxonID,
        accepted_name_w = name,
        accepted_family = family,
        accepted_genus = genus
      ) %>%
      dplyr::distinct(accepted_id, .keep_all = TRUE)
  )

  for (col in c("genus", "species", "family")) {
    if (!col %in% names(lin_wfo)) lin_wfo[[col]] <- NA_character_
  }

  lin_wfo <- as_base_df(
    lin_wfo %>%
      dplyr::mutate(
        genus = as.character(genus),
        species = as.character(species),
        family = as.character(family)
      )
  )

  lin_applied <- as_base_df(
    lin_wfo %>%
      dplyr::mutate(
        .corr_scientific = dplyr::coalesce(accepted_name, .original_name),
        .corr_scientific = stringr::str_squish(.corr_scientific)
      ) %>%
      dplyr::left_join(accepted_meta, by = "accepted_id")
  )

  lin_applied <- as_base_df(cbind(lin_applied, split_scientific2(lin_applied$.corr_scientific)))

  lin_applied <- as_base_df(
    lin_applied %>%
      dplyr::mutate(
        corrected_family = dplyr::coalesce(accepted_family, family),
        taxonomy_action = dplyr::case_when(
          is.na(accepted_name) ~ "unresolved",
          canon_name(.original_name) == canon_name(accepted_name) ~ "accepted_confirmed",
          TRUE ~ "synonym_replaced"
        ),
        tax_checked_at = TODAY_STR
      ) %>%
      dplyr::mutate(
        species_original = species,
        genus_original = genus,
        family_original = family,
        species = dplyr::coalesce(.corr_scientific, species),
        genus = dplyr::coalesce(corrected_genus, genus),
        family = dplyr::coalesce(corrected_family, family)
      )
  )

  lin_enriched <- lin_applied
  cat("[WFO] Taxonomic normalization completed successfully.\n")

  genus_or_lower <- function(x) tolower(tidy_space(as.character(x)))
  if (nrow(lin_enriched) > 0) {
    genus_original <- genus_or_lower(dplyr::coalesce(lin_enriched$genus_original, lin_enriched$genus))
    genus_final <- genus_or_lower(dplyr::coalesce(
      lin_enriched$accepted_genus,
      lin_enriched$corrected_genus,
      lin_enriched$genus
    ))
    target_genera_norm <- genus_or_lower(cfg$taxon_values)

    unresolved_flag <- grepl(
      "unresolved",
      tolower(dplyr::coalesce(lin_enriched$tax_status, lin_enriched$taxonomy_action, "")),
      fixed = TRUE
    )
    reclassified_flag <- nzchar(genus_original) & nzchar(genus_final) & (genus_original != genus_final)
    outside_target <- !(genus_final %in% target_genera_norm)

    remove_vec <- rep(FALSE, nrow(lin_enriched))
    if (identical(tolower(cfg$taxon_mode), "genus")) {
      remove_vec <- outside_target | reclassified_flag
    }
    keep_row <- (!unresolved_flag) & (!remove_vec)

    lin_enriched <- as_base_df(lin_enriched[keep_row, , drop = FALSE])
    cat(sprintf("[WFO] Canonical genus filtering retained %d rows.\n", nrow(lin_enriched)))
  }
}

lin_enriched <- as_base_df(
  lin_enriched %>%
    dplyr::mutate(
      genus = canon_genus(genus),
      species = gsub("\\bsp\\.?\\b", "", species, ignore.case = TRUE),
      species = stringr::str_squish(species),
      species = fix_glued_species(genus, species)
    )
)

dedup_before <- nrow(lin_enriched)

lin_enriched <- as_base_df(
  lin_enriched %>%
    dplyr::arrange(inchikey, family, genus, species, ref_id, source) %>%
    dplyr::distinct(inchikey, family, genus, species, ref_id, .keep_all = TRUE)
)

dedup_after <- nrow(lin_enriched)
cat(sprintf(
  "[DEDUP] Compound x Species x Reference: %d -> %d rows (%d duplicates removed).\n",
  dedup_before,
  dedup_after,
  dedup_before - dedup_after
))

build_map_tax_inchi <- function(lin_tbl, tax_col = c("family", "genus", "species")) {
  tax_col <- match.arg(tax_col)

  out <- as_base_df(
    lin_tbl %>%
      dplyr::transmute(
        inchikey = as.character(inchikey),
        family = as.character(family),
        genus = canon_genus(genus),
        species = species
      ) %>%
      dplyr::mutate(
        species = gsub("\\bsp\\.?\\b", "", species, ignore.case = TRUE),
        species = stringr::str_squish(species),
        species = fix_glued_species(genus, species)
      ) %>%
      dplyr::filter(!is.na(inchikey), nzchar(inchikey))
  )

  if (identical(tax_col, "genus")) {
    out <- as_base_df(out %>% dplyr::filter(!is.na(genus), nzchar(genus)))
  } else if (identical(tax_col, "species")) {
    out <- as_base_df(out %>% dplyr::filter(is_binomial(species)))
  } else if (identical(tax_col, "family")) {
    out <- as_base_df(out %>% dplyr::filter(!is.na(family), nzchar(family)))
  }

  out <- as_base_df(
    out %>%
      dplyr::mutate(
        taxon = dplyr::case_when(
          identical(tax_col, "family") ~ family,
          identical(tax_col, "genus") ~ genus,
          identical(tax_col, "species") ~ species,
          TRUE ~ species
        )
      ) %>%
      dplyr::filter(!is.na(taxon), nzchar(taxon)) %>%
      dplyr::distinct(inchikey, taxon, family, genus, species)
  )

  if (!tax_col %in% names(out)) out[[tax_col]] <- out$taxon
  as_base_df(out)
}

tax_col_for_map <- cfg$analysis_tax_level %||% "genus"
map_tax_inchi <- as_base_df(build_map_tax_inchi(lin_enriched, tax_col = tax_col_for_map))
cat(sprintf("[INFO] Taxon-to-InChIKey map generated at the %s level. Rows: %d\n", tax_col_for_map, nrow(map_tax_inchi)))

uni <- as_base_df(
  lin_enriched %>%
    dplyr::group_by(inchikey) %>%
    dplyr::summarise(
      lotus_id = safe_first(lotus_id),
      smiles = safe_first(smiles),
      iupac_name = safe_first(iupac_name),
      molecular_formula = safe_first(molecular_formula),
      genus = collapse(genus),
      family = collapse(family),
      species = collapse(species),
      ref_ids = collapse(ref_id),
      .groups = "drop"
    )
)

uni_enriched <- as_base_df(
  uni %>%
    dplyr::left_join(props_core, by = "inchikey") %>%
    unify_dupes() %>%
    dplyr::mutate(
      dplyr::across(all_of(NUMERIC_PROPS), ~ suppressWarnings(as.numeric(.))),
      dplyr::across(all_of(LOGICAL_PROPS), ~ {
        if (is.logical(.)) . else tolower(as.character(.)) %in% c("true", "t", "1")
      })
    )
)

if (isTRUE(cfg$verbose %||% TRUE)) {
  cat("[OK] UNI tables generated. UNI:", nrow(uni), "rows | UNI_ENRICHED:", nrow(uni_enriched), "rows.\n")
}

dedup_compound_species <- function(df) {
  num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  logi_cols <- names(df)[vapply(df, is.logical, logical(1))]
  chr_cols <- names(df)[vapply(df, is.character, logical(1))]
  chr_cols <- setdiff(chr_cols, c("inchikey", "family", "genus", "species", "ref_id", "source"))

  as_base_df(
    df %>%
      dplyr::group_by(inchikey, species) %>%
      dplyr::summarise(
        family = collapse(family),
        genus = collapse(genus),
        ref_id = collapse(ref_id),
        source = collapse(source),
        dplyr::across(dplyr::all_of(num_cols), ~ {
          v <- .[!is.na(.)]
          if (length(v)) v[1] else NA_real_
        }, .names = "{.col}"),
        dplyr::across(dplyr::all_of(logi_cols), ~ any(., na.rm = TRUE), .names = "{.col}"),
        dplyr::across(dplyr::all_of(chr_cols), ~ safe_first(.), .names = "{.col}"),
        .groups = "drop"
      )
  )
}

lin_cs <- dedup_compound_species(lin_enriched)
cat(sprintf("[INFO] Compound x Species table generated. Rows: %d\n", nrow(lin_cs)))

EXPORT_EXCEL <- isTRUE(cfg$export_excel %||% TRUE)
EXPORT_PARQUET <- isTRUE(cfg$export_parquet %||% TRUE)

if (EXPORT_EXCEL) {
  if (!requireNamespace("writexl", quietly = TRUE)) {
    stop("Missing package: writexl. Run renv::restore() or install.packages('writexl').")
  }

  MAX_XLSX <- cfg$max_xlsx_cell_chars %||% 32767L
  BIG_COLS <- cfg$big_cols %||% c(
    "allWikidataIds", "xrefs", "pubchemBitsString", "pubchemBits",
    "circularFingerprint", "extendedFingerprint", "pfCounts",
    "ertlFunctionalFragmentsPseudoSmiles"
  )

  collapse_list <- function(x) {
    sapply(x, function(v) {
      if (is.null(v)) return(NA_character_)
      v <- tryCatch(unlist(v, use.names = FALSE), error = function(e) v)
      v <- v[!is.na(v)]
      if (!length(v)) return(NA_character_)
      paste(unique(as.character(v)), collapse = ";")
    })
  }

  trim_cell <- function(x, max_len = MAX_XLSX) {
    if (is.null(x)) return(x)
    if (!is.character(x)) x <- as.character(x)
    n <- nchar(x, allowNA = TRUE)
    too_long <- !is.na(n) & n > max_len
    x[too_long] <- paste0(substr(x[too_long], 1, max_len - 3), "...")
    x
  }

  sanitize_for_excel <- function(df, drop_cols = NULL) {
    df <- as_base_df(df)
    df <- as_base_df(dplyr::mutate(df, dplyr::across(where(is.factor), as.character)))
    df <- as_base_df(dplyr::mutate(df, dplyr::across(where(is.list), collapse_list)))
    if (!is.null(drop_cols)) {
      keep <- setdiff(names(df), drop_cols)
      df <- as_base_df(df[keep])
    }
    as_base_df(dplyr::mutate(df, dplyr::across(where(is.character), trim_cell)))
  }

  lin_x <- sanitize_for_excel(lin, drop_cols = intersect(BIG_COLS, names(lin)))
  lin_enriched_x <- sanitize_for_excel(lin_enriched, drop_cols = intersect(BIG_COLS, names(lin_enriched)))
  uni_x <- sanitize_for_excel(uni, drop_cols = intersect(BIG_COLS, names(uni)))
  uni_enriched_x <- sanitize_for_excel(uni_enriched, drop_cols = intersect(BIG_COLS, names(uni_enriched)))
  lin_cs_x <- sanitize_for_excel(lin_cs, drop_cols = intersect(BIG_COLS, names(lin_cs)))

  xlsx_path <- safe_file(base_tag, ".xlsx")
  writexl::write_xlsx(
    list(
      lin = lin_x,
      lin_enriched = lin_enriched_x,
      lin_compound_species = lin_cs_x,
      uni = uni_x,
      uni_enriched = uni_enriched_x
    ),
    path = xlsx_path
  )
  cat("[OK] Excel workbook written to:", normalizePath(xlsx_path), "\n")
}

if (EXPORT_PARQUET) {
  if (!requireNamespace("arrow", quietly = TRUE)) {
    if (isTRUE(cfg$auto_install_missing_packages %||% TRUE)) install.packages("arrow")
  }
  if (requireNamespace("arrow", quietly = TRUE)) {
    arrow::write_parquet(as_base_df(lin_enriched), safe_file(paste0(base_tag, "_lin_enriched"), ".parquet"))
    arrow::write_parquet(as_base_df(uni_enriched), safe_file(paste0(base_tag, "_uni_enriched"), ".parquet"))
    arrow::write_parquet(as_base_df(lin_cs), safe_file(paste0(base_tag, "_lin_compound_species"), ".parquet"))
    cat("[OK] Parquet files written to:", normalizePath(OUT_DIR), "\n")
  } else {
    warning("Package 'arrow' is not available; Parquet export was skipped.")
  }
}
