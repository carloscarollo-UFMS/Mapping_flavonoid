# Part I v2.2 — extraction, WFO hierarchical taxonomic resolution, and audit
# Hierarchy: exact accepted name > exact synonym consensus > canonical consensus.
# Conflicting taxonomies remain unresolved and are exported for manual review.

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

if (!exists("cfg")) {
  stop("Part I requires object 'cfg' created by Main_Pipeline_end.R.")
}

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
    startsWith(tolower(species), tolower(genus)) &
    !startsWith(tolower(species), paste0(tolower(genus), " "))

  idx <- which(needs)
  if (length(idx)) {
    species[idx] <- vapply(idx, function(i) {
      sub(
        paste0("^(", regex_escape(genus[i]), ")([A-Za-z])"),
        "\\1 \\2",
        species[i],
        perl = TRUE
      )
    }, FUN.VALUE = character(1))
  }
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


count_distinct_nonempty <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(trimws(x))]
  length(unique(x))
}

split_multivalue_tokens <- function(x, sep = ";") {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(trimws(x))]
  if (!length(x)) return(character(0))
  tokens <- unlist(strsplit(x, split = sep, fixed = TRUE), use.names = FALSE)
  tokens <- trimws(tokens)
  unique(tokens[nzchar(tokens)])
}

count_distinct_tokens <- function(x, sep = ";") {
  length(split_multivalue_tokens(x, sep = sep))
}

snapshot_counts <- function(stage, df = NULL, notes = NA_character_, n_rows_override = NULL) {
  if (is.null(df)) {
    return(data.frame(
      stage = stage,
      n_rows = if (is.null(n_rows_override)) NA_real_ else as.numeric(n_rows_override),
      n_compounds = NA_real_,
      n_families = NA_real_,
      n_genera = NA_real_,
      n_species = NA_real_,
      n_species_all = NA_real_,
      n_references = NA_real_,
      notes = notes,
      stringsAsFactors = FALSE
    ))
  }

  get_n <- function(nm, split_tokens = FALSE) {
    if (!nm %in% names(df)) return(NA_real_)
    if (isTRUE(split_tokens)) {
      return(as.numeric(count_distinct_tokens(df[[nm]])))
    }
    as.numeric(count_distinct_nonempty(df[[nm]]))
  }

  species_tokens <- if ("species" %in% names(df)) {
    split_multivalue_tokens(df$species)
  } else {
    character(0)
  }

  data.frame(
    stage = stage,
    n_rows = as.numeric(nrow(df)),
    n_compounds = get_n("inchikey"),
    n_families = get_n("family", split_tokens = TRUE),
    n_genera = get_n("genus", split_tokens = TRUE),
    n_species = as.numeric(count_distinct_nonempty(species_tokens[is_binomial(species_tokens)])),
    n_species_all = as.numeric(length(species_tokens)),
    n_references = get_n("ref_id", split_tokens = TRUE),
    notes = notes,
    stringsAsFactors = FALSE
  )
}

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

props_core_raw <- as_base_df(dplyr::bind_rows(props_list))


if (nrow(props_core_raw) > 0) {
  props_core <- as_base_df(
    props_core_raw %>%
      dplyr::group_by(inchikey) %>%
      dplyr::summarise(
        dplyr::across(
          dplyr::everything(),
          ~ safe_first(as.character(.)),
          .names = "{.col}"
        ),
        .groups = "drop"
      )
  )
} else {
  props_core <- as_base_df(
    stats::setNames(
      replicate(length(PROPS_CORE_FIELDS), character(0), simplify = FALSE),
      PROPS_CORE_FIELDS
    )
  )
}

if (isTRUE(cfg$verbose %||% TRUE)) {
  cat("[OK] Chemical properties retrieved for unique InChIKeys:", nrow(props_core), "\n")
}

property_match_summary <- data.frame(
  metric = c(
    "unique_inchikeys_requested",
    "unique_inchikeys_with_properties",
    "unique_inchikeys_without_properties",
    "property_match_percent"
  ),
  value = c(
    length(inchis),
    count_distinct_nonempty(props_core$inchikey),
    length(setdiff(inchis, props_core$inchikey)),
    if (length(inchis) > 0) {
      100 * count_distinct_nonempty(props_core$inchikey) / length(inchis)
    } else {
      NA_real_
    }
  ),
  stringsAsFactors = FALSE
)

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

lin_enriched_pre_wfo <- as_base_df(
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

lin_enriched_pre_wfo <- as_base_df(
  lin_enriched_pre_wfo %>%
    dplyr::mutate(dplyr::across(dplyr::any_of(NUMERIC_PROPS), ~ suppressWarnings(as.numeric(.)))) %>%
    dplyr::mutate(dplyr::across(dplyr::any_of(LOGICAL_PROPS), ~ {
      if (is.logical(.)) . else tolower(as.character(.)) %in% c("true", "t", "1")
    }))
)

USE_WFO <- isTRUE(cfg$use_WFO_normalization %||% TRUE)
WFO_CORE_POLICY <- tolower(as.character(cfg$wfo_core_policy %||% "resolved_only"))
valid_wfo_policies <- c("resolved_only", "family_available", "all")
if (!WFO_CORE_POLICY %in% valid_wfo_policies) {
  stop(
    "Invalid cfg$wfo_core_policy: ", WFO_CORE_POLICY,
    ". Use one of: ", paste(valid_wfo_policies, collapse = ", "), "."
  )
}

crosswalk_wfo <- data.frame(
  original_name = character(0), accepted_name = character(0),
  accepted_id = character(0), tax_status = character(0),
  tax_provider = character(0), tax_checked_at = character(0),
  stringsAsFactors = FALSE
)
wfo_ambiguous_keys <- data.frame(stringsAsFactors = FALSE)
wfo_unresolved_names <- data.frame(stringsAsFactors = FALSE)
wfo_unresolved_priority <- data.frame(stringsAsFactors = FALSE)
wfo_resolution_rule_summary <- data.frame(stringsAsFactors = FALSE)
wfo_resolution_comparison <- data.frame(stringsAsFactors = FALSE)
wfo_resolution_candidates <- data.frame(stringsAsFactors = FALSE)
wfo_family_species_totals <- data.frame(
  family = character(0), n_accepted_species_wfo = integer(0),
  stringsAsFactors = FALSE
)

build_original_scientific_name <- function(genus, species) {
  genus <- tidy_space(as.character(genus))
  species <- tidy_space(as.character(species))
  genus[!nzchar(genus)] <- NA_character_
  species[!nzchar(species)] <- NA_character_

  out <- species
  species_is_epithet <- !is.na(species) & !grepl("\\s", species)
  out[species_is_epithet & !is.na(genus)] <- paste(genus[species_is_epithet & !is.na(genus)], species[species_is_epithet & !is.na(genus)])
  out[is.na(out)] <- genus[is.na(out)]
  stringr::str_squish(out)
}

canon_name <- function(x) {
  x <- norm_ascii(as.character(x))
  x <- tolower(stringr::str_squish(x))
  x <- gsub("[×]", "x", x)
  x <- gsub("\\s+\\(.*?\\)", "", x)
  x <- gsub("\\b(ex|sensu|auct\\.|non)\\b.*$", "", x)
  x <- gsub("[,;]", " ", x)
  x <- stringr::str_squish(x)

  rank_tokens <- c("subsp.", "subsp", "ssp.", "ssp", "var.", "var", "subvar.", "subvar", "f.", "f", "forma")
  parts <- strsplit(x, "\\s+")

  vapply(parts, function(tt) {
    tt <- tt[nzchar(tt)]
    if (!length(tt)) return(NA_character_)
    if (length(tt) == 1L) return(tt[1])

    
    
    rank_pos <- which(tt %in% rank_tokens)
    if (length(rank_pos) && rank_pos[1] >= 3L && length(tt) >= rank_pos[1] + 1L) {
      return(paste(tt[c(1L, 2L, rank_pos[1], rank_pos[1] + 1L)], collapse = " "))
    }
    paste(tt[1:2], collapse = " ")
  }, FUN.VALUE = character(1))
}

first_non_na_chr <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x)) x[1] else NA_character_
}

split_scientific_name <- function(nm) {
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

if (USE_WFO) {
  cat("[WFO] Starting hierarchical taxonomic normalization.\n")

  WFO_CSV_PATH <- cfg$wfo_csv_path %||% here::here("data", "wfo", "classification.tsv")
  if (!file.exists(WFO_CSV_PATH)) {
    stop("WFO classification file not found: ", WFO_CSV_PATH)
  }

  TODAY_STR <- format(Sys.Date())
  wfo_cols <- cfg$wfo_cols %||% list()
  col_taxonID <- wfo_cols$taxonID %||% "taxonID"
  col_sciName <- wfo_cols$scientificName %||% "scientificName"
  col_status <- wfo_cols$taxonomicStatus %||% "taxonomicStatus"
  col_accID <- wfo_cols$acceptedNameUsageID %||% "acceptedNameUsageID"
  col_family <- wfo_cols$family %||% "family"
  col_genus <- wfo_cols$genus %||% "genus"
  col_rank <- wfo_cols$taxonRank %||% "taxonRank"

  normalize_exact_name <- function(x) {
    x <- norm_ascii(as.character(x))
    x <- tolower(stringr::str_squish(x))
    x[is.na(x) | !nzchar(x)] <- NA_character_
    x
  }

  make_taxonomy_signature <- function(name, family, genus, rank) {
    paste(
      dplyr::coalesce(as.character(name), "<NA>"),
      dplyr::coalesce(as.character(family), "<NA>"),
      dplyr::coalesce(as.character(genus), "<NA>"),
      dplyr::coalesce(as.character(rank), "<NA>"),
      sep = "||"
    )
  }

  collapse_sorted_chr <- function(x, sep = ";") {
    x <- sort(unique(as.character(x[!is.na(x) & nzchar(as.character(x))])))
    if (!length(x)) return(NA_character_)
    paste(x, collapse = sep)
  }

  wfo_header <- names(readr::read_tsv(
    file = WFO_CSV_PATH,
    n_max = 0,
    progress = FALSE,
    locale = readr::locale(encoding = "UTF-8"),
    name_repair = "minimal"
  ))

  required_wfo_cols <- c(
    col_taxonID, col_sciName, col_status, col_accID, col_family, col_genus
  )
  missing_wfo_cols <- setdiff(required_wfo_cols, wfo_header)
  if (length(missing_wfo_cols)) {
    stop("WFO file is missing required columns: ", paste(missing_wfo_cols, collapse = ", "))
  }

  rank_available <- col_rank %in% wfo_header
  wfo_cols_to_read <- unique(c(required_wfo_cols, if (rank_available) col_rank))

  wfo_raw <- as_base_df(readr::read_tsv(
    file = WFO_CSV_PATH,
    col_select = dplyr::all_of(wfo_cols_to_read),
    col_types = readr::cols(.default = readr::col_character()),
    progress = TRUE,
    locale = readr::locale(encoding = "UTF-8"),
    na = c("", "NA", "NULL"),
    name_repair = "minimal"
  ))
  if (ncol(wfo_raw) == 1L) stop("The WFO file appears invalid. Verify that it is tab-delimited.")

  wfo <- as_base_df(
    wfo_raw %>%
      dplyr::transmute(
        taxonID = as.character(.data[[col_taxonID]]),
        name = stringr::str_squish(as.character(.data[[col_sciName]])),
        status = tolower(stringr::str_squish(as.character(.data[[col_status]]))),
        accID = as.character(.data[[col_accID]]),
        family = stringr::str_squish(as.character(.data[[col_family]])),
        genus = stringr::str_squish(as.character(.data[[col_genus]])),
        taxon_rank = if (rank_available) {
          tolower(stringr::str_squish(as.character(.data[[col_rank]])))
        } else {
          NA_character_
        }
      ) %>%
      dplyr::mutate(
        exact_name = normalize_exact_name(name),
        canonical_key = canon_name(name)
      )
  )

  accepted <- as_base_df(
    wfo %>%
      dplyr::filter(status == "accepted", !is.na(taxonID), nzchar(taxonID)) %>%
      dplyr::transmute(
        accepted_id = taxonID,
        accepted_name = name,
        accepted_family = family,
        accepted_genus = genus,
        accepted_rank = taxon_rank
      ) %>%
      dplyr::distinct()
  )

  accepted_candidates <- as_base_df(
    accepted %>%
      dplyr::transmute(
        source_name = accepted_name,
        source_id = accepted_id,
        source_status = "accepted",
        accepted_id,
        accepted_name,
        accepted_family,
        accepted_genus,
        accepted_rank
      )
  )

  synonym_candidates <- as_base_df(
    wfo %>%
      dplyr::filter(status != "accepted", !is.na(accID), nzchar(accID)) %>%
      dplyr::transmute(
        source_name = name,
        source_id = taxonID,
        source_status = "synonym",
        accepted_id = accID
      ) %>%
      dplyr::left_join(accepted, by = "accepted_id") %>%
      dplyr::filter(!is.na(accepted_name), nzchar(accepted_name))
  )

  dict_candidates <- as_base_df(
    dplyr::bind_rows(accepted_candidates, synonym_candidates) %>%
      dplyr::mutate(
        candidate_exact_name = normalize_exact_name(source_name),
        candidate_key = canon_name(source_name),
        taxonomy_signature = make_taxonomy_signature(
          accepted_name, accepted_family, accepted_genus, accepted_rank
        )
      ) %>%
      dplyr::filter(
        !is.na(accepted_id), nzchar(accepted_id),
        !is.na(accepted_name), nzchar(accepted_name)
      ) %>%
      dplyr::distinct()
  )

  rm(wfo_raw, wfo, accepted_candidates, synonym_candidates)
  invisible(gc())

  original_names <- build_original_scientific_name(
    lin_enriched_pre_wfo$genus,
    lin_enriched_pre_wfo$species
  )

  target_names <- as_base_df(
    data.frame(
      original_name = original_names,
      exact_name = normalize_exact_name(original_names),
      key = canon_name(original_names),
      stringsAsFactors = FALSE
    ) %>%
      dplyr::distinct(original_name, exact_name, key) %>%
      dplyr::mutate(name_id = dplyr::row_number()) %>%
      dplyr::select(name_id, dplyr::everything())
  )

  exact_hits <- as_base_df(
    target_names %>%
      dplyr::filter(!is.na(exact_name), nzchar(exact_name)) %>%
      dplyr::inner_join(
        dict_candidates,
        by = c("exact_name" = "candidate_exact_name")
      )
  )

  canonical_hits <- as_base_df(
    target_names %>%
      dplyr::filter(!is.na(key), nzchar(key)) %>%
      dplyr::inner_join(
        dict_candidates,
        by = c("key" = "candidate_key")
      )
  )

  exact_split <- split(exact_hits, exact_hits$name_id)
  canonical_split <- split(canonical_hits, canonical_hits$name_id)

  make_resolution_row <- function(target_row) {
    id <- as.character(target_row$name_id)
    exact_group <- exact_split[[id]]
    canonical_group <- canonical_split[[id]]

    empty_result <- list(
      accepted_name = NA_character_,
      accepted_id = NA_character_,
      accepted_ids = NA_character_,
      accepted_family = NA_character_,
      accepted_genus = NA_character_,
      accepted_rank = NA_character_,
      source_status = NA_character_,
      tax_status = "unresolved",
      match_rule = "unresolved",
      resolution_status = "no_match",
      n_candidate_rows = 0L,
      n_candidate_accepted_ids = 0L,
      n_candidate_taxonomies = 0L
    )

    resolve_consensus <- function(candidates, prefix, source_preference = NULL) {
      if (is.null(candidates) || !nrow(candidates)) return(NULL)

      candidates <- as_base_df(candidates)
      candidates <- candidates %>%
        dplyr::filter(
          !is.na(accepted_id), nzchar(accepted_id),
          !is.na(accepted_name), nzchar(accepted_name)
        ) %>%
        dplyr::distinct()
      if (!nrow(candidates)) return(NULL)

      signatures <- sort(unique(candidates$taxonomy_signature))
      n_ids <- dplyr::n_distinct(candidates$accepted_id)
      n_taxonomies <- length(signatures)

      if (n_taxonomies != 1L) {
        return(list(
          resolved = FALSE,
          resolution_status = paste0(prefix, "_conflict"),
          n_candidate_rows = nrow(candidates),
          n_candidate_accepted_ids = n_ids,
          n_candidate_taxonomies = n_taxonomies
        ))
      }

      if (!is.null(source_preference)) {
        preferred <- candidates[candidates$source_status == source_preference, , drop = FALSE]
        if (nrow(preferred)) candidates <- preferred
      }

      candidates <- as_base_df(
        candidates %>%
          dplyr::arrange(
            dplyr::desc(source_status == "accepted"),
            accepted_id,
            source_id
          )
      )
      selected <- candidates[1, , drop = FALSE]
      accepted_ids <- collapse_sorted_chr(candidates$accepted_id)
      status_suffix <- if (n_ids == 1L) "unique" else "consensus_multiple_ids"

      list(
        resolved = TRUE,
        accepted_name = selected$accepted_name[1],
        accepted_id = selected$accepted_id[1],
        accepted_ids = accepted_ids,
        accepted_family = selected$accepted_family[1],
        accepted_genus = selected$accepted_genus[1],
        accepted_rank = selected$accepted_rank[1],
        source_status = selected$source_status[1],
        tax_status = ifelse(selected$source_status[1] == "accepted", "accepted", "synonym"),
        match_rule = prefix,
        resolution_status = paste0(prefix, "_", status_suffix),
        n_candidate_rows = nrow(candidates),
        n_candidate_accepted_ids = n_ids,
        n_candidate_taxonomies = n_taxonomies
      )
    }

    result <- NULL

    if (!is.null(exact_group) && nrow(exact_group)) {
      exact_accepted <- exact_group[
        exact_group$source_status == "accepted",
        ,
        drop = FALSE
      ]

      if (nrow(exact_accepted)) {
        result <- resolve_consensus(
          exact_accepted,
          prefix = "exact_accepted",
          source_preference = "accepted"
        )
      } else {
        result <- resolve_consensus(
          exact_group,
          prefix = "exact_synonym",
          source_preference = "synonym"
        )
      }
    } else if (!is.null(canonical_group) && nrow(canonical_group)) {
      result <- resolve_consensus(
        canonical_group,
        prefix = "canonical",
        source_preference = NULL
      )
    }

    if (is.null(result)) result <- empty_result

    if (!isTRUE(result$resolved %||% FALSE)) {
      empty_result$resolution_status <- result$resolution_status %||% "no_match"
      empty_result$match_rule <- ifelse(
        grepl("^exact", empty_result$resolution_status),
        "exact_conflict",
        ifelse(grepl("^canonical", empty_result$resolution_status), "canonical_conflict", "unresolved")
      )
      empty_result$n_candidate_rows <- result$n_candidate_rows %||% 0L
      empty_result$n_candidate_accepted_ids <- result$n_candidate_accepted_ids %||% 0L
      empty_result$n_candidate_taxonomies <- result$n_candidate_taxonomies %||% 0L
      result <- empty_result
    }

    data.frame(
      name_id = target_row$name_id,
      original_name = target_row$original_name,
      exact_name = target_row$exact_name,
      key = target_row$key,
      accepted_name = result$accepted_name,
      accepted_id = result$accepted_id,
      accepted_ids = result$accepted_ids,
      accepted_family = result$accepted_family,
      accepted_genus = result$accepted_genus,
      accepted_rank = result$accepted_rank,
      source_status = result$source_status,
      tax_status = result$tax_status,
      match_rule = result$match_rule,
      resolution_status = result$resolution_status,
      n_candidate_rows = as.integer(result$n_candidate_rows),
      n_candidate_accepted_ids = as.integer(result$n_candidate_accepted_ids),
      n_candidate_taxonomies = as.integer(result$n_candidate_taxonomies),
      stringsAsFactors = FALSE
    )
  }

  resolution_rows <- lapply(
    seq_len(nrow(target_names)),
    function(i) make_resolution_row(target_names[i, , drop = FALSE])
  )

  crosswalk_wfo <- as_base_df(
    dplyr::bind_rows(resolution_rows) %>%
      dplyr::mutate(
        tax_provider = ifelse(
          !is.na(accepted_name),
          "WFO_offline_hierarchical",
          NA_character_
        ),
        tax_checked_at = TODAY_STR
      ) %>%
      dplyr::select(
        original_name, exact_name, key,
        accepted_name, accepted_id, accepted_ids,
        accepted_family, accepted_genus, accepted_rank,
        source_status, tax_status, match_rule, resolution_status,
        n_candidate_rows, n_candidate_accepted_ids,
        n_candidate_taxonomies,
        tax_provider, tax_checked_at, name_id
      )
  )

  legacy_key_summary <- as_base_df(
    dict_candidates %>%
      dplyr::filter(!is.na(candidate_key), nzchar(candidate_key)) %>%
      dplyr::group_by(key = candidate_key) %>%
      dplyr::summarise(
        legacy_n_accepted_ids = dplyr::n_distinct(accepted_id),
        .groups = "drop"
      )
  )

  wfo_resolution_comparison <- as_base_df(
    crosswalk_wfo %>%
      dplyr::left_join(legacy_key_summary, by = "key") %>%
      dplyr::mutate(
        legacy_n_accepted_ids = dplyr::coalesce(legacy_n_accepted_ids, 0L),
        legacy_resolved = legacy_n_accepted_ids == 1L,
        hierarchical_resolved = !is.na(accepted_name),
        comparison_status = dplyr::case_when(
          legacy_resolved & hierarchical_resolved ~ "resolved_by_both",
          !legacy_resolved & hierarchical_resolved ~ "recovered_by_hierarchical_rules",
          legacy_resolved & !hierarchical_resolved ~ "lost_relative_to_legacy",
          TRUE ~ "unresolved_by_both"
        )
      ) %>%
      dplyr::select(
        original_name, exact_name, key,
        legacy_n_accepted_ids, legacy_resolved,
        hierarchical_resolved, comparison_status,
        match_rule, resolution_status,
        accepted_name, accepted_id, accepted_family,
        accepted_genus, accepted_rank
      )
  )

  wfo_resolution_rule_summary <- as_base_df(
    crosswalk_wfo %>%
      dplyr::count(
        match_rule,
        resolution_status,
        tax_status,
        name = "n_original_names",
        sort = TRUE
      )
  )

  exact_conflict_ids <- crosswalk_wfo$name_id[
    grepl("^exact_.*_conflict$", crosswalk_wfo$resolution_status)
  ]
  canonical_conflict_ids <- crosswalk_wfo$name_id[
    crosswalk_wfo$resolution_status == "canonical_conflict"
  ]

  conflict_exact_candidates <- if (length(exact_conflict_ids)) {
    exact_hits %>%
      dplyr::filter(name_id %in% exact_conflict_ids) %>%
      dplyr::mutate(candidate_match_stage = "exact")
  } else {
    data.frame(stringsAsFactors = FALSE)
  }

  conflict_canonical_candidates <- if (length(canonical_conflict_ids)) {
    canonical_hits %>%
      dplyr::filter(name_id %in% canonical_conflict_ids) %>%
      dplyr::mutate(candidate_match_stage = "canonical")
  } else {
    data.frame(stringsAsFactors = FALSE)
  }

  wfo_ambiguous_keys <- as_base_df(
    dplyr::bind_rows(conflict_exact_candidates, conflict_canonical_candidates) %>%
      dplyr::select(
        dplyr::any_of(c(
          "name_id", "original_name", "exact_name", "key",
          "candidate_match_stage", "source_name", "source_id",
          "source_status", "accepted_id", "accepted_name",
          "accepted_family", "accepted_genus", "accepted_rank",
          "taxonomy_signature"
        ))
      ) %>%
      dplyr::distinct() %>%
      dplyr::arrange(original_name, candidate_match_stage, accepted_name, accepted_id)
  )
  wfo_resolution_candidates <- wfo_ambiguous_keys

  if (isTRUE(cfg$verbose %||% TRUE)) {
    n_recovered <- sum(
      wfo_resolution_comparison$comparison_status == "recovered_by_hierarchical_rules",
      na.rm = TRUE
    )
    cat("[WFO] Original names evaluated:", nrow(crosswalk_wfo), "\n")
    cat("[WFO] Names resolved by hierarchical rules:", sum(!is.na(crosswalk_wfo$accepted_name)), "\n")
    cat("[WFO] Names recovered relative to the legacy key-only rule:", n_recovered, "\n")
    cat("[WFO] Names remaining unresolved:", sum(is.na(crosswalk_wfo$accepted_name)), "\n")
  }

  lin_wfo <- as_base_df(
    lin_enriched_pre_wfo %>%
      dplyr::mutate(
        .original_name = build_original_scientific_name(genus, species)
      ) %>%
      dplyr::left_join(
        crosswalk_wfo,
        by = c(".original_name" = "original_name")
      )
  )

  for (col in c("genus", "species", "family")) {
    if (!col %in% names(lin_wfo)) lin_wfo[[col]] <- NA_character_
  }

  compound_resolution_status <- as_base_df(
    lin_wfo %>%
      dplyr::group_by(inchikey) %>%
      dplyr::summarise(
        compound_has_resolved_occurrence = any(
          !is.na(accepted_name) & nzchar(accepted_name)
        ),
        .groups = "drop"
      )
  )

  unresolved_occurrence_audit <- as_base_df(
    lin_wfo %>%
      dplyr::filter(is.na(accepted_name) | !nzchar(accepted_name)) %>%
      dplyr::left_join(compound_resolution_status, by = "inchikey")
  )

  wfo_unresolved_names <- as_base_df(
    unresolved_occurrence_audit %>%
      dplyr::group_by(
        original_name = .original_name,
        exact_name,
        key,
        match_rule,
        resolution_status,
        n_candidate_rows,
        n_candidate_accepted_ids,
        n_candidate_taxonomies
      ) %>%
      dplyr::summarise(
        n_occurrence_rows = dplyr::n(),
        n_compounds = dplyr::n_distinct(inchikey),
        n_compounds_exclusive_to_unresolved = dplyr::n_distinct(
          inchikey[!compound_has_resolved_occurrence]
        ),
        n_families_original = dplyr::n_distinct(family),
        n_genera_original = dplyr::n_distinct(genus),
        n_references = dplyr::n_distinct(ref_id),
        original_families = collapse_sorted_chr(family),
        original_genera = collapse_sorted_chr(genus),
        .groups = "drop"
      ) %>%
      dplyr::arrange(
        dplyr::desc(n_compounds_exclusive_to_unresolved),
        dplyr::desc(n_compounds),
        original_name
      )
  )

  wfo_unresolved_priority <- as_base_df(
    wfo_unresolved_names %>%
      dplyr::filter(n_compounds_exclusive_to_unresolved > 0) %>%
      dplyr::arrange(
        dplyr::desc(n_compounds_exclusive_to_unresolved),
        dplyr::desc(n_occurrence_rows),
        original_name
      )
  )

  lin_wfo <- as_base_df(
    lin_wfo %>%
      dplyr::mutate(
        genus = as.character(genus),
        species = as.character(species),
        family = as.character(family),
        original_name = .original_name,
        .corrected_scientific = dplyr::coalesce(accepted_name, .original_name),
        .corrected_scientific = stringr::str_squish(.corrected_scientific)
      )
  )

  name_parts <- split_scientific_name(lin_wfo$.corrected_scientific)
  lin_applied <- as_base_df(cbind(lin_wfo, name_parts))

  lin_applied <- as_base_df(
    lin_applied %>%
      dplyr::mutate(
        species_original = species,
        genus_original = genus,
        family_original = family,
        taxonomy_action = dplyr::case_when(
          is.na(accepted_name) ~ "unresolved",
          match_rule == "exact_synonym" ~ "synonym_replaced",
          match_rule == "canonical" & tax_status == "synonym" ~ "synonym_canonical_match",
          match_rule == "canonical" ~ "accepted_canonical_match",
          TRUE ~ "accepted_confirmed"
        ),
        species = dplyr::coalesce(accepted_name, species),
        genus = dplyr::coalesce(accepted_genus, corrected_genus, genus),
        family = dplyr::coalesce(accepted_family, family),
        tax_checked_at = TODAY_STR
      )
  )

  genus_final_norm <- tolower(tidy_space(lin_applied$genus))
  target_genera_norm <- tolower(tidy_space(cfg$taxon_values %||% character(0)))

  genus_scope_ok <- rep(TRUE, nrow(lin_applied))
  if (identical(TAXON_MODE, "genus")) {
    genus_scope_ok <- !is.na(genus_final_norm) & genus_final_norm %in% target_genera_norm
  }

  resolved_ok <- !is.na(lin_applied$accepted_name) & nzchar(lin_applied$accepted_name)
  family_available <- !is.na(lin_applied$family) & nzchar(trimws(lin_applied$family))

  policy_ok <- switch(
    WFO_CORE_POLICY,
    resolved_only = resolved_ok,
    family_available = resolved_ok | family_available,
    all = rep(TRUE, nrow(lin_applied))
  )

  lin_enriched_all <- as_base_df(
    lin_applied %>%
      dplyr::mutate(
        taxonomy_core_eligible = policy_ok & genus_scope_ok,
        taxonomy_exclusion_reason = dplyr::case_when(
          !genus_scope_ok ~ "outside_requested_genus_after_WFO",
          !policy_ok & taxonomy_action == "unresolved" ~ "unresolved_name",
          !policy_ok ~ "excluded_by_wfo_policy",
          TRUE ~ NA_character_
        )
      )
  )

  lin_enriched <- as_base_df(
    lin_enriched_all %>%
      dplyr::filter(taxonomy_core_eligible)
  )

  if (rank_available) {
    wfo_family_species_totals <- as_base_df(
      accepted %>%
        dplyr::filter(
          accepted_rank == "species",
          !is.na(accepted_family), nzchar(accepted_family)
        ) %>%
        dplyr::group_by(family = accepted_family) %>%
        dplyr::summarise(
          n_accepted_species_wfo = dplyr::n_distinct(accepted_id),
          .groups = "drop"
        )
    )
  }

  cat(sprintf(
    "[WFO] Hierarchical normalization completed. Core policy='%s': %d of %d rows retained.\n",
    WFO_CORE_POLICY, nrow(lin_enriched), nrow(lin_enriched_all)
  ))
} else {
  lin_enriched_all <- as_base_df(
    lin_enriched_pre_wfo %>%
      dplyr::mutate(
        species_original = species,
        genus_original = genus,
        family_original = family,
        original_name = build_original_scientific_name(genus, species),
        accepted_name = NA_character_,
        accepted_id = NA_character_,
        accepted_ids = NA_character_,
        accepted_family = NA_character_,
        accepted_genus = NA_character_,
        accepted_rank = NA_character_,
        source_status = NA_character_,
        tax_status = "not_checked",
        match_rule = "not_checked",
        resolution_status = "not_checked",
        n_candidate_rows = NA_integer_,
        n_candidate_accepted_ids = NA_integer_,
        n_candidate_taxonomies = NA_integer_,
        tax_provider = NA_character_,
        tax_checked_at = NA_character_,
        taxonomy_action = "not_checked",
        taxonomy_core_eligible = TRUE,
        taxonomy_exclusion_reason = NA_character_
      )
  )
  lin_enriched <- lin_enriched_all
}

clean_taxonomic_fields <- function(df) {
  df <- as_base_df(df)
  df <- as_base_df(
    df %>%
      dplyr::mutate(
        genus = canon_genus(genus),
        species = as.character(species),
        species = gsub("\\bsp\\.?\\b", "", species, ignore.case = TRUE),
        species = stringr::str_squish(species),
        species = fix_glued_species(genus, species),
        species = dplyr::na_if(species, ""),
        family = stringr::str_squish(as.character(family)),
        family = dplyr::na_if(family, "")
      )
  )
  df
}

lin_enriched_all <- clean_taxonomic_fields(lin_enriched_all)
lin_enriched <- clean_taxonomic_fields(lin_enriched)

if (nrow(lin_enriched) == 0L) {
  stop(
    "The WFO/core filtering policy retained zero rows. ",
    "Inspect PartI_ALL/wfo_taxonomy_crosswalk.csv or use a less restrictive cfg$wfo_core_policy."
  )
}

lin_enriched_all_pre_dedup <- as_base_df(lin_enriched_all)
lin_enriched_pre_dedup <- as_base_df(lin_enriched)
taxonomy_excluded_pre_dedup <- as_base_df(
  lin_enriched_all_pre_dedup %>%
    dplyr::filter(!taxonomy_core_eligible)
)

deduplicate_occurrences <- function(df) {
  as_base_df(
    df %>%
      dplyr::arrange(inchikey, family, genus, species, ref_id, source) %>%
      dplyr::distinct(inchikey, family, genus, species, ref_id, .keep_all = TRUE)
  )
}

dedup_before_core <- nrow(lin_enriched_pre_dedup)
lin_enriched <- deduplicate_occurrences(lin_enriched_pre_dedup)
dedup_after_core <- nrow(lin_enriched)

taxonomy_excluded <- deduplicate_occurrences(taxonomy_excluded_pre_dedup)

lin_enriched_all <- as_base_df(
  dplyr::bind_rows(lin_enriched, taxonomy_excluded) %>%
    dplyr::arrange(
      dplyr::desc(taxonomy_core_eligible),
      inchikey, family, genus, species, ref_id, source
    )
)

cat(sprintf(
  "[DEDUP] Core compound x species x reference: %d -> %d rows (%d duplicates removed).\n",
  dedup_before_core,
  dedup_after_core,
  dedup_before_core - dedup_after_core
))
cat(sprintf(
  "[DEDUP] Excluded-taxonomy records: %d -> %d rows (%d duplicates removed).\n",
  nrow(taxonomy_excluded_pre_dedup),
  nrow(taxonomy_excluded),
  nrow(taxonomy_excluded_pre_dedup) - nrow(taxonomy_excluded)
))

build_map_tax_inchi <- function(lin_tbl, tax_col = c("family", "genus", "species")) {
  tax_col <- match.arg(tax_col)

  out <- as_base_df(
    lin_tbl %>%
      dplyr::transmute(
        inchikey = as.character(inchikey),
        family = as.character(family),
        genus = canon_genus(genus),
        species = as.character(species)
      ) %>%
      dplyr::filter(!is.na(inchikey), nzchar(inchikey))
  )

  if (identical(tax_col, "genus")) {
    out <- as_base_df(out %>% dplyr::filter(!is.na(genus), nzchar(genus)))
  } else if (identical(tax_col, "species")) {
    out <- as_base_df(out %>% dplyr::filter(is_binomial(species)))
  } else {
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

tax_col_for_map <- cfg$analysis_tax_level %||% "family"
map_tax_inchi <- as_base_df(build_map_tax_inchi(lin_enriched, tax_col = tax_col_for_map))
cat(sprintf(
  "[INFO] Taxon-to-InChIKey map generated at the %s level. Rows: %d\n",
  tax_col_for_map, nrow(map_tax_inchi)
))

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
      dplyr::across(dplyr::any_of(NUMERIC_PROPS), ~ suppressWarnings(as.numeric(.))),
      dplyr::across(dplyr::any_of(LOGICAL_PROPS), ~ {
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
      dplyr::filter(!is.na(species), nzchar(species)) %>%
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

build_family_coverage <- function(df, dataset_label) {
  if (!"murko_framework" %in% names(df)) df$murko_framework <- NA_character_
  if (!"taxonomy_action" %in% names(df)) df$taxonomy_action <- NA_character_

  out <- as_base_df(
    df %>%
      dplyr::filter(!is.na(family), nzchar(family)) %>%
      dplyr::group_by(family) %>%
      dplyr::summarise(
        n_occurrence_rows = dplyr::n(),
        n_compounds = count_distinct_nonempty(inchikey),
        n_scaffolds_murcko = count_distinct_nonempty(murko_framework),
        n_compounds_with_scaffold = count_distinct_nonempty(
          inchikey[!is.na(murko_framework) & nzchar(as.character(murko_framework))]
        ),
        n_genera = count_distinct_nonempty(genus),
        n_species = count_distinct_nonempty(species[is_binomial(species)]),
        n_references = count_distinct_nonempty(ref_id),
        n_sources = count_distinct_nonempty(source),
        n_unresolved_rows = sum(taxonomy_action == "unresolved", na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        dataset = dataset_label,
        scaffold_per_compound = dplyr::if_else(
          n_compounds > 0,
          n_scaffolds_murcko / n_compounds,
          NA_real_
        ),
        pct_compounds_with_scaffold = dplyr::if_else(
          n_compounds > 0,
          100 * n_compounds_with_scaffold / n_compounds,
          NA_real_
        ),
        compounds_per_species = dplyr::if_else(
          n_species > 0,
          n_compounds / n_species,
          NA_real_
        ),
        references_per_compound = dplyr::if_else(
          n_compounds > 0,
          n_references / n_compounds,
          NA_real_
        )
      )
  )

  if (nrow(wfo_family_species_totals) > 0) {
    out <- as_base_df(
      out %>%
        dplyr::left_join(wfo_family_species_totals, by = "family") %>%
        dplyr::mutate(
          pct_wfo_species_represented = dplyr::if_else(
            !is.na(n_accepted_species_wfo) & n_accepted_species_wfo > 0,
            100 * n_species / n_accepted_species_wfo,
            NA_real_
          )
        )
    )
  } else {
    out$n_accepted_species_wfo <- NA_integer_
    out$pct_wfo_species_represented <- NA_real_
  }

  thresholds <- sort(unique(as.integer(cfg$analysis_compound_thresholds %||% c(1L, 5L, 10L, 20L))))
  thresholds <- thresholds[is.finite(thresholds) & thresholds >= 1L]
  for (thr in thresholds) {
    out[[paste0("ncomp_ge_", thr)]] <- out$n_compounds >= thr
  }

  min_comp <- as.integer(cfg$analysis_min_compounds_per_taxon %||% 5L)
  min_species <- as.integer(cfg$analysis_min_species_per_taxon %||% 1L)
  out$eligible_primary <- out$n_compounds >= min_comp & out$n_species >= min_species

  as_base_df(out %>% dplyr::arrange(dplyr::desc(n_compounds), family))
}

family_coverage <- build_family_coverage(lin_enriched, "core")
family_coverage_all <- build_family_coverage(lin_enriched_all, "all_taxonomic_records")

build_threshold_summary <- function(coverage_tbl, occurrence_tbl) {
  thresholds <- sort(unique(as.integer(cfg$analysis_compound_thresholds %||% c(1L, 5L, 10L, 20L))))
  thresholds <- thresholds[is.finite(thresholds) & thresholds >= 1L]
  total_families <- nrow(coverage_tbl)
  total_compounds <- count_distinct_nonempty(occurrence_tbl$inchikey)

  as_base_df(dplyr::bind_rows(lapply(thresholds, function(thr) {
    retained_families <- coverage_tbl$family[coverage_tbl$n_compounds >= thr]
    retained_compounds <- occurrence_tbl %>%
      dplyr::filter(family %in% retained_families) %>%
      dplyr::pull(inchikey)

    data.frame(
      min_compounds = thr,
      n_families_retained = length(retained_families),
      pct_families_retained = if (total_families > 0) 100 * length(retained_families) / total_families else NA_real_,
      n_compounds_represented = count_distinct_nonempty(retained_compounds),
      pct_compounds_represented = if (total_compounds > 0) 100 * count_distinct_nonempty(retained_compounds) / total_compounds else NA_real_,
      stringsAsFactors = FALSE
    )
  })))
}

coverage_threshold_summary <- build_threshold_summary(family_coverage, lin_enriched)

taxonomy_status_summary <- as_base_df(
  lin_enriched_all %>%
    dplyr::group_by(taxonomy_action, taxonomy_core_eligible, taxonomy_exclusion_reason) %>%
    dplyr::summarise(
      n_rows = dplyr::n(),
      n_compounds = count_distinct_nonempty(inchikey),
      n_families = count_distinct_nonempty(family),
      n_genera = count_distinct_nonempty(genus),
      n_species = count_distinct_nonempty(species),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(n_rows))
)

taxonomic_rank_summary <- as_base_df(
  lin_enriched %>%
    dplyr::mutate(
      accepted_rank = dplyr::coalesce(as.character(accepted_rank), "unknown"),
      valid_binomial = is_binomial(species)
    ) %>%
    dplyr::group_by(accepted_rank, valid_binomial) %>%
    dplyr::summarise(
      n_rows = dplyr::n(),
      n_compounds = count_distinct_nonempty(inchikey),
      n_families = count_distinct_nonempty(family),
      n_genera = count_distinct_nonempty(genus),
      n_species_labels = count_distinct_nonempty(species),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(n_rows))
)

non_binomial_core_records <- as_base_df(
  lin_enriched %>%
    dplyr::filter(
      !is.na(species),
      nzchar(as.character(species)),
      !is_binomial(species)
    ) %>%
    dplyr::select(
      dplyr::any_of(c(
        "inchikey", "lotus_id", "family", "genus", "species",
        "accepted_name", "accepted_rank", "taxonomy_action",
        "ref_id", "source"
      ))
    ) %>%
    dplyr::arrange(family, genus, species, inchikey, ref_id)
)

quality_metrics <- list(
  missing_inchikey = sum(is.na(lin_enriched$inchikey) | !nzchar(as.character(lin_enriched$inchikey))),
  missing_smiles = sum(is.na(lin_enriched$smiles) | !nzchar(as.character(lin_enriched$smiles))),
  missing_family = sum(is.na(lin_enriched$family) | !nzchar(as.character(lin_enriched$family))),
  missing_genus = sum(is.na(lin_enriched$genus) | !nzchar(as.character(lin_enriched$genus))),
  missing_species = sum(is.na(lin_enriched$species) | !nzchar(as.character(lin_enriched$species))),
  non_binomial_species = sum(!is.na(lin_enriched$species) & !is_binomial(lin_enriched$species)),
  missing_murcko_framework = if ("murko_framework" %in% names(lin_enriched)) {
    sum(is.na(lin_enriched$murko_framework) | !nzchar(as.character(lin_enriched$murko_framework)))
  } else {
    nrow(lin_enriched)
  },
  missing_npclassifier_class = if ("chemicalTaxonomyNPclassifierClass" %in% names(lin_enriched)) {
    sum(is.na(lin_enriched$chemicalTaxonomyNPclassifierClass) |
          !nzchar(as.character(lin_enriched$chemicalTaxonomyNPclassifierClass)))
  } else {
    nrow(lin_enriched)
  }
)

data_quality_summary <- data.frame(
  metric = names(quality_metrics),
  n_rows = as.numeric(unlist(quality_metrics, use.names = FALSE)),
  denominator_rows = nrow(lin_enriched),
  percent = if (nrow(lin_enriched) > 0) {
    100 * as.numeric(unlist(quality_metrics, use.names = FALSE)) / nrow(lin_enriched)
  } else {
    NA_real_
  },
  stringsAsFactors = FALSE
)

pipeline_audit <- as_base_df(dplyr::bind_rows(
  snapshot_counts(
    "01_mongodb_matched",
    notes = "Rows matched before local distinct filtering.",
    n_rows_override = total_lines
  ),
  snapshot_counts(
    "02_extracted_distinct",
    lin,
    notes = "Distinct compound x reference x source x taxon rows after extraction."
  ),
  snapshot_counts(
    "03_enriched_pre_wfo",
    lin_enriched_pre_wfo,
    notes = "Rows after joining chemical properties, before WFO normalization."
  ),
  snapshot_counts(
    "04_wfo_all_before_dedup",
    lin_enriched_all_pre_dedup,
    notes = paste0("All rows after WFO processing, before duplicate removal; core policy=", WFO_CORE_POLICY, ".")
  ),
  snapshot_counts(
    "05_wfo_core_eligible_before_dedup",
    lin_enriched_pre_dedup,
    notes = "Rows eligible for the primary dataset before duplicate removal."
  ),
  snapshot_counts(
    "06_wfo_excluded_before_dedup",
    taxonomy_excluded_pre_dedup,
    notes = "Rows excluded by the taxonomic policy before duplicate removal."
  ),
  snapshot_counts(
    "07_core_after_dedup",
    lin_enriched,
    notes = "Primary occurrence-level analytical dataset after duplicate removal."
  ),
  snapshot_counts(
    "08_taxonomy_excluded_after_dedup",
    taxonomy_excluded,
    notes = "Taxonomically excluded records retained for sensitivity/audit after duplicate removal."
  ),
  snapshot_counts(
    "09_all_taxonomic_records_after_dedup",
    lin_enriched_all,
    notes = "Union of separately deduplicated core and excluded taxonomic records."
  ),
  snapshot_counts(
    "10_compound_species",
    lin_cs,
    notes = "One row per compound x species."
  ),
  snapshot_counts(
    "11_unique_compounds",
    uni_enriched,
    notes = "One row per InChIKey; semicolon-delimited taxonomic fields are tokenized for audit counts."
  )
))

readr::write_csv(family_coverage, file.path(PART1_DIR, "family_coverage_core.csv"))
readr::write_csv(family_coverage_all, file.path(PART1_DIR, "family_coverage_all_taxonomic_records.csv"))
readr::write_csv(coverage_threshold_summary, file.path(PART1_DIR, "coverage_threshold_summary.csv"))
readr::write_csv(taxonomy_status_summary, file.path(PART1_DIR, "taxonomy_status_summary.csv"))
readr::write_csv(taxonomic_rank_summary, file.path(PART1_DIR, "taxonomic_rank_summary.csv"))
readr::write_csv(non_binomial_core_records, file.path(PART1_DIR, "non_binomial_core_records.csv"))
readr::write_csv(crosswalk_wfo, file.path(PART1_DIR, "wfo_taxonomy_crosswalk.csv"))
readr::write_csv(taxonomy_excluded, file.path(PART1_DIR, "taxonomy_excluded_records.csv"))
readr::write_csv(property_match_summary, file.path(PART1_DIR, "property_match_summary.csv"))
readr::write_csv(data_quality_summary, file.path(PART1_DIR, "data_quality_summary.csv"))
readr::write_csv(pipeline_audit, file.path(PART1_DIR, "pipeline_audit.csv"))
readr::write_csv(wfo_resolution_rule_summary, file.path(PART1_DIR, "wfo_resolution_rule_summary.csv"))
readr::write_csv(wfo_resolution_comparison, file.path(PART1_DIR, "wfo_resolution_legacy_vs_hierarchical.csv"))
readr::write_csv(wfo_unresolved_names, file.path(PART1_DIR, "wfo_unresolved_names.csv"))
readr::write_csv(wfo_unresolved_priority, file.path(PART1_DIR, "wfo_unresolved_priority_by_chemical_loss.csv"))
if (nrow(wfo_ambiguous_keys) > 0) {
  readr::write_csv(wfo_ambiguous_keys, file.path(PART1_DIR, "wfo_ambiguous_keys.csv"))
  readr::write_csv(wfo_resolution_candidates, file.path(PART1_DIR, "wfo_unresolved_conflict_candidates.csv"))
}

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
      uni_enriched = uni_enriched_x,
      family_coverage = sanitize_for_excel(family_coverage),
      family_cov_all = sanitize_for_excel(family_coverage_all),
      coverage_thresholds = sanitize_for_excel(coverage_threshold_summary),
      taxonomy_status = sanitize_for_excel(taxonomy_status_summary),
      taxonomic_ranks = sanitize_for_excel(taxonomic_rank_summary),
      non_binomial_core = sanitize_for_excel(non_binomial_core_records),
      pipeline_audit = sanitize_for_excel(pipeline_audit),
      data_quality = sanitize_for_excel(data_quality_summary),
      property_match = sanitize_for_excel(property_match_summary),
      wfo_resolution_rules = sanitize_for_excel(wfo_resolution_rule_summary),
      wfo_legacy_comparison = sanitize_for_excel(wfo_resolution_comparison),
      wfo_unresolved_names = sanitize_for_excel(wfo_unresolved_names),
      wfo_unresolved_priority = sanitize_for_excel(wfo_unresolved_priority)
    ),
    path = xlsx_path
  )
  cat("[OK] Excel workbook written to:", normalizePath(xlsx_path), "\n")
}

if (EXPORT_PARQUET) {
  if (!requireNamespace("arrow", quietly = TRUE) && isTRUE(cfg$auto_install_missing_packages %||% FALSE)) {
    install.packages("arrow")
  }

  if (requireNamespace("arrow", quietly = TRUE)) {
    arrow::write_parquet(as_base_df(lin_enriched), safe_file(paste0(base_tag, "_lin_enriched"), ".parquet"))
    arrow::write_parquet(as_base_df(uni_enriched), safe_file(paste0(base_tag, "_uni_enriched"), ".parquet"))
    arrow::write_parquet(as_base_df(lin_cs), safe_file(paste0(base_tag, "_lin_compound_species"), ".parquet"))
    arrow::write_parquet(as_base_df(lin_enriched_all), file.path(PART1_DIR, "lin_enriched_all_taxonomic_records.parquet"))
    arrow::write_parquet(as_base_df(taxonomy_excluded), file.path(PART1_DIR, "taxonomy_excluded_records.parquet"))
    arrow::write_parquet(as_base_df(family_coverage), file.path(PART1_DIR, "family_coverage_core.parquet"))
    cat("[OK] Parquet files written to:", normalizePath(OUT_DIR), "\n")
  } else {
    warning("Package 'arrow' is not available; Parquet export was skipped.")
  }
}

cat("-----------------------------------------------------\n")
cat("[PART I SUMMARY]\n")
cat(" Core occurrence rows: ", nrow(lin_enriched), "\n", sep = "")
cat(" Unique compounds: ", count_distinct_nonempty(lin_enriched$inchikey), "\n", sep = "")
cat(" Families: ", count_distinct_nonempty(lin_enriched$family), "\n", sep = "")
cat(" Species: ", count_distinct_nonempty(lin_enriched$species[is_binomial(lin_enriched$species)]), "\n", sep = "")
cat(" Excluded taxonomy rows retained for audit: ", nrow(taxonomy_excluded), "\n", sep = "")
cat(" Audit directory: ", normalizePath(PART1_DIR), "\n", sep = "")
cat("-----------------------------------------------------\n")
