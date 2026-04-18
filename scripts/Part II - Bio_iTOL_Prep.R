# Part II: bioactivity integration, iTOL preparation, RDKit export, and atlas generation.

`%||%` <- function(a, b) if (is.null(a)) b else a

# Keep named analytical tables as base data.frame objects.
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

pkgs <- c(
  "dplyr", "tidyr", "readr", "stringr", "readxl", "writexl",
  "httr", "jsonlite", "progress", "mongolite",
  "ggplot2", "grid", "gridExtra", "scales"
)
require_pkgs(pkgs)

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
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(scales)
})

if (!exists("cfg")) stop("Part II requires object 'cfg' created by the main pipeline.")
if (!exists("runtime")) stop("Part II requires object 'runtime' created by the main pipeline.")

OUT_DIR <- if (exists("OUT_DIR")) OUT_DIR else runtime$OUT_DIR
base_tag <- if (exists("base_tag")) base_tag else runtime$base_tag

if (is.null(OUT_DIR) || !dir.exists(OUT_DIR)) stop("OUT_DIR does not exist: ", OUT_DIR)
if (is.null(base_tag) || !nzchar(base_tag)) stop("base_tag is missing or empty.")

if (exists("runtime")) OUT_DIR <- runtime$OUT_DIR
stopifnot(dir.exists(OUT_DIR))

out2_dir <- file.path(OUT_DIR, "PartII_ALL")
dir.create(out2_dir, showWarnings = FALSE, recursive = TRUE)

out2_path <- function(x) file.path(out2_dir, x)

path_tax_excel <- file.path(OUT_DIR, paste0(base_tag, ".xlsx"))
if (!file.exists(path_tax_excel)) {
  message("[WARNING] LOTUS Excel workbook not found at: ", path_tax_excel)
  message("This is acceptable when 'lin_enriched' and 'uni_enriched' are already available in memory.")
}

message("----------------------------------------------------------------")
message(">>> STEP A: RDKit input CSV (flavonoid subset)")
message("----------------------------------------------------------------")

path_out_csv <- file.path(out2_dir, paste0(base_tag, "__flavonoids_for_rdkit.csv"))

if (exists("lin_enriched")) {
  lotus_raw <- as_base_df(lin_enriched)
} else {
  if (!file.exists(path_tax_excel)) {
    stop("No 'lin_enriched' object is available in memory and the LOTUS Excel workbook was not found: ", path_tax_excel)
  }
  sheets <- readxl::excel_sheets(path_tax_excel)
  if (!"lin_enriched" %in% sheets) {
    stop("Sheet 'lin_enriched' was not found in: ", path_tax_excel)
  }
  lotus_raw <- as_base_df(
    readxl::read_excel(path_tax_excel, sheet = "lin_enriched", .name_repair = "minimal")
  )
}

req_cols <- c("inchikey", "smiles", "chemicalTaxonomyNPclassifierClass", "family", "genus", "species")
missing_cols <- setdiff(req_cols, names(lotus_raw))
if (length(missing_cols) > 0) {
  stop("The LOTUS source table is missing required columns: ", paste(missing_cols, collapse = ", "))
}

lotus_raw <- as_base_df(
  lotus_raw %>%
    dplyr::rename(class_np = chemicalTaxonomyNPclassifierClass) %>%
    dplyr::mutate(
      class_np = stringr::str_trim(class_np),
      class_np = dplyr::recode(class_np, "Flavandiols (Leucoanthocyanidins)" = "Flavandiols")
    )
)

flav_classes <- cfg$flav_classes %||% c(
  "Flavanones", "Isoflavones", "Flavan-3-ols", "Chalcones", "Dihydroflavonols",
  "Pterocarpan", "Proanthocyanins", "Isoflavanones", "2-arylbenzofurans",
  "Flavandiols", "Rotenoids", "Aurones", "Flavans", "Coumestan", "Anthocyanidins",
  "Neoflavonoids", "Flavonolignans"
)

flav_only <- as_base_df(
  lotus_raw %>%
    dplyr::filter(!is.na(class_np), class_np %in% flav_classes)
)

if (nrow(flav_only) == 0) {
  stop("The flavonoid filter returned zero rows. Review the class labels in 'cfg$flav_classes'.")
}

flav_for_rdkit <- as_base_df(
  flav_only %>%
    dplyr::mutate(core14 = stringr::str_sub(inchikey, 1, 14)) %>%
    dplyr::distinct(inchikey, core14, smiles, class_np, family, genus, species) %>%
    dplyr::arrange(class_np, family, genus, species, inchikey)
)

readr::write_csv(x = flav_for_rdkit, file = path_out_csv)
message("[OK] RDKit input table written: ", basename(path_out_csv))
message("Rows exported: ", nrow(flav_for_rdkit))

message("----------------------------------------------------------------")
message(">>> STEP B: Bioactivity, reference integration, and global context")
message("----------------------------------------------------------------")

if (exists("uni_enriched")) {
  uni_df <- as_base_df(uni_enriched)
} else {
  if (!file.exists(path_tax_excel)) {
    stop("The object 'uni_enriched' is not available in memory and the LOTUS Excel workbook was not found: ", path_tax_excel)
  }
  sheets <- readxl::excel_sheets(path_tax_excel)
  if (!"uni_enriched" %in% sheets) {
    stop("Sheet 'uni_enriched' was not found in: ", path_tax_excel)
  }
  uni_df <- as_base_df(
    readxl::read_excel(path_tax_excel, sheet = "uni_enriched", .name_repair = "minimal")
  )
}

if (!"inchikey" %in% names(uni_df)) stop("Column 'inchikey' is missing in 'uni_enriched'.")

if ("chemicalTaxonomyNPclassifierClass" %in% names(uni_df)) {
  uni_df <- as_base_df(
    uni_df %>%
      dplyr::mutate(hybrid_class = chemicalTaxonomyNPclassifierClass)
  )
} else {
  stop("Column 'chemicalTaxonomyNPclassifierClass' is missing in 'uni_enriched'.")
}

if (!"ref_ids" %in% names(uni_df)) {
  uni_df$ref_ids <- NA_character_
}

ACTIVITY_CUTOFF_NM <- cfg$activity_cutoff_nm %||% 10000
VALID_TYPES <- cfg$chembl_valid_types %||% c("IC50", "Ki", "EC50", "Kd", "MIC", "AC50", "GI50")

inchis_to_map <- uni_df %>%
  dplyr::filter(!is.na(inchikey), nzchar(inchikey), inchikey != "NA") %>%
  dplyr::distinct(inchikey) %>%
  dplyr::pull(inchikey)

message("[INFO] Unique compounds to query: ", length(inchis_to_map))

get_chembl_id_unichem <- function(ik) {
  url <- paste0("https://www.ebi.ac.uk/unichem/rest/inchikey/", ik)
  resp <- tryCatch(
    httr::GET(url, httr::timeout(10)),
    error = function(e) NULL
  )

  if (is.null(resp) || httr::status_code(resp) != 200) return(NULL)

  cont <- httr::content(resp, as = "text", encoding = "UTF-8")
  json <- tryCatch(jsonlite::fromJSON(cont, flatten = TRUE), error = function(e) NULL)

  if (is.data.frame(json)) {
    hit <- as_base_df(json) %>%
      dplyr::filter(as.character(src_id) == "1") %>%
      head(1)
    if (nrow(hit) > 0) {
      return(data.frame(
        inchikey = ik,
        chembl_id = hit$src_compound_id,
        stringsAsFactors = FALSE
      ))
    }
  }

  NULL
}

message("[INFO] Mapping InChIKeys through UniChem.")
res_list <- vector("list", 0)
for (i in seq_along(inchis_to_map)) {
  out <- get_chembl_id_unichem(inchis_to_map[i])
  if (!is.null(out)) res_list[[length(res_list) + 1]] <- out
  if (i %% 50 == 0) Sys.sleep(0.1)
}

if (length(res_list) > 0) {
  chembl_map <- as_base_df(
    dplyr::bind_rows(res_list) %>%
      dplyr::distinct(inchikey, chembl_id)
  )
  message("[INFO] ChEMBL identifiers mapped: ", nrow(chembl_map))
} else {
  chembl_map <- data.frame(
    inchikey = character(0),
    chembl_id = character(0),
    stringsAsFactors = FALSE
  )
  message("[INFO] No ChEMBL identifiers were retrieved. The workflow will continue without the bioactivity layer.")
}

get_chembl_bioactivity <- function(chembl_ids) {
  url <- "https://www.ebi.ac.uk/chembl/api/data/activity.json"
  ids_str <- paste(chembl_ids, collapse = ",")
  resp <- tryCatch(
    httr::GET(
      url,
      query = list(
        molecule_chembl_id__in = ids_str,
        standard_type__in = paste(VALID_TYPES, collapse = ","),
        limit = 1000
      )
    ),
    error = function(e) NULL
  )

  if (is.null(resp) || httr::status_code(resp) != 200) return(NULL)

  cont <- httr::content(resp, as = "text", encoding = "UTF-8")
  json <- tryCatch(jsonlite::fromJSON(cont, flatten = TRUE), error = function(e) NULL)
  if (is.null(json) || length(json$activities) == 0) return(NULL)

  as_base_df(
    json$activities %>%
      dplyr::select(
        molecule_chembl_id,
        standard_type,
        standard_relation,
        standard_value,
        standard_units,
        target_pref_name,
        target_organism,
        document_chembl_id
      )
  )
}

activities_list <- list()
if (nrow(chembl_map) > 0) {
  ids_vec <- unique(chembl_map$chembl_id)
  batches <- split(ids_vec, ceiling(seq_along(ids_vec) / 50))
  pb_act <- txtProgressBar(min = 0, max = length(batches), style = 3)
  for (i in seq_along(batches)) {
    setTxtProgressBar(pb_act, i)
    act_df <- get_chembl_bioactivity(batches[[i]])
    if (!is.null(act_df)) activities_list[[length(activities_list) + 1]] <- as_base_df(act_df)
    Sys.sleep(0.5)
  }
  close(pb_act)
}

full_docs_df <- data.frame(
  document_chembl_id = character(0),
  doi = character(0),
  pubmed_id = character(0),
  year = character(0),
  title = character(0),
  journal = character(0),
  stringsAsFactors = FALSE
)

if (length(activities_list) > 0) {
  raw_activities <- as_base_df(
    dplyr::bind_rows(activities_list) %>%
      dplyr::rename(chembl_id = molecule_chembl_id)
  )

  doc_ids <- unique(na.omit(raw_activities$document_chembl_id))
  doc_batches <- split(doc_ids, ceiling(seq_along(doc_ids) / 50))
  docs_list <- list()

  pb_doc <- txtProgressBar(min = 0, max = length(doc_batches), style = 3)
  for (i in seq_along(doc_batches)) {
    setTxtProgressBar(pb_doc, i)
    ids_str <- paste(doc_batches[[i]], collapse = ",")
    url_doc <- paste0(
      "https://www.ebi.ac.uk/chembl/api/data/document.json?document_chembl_id__in=",
      ids_str
    )

    resp <- tryCatch(httr::GET(url_doc), error = function(e) NULL)
    if (!is.null(resp) && httr::status_code(resp) == 200) {
      cont <- httr::content(resp, as = "text", encoding = "UTF-8")
      json <- tryCatch(jsonlite::fromJSON(cont, flatten = TRUE), error = function(e) NULL)
      if (!is.null(json) && length(json$documents) > 0) {
        docs_df <- as_base_df(
          json$documents %>%
            dplyr::select(any_of(c("document_chembl_id", "doi", "pubmed_id", "year", "title", "journal")))
        )
        docs_list[[i]] <- docs_df
      }
    }
    Sys.sleep(0.2)
  }
  close(pb_doc)

  if (length(docs_list) > 0) {
    full_docs_df <- as_base_df(dplyr::bind_rows(docs_list))
  }
}

bio_evidence_A <- data.frame(inchikey = character(0), stringsAsFactors = FALSE)
outfile_sum <- out2_path(paste0(base_tag, "_BIO_A_Summary.xlsx"))

if (exists("raw_activities") && nrow(raw_activities) > 0) {
  clean_activities <- as_base_df(
    raw_activities %>%
      dplyr::mutate(standard_value = as.numeric(standard_value)) %>%
      dplyr::filter(
        !is.na(standard_value),
        standard_units == "nM",
        standard_value <= ACTIVITY_CUTOFF_NM
      ) %>%
      dplyr::left_join(chembl_map, by = "chembl_id") %>%
      dplyr::left_join(full_docs_df, by = "document_chembl_id") %>%
      dplyr::mutate(
        ref_ids = ifelse(
          !is.na(doi),
          doi,
          ifelse(!is.na(pubmed_id), paste0("PMID:", pubmed_id), NA_character_)
        )
      )
  )

  bio_compound_summary <- as_base_df(
    clean_activities %>%
      dplyr::group_by(inchikey) %>%
      dplyr::arrange(standard_value) %>%
      dplyr::summarise(
        Best_Potency_nM = first(standard_value),
        Best_Target = first(target_pref_name),
        Best_Organism = first(target_organism),
        N_Assays = n(),
        Evidence_Flag_A = "E1 (Experimental Activity)",
        Best_Ref_DOI = first(ref_ids),
        Best_Ref_Year = first(year),
        .groups = "drop"
      )
  )

  writexl::write_xlsx(
    list(
      Summary = as_base_df(bio_compound_summary),
      Raw = as_base_df(clean_activities)
    ),
    path = outfile_sum
  )
  bio_evidence_A <- as_base_df(bio_compound_summary)
  message("[OK] Bioactivity summary written: ", basename(outfile_sum))
} else {
  message("[INFO] No qualifying bioactivity records were identified under the selected cutoff.")
}

MONGO_URL <- cfg$mongo_url %||% "mongodb://127.0.0.1:27017"
lotus_db <- mongolite::mongo(collection = "lotusUniqueNaturalProduct", db = "lotusdb", url = MONGO_URL)

target_inchikeys <- unique(na.omit(uni_df$inchikey))
message("[INFO] Querying global occurrence data for ", length(target_inchikeys), " compounds.")

extract_family_nuclear <- function(tax_obj) {
  if (is.null(tax_obj)) return(character(0))
  flat <- unlist(tax_obj)
  idx <- grep("family$", names(flat), ignore.case = TRUE)
  if (length(idx) > 0) {
    fams <- as.character(flat[idx])
    return(unique(fams[nzchar(fams) & fams != "NA"]))
  }
  character(0)
}

BATCH_SIZE <- as.integer(cfg$mongo_batch_size %||% 2000)
chunks <- split(target_inchikeys, ceiling(seq_along(target_inchikeys) / BATCH_SIZE))
results_list <- list()

pb <- progress::progress_bar$new(
  format = "   Processing [:bar] :percent | Batch :current/:total | ETA: :eta",
  total = length(chunks),
  width = 70
)

for (i in seq_along(chunks)) {
  pb$tick()
  ids_vec <- chunks[[i]]
  ids_json <- if (length(ids_vec) == 1) {
    paste0("[\"", ids_vec, "\"]")
  } else {
    jsonlite::toJSON(ids_vec, auto_unbox = TRUE)
  }
  q_json <- sprintf('{"inchikey": {"$in": %s}}', ids_json)

  tryCatch({
    b_data <- as_base_df(
      lotus_db$find(
        query = q_json,
        fields = '{"inchikey": 1, "taxonomyReferenceObjects": 1, "_id": 0}'
      )
    )

    if (nrow(b_data) > 0) {
      row_indices <- seq_len(nrow(b_data))
      families_list <- lapply(row_indices, function(idx) {
        if (!is.null(b_data$taxonomyReferenceObjects)) {
          if (is.data.frame(b_data$taxonomyReferenceObjects)) {
            return(extract_family_nuclear(b_data$taxonomyReferenceObjects[idx, , drop = FALSE]))
          }
          if (is.list(b_data$taxonomyReferenceObjects)) {
            return(extract_family_nuclear(b_data$taxonomyReferenceObjects[[idx]]))
          }
        }
        character(0)
      })

      counts <- vapply(families_list, length, integer(1))
      strings <- vapply(families_list, function(x) {
        if (length(x) == 0) return("Unknown")
        paste(head(as.character(x), 5), collapse = "; ")
      }, character(1))

      results_list[[i]] <- data.frame(
        inchikey = b_data$inchikey,
        Global_Family_Count = counts,
        Family_List_String = strings,
        stringsAsFactors = FALSE
      )
    }
  }, error = function(e) {
    warning("MongoDB batch ", i, " returned an error: ", e$message)
  })

  if (i %% 10 == 0) gc()
}

if (length(results_list) > 0) {
  df_pass1 <- as_base_df(dplyr::bind_rows(results_list))
} else {
  df_pass1 <- data.frame(
    inchikey = character(0),
    Global_Family_Count = integer(0),
    Family_List_String = character(0),
    stringsAsFactors = FALSE
  )
}

missing_keys <- setdiff(target_inchikeys, df_pass1$inchikey)

if (length(missing_keys) > 0) {
  df_missing <- data.frame(
    inchikey = missing_keys,
    Global_Family_Count = 0,
    Family_List_String = "NOT_IN_DB_OR_NO_FAMILY",
    stringsAsFactors = FALSE
  )
  df_pass1 <- as_base_df(dplyr::bind_rows(df_pass1, df_missing))
}

global_context_C <- as_base_df(
  df_pass1 %>%
    dplyr::distinct(inchikey, .keep_all = TRUE) %>%
    dplyr::mutate(
      Biogeography_Status = dplyr::case_when(
        Global_Family_Count == 0 ~ "Not Found / Error",
        Global_Family_Count == 1 ~ "Exclusive (1 Family)",
        Global_Family_Count <= 10 ~ "Restricted",
        TRUE ~ "Ubiquitous"
      )
    )
)

outfile_gc <- out2_path(paste0(base_tag, "_BIO_C_Global_Context.xlsx"))
writexl::write_xlsx(as_base_df(global_context_C), outfile_gc)
message("[OK] Global context table written: ", basename(outfile_gc))

blacklist_terms <- cfg$blacklist_terms %||% c(
  "oleic", "linoleic", "palmitic", "stearic", "myristic",
  "ergosterol", "sitosterol", "stigmasterol", "campesterol", "cholesterol",
  "quercetin", "rutin", "kaempferol", "catechin", "epicatechin",
  "gallic acid", "ellagic acid", "chlorogenic", "caffeic",
  "lupeol", "amyrin", "friedelin", "betulin",
  "chlorophyll", "carotene", "squalene", "tocopherol", "sucrose", "glucose"
)

if (!"murko_framework" %in% names(uni_df)) uni_df$murko_framework <- NA_character_
if (!"iupac_name" %in% names(uni_df)) uni_df$iupac_name <- NA_character_
if (!"smiles" %in% names(uni_df)) uni_df$smiles <- NA_character_
if (!"molecular_formula" %in% names(uni_df)) uni_df$molecular_formula <- NA_character_
if (!"molecular_weight" %in% names(uni_df)) uni_df$molecular_weight <- NA_real_
if (!"xlogp" %in% names(uni_df)) uni_df$xlogp <- NA_real_
if (!"heavy_atom_number" %in% names(uni_df)) uni_df$heavy_atom_number <- NA_integer_
if (!"family" %in% names(uni_df)) uni_df$family <- NA_character_
if (!"genus" %in% names(uni_df)) uni_df$genus <- NA_character_
if (!"species" %in% names(uni_df)) uni_df$species <- NA_character_

scaffold_counts <- as_base_df(
  uni_df %>%
    dplyr::filter(!is.na(murko_framework), murko_framework != "") %>%
    dplyr::count(murko_framework, name = "Local_Scaffold_Freq")
)

master <- as_base_df(
  uni_df %>%
    dplyr::transmute(
      inchikey,
      smiles,
      iupac_name,
      molecular_formula,
      molecular_weight,
      xlogp,
      heavy_atom_number,
      Primary_Class = hybrid_class,
      Scaffold = murko_framework,
      Taxon_Families = family,
      Taxon_Genera = genus,
      Taxon_Species = species,
      ref_ids
    ) %>%
    dplyr::left_join(scaffold_counts, by = c("Scaffold" = "murko_framework")) %>%
    dplyr::left_join(bio_evidence_A, by = "inchikey") %>%
    dplyr::left_join(global_context_C, by = "inchikey") %>%
    dplyr::mutate(
      Global_Family_Count = ifelse(is.na(Global_Family_Count), 0, Global_Family_Count),
      Best_Potency_nM = suppressWarnings(as.numeric(Best_Potency_nM)),
      Has_Activity = !is.na(Best_Potency_nM),
      pIC50 = ifelse(Has_Activity & Best_Potency_nM > 0, -log10(Best_Potency_nM * 1e-9), 0),
      xlogp_val = suppressWarnings(as.numeric(xlogp)),
      LLE_Score = ifelse(!is.na(xlogp_val) & Has_Activity, pIC50 - xlogp_val, -99),
      Is_Primary = grepl(paste(blacklist_terms, collapse = "|"), iupac_name, ignore.case = TRUE),
      Score_Local = dplyr::case_when(
        tidyr::replace_na(Local_Scaffold_Freq, 1) == 1 ~ 3,
        Local_Scaffold_Freq <= 5 ~ 1,
        TRUE ~ 0
      ),
      Score_Global = dplyr::case_when(
        Global_Family_Count == 0 ~ -10,
        Global_Family_Count == 1 ~ 3,
        Global_Family_Count <= 3 ~ 1,
        TRUE ~ -5
      ),
      Score_Potency = dplyr::case_when(
        Has_Activity & (LLE_Score > 2 | Best_Potency_nM < 200) ~ 3,
        Has_Activity ~ 1,
        TRUE ~ 0
      ),
      PRIORITY_SCORE = Score_Local + Score_Global + (Score_Potency * 2),
      Class = dplyr::case_when(
        Is_Primary ~ "DISCARD (Primary Metabolite)",
        Global_Family_Count == 0 ~ "DISCARD (Unverified Data)",
        Global_Family_Count > 10 ~ "DISCARD (Ubiquitous)",
        Has_Activity & Score_Global > 0 & (LLE_Score > 0 | Best_Potency_nM < 500) ~ "STAR (Exclusive and Active)",
        Has_Activity & Score_Global > 0 ~ "HIT (Low Efficiency)",
        Score_Global > 0 ~ "HIDDEN GEM (High Novelty)",
        TRUE ~ "BASELINE"
      ),
      Combined_Refs = paste(na.omit(ref_ids), na.omit(Best_Ref_DOI), sep = "; ")
    ) %>%
    dplyr::arrange(desc(PRIORITY_SCORE))
)

outfile_master <- out2_path(paste0(base_tag, "_MASTER_LIST_GLOBAL.xlsx"))
writexl::write_xlsx(
  list(
    All_Ranked = as_base_df(master),
    Top_STARS = as_base_df(dplyr::filter(master, grepl("STAR", Class)) %>% head(50)),
    Top_GEMS = as_base_df(dplyr::filter(master, grepl("GEM", Class)) %>% head(50))
  ),
  path = outfile_master
)

message("[OK] Master prioritization table written: ", basename(outfile_master))
message("STAR-class entries identified: ", sum(grepl("STAR", master$Class)))

message("----------------------------------------------------------------")
message(">>> STEP C: target categorization and atlas rendering")
message("----------------------------------------------------------------")

path_bio <- outfile_sum
if (!file.exists(path_bio)) {
  stop(
    "The bioactivity summary file required for atlas rendering was not found: ",
    path_bio,
    "\nRun STEP B successfully before STEP C."
  )
}

file_xls_out <- out2_path(paste0(base_tag, "_Lotus_Final_Database.xlsx"))
file_pdf_macro <- out2_path(paste0(base_tag, "_Lotus_Atlas_MACRO.pdf"))
file_pdf_detail <- out2_path(paste0(base_tag, "_Lotus_Atlas_DETAILED.pdf"))

normalize_name <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("\\s+", " ") %>%
    str_squish()
}

rx <- function(pattern) regex(pattern, ignore_case = TRUE)

categorize_targets_v9 <- function(df_targets) {
  as_base_df(
    df_targets %>%
      dplyr::mutate(
        target_pref_name = normalize_name(target_pref_name),
        Target_Category_L3 = dplyr::case_when(
          str_detect(target_pref_name, rx("Non-protein|Unchecked|No relevant|Molecular identity|Lytechinus|Rattus|Homo sapiens|Mus musculus|Assay|Tumor Cell")) ~ "Assay: Phenotypic/General",
          str_detect(target_pref_name, rx("Antioxidant|Radical scavenging")) ~ "Assay: Antioxidant",
          str_detect(target_pref_name, rx("ADMET|hERG|P-glycoprotein|MDR|ABC transporter|Liver carboxylesterase")) ~ "ADMET/Safety",
          str_detect(target_pref_name, rx("A549|NCI-H|Lewis|HOP|DMS|LXFL|Calu|SF-539|MSTO|EKVX")) ~ "Cells: Lung Cancer",
          str_detect(target_pref_name, rx("MCF7|MDA|T47D|SK-BR|Hs-578|BT-549|CAL-51|HBL-100")) ~ "Cells: Breast Cancer",
          str_detect(target_pref_name, rx("HepG2|Hep 3B|Liver")) ~ "Cells: Liver Cancer",
          str_detect(target_pref_name, rx("HCT|Caco|HT-29|SW480|SW-620|COLO|LoVo|DLD|KM12|HCC")) ~ "Cells: Colorectal Cancer",
          str_detect(target_pref_name, rx("K562|HL-60|Jurkat|CCRF|MOLT|RPMI|P388|L1210|THP|U-937|L5178Y|RL|U-266|Daudi|Raji|697|MV4|SR|NCI/ADR")) ~ "Cells: Leukemia/Lymphoma",
          str_detect(target_pref_name, rx("HeLa|SiHa|KB")) ~ "Cells: Cervical Cancer",
          str_detect(target_pref_name, rx("PC-3|DU-145|LNCaP")) ~ "Cells: Prostate Cancer",
          str_detect(target_pref_name, rx("Melanoma|SK-MEL|A-375|B16|M14|Malme|HaCaT|UACC|M19|A2058|LOX IMVI")) ~ "Cells: Melanoma",
          str_detect(target_pref_name, rx("786-0|A498|ACHN|CAKI|RXF|SN12|MDCK|TK-10|UO-31|XF498")) ~ "Cells: Renal Cancer",
          str_detect(target_pref_name, rx("OVCAR|IGROV|SK-OV|A2780")) ~ "Cells: Ovarian Cancer",
          str_detect(target_pref_name, rx("BGC-823|SGC-7901|MGC-803|AGS|MIA PaCa|PANC-1")) ~ "Cells: GI/Pancreatic Cancer",
          str_detect(target_pref_name, rx("SF-268|SF-295|SNB|U-251|SH-SY5Y|BV-2|HT-22|N9|PC-12")) ~ "Cells: CNS/Neuro",
          str_detect(target_pref_name, rx("HUVEC|Aorta|VSMC|Endothelial")) ~ "Cells: Endothelial/Vascular",
          str_detect(target_pref_name, rx("RAW264|J774|Macrophage|PBMC|Splenocyte|Neutrophil|HMC1")) ~ "Cells: Immune/Macrophage",
          str_detect(target_pref_name, rx("Vero|HEK|CHO|NIH3T3|Fibroblast|L929|BJ|ANN-1|Balb|HT-1080|S1|MEF|Ishikawa|Endometrial")) ~ "Cells: Non-Cancer/General",
          str_detect(target_pref_name, rx("Cell line")) ~ "Cells: Other",
          str_detect(target_pref_name, rx("Virus|HIV|Influenza|Dengue|SARS|Hepatitis|Chikungunya|NS3|Replicase|Large T")) ~ "Pathogen: Virus",
          str_detect(target_pref_name, rx("Bacteria|Staphylococcus|Escherichia|Mycobacterium|Pseudomonas|Streptococcus|Bacillus|Enterococcus|Sortase|Lethal factor")) ~ "Pathogen: Bacteria",
          str_detect(target_pref_name, rx("Fungus|Candida|Aspergillus|Cryptococcus|Trichophyton")) ~ "Pathogen: Fungus",
          str_detect(target_pref_name, rx("Parasite|Plasmodium|Trypanosoma|Leishmania|Entamoeba|Taenia|Cruzipain")) ~ "Pathogen: Parasite",
          str_detect(target_pref_name, rx("COX|Cyclooxygenase|Prostaglandin|Thromboxane")) ~ "Enzymes: Inflammatory (COX/PG)",
          str_detect(target_pref_name, rx("LOX|Lipoxygenase|ALOX")) ~ "Enzymes: Lipoxygenase",
          str_detect(target_pref_name, rx("NOS|Nitric oxide")) ~ "Enzymes: Nitric Oxide Synthase",
          str_detect(target_pref_name, rx("NF-kB|NF-kappa|Tumor necrosis factor|TNF")) ~ "Proteins: NF-kB/TNF Signaling",
          str_detect(target_pref_name, rx("AChE|Cholinesterase|Butyrylcholinesterase")) ~ "Enzymes: Cholinesterase",
          str_detect(target_pref_name, rx("MAO|Monoamine oxidase")) ~ "Enzymes: Monoamine Oxidase",
          str_detect(target_pref_name, rx("Alzheimer|Amyloid|BACE|Tau|Secretase|Synaptojanin|Alpha-synuclein")) ~ "Proteins: Alzheimer/Neuro",
          str_detect(target_pref_name, rx("Serotonin|5-Hydroxytryptamine|5-HT")) ~ "Receptors: Serotonin (5-HT)",
          str_detect(target_pref_name, rx("Dopamine|Adrenergic|Histamine|Muscarinic|GABA|Opioid|Cannabinoid|Adenosine|Glutamate|Taste")) ~ "Receptors: GPCR/Ion (Neuro)",
          str_detect(target_pref_name, rx("Tyrosinase")) ~ "Enzymes: Tyrosinase",
          str_detect(target_pref_name, rx("Sialidase|Neuraminidase")) ~ "Enzymes: Sialidase",
          str_detect(target_pref_name, rx("Glucosidase|Amylase|Maltase|Aldose reductase|Diabetes")) ~ "Enzymes: Diabetes Related",
          str_detect(target_pref_name, rx("Lipase|Fatty acid synthase|Perilipin|ABHD5|Squalene|HMG-CoA")) ~ "Enzymes: Lipid Metabolism",
          str_detect(target_pref_name, rx("Steroid|Aromatase|17-beta|Estrogen synthase|Hydroxysteroid")) ~ "Enzymes: Steroid Metabolism",
          str_detect(target_pref_name, rx("CYP|Cytochrome P450")) ~ "Enzymes: CYP450",
          str_detect(target_pref_name, rx("UGT|Sulfotransferase|Phase II|GST")) ~ "Enzymes: Phase II Metabolism",
          str_detect(target_pref_name, rx("Carbonic anhydrase")) ~ "Enzymes: Carbonic Anhydrase",
          str_detect(target_pref_name, rx("PDE|Phosphodiesterase")) ~ "Enzymes: PDE",
          str_detect(target_pref_name, rx("HDAC|Sirtuin|Histone|Methyltransferase|Polycomb|Demethylase|Epigenetic")) ~ "Epigenetics: HDAC/Chromatin",
          str_detect(target_pref_name, rx("Protease|Peptidase|Cathepsin|Calpain|Thrombin|Trypsin|Kallikrein|Caspase|Renin|Proteasome")) ~ "Enzymes: Proteases (General)",
          str_detect(target_pref_name, rx("MMP|Matrix metallo|Collagenase|Stromelysin|Matrilysin|Gelatinase|Disintegrin")) ~ "Enzymes: Proteases (MMP)",
          str_detect(target_pref_name, rx("Kinase|CDK|Aurora|MAPK|PI3K|Akt|mTOR|GSK|PKC|EGFR|VEGF|FLT3|Src|Abl|JAK|Tyrosine")) ~ "Enzymes: Kinases",
          str_detect(target_pref_name, rx("Phosphatase")) ~ "Enzymes: Phosphatases",
          str_detect(target_pref_name, rx("Estrogen|Androgen|PPAR|Glucocorticoid|Progesterone|Retinoic|Vitamin D|RXR|ERR|Nuclear receptor")) ~ "Receptors: Nuclear",
          str_detect(target_pref_name, rx("Tubulin|Microtubule|Cytoskeleton")) ~ "Proteins: Tubulin/Cytoskeleton",
          str_detect(target_pref_name, rx("Apoptosis|Bcl-2|Bcl-xL|Mcl-1|p53|MDM2|Survivin|BAX|BID")) ~ "Proteins: Apoptosis",
          str_detect(target_pref_name, rx("HSP90|Heat shock|Chaperone|GroEL")) ~ "Proteins: Heat Shock (HSP)",
          str_detect(target_pref_name, rx("HIF|Hypoxia")) ~ "Proteins: Hypoxia (HIF)",
          str_detect(target_pref_name, rx("Channel|Transporter|Pump|CFTR|SGLT|Sodium|Potassium|Calcium")) ~ "Transporters & Channels",
          str_detect(target_pref_name, rx("DNA|RNA|Polymerase|Helicase|Telomerase|Zinc finger|Transcription")) ~ "Proteins: DNA/RNA Regulation",
          str_detect(target_pref_name, rx("Growth factor|Angiopoietin|EGF|HGF|IGF")) ~ "Proteins: Growth Factors",
          str_detect(target_pref_name, rx("Oxidoreductase|Dehydrogenase|Reductase|Oxidase")) ~ "Enzymes: Oxidoreductases (Other)",
          TRUE ~ "Other Targets (Miscellaneous)"
        ),
        Target_Category_Macro = dplyr::case_when(
          str_detect(Target_Category_L3, "Cells:") ~ "Phenotypic: Cancer & Cells",
          str_detect(Target_Category_L3, "Pathogen:") ~ "Infectious Diseases",
          str_detect(Target_Category_L3, "Kinases|Growth|NF-kB|Apoptosis|HSP") ~ "Cell Signaling & Survival",
          str_detect(Target_Category_L3, "Inflammatory|LOX|Nitric|Proteases|MMP") ~ "Inflammation & Proteolysis",
          str_detect(Target_Category_L3, "Neuro|Cholinesterase|Monoamine|GPCR|Serotonin") ~ "Neuroscience Targets",
          str_detect(Target_Category_L3, "Lipid|Diabetes|Steroid|CYP|Phase II|Carbonic|PDE|Tyrosinase|Sialidase") ~ "Metabolic Enzymes",
          str_detect(Target_Category_L3, "Nuclear|Epigenetics|DNA") ~ "Gene Regulation",
          str_detect(Target_Category_L3, "Transporters") ~ "Transporters & Ion Channels",
          str_detect(Target_Category_L3, "Assay|Safety|ADMET") ~ "General Assays & Safety",
          TRUE ~ "Miscellaneous"
        )
      )
  )
}

message("[INFO] Loading bioactivity and taxonomy tables for atlas assembly.")
raw_bio <- as_base_df(readxl::read_excel(path_bio, sheet = "Raw"))
raw_tax <- as_base_df(readxl::read_excel(path_tax_excel, sheet = "uni_enriched"))

col_class <- grep("chemicalTaxonomyNPclassifierClass", colnames(raw_tax), ignore.case = TRUE, value = TRUE)[1]
col_id <- grep("^inchikey$", colnames(raw_tax), ignore.case = TRUE, value = TRUE)[1]

if (is.na(col_class) || is.na(col_id)) {
  stop("Required columns could not be identified in the 'uni_enriched' sheet.")
}

clean_numeric <- function(x) as.numeric(gsub("[^0-9.]", "", as.character(x)))

assign_tier <- function(pActivity, is_censored_gt) {
  dplyr::case_when(
    is_censored_gt ~ "Tier 5 (Limit/Inactive)",
    pActivity > 8.0 ~ "Tier 1 (Elite)",
    pActivity >= 7.0 ~ "Tier 2 (Potent)",
    pActivity >= 6.0 ~ "Tier 3 (Moderate)",
    pActivity > 5.0 ~ "Tier 4 (Weak)",
    TRUE ~ "Tier 5 (Limit/Inactive)"
  )
}

tier_colors <- c(
  "Tier 1 (Elite)" = "#DC0000",
  "Tier 2 (Potent)" = "#FF681F",
  "Tier 3 (Moderate)" = "#FFC500",
  "Tier 4 (Weak)" = "#3B9AB2",
  "Tier 5 (Limit/Inactive)" = "#E0E0E0"
)

raw_tax_index <- as_base_df(
  raw_tax %>%
    dplyr::rename(inchikey = all_of(col_id)) %>%
    dplyr::distinct(inchikey, .keep_all = TRUE)
)

df_master <- as_base_df(
  raw_bio %>%
    dplyr::mutate(
      inchikey = as.character(inchikey),
      val_num = clean_numeric(standard_value),
      is_censored_gt = standard_relation %in% c(">", ">=")
    ) %>%
    dplyr::filter(!is.na(val_num), val_num > 0) %>%
    dplyr::filter(standard_type %in% c("IC50", "EC50", "Ki", "Kd", "MIC", "AC50", "GI50")) %>%
    categorize_targets_v9() %>%
    dplyr::inner_join(raw_tax_index, by = "inchikey") %>%
    dplyr::group_by(inchikey, Target_Category_L3) %>%
    dplyr::arrange(is_censored_gt, val_num) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      pActivity = -log10(val_num * 1e-9),
      Plot_Class = .data[[col_class]],
      Tier = assign_tier(pActivity, is_censored_gt),
      Tier = factor(Tier, levels = names(tier_colors))
    ) %>%
    dplyr::filter(!is.na(Plot_Class), Plot_Class != "Unknown")
)

writexl::write_xlsx(list(MASTER_DATA = as_base_df(df_master)), path = file_xls_out)
message("[OK] Final Excel workbook written: ", basename(file_xls_out))

theme_science <- theme_bw(base_size = 14) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

plot_category_loop <- function(cat_list, cat_column, filename) {
  pdf(filename, width = 14, height = 8)
  for (cat in cat_list) {
    data_cat <- as_base_df(df_master %>% dplyr::filter(.data[[cat_column]] == cat))
    if (nrow(data_cat) < 3) next

    top_classes <- as_base_df(
      data_cat %>%
        dplyr::count(Plot_Class, sort = TRUE) %>%
        head(15)
    )
    data_viz <- as_base_df(data_cat %>% dplyr::filter(Plot_Class %in% top_classes$Plot_Class))

    p1 <- ggplot(data_viz, aes(x = factor(Plot_Class, levels = rev(top_classes$Plot_Class)))) +
      geom_bar(fill = "#2C3E50", width = 0.7) +
      geom_text(stat = "count", aes(label = after_stat(count)), hjust = -0.2, size = 3) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
      coord_flip() +
      labs(title = "Diversity", subtitle = paste0("Total: ", nrow(data_cat)), x = "", y = "Count") +
      theme_science

    p2 <- ggplot(data_viz, aes(x = factor(Plot_Class, levels = rev(top_classes$Plot_Class)), y = pActivity)) +
      geom_boxplot(fill = "gray95", outlier.shape = NA) +
      geom_jitter(aes(color = Tier), width = 0.2, size = 1.8, alpha = 0.7) +
      scale_color_manual(values = tier_colors, drop = FALSE) +
      geom_hline(yintercept = c(5, 6, 7, 8), linetype = "dotted") +
      coord_flip() +
      labs(title = "Potency", x = "", y = "pActivity (-log M)") +
      theme_science +
      theme(axis.text.y = element_blank())

    gridExtra::grid.arrange(
      grid::textGrob(as.character(cat), gp = grid::gpar(fontsize = 16, fontface = "bold"), x = 0, hjust = 0),
      gridExtra::arrangeGrob(p1, p2, ncol = 2, widths = c(1, 1.5)),
      heights = c(0.1, 0.9)
    )
    message("[INFO] Atlas page rendered for category: ", cat)
  }
  dev.off()
}

plot_category_loop(sort(unique(df_master$Target_Category_Macro)), "Target_Category_Macro", file_pdf_macro)
plot_category_loop(sort(unique(df_master$Target_Category_L3)), "Target_Category_L3", file_pdf_detail)

message("[OK] Atlas PDFs written:")
message(" - ", basename(file_pdf_macro))
message(" - ", basename(file_pdf_detail))
message("----------------------------------------------------------------")
message(">>> PART II COMPLETED SUCCESSFULLY")
message("OUT_DIR: ", normalizePath(OUT_DIR))
message("----------------------------------------------------------------")
