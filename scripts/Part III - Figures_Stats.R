# Part III: statistical analyses and figure generation for macro-group comparisons.

`%||%` <- function(a, b) {
  if (is.null(a)) return(b)
  if (length(a) == 0) return(b)
  if (is.character(a) && length(a) == 1 && !nzchar(a)) return(b)
  a
}

as_base_df <- function(x) as.data.frame(x, stringsAsFactors = FALSE)

bind_rows_df <- function(...) {
  as_base_df(dplyr::bind_rows(...))
}

export_df <- function(x) {
  x <- as_base_df(x)
  if (ncol(x) == 0) {
    return(data.frame(note = character(0), stringsAsFactors = FALSE))
  }
  x
}

safe_first_chr <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x)) x[1] else NA_character_
}

safe_num <- function(x) {
  suppressWarnings(as.numeric(gsub("[^0-9.+-eE]", "", as.character(x))))
}

safe_median_num <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x)) median(x, na.rm = TRUE) else NA_real_
}

safe_quantile_num <- function(x, prob) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x)) as.numeric(stats::quantile(x, prob = prob, na.rm = TRUE, names = FALSE)) else NA_real_
}

safe_iqr_num <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x)) stats::IQR(x, na.rm = TRUE) else NA_real_
}

safe_max_num <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x)) max(x, na.rm = TRUE) else NA_real_
}

safe_prop_true <- function(x) {
  x <- as.logical(x %in% TRUE)
  if (length(x)) mean(x, na.rm = TRUE) else NA_real_
}

normalize_unit <- function(x) {
  x <- iconv(as.character(x), from = "", to = "ASCII//TRANSLIT")
  x <- tolower(gsub("[^a-z]", "", x))
  x
}

wilson_ci <- function(k, n, z = 1.96) {
  if (is.na(k) || is.na(n) || n <= 0) return(c(NA_real_, NA_real_))
  phat <- k / n
  denom <- 1 + z^2 / n
  center <- (phat + z^2 / (2 * n)) / denom
  half <- (z * sqrt((phat * (1 - phat) + z^2 / (4 * n)) / n)) / denom
  c(max(0, center - half), min(1, center + half))
}

theme_pub <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = base_size - 1),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )
}

placeholder_plot <- function(title, subtitle) {
  ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::labs(title = title, subtitle = subtitle)
}

add_group_hulls <- function(p, df, x, y, grp, alpha = 0.10) {
  gvals <- unique(df[[grp]])
  hull_df <- list()

  for (g in gvals) {
    sub <- df[df[[grp]] == g & is.finite(df[[x]]) & is.finite(df[[y]]), , drop = FALSE]
    if (nrow(sub) < 3) next
    idx <- chull(sub[[x]], sub[[y]])
    hull_df[[as.character(g)]] <- sub[idx, , drop = FALSE]
  }

  hull <- bind_rows_df(hull_df)
  if (nrow(hull) == 0) return(p)

  p + ggplot2::geom_polygon(
    data = hull,
    ggplot2::aes(x = .data[[x]], y = .data[[y]], fill = .data[[grp]]),
    alpha = alpha,
    color = NA,
    inherit.aes = FALSE,
    show.legend = FALSE
  )
}

cld_from_dunn <- function(dunn_df) {
  comps <- paste(dunn_df$group1, dunn_df$group2, sep = "-")
  pvals <- dunn_df$p.adj
  names(pvals) <- comps
  multcompView::multcompLetters(pvals)$Letters
}

find_first_existing <- function(x) {
  x <- unique(x[!is.na(x) & nzchar(x)])
  hit <- x[file.exists(x)]
  if (length(hit)) hit[1] else NA_character_
}

resolve_file <- function(label, candidates, required = TRUE) {
  hit <- find_first_existing(candidates)
  if (is.na(hit) && required) {
    stop("Required file not found [", label, "]. Tried:\n- ", paste(candidates, collapse = "\n- "))
  }
  hit
}

if (!exists("runtime")) {
  stop("Object 'runtime' was not found. This script must be called by the main pipeline.")
}

OUT_DIR <- runtime$OUT_DIR %||% NA_character_
if (is.na(OUT_DIR) || !nzchar(OUT_DIR) || !dir.exists(OUT_DIR)) {
  stop("OUT_DIR does not exist: ", OUT_DIR)
}

cfg <- if (exists("cfg")) cfg else list()
cfg$auto_install_pkgs <- cfg$auto_install_pkgs %||% FALSE

BASE_TAG <- runtime$base_tag %||% if (exists("base_tag")) base_tag else "lotus_run"

out_dir <- file.path(OUT_DIR, "PartIII_ALL")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_xlsx <- file.path(out_dir, "macro_groups_STATS_MASTER.xlsx")
out_fig1 <- file.path(out_dir, "Fig1_ChemistryStats.pdf")
out_fig2 <- file.path(out_dir, "Fig2_OrdinationStats.pdf")
out_fig3 <- file.path(out_dir, "Fig3_BioprospectingLandscape.pdf")

set.seed(1)

pkgs <- c(
  "readr", "readxl", "openxlsx", "janitor", "dplyr", "tidyr", "stringr",
  "ggplot2", "ggrepel", "scales", "vegan", "rstatix", "multcompView", "broom"
)

missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_pkgs)) {
  if (isTRUE(cfg$auto_install_pkgs)) {
    install.packages(missing_pkgs, dependencies = TRUE)
  } else {
    stop(
      "Missing packages: ", paste(missing_pkgs, collapse = ", "),
      "\nSet cfg$auto_install_pkgs <- TRUE if automatic installation is desired."
    )
  }
}

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(openxlsx)
  library(janitor)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(vegan)
  library(rstatix)
  library(multcompView)
  library(broom)
})

project_root <- getwd()

normalize_macro_group_name <- function(x) {
  x <- stringr::str_trim(as.character(x))
  dplyr::case_when(
    x == "Basal Angiosperms" ~ "ANA",
    x %in% c("Caryophyllales", "Santalales") ~ "Caryophy_Santalales",
    TRUE ~ x
  )
}

auto_generate_groups <- function(project_root, OUT_DIR, cfg = list()) {
  groups_candidates <- c(
    cfg$path_groups %||% NA_character_,
    file.path(project_root, "inputs", "groups.xlsx"),
    file.path(project_root, "groups.xlsx"),
    file.path(OUT_DIR, "inputs", "groups.xlsx"),
    file.path(OUT_DIR, "groups.xlsx")
  )
  
  existing_groups <- find_first_existing(groups_candidates)
  if (!is.na(existing_groups)) {
    message("[INFO] Using existing groups.xlsx: ", existing_groups)
    return(existing_groups)
  }
  
  apg_candidates <- c(
    cfg$path_apg_assignments %||% NA_character_,
    file.path(project_root, "phylo_outputs_FINAL", "APG_clade_assignments.csv"),
    file.path(project_root, "outputs", "phylo_outputs_FINAL", "APG_clade_assignments.csv"),
    file.path(OUT_DIR, "phylo_outputs_FINAL", "APG_clade_assignments.csv"),
    file.path(OUT_DIR, "APG_clade_assignments.csv"),
    Sys.glob(file.path(project_root, "**", "APG_clade_assignments.csv")),
    Sys.glob(file.path(OUT_DIR, "**", "APG_clade_assignments.csv"))
  )
  
  path_apg_assign <- find_first_existing(apg_candidates)
  
  if (is.na(path_apg_assign)) {
    stop(
      "groups.xlsx was not found, and APG_clade_assignments.csv was also not found.\n",
      "Provide an existing groups.xlsx or make APG_clade_assignments.csv available."
    )
  }
  
  apg_df <- readr::read_csv(path_apg_assign, show_col_types = FALSE)
  
  required_cols <- c("ID", "Color", "FinalClade")
  missing_cols <- setdiff(required_cols, names(apg_df))
  if (length(missing_cols) > 0) {
    stop(
      "APG_clade_assignments.csv is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  groups_auto <- apg_df %>%
    dplyr::transmute(
      family = stringr::str_trim(as.character(ID)),
      color = stringr::str_trim(as.character(Color)),
      macro_group = normalize_macro_group_name(FinalClade)
    ) %>%
    dplyr::filter(!is.na(family), family != "", !is.na(macro_group), macro_group != "") %>%
    dplyr::mutate(
      color = ifelse(stringr::str_detect(color, "^#"), color, paste0("#", color)),
      color = stringr::str_replace_all(color, "\\s+", "")
    ) %>%
    dplyr::distinct(family, .keep_all = TRUE)
  
  if (nrow(groups_auto) == 0) {
    stop("Automatic generation of groups.xlsx failed: no valid rows were produced.")
  }
  
  out_groups_dir <- file.path(project_root, "inputs")
  dir.create(out_groups_dir, showWarnings = FALSE, recursive = TRUE)
  
  out_groups <- file.path(out_groups_dir, "groups.xlsx")
  openxlsx::write.xlsx(groups_auto, out_groups, overwrite = TRUE)
  
  message("[OK] groups.xlsx was generated automatically at: ", out_groups)
  out_groups
}

path_groups <- auto_generate_groups(project_root = project_root, OUT_DIR = OUT_DIR, cfg = cfg)

path_rdkit <- resolve_file(
  "RDKit annotations",
  c(
    cfg$path_rdkit %||% NA_character_,
    file.path(project_root, "inputs", "lotus_flavonoids_rdkit_annotations.csv"),
    file.path(project_root, "lotus_flavonoids_rdkit_annotations.csv"),
    file.path(OUT_DIR, "PartII_ALL", paste0(BASE_TAG, "__flavonoids_for_rdkit.csv")),
    file.path(OUT_DIR, "inputs", "lotus_flavonoids_rdkit_annotations.csv"),
    file.path(OUT_DIR, "lotus_flavonoids_rdkit_annotations.csv"),
    Sys.glob(file.path(OUT_DIR, "PartII_ALL", "*__flavonoids_for_rdkit.csv"))
  ),
  required = TRUE
)

path_chembl <- resolve_file(
  "final bioactivity workbook",
  c(
    cfg$path_chembl %||% NA_character_,
    file.path(project_root, "inputs", "Lotus_Final_Database.xlsx"),
    file.path(project_root, "inputs", "Lotus_Final_Database_v9_Fixed.xlsx"),
    file.path(project_root, "Lotus_Final_Database.xlsx"),
    file.path(project_root, "Lotus_Final_Database_v9_Fixed.xlsx"),
    file.path(OUT_DIR, "PartII_ALL", paste0(BASE_TAG, "_Lotus_Final_Database.xlsx")),
    file.path(OUT_DIR, "PartII_ALL", paste0(BASE_TAG, "_Lotus_Final_Database_v9_Fixed.xlsx")),
    file.path(OUT_DIR, "inputs", "Lotus_Final_Database.xlsx"),
    file.path(OUT_DIR, "inputs", "Lotus_Final_Database_v9_Fixed.xlsx"),
    Sys.glob(file.path(OUT_DIR, "PartII_ALL", "*Lotus_Final_Database*.xlsx"))
  ),
  required = TRUE
)

path_lotus <- NA_character_
if (!exists("lin_enriched", inherits = TRUE)) {
  path_lotus <- resolve_file(
    "LOTUS Excel workbook",
    c(
      cfg$path_lotus %||% NA_character_,
      file.path(OUT_DIR, paste0(BASE_TAG, ".xlsx")),
      Sys.glob(file.path(OUT_DIR, "lotus_*.xlsx"))
    ),
    required = FALSE
  )
}

groups_raw <- as_base_df(readxl::read_excel(path_groups, col_names = FALSE, .name_repair = "minimal"))
if (ncol(groups_raw) < 3) {
  stop("The macro-group mapping file must contain at least three columns: family, color, and macro_group.")
}

groups <- as_base_df(
  data.frame(
    family = stringr::str_trim(as.character(groups_raw[[1]])),
    color = stringr::str_trim(as.character(groups_raw[[2]])),
    macro_group = stringr::str_trim(as.character(groups_raw[[3]])),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(!is.na(family), family != "", !is.na(macro_group), macro_group != "") %>%
    dplyr::filter(tolower(family) != "family") %>%
    dplyr::mutate(
      color = ifelse(stringr::str_detect(color, "^#"), color, paste0("#", color)),
      color = stringr::str_replace_all(color, "\\s+", "")
    ) %>%
    dplyr::distinct(family, .keep_all = TRUE)
)

if (nrow(groups) == 0) {
  stop("No valid macro-group mapping entries were found in groups.xlsx.")
}

apg_order <- c(
  "Lycophytes", "Ferns", "Gymnosperms", "ANA", "Magnoliids",
  "Monocots", "Monocots (Commelinids)", "Basal Eudicots", "Saxifragales",
  "Rosids (Fabids)", "Rosids (Malvids)", "Caryophy_Santalales",
  "Asterids (Basal)", "Asterids (Lamiids)", "Asterids (Campanulids)"
)

present_groups <- unique(groups$macro_group)
macro_levels <- c(apg_order[apg_order %in% present_groups], setdiff(sort(present_groups), apg_order))

palette_df <- as_base_df(groups[!duplicated(groups$macro_group), c("macro_group", "color"), drop = FALSE])
pal_macro <- setNames(palette_df$color, palette_df$macro_group)

rdkit <- as_base_df(
  janitor::clean_names(
    as_base_df(
      readr::read_csv(path_rdkit, show_col_types = FALSE, progress = FALSE)
    )
  )
)

if ("chemical_taxonomy_npclassifier_class" %in% names(rdkit) && !"class_np" %in% names(rdkit)) {
  rdkit$class_np <- rdkit$chemical_taxonomy_npclassifier_class
}

for (nm in c("family", "genus", "species", "inchikey", "smiles", "class_np", "core14")) {
  if (nm %in% names(rdkit)) rdkit[[nm]] <- stringr::str_trim(as.character(rdkit[[nm]]))
}

if (!"core14" %in% names(rdkit) && "inchikey" %in% names(rdkit)) {
  rdkit$core14 <- substr(as.character(rdkit$inchikey), 1, 14)
}

required_rdkit_cols <- c("inchikey", "family", "class_np", "core14")
missing_rdkit_cols <- setdiff(required_rdkit_cols, names(rdkit))
if (length(missing_rdkit_cols) > 0) {
  stop("The RDKit annotation table is missing required columns: ", paste(missing_rdkit_cols, collapse = ", "))
}

if (!"smiles" %in% names(rdkit)) rdkit$smiles <- NA_character_

if ("hba" %in% names(rdkit) && !"num_hba" %in% names(rdkit)) rdkit$num_hba <- rdkit$hba
if ("hbd" %in% names(rdkit) && !"num_hbd" %in% names(rdkit)) rdkit$num_hbd <- rdkit$hbd
if ("rings" %in% names(rdkit) && !"num_rings" %in% names(rdkit)) rdkit$num_rings <- rdkit$rings

for (nm in c("mw", "logp", "tpsa", "num_hba", "num_hbd", "num_rings", "heavy_atom_number")) {
  if (nm %in% names(rdkit)) rdkit[[nm]] <- safe_num(rdkit[[nm]])
}

for (nm in c("has_phenolic_oh", "has_methoxy_aryl", "has_prenyl_like", "has_conj_carbonyl", "has_probable_sugar", "contains_sugar")) {
  if (nm %in% names(rdkit)) {
    rdkit[[nm]] <- tolower(as.character(rdkit[[nm]])) %in% c("true", "t", "1", "yes")
  }
}

rdkit <- as_base_df(
  rdkit %>%
    dplyr::filter(!is.na(family), family != "", !is.na(inchikey), inchikey != "") %>%
    dplyr::left_join(groups[, c("family", "macro_group"), drop = FALSE], by = "family") %>%
    dplyr::filter(!is.na(macro_group), macro_group != "") %>%
    dplyr::mutate(macro_group = factor(macro_group, levels = macro_levels))
)

lotus_available <- FALSE
lotus <- data.frame(
  inchikey = character(0),
  family = character(0),
  macro_group = character(0),
  stringsAsFactors = FALSE
)

if (exists("lin_enriched", inherits = TRUE)) {
  lotus <- as_base_df(janitor::clean_names(as_base_df(get("lin_enriched", inherits = TRUE))))
  lotus_available <- TRUE
} else if (!is.na(path_lotus) && nzchar(path_lotus) && file.exists(path_lotus)) {
  lotus <- as_base_df(
    janitor::clean_names(
      as_base_df(
        readxl::read_excel(path_lotus, sheet = "lin_enriched", .name_repair = "minimal")
      )
    )
  )
  lotus_available <- TRUE
}

if (lotus_available) {
  for (nm in c("family", "genus", "species", "inchikey", "smiles", "lotus_id")) {
    if (nm %in% names(lotus)) lotus[[nm]] <- stringr::str_trim(as.character(lotus[[nm]]))
  }

  lotus <- as_base_df(
    lotus %>%
      dplyr::filter(!is.na(family), family != "", !is.na(inchikey), inchikey != "") %>%
      dplyr::left_join(groups[, c("family", "macro_group"), drop = FALSE], by = "family") %>%
      dplyr::filter(!is.na(macro_group), macro_group != "") %>%
      dplyr::mutate(macro_group = factor(macro_group, levels = macro_levels))
  )
}

if (lotus_available && nrow(lotus) > 0) {
  if ("topopsa" %in% names(lotus) && !"topo_psa" %in% names(lotus)) lotus$topo_psa <- lotus$topopsa

  lotus_desc_lookup <- as_base_df(
    lotus %>%
      dplyr::select(
        inchikey,
        dplyr::any_of(c("molecular_weight", "xlogp", "topo_psa", "heavy_atom_number"))
      ) %>%
      dplyr::mutate(
        mw_lotus = if ("molecular_weight" %in% names(.)) safe_num(molecular_weight) else NA_real_,
        logp_lotus = if ("xlogp" %in% names(.)) safe_num(xlogp) else NA_real_,
        tpsa_lotus = if ("topo_psa" %in% names(.)) safe_num(topo_psa) else NA_real_,
        heavy_lotus = if ("heavy_atom_number" %in% names(.)) safe_num(heavy_atom_number) else NA_real_
      ) %>%
      dplyr::group_by(inchikey) %>%
      dplyr::summarise(
        mw_lotus = safe_median_num(mw_lotus),
        logp_lotus = safe_median_num(logp_lotus),
        tpsa_lotus = safe_median_num(tpsa_lotus),
        heavy_lotus = safe_median_num(heavy_lotus),
        .groups = "drop"
      )
  )
} else {
  lotus_desc_lookup <- data.frame(
    inchikey = character(0),
    mw_lotus = numeric(0),
    logp_lotus = numeric(0),
    tpsa_lotus = numeric(0),
    heavy_lotus = numeric(0),
    stringsAsFactors = FALSE
  )
}

rdkit_aug <- as_base_df(rdkit %>% dplyr::left_join(lotus_desc_lookup, by = "inchikey"))

rdkit_has_mw <- "mw" %in% names(rdkit_aug)
rdkit_has_logp <- "logp" %in% names(rdkit_aug)
rdkit_has_tpsa <- "tpsa" %in% names(rdkit_aug)
rdkit_has_heavy <- "heavy_atom_number" %in% names(rdkit_aug)
rdkit_has_hba <- "num_hba" %in% names(rdkit_aug)
rdkit_has_hbd <- "num_hbd" %in% names(rdkit_aug)
rdkit_has_rings <- "num_rings" %in% names(rdkit_aug)
rdkit_has_phoh <- "has_phenolic_oh" %in% names(rdkit_aug)
rdkit_has_methoxy <- "has_methoxy_aryl" %in% names(rdkit_aug)
rdkit_has_prenyl <- "has_prenyl_like" %in% names(rdkit_aug)
rdkit_has_carbonyl <- "has_conj_carbonyl" %in% names(rdkit_aug)
rdkit_has_sugar <- "has_probable_sugar" %in% names(rdkit_aug)

chem_compounds <- as_base_df(
  rdkit_aug %>%
    dplyr::group_by(inchikey) %>%
    dplyr::summarise(
      smiles = safe_first_chr(smiles),
      class_np = safe_first_chr(class_np),
      core14 = safe_first_chr(core14),
      mw = if (rdkit_has_mw) safe_median_num(mw) else safe_median_num(mw_lotus),
      logp = if (rdkit_has_logp) safe_median_num(logp) else safe_median_num(logp_lotus),
      tpsa = if (rdkit_has_tpsa) safe_median_num(tpsa) else safe_median_num(tpsa_lotus),
      num_hba = if (rdkit_has_hba) safe_median_num(num_hba) else NA_real_,
      num_hbd = if (rdkit_has_hbd) safe_median_num(num_hbd) else NA_real_,
      num_rings = if (rdkit_has_rings) safe_median_num(num_rings) else NA_real_,
      heavy_atom_number = if (rdkit_has_heavy) safe_median_num(heavy_atom_number) else safe_median_num(heavy_lotus),
      has_phenolic_oh = if (rdkit_has_phoh) any(has_phenolic_oh %in% TRUE, na.rm = TRUE) else NA,
      has_methoxy_aryl = if (rdkit_has_methoxy) any(has_methoxy_aryl %in% TRUE, na.rm = TRUE) else NA,
      has_prenyl_like = if (rdkit_has_prenyl) any(has_prenyl_like %in% TRUE, na.rm = TRUE) else NA,
      has_conj_carbonyl = if (rdkit_has_carbonyl) any(has_conj_carbonyl %in% TRUE, na.rm = TRUE) else NA,
      has_probable_sugar = if (rdkit_has_sugar) any(has_probable_sugar %in% TRUE, na.rm = TRUE) else NA,
      .groups = "drop"
    ) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ ifelse(is.infinite(.x), NA_real_, .x)))
)

fam_comp <- as_base_df(rdkit %>% dplyr::distinct(family, macro_group, inchikey))

fam_scaf <- as_base_df(
  rdkit %>%
    dplyr::filter(!is.na(core14), core14 != "") %>%
    dplyr::distinct(family, macro_group, core14, inchikey)
)

fam_class <- as_base_df(
  rdkit %>%
    dplyr::filter(!is.na(class_np), class_np != "") %>%
    dplyr::distinct(family, macro_group, class_np, inchikey)
)

family_metrics <- as_base_df(
  fam_scaf %>%
    dplyr::group_by(family, macro_group) %>%
    dplyr::summarise(
      n_compounds = dplyr::n_distinct(inchikey),
      n_scaffolds = dplyr::n_distinct(core14),
      novelty_ratio = n_scaffolds / pmax(n_compounds, 1),
      .groups = "drop"
    )
)

if (nrow(family_metrics) == 0) {
  stop("No family-level scaffold records were available after filtering and mapping.")
}

if (any(c("has_prenyl_like", "has_probable_sugar", "has_methoxy_aryl", "has_phenolic_oh", "has_conj_carbonyl") %in% names(chem_compounds))) {
  deco_by_family <- as_base_df(
    fam_comp %>%
      dplyr::left_join(
        chem_compounds[, intersect(c("inchikey", "has_prenyl_like", "has_probable_sugar", "has_methoxy_aryl", "has_phenolic_oh", "has_conj_carbonyl"), names(chem_compounds)), drop = FALSE],
        by = "inchikey"
      ) %>%
      dplyr::group_by(family) %>%
      dplyr::summarise(
        prop_prenyl = safe_prop_true(has_prenyl_like),
        prop_sugar = safe_prop_true(has_probable_sugar),
        prop_methoxy = safe_prop_true(has_methoxy_aryl),
        prop_phenolic_oh = safe_prop_true(has_phenolic_oh),
        prop_conj_carbonyl = safe_prop_true(has_conj_carbonyl),
        .groups = "drop"
      )
  )
} else {
  deco_by_family <- data.frame(family = unique(family_metrics$family), stringsAsFactors = FALSE)
}

desc_present <- intersect(c("mw", "logp", "tpsa", "num_hba", "num_hbd", "num_rings", "heavy_atom_number"), names(chem_compounds))

if (length(desc_present) > 0) {
  desc_by_family <- as_base_df(
    fam_comp %>%
      dplyr::left_join(chem_compounds[, c("inchikey", desc_present), drop = FALSE], by = "inchikey") %>%
      dplyr::group_by(family) %>%
      dplyr::summarise(
        dplyr::across(dplyr::all_of(desc_present), ~ safe_median_num(.x), .names = "med_{.col}"),
        .groups = "drop"
      ) %>%
      dplyr::mutate(dplyr::across(where(is.numeric), ~ ifelse(is.infinite(.x), NA_real_, .x)))
  )
} else {
  desc_by_family <- data.frame(family = unique(family_metrics$family), stringsAsFactors = FALSE)
}

family_metrics <- as_base_df(
  family_metrics %>%
    dplyr::left_join(deco_by_family, by = "family") %>%
    dplyr::left_join(desc_by_family, by = "family") %>%
    dplyr::mutate(macro_group = factor(macro_group, levels = macro_levels))
)

macro_metrics <- as_base_df(
  family_metrics %>%
    dplyr::group_by(macro_group) %>%
    dplyr::summarise(
      n_families = dplyr::n_distinct(family),
      n_compounds_median = safe_median_num(n_compounds),
      n_compounds_iqr = safe_iqr_num(n_compounds),
      n_scaffolds_median = safe_median_num(n_scaffolds),
      novelty_ratio_median = safe_median_num(novelty_ratio),
      .groups = "drop"
    )
)

message("[INFO] Running univariate tests across macro-groups.")

univar_targets <- c(
  "n_compounds", "n_scaffolds", "novelty_ratio",
  "prop_prenyl", "prop_sugar", "prop_methoxy", "prop_phenolic_oh", "prop_conj_carbonyl",
  paste0("med_", desc_present)
)
univar_targets <- univar_targets[univar_targets %in% names(family_metrics)]

univar_kw <- list()
univar_dunn <- list()
univar_cld <- list()

for (m in univar_targets) {
  dfm <- as_base_df(
    family_metrics %>%
      dplyr::filter(!is.na(macro_group), macro_group != "") %>%
      dplyr::select(macro_group, value = all_of(m)) %>%
      dplyr::filter(is.finite(value))
  )

  if (nrow(dfm) < 12 || dplyr::n_distinct(dfm$macro_group) < 2) next

  kw <- as_base_df(rstatix::kruskal_test(dfm, value ~ macro_group))
  eff <- as_base_df(rstatix::kruskal_effsize(dfm, value ~ macro_group))
  kw$metric <- m
  kw$eps2 <- eff$effsize[1]
  univar_kw[[m]] <- kw

  dunn <- as_base_df(
    rstatix::dunn_test(dfm, value ~ macro_group, p.adjust.method = "BH") %>%
      dplyr::mutate(metric = m)
  )
  univar_dunn[[m]] <- dunn

  letters <- tryCatch(cld_from_dunn(dunn), error = function(e) NULL)
  if (!is.null(letters)) {
    univar_cld[[m]] <- data.frame(
      metric = m,
      macro_group = names(letters),
      cld = unname(letters),
      stringsAsFactors = FALSE
    )
  }
}

univar_kw_tbl <- bind_rows_df(univar_kw)
univar_dunn_tbl <- bind_rows_df(univar_dunn)
univar_cld_tbl <- bind_rows_df(univar_cld)

message("[INFO] Building multivariate incidence matrices.")

min_scaffold_family_support <- 3L
max_scaffold_features <- 5000L

scaf_freq <- as_base_df(
  fam_scaf %>%
    dplyr::distinct(family, core14) %>%
    dplyr::count(core14, name = "n_families") %>%
    dplyr::arrange(desc(n_families))
)

keep_scaf <- scaf_freq %>%
  dplyr::filter(n_families >= min_scaffold_family_support) %>%
  dplyr::slice_head(n = max_scaffold_features) %>%
  dplyr::pull(core14)

fam_scaf_reduced <- as_base_df(fam_scaf %>% dplyr::filter(core14 %in% keep_scaf))
message("[INFO] Scaffold reduction retained ", length(keep_scaf), " scaffold features.")

mat_scaf <- as_base_df(
  fam_scaf_reduced %>%
    dplyr::mutate(val = 1L) %>%
    dplyr::select(family, macro_group, core14, val) %>%
    dplyr::distinct() %>%
    tidyr::pivot_wider(names_from = core14, values_from = val, values_fill = 0L)
)

meta_fam <- as_base_df(mat_scaf[, c("family", "macro_group"), drop = FALSE])
X_scaf <- as_base_df(mat_scaf[, setdiff(names(mat_scaf), c("family", "macro_group")), drop = FALSE])

mat_class <- as_base_df(
  fam_class %>%
    dplyr::mutate(val = 1L) %>%
    dplyr::select(family, macro_group, class_np, val) %>%
    dplyr::distinct() %>%
    tidyr::pivot_wider(names_from = class_np, values_from = val, values_fill = 0L)
)

meta_class <- as_base_df(mat_class[, c("family", "macro_group"), drop = FALSE])
X_class <- as_base_df(mat_class[, setdiff(names(mat_class), c("family", "macro_group")), drop = FALSE])

med_cols <- grep("^med_", names(family_metrics), value = TRUE)
X_desc <- as_base_df(family_metrics[, c("family", "macro_group", med_cols), drop = FALSE])

empty_tbl <- function() data.frame(stringsAsFactors = FALSE)

run_permanova_block <- function(X, meta, dist_method, binary = FALSE, tag = "X", permutations = 999, run_nmds = TRUE) {
  keep <- complete.cases(meta$macro_group) & !is.na(meta$macro_group)
  X2 <- as_base_df(X[keep, , drop = FALSE])
  meta2 <- as_base_df(meta[keep, , drop = FALSE])
  meta2$macro_group <- droplevels(factor(meta2$macro_group))

  if (!binary) {
    for (nm in names(X2)) X2[[nm]] <- safe_num(X2[[nm]])
    ok_cols <- vapply(X2, function(v) !all(is.na(v)), logical(1))
    X2 <- as_base_df(X2[, ok_cols, drop = FALSE])
    for (nm in names(X2)) {
      v <- X2[[nm]]
      if (any(is.na(v))) X2[[nm]][is.na(v)] <- safe_median_num(v)
    }
    X2 <- as.matrix(X2)
  } else {
    X2 <- as.matrix(X2)
  }

  if (binary) {
    row_support <- rowSums(X2 > 0, na.rm = TRUE)
    keep_rows <- row_support > 0
    if (sum(!keep_rows) > 0) {
      message("[INFO] Block ", tag, ": dropping ", sum(!keep_rows), " empty rows before distance calculation.")
      X2 <- X2[keep_rows, , drop = FALSE]
      meta2 <- as_base_df(meta2[keep_rows, , drop = FALSE])
      meta2$macro_group <- droplevels(factor(meta2$macro_group))
    }
  }

  if (nrow(X2) < 10 || ncol(X2) < 2 || nlevels(meta2$macro_group) < 2) {
    return(list(tag = tag, permanova = empty_tbl(), betadisper = empty_tbl(), pairwise = empty_tbl(), ord_pcoa = empty_tbl(), ord_nmds = empty_tbl()))
  }

  dist_obj <- tryCatch(
    vegan::vegdist(X2, method = dist_method, binary = binary),
    error = function(e) NULL
  )

  if (is.null(dist_obj) || any(is.na(dist_obj))) {
    return(list(tag = tag, permanova = empty_tbl(), betadisper = empty_tbl(), pairwise = empty_tbl(), ord_pcoa = empty_tbl(), ord_nmds = empty_tbl()))
  }

  perma_tbl <- empty_tbl()
  tryCatch({
    perma <- vegan::adonis2(dist_obj ~ macro_group, data = meta2, permutations = permutations)
    perma_tbl <- as_base_df(broom::tidy(perma))
    perma_tbl$tag <- tag
  }, error = function(e) NULL)

  bd_tbl <- empty_tbl()
  tryCatch({
    bd <- vegan::betadisper(dist_obj, meta2$macro_group)
    bd_tbl <- as_base_df(broom::tidy(stats::anova(bd)))
    bd_tbl$tag <- tag
  }, error = function(e) NULL)

  pw_tbl <- empty_tbl()
  tryCatch({
    levs <- levels(meta2$macro_group)
    if (length(levs) >= 2) {
      pw_list <- list()
      combos <- combn(levs, 2, simplify = FALSE)
      for (pair in combos) {
        g1 <- pair[1]
        g2 <- pair[2]
        idx <- which(meta2$macro_group %in% c(g1, g2))
        if (length(idx) < 6) next

        d_sub <- as.dist(as.matrix(dist_obj)[idx, idx])
        m_sub <- as_base_df(meta2[idx, , drop = FALSE])
        m_sub$macro_group <- droplevels(factor(m_sub$macro_group))
        if (nlevels(m_sub$macro_group) < 2) next

        try({
          p_sub <- vegan::adonis2(d_sub ~ macro_group, data = m_sub, permutations = permutations)
          t_sub <- as_base_df(broom::tidy(p_sub))
          t_sub <- as_base_df(t_sub[t_sub$term == "macro_group", , drop = FALSE])
          if (nrow(t_sub) > 0) {
            t_sub$tag <- tag
            t_sub$group1 <- g1
            t_sub$group2 <- g2
            pw_list[[paste(g1, g2)]] <- t_sub
          }
        }, silent = TRUE)
      }

      pw_tbl <- bind_rows_df(pw_list)
      if (nrow(pw_tbl) > 0 && "p.value" %in% names(pw_tbl)) {
        pw_tbl$p_adj <- stats::p.adjust(pw_tbl$p.value, method = "BH")
      }
    }
  }, error = function(e) NULL)

  ord_pcoa <- empty_tbl()
  tryCatch({
    pcoa <- stats::cmdscale(dist_obj, k = 2, eig = TRUE)
    pts <- as_base_df(pcoa$points)
    if (ncol(pts) >= 2) {
      names(pts)[1:2] <- c("Axis1", "Axis2")
      ord_pcoa <- bind_rows_df(cbind(meta2, pts))
      ord_pcoa$tag <- tag
    }
  }, error = function(e) NULL)

  ord_nmds <- empty_tbl()
  if (isTRUE(run_nmds) && ncol(X2) < 3000) {
    tryCatch({
      nmds <- vegan::metaMDS(dist_obj, k = 2, trymax = 20, autotransform = FALSE, trace = FALSE)
      pts <- as_base_df(nmds$points)
      if (ncol(pts) >= 2) {
        names(pts)[1:2] <- c("NMDS1", "NMDS2")
        ord_nmds <- bind_rows_df(cbind(meta2, pts))
        ord_nmds$tag <- tag
        ord_nmds$stress <- nmds$stress
      }
    }, error = function(e) NULL)
  }

  list(tag = tag, permanova = perma_tbl, betadisper = bd_tbl, pairwise = pw_tbl, ord_pcoa = ord_pcoa, ord_nmds = ord_nmds)
}

message("[INFO] Running scaffold composition analysis.")
res_scaf <- run_permanova_block(X_scaf, meta_fam, dist_method = "jaccard", binary = TRUE, tag = "Scaffold_Jaccard")

message("[INFO] Running class composition analysis.")
res_class <- run_permanova_block(X_class, meta_class, dist_method = "jaccard", binary = TRUE, tag = "Class_Jaccard")

message("[INFO] Running descriptor-space analysis.")
meta_desc <- as_base_df(X_desc[, c("family", "macro_group"), drop = FALSE])
Xd2 <- as_base_df(X_desc[, setdiff(names(X_desc), c("family", "macro_group")), drop = FALSE])
for (nm in names(Xd2)) Xd2[[nm]] <- safe_num(Xd2[[nm]])

ok_desc_cols <- vapply(Xd2, function(v) !all(is.na(v)), logical(1))
Xd2 <- as_base_df(Xd2[, ok_desc_cols, drop = FALSE])
for (nm in names(Xd2)) {
  if (any(is.na(Xd2[[nm]]))) Xd2[[nm]][is.na(Xd2[[nm]])] <- safe_median_num(Xd2[[nm]])
}

vars <- if (ncol(Xd2) > 0) vapply(Xd2, function(x) stats::var(x, na.rm = TRUE), numeric(1)) else numeric(0)
keep_var <- !is.na(vars) & vars > 1e-12
Xd2 <- as_base_df(Xd2[, keep_var, drop = FALSE])

if (ncol(Xd2) < 1) {
  res_desc <- list(
    tag = "Descriptors_Euclidean",
    permanova = data.frame(term = "SKIPPED_NO_DATA", p.value = NA_real_, tag = "Descriptors_Euclidean", stringsAsFactors = FALSE),
    betadisper = empty_tbl(),
    pairwise = empty_tbl(),
    ord_pcoa = empty_tbl(),
    ord_nmds = empty_tbl()
  )
} else {
  Xd_final <- scale(as.matrix(Xd2))
  Xd_final[!is.finite(Xd_final)] <- 0
  keep_rows <- !is.na(meta_desc$macro_group)
  Xd_final <- Xd_final[keep_rows, , drop = FALSE]
  meta_final <- as_base_df(meta_desc[keep_rows, , drop = FALSE])
  meta_final$macro_group <- droplevels(factor(meta_final$macro_group))

  if (nrow(Xd_final) < 5 || nlevels(meta_final$macro_group) < 2) {
    res_desc <- list(tag = "Descriptors_Euclidean", permanova = empty_tbl(), betadisper = empty_tbl(), pairwise = empty_tbl(), ord_pcoa = empty_tbl(), ord_nmds = empty_tbl())
  } else {
    dist_desc <- tryCatch(
      vegan::vegdist(Xd_final, method = "euclidean"),
      error = function(e) NULL
    )

    if (is.null(dist_desc) || any(!is.finite(dist_desc))) {
      res_desc <- list(tag = "Descriptors_Euclidean", permanova = empty_tbl(), betadisper = empty_tbl(), pairwise = empty_tbl(), ord_pcoa = empty_tbl(), ord_nmds = empty_tbl())
    } else {
      perma_tbl_d <- empty_tbl()
      tryCatch({
        perma_desc <- vegan::adonis2(dist_desc ~ macro_group, data = meta_final, permutations = 999)
        perma_tbl_d <- as_base_df(broom::tidy(perma_desc))
        perma_tbl_d$tag <- "Descriptors_Euclidean"
      }, error = function(e) NULL)

      bd_tbl_d <- empty_tbl()
      tryCatch({
        bd_desc <- vegan::betadisper(dist_desc, meta_final$macro_group)
        bd_tbl_d <- as_base_df(broom::tidy(stats::anova(bd_desc)))
        bd_tbl_d$tag <- "Descriptors_Euclidean"
      }, error = function(e) NULL)

      ord_pcoa_d <- empty_tbl()
      tryCatch({
        pcoa_desc <- stats::cmdscale(dist_desc, k = 2, eig = TRUE)
        pts <- as_base_df(pcoa_desc$points)
        if (ncol(pts) >= 2) {
          names(pts)[1:2] <- c("Axis1", "Axis2")
          ord_pcoa_d <- bind_rows_df(cbind(meta_final, pts))
          ord_pcoa_d$tag <- "Descriptors_Euclidean"
        }
      }, error = function(e) NULL)

      pw_tbl_d <- empty_tbl()
      tryCatch({
        levs <- levels(meta_final$macro_group)
        if (length(levs) >= 2) {
          pw_list <- list()
          combos <- combn(levs, 2, simplify = FALSE)
          for (pair in combos) {
            g1 <- pair[1]
            g2 <- pair[2]
            idx <- which(meta_final$macro_group %in% c(g1, g2))
            if (length(idx) < 6) next
            d_sub <- as.dist(as.matrix(dist_desc)[idx, idx])
            m_sub <- as_base_df(meta_final[idx, , drop = FALSE])
            m_sub$macro_group <- droplevels(factor(m_sub$macro_group))
            if (nlevels(m_sub$macro_group) < 2) next

            try({
              p_sub <- vegan::adonis2(d_sub ~ macro_group, data = m_sub, permutations = 999)
              t_sub <- as_base_df(broom::tidy(p_sub))
              t_sub <- as_base_df(t_sub[t_sub$term == "macro_group", , drop = FALSE])
              if (nrow(t_sub) > 0) {
                t_sub$tag <- "Descriptors_Euclidean"
                t_sub$group1 <- g1
                t_sub$group2 <- g2
                pw_list[[paste(g1, g2)]] <- t_sub
              }
            }, silent = TRUE)
          }
          pw_tbl_d <- bind_rows_df(pw_list)
          if (nrow(pw_tbl_d) > 0 && "p.value" %in% names(pw_tbl_d)) {
            pw_tbl_d$p_adj <- stats::p.adjust(pw_tbl_d$p.value, method = "BH")
          }
        }
      }, error = function(e) NULL)

      res_desc <- list(
        tag = "Descriptors_Euclidean",
        permanova = perma_tbl_d,
        betadisper = bd_tbl_d,
        pairwise = pw_tbl_d,
        ord_pcoa = ord_pcoa_d,
        ord_nmds = empty_tbl()
      )
    }
  }
}

perma_all <- bind_rows_df(res_scaf$permanova, res_class$permanova, res_desc$permanova)
bd_all <- bind_rows_df(res_scaf$betadisper, res_class$betadisper, res_desc$betadisper)
pw_all <- bind_rows_df(res_scaf$pairwise, res_class$pairwise, res_desc$pairwise)

ind_scaf_tbl <- empty_tbl()
ind_class_tbl <- empty_tbl()

if (requireNamespace("indicspecies", quietly = TRUE) && requireNamespace("permute", quietly = TRUE)) {
  if (ncol(X_scaf) > 0) {
    scaf_support <- colSums(as.matrix(X_scaf) > 0, na.rm = TRUE)
    keep_scaf_ind <- names(scaf_support)[scaf_support >= 5]
    if (length(keep_scaf_ind) > 200) {
      keep_scaf_ind <- names(sort(scaf_support[keep_scaf_ind], decreasing = TRUE))[1:200]
    }

    X_scaf_f <- as.matrix(X_scaf[, keep_scaf_ind, drop = FALSE])
    grp_scaf <- meta_fam$macro_group
    if (nrow(X_scaf_f) >= 10 && ncol(X_scaf_f) >= 5 && dplyr::n_distinct(grp_scaf) >= 2) {
      message("[INFO] Running scaffold indicator analysis.")
      ind_scaf <- tryCatch(
        indicspecies::multipatt(X_scaf_f, grp_scaf, func = "r.g", control = permute::how(nperm = 199)),
        error = function(e) NULL
      )

      if (!is.null(ind_scaf)) {
        res <- as_base_df(ind_scaf$sign)
        if (nrow(res) > 0) {
          res$feature <- rownames(res)
          res$tag <- "indicator_scaffold"
          ind_scaf_tbl <- as_base_df(res[order(res$p.value), , drop = FALSE])
        }
      }
    }
  }

  if (ncol(X_class) > 0) {
    class_support <- colSums(as.matrix(X_class) > 0, na.rm = TRUE)
    keep_class_ind <- names(class_support)[class_support >= 3]
    X_class_f <- as.matrix(X_class[, keep_class_ind, drop = FALSE])
    grp_class <- meta_class$macro_group
    if (nrow(X_class_f) >= 10 && ncol(X_class_f) >= 3 && dplyr::n_distinct(grp_class) >= 2) {
      message("[INFO] Running class indicator analysis.")
      ind_class <- tryCatch(
        indicspecies::multipatt(X_class_f, grp_class, func = "r.g", control = permute::how(nperm = 199)),
        error = function(e) NULL
      )

      if (!is.null(ind_class)) {
        res2 <- as_base_df(ind_class$sign)
        if (nrow(res2) > 0) {
          res2$feature <- rownames(res2)
          res2$tag <- "indicator_class"
          ind_class_tbl <- as_base_df(res2[order(res2$p.value), , drop = FALSE])
        }
      }
    }
  }
} else {
  message("[INFO] Optional package 'indicspecies' or 'permute' is not available. Indicator analyses were skipped.")
}

message("[INFO] Loading the final bioactivity workbook.")
chembl <- as_base_df(
  janitor::clean_names(
    as_base_df(
      readxl::read_excel(path_chembl, sheet = "MASTER_DATA", .name_repair = "minimal")
    )
  )
)

if (!"family" %in% names(chembl)) {
  stop("The final bioactivity workbook is missing the required 'family' column in sheet 'MASTER_DATA'.")
}

if (!"val_num" %in% names(chembl)) chembl$val_num <- NA_real_
if (!"standard_value" %in% names(chembl)) chembl$standard_value <- NA_real_
if (!"standard_units" %in% names(chembl)) chembl$standard_units <- NA_character_
if (!"target_category_macro" %in% names(chembl)) chembl$target_category_macro <- NA_character_
if (!"target_category_l3" %in% names(chembl)) chembl$target_category_l3 <- NA_character_
if (!"inchikey" %in% names(chembl)) chembl$inchikey <- NA_character_

if (is.list(chembl$family)) {
  chembl$family <- vapply(chembl$family, function(x) {
    if (length(x) == 0) return(NA_character_)
    paste(as.character(x), collapse = ";")
  }, character(1))
} else {
  chembl$family <- as.character(chembl$family)
}

chembl$family <- stringr::str_trim(chembl$family)
chembl$inchikey <- stringr::str_trim(as.character(chembl$inchikey))
chembl$target_category_macro <- stringr::str_trim(as.character(chembl$target_category_macro))
chembl$target_category_l3 <- stringr::str_trim(as.character(chembl$target_category_l3))
chembl$standard_units <- normalize_unit(chembl$standard_units)

chembl2 <- as_base_df(
  chembl %>%
    dplyr::filter(!is.na(inchikey), inchikey != "", !is.na(family), family != "") %>%
    tidyr::separate_rows(family, sep = ";") %>%
    dplyr::mutate(family = stringr::str_trim(family)) %>%
    dplyr::filter(!is.na(family), family != "") %>%
    dplyr::left_join(groups[, c("family", "macro_group"), drop = FALSE], by = "family") %>%
    dplyr::filter(!is.na(macro_group), macro_group != "") %>%
    dplyr::mutate(
      macro_group = factor(macro_group, levels = macro_levels),
      val_num = safe_num(val_num),
      standard_value = safe_num(standard_value),
      value_nm = dplyr::case_when(
        !is.na(val_num) ~ val_num,
        standard_units == "nm" ~ standard_value,
        standard_units == "um" ~ standard_value * 1000,
        TRUE ~ NA_real_
      ),
      p_activity = ifelse(!is.na(value_nm) & value_nm > 0, -log10(value_nm * 1e-9), NA_real_),
      is_active_10uM = ifelse(!is.na(value_nm) & value_nm <= 10000, TRUE, FALSE)
    )
)

bio_family <- as_base_df(
  chembl2 %>%
    dplyr::group_by(family, macro_group) %>%
    dplyr::summarise(
      n_records = dplyr::n(),
      n_compounds_tested = dplyr::n_distinct(inchikey),
      n_targets_macro = dplyr::n_distinct(target_category_macro),
      n_targets_l3 = dplyr::n_distinct(target_category_l3),
      n_targets = dplyr::n_distinct(target_category_macro),
      n_active = sum(is_active_10uM %in% TRUE, na.rm = TRUE),
      active_rate = mean(is_active_10uM, na.rm = TRUE),
      active_rate_cc = (n_active + 0.5) / (n_records + 1),
      median_p_activity = safe_median_num(p_activity),
      q25_p_activity = safe_quantile_num(p_activity, 0.25),
      q75_p_activity = safe_quantile_num(p_activity, 0.75),
      best_p_activity = safe_max_num(p_activity),
      .groups = "drop"
    )
)

land <- as_base_df(
  family_metrics %>%
    dplyr::left_join(bio_family, by = c("family", "macro_group"))
)

ci <- t(mapply(wilson_ci, land$n_active, land$n_records))
if (nrow(ci) == nrow(land)) {
  land$active_lo <- ci[, 1]
  land$active_hi <- ci[, 2]
} else {
  land$active_lo <- NA_real_
  land$active_hi <- NA_real_
}

df_cor <- as_base_df(
  land %>%
    dplyr::filter(
      !is.na(novelty_ratio),
      !is.na(median_p_activity),
      !is.na(n_records),
      !is.na(n_targets),
      n_records >= 20,
      n_targets >= 3
    )
)

rho <- NA_real_
p_rho <- NA_real_
if (nrow(df_cor) >= 10) {
  rho <- suppressWarnings(stats::cor(df_cor$novelty_ratio, df_cor$median_p_activity, method = "spearman"))
  p_rho <- tryCatch(
    suppressWarnings(stats::cor.test(df_cor$novelty_ratio, df_cor$median_p_activity, method = "spearman")$p.value),
    error = function(e) NA_real_
  )
}

df_pc <- as_base_df(
  land %>%
    dplyr::filter(!is.na(novelty_ratio), !is.na(median_p_activity), !is.na(n_records)) %>%
    dplyr::mutate(log_effort = log10(n_records + 1))
)

rho_partial <- NA_real_
p_partial <- NA_real_
if (nrow(df_pc) >= 20) {
  resid_novelty <- stats::residuals(stats::lm(novelty_ratio ~ log_effort, data = df_pc))
  resid_potency <- stats::residuals(stats::lm(median_p_activity ~ log_effort, data = df_pc))
  rho_partial <- suppressWarnings(stats::cor(resid_novelty, resid_potency, method = "spearman"))
  p_partial <- tryCatch(
    suppressWarnings(stats::cor.test(resid_novelty, resid_potency, method = "spearman")$p.value),
    error = function(e) NA_real_
  )
}

land2 <- as_base_df(
  land %>%
    dplyr::mutate(
      z_novelty = as.numeric(scale(novelty_ratio)),
      z_potency = as.numeric(scale(median_p_activity)),
      tier = dplyr::case_when(
        z_novelty >= 1 & z_potency >= 1 ~ "Priority I: High novelty + High potency",
        z_novelty >= 1 & z_potency < 1 ~ "Hidden Gems: High novelty, lower evidence",
        z_novelty < 1 & z_potency >= 1 ~ "Stars: High evidence, lower novelty",
        TRUE ~ "Lower priority / baseline"
      ),
      score = scales::rescale(novelty_ratio, to = c(0, 1), na.rm = TRUE) +
        scales::rescale(median_p_activity, to = c(0, 1), na.rm = TRUE)
    ) %>%
    dplyr::arrange(desc(score))
)

top_priority <- as_base_df(
  land2 %>%
    dplyr::filter(!is.na(novelty_ratio), !is.na(median_p_activity)) %>%
    dplyr::slice_head(n = 40)
)

plot_metric <- function(df, ycol, ylab, title, log_y = FALSE) {
  df <- as_base_df(
    df %>%
      dplyr::filter(!is.na(macro_group), macro_group != "") %>%
      dplyr::filter(is.finite(.data[[ycol]])) %>%
      dplyr::mutate(macro_group = factor(macro_group, levels = macro_levels))
  )

  if (nrow(df) == 0) {
    return(placeholder_plot(title, "No data were available for this metric."))
  }

  n_df <- as_base_df(df %>% dplyr::group_by(macro_group) %>% dplyr::summarise(n = dplyr::n(), .groups = "drop"))
  lbl_map <- setNames(
    paste0(as.character(n_df$macro_group), "\n(n=", n_df$n, ")"),
    as.character(n_df$macro_group)
  )

  subtitle_text <- "Kruskal-Wallis statistics were not computed because the available data were insufficient."
  cld_pos <- data.frame(macro_group = character(0), y = numeric(0), cld = character(0), stringsAsFactors = FALSE)

  if (nrow(df) >= 12 && dplyr::n_distinct(df$macro_group) >= 2) {
    kw <- as_base_df(rstatix::kruskal_test(df, reformulate("macro_group", ycol)))
    eff <- as_base_df(rstatix::kruskal_effsize(df, reformulate("macro_group", ycol)))
    subtitle_text <- paste0(
      "Kruskal-Wallis p = ", format.pval(kw$p[1], digits = 2),
      " | effect size eps2 = ", round(eff$effsize[1], 3),
      " | Dunn post hoc letters are shown when available."
    )

    dunn <- tryCatch(
      as_base_df(rstatix::dunn_test(df, reformulate("macro_group", ycol), p.adjust.method = "BH")),
      error = function(e) NULL
    )

    if (!is.null(dunn) && nrow(dunn) > 0) {
      letters <- tryCatch(cld_from_dunn(dunn), error = function(e) NULL)
      if (!is.null(letters)) {
        cld_tbl <- data.frame(
          macro_group = factor(names(letters), levels = levels(df$macro_group)),
          cld = unname(letters),
          stringsAsFactors = FALSE
        )
        y_max <- as_base_df(
          df %>%
            dplyr::group_by(macro_group) %>%
            dplyr::summarise(y = max(.data[[ycol]], na.rm = TRUE), .groups = "drop")
        )
        cld_pos <- as_base_df(
          y_max %>%
            dplyr::left_join(cld_tbl, by = "macro_group") %>%
            dplyr::mutate(y = y * 1.08)
        )
      }
    }
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = macro_group, y = .data[[ycol]], fill = macro_group)) +
    ggplot2::geom_violin(width = 0.9, alpha = 0.55, color = NA) +
    ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.9) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.12), size = 0.75, alpha = 0.45) +
    ggplot2::scale_x_discrete(labels = lbl_map, drop = FALSE) +
    ggplot2::scale_fill_manual(values = pal_macro, guide = "none") +
    ggplot2::labs(title = title, subtitle = subtitle_text, x = NULL, y = ylab) +
    theme_pub(10)

  if (nrow(cld_pos) > 0) {
    p <- p + ggplot2::geom_text(
      data = cld_pos,
      ggplot2::aes(x = macro_group, y = y, label = cld),
      inherit.aes = FALSE,
      size = 3.2,
      fontface = "bold"
    )
  }

  if (log_y) {
    p <- p + ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(base = 10))
  }

  p
}

pA <- plot_metric(
  family_metrics,
  "n_compounds",
  "Unique compounds per family (pseudo-log10)",
  "A  Compound richness across macro-groups",
  log_y = TRUE
)

pB <- plot_metric(
  family_metrics,
  "n_scaffolds",
  "Unique scaffolds per family (pseudo-log10)",
  "B  Scaffold richness across macro-groups",
  log_y = TRUE
)

pC <- plot_metric(
  family_metrics,
  "novelty_ratio",
  "Novelty ratio (n_scaffolds / n_compounds)",
  "C  Structural novelty across macro-groups",
  log_y = FALSE
)

pdf(out_fig1, width = 13, height = 10, useDingbats = FALSE)
print(pA)
print(pB)
print(pC)
dev.off()

pcoa_df <- as_base_df(res_scaf$ord_pcoa)
if (nrow(pcoa_df) > 0) {
  pcoa_df$macro_group <- factor(pcoa_df$macro_group, levels = macro_levels)

  perma_line <- as_base_df(
    perma_all %>%
      dplyr::filter(tag == "Scaffold_Jaccard", term == "macro_group") %>%
      dplyr::slice_head(n = 1)
  )

  bd_line <- as_base_df(
    bd_all %>%
      dplyr::filter(tag == "Scaffold_Jaccard") %>%
      dplyr::slice_head(n = 1)
  )

  subtitle2 <- paste0(
    "Distance = Jaccard (binary core14). PERMANOVA R2 = ",
    ifelse(nrow(perma_line) > 0, round(perma_line$r.squared[1], 2), NA),
    ", p = ",
    ifelse(nrow(perma_line) > 0, format.pval(perma_line$p.value[1], digits = 2), NA),
    " | betadisper p = ",
    ifelse(nrow(bd_line) > 0, format.pval(bd_line$p.value[1], digits = 2), NA)
  )

  p2 <- ggplot2::ggplot(pcoa_df, ggplot2::aes(Axis1, Axis2, color = macro_group)) +
    ggplot2::geom_point(size = 2.0, alpha = 0.85) +
    ggplot2::scale_color_manual(values = pal_macro) +
    ggplot2::labs(
      title = "Family-level scaffold composition (PCoA on Jaccard distance)",
      subtitle = subtitle2,
      x = "PCoA1",
      y = "PCoA2",
      color = "Macro-group"
    ) +
    theme_pub(10)

  p2 <- add_group_hulls(p2, pcoa_df, "Axis1", "Axis2", "macro_group", alpha = 0.10)
} else {
  p2 <- placeholder_plot(
    "Family-level scaffold composition (PCoA on Jaccard distance)",
    "Ordination could not be generated from the available scaffold matrix."
  )
}

nmds_df <- as_base_df(res_scaf$ord_nmds)
p2b <- NULL
if (nrow(nmds_df) > 0) {
  nmds_df$macro_group <- factor(nmds_df$macro_group, levels = macro_levels)
  stress_value <- unique(nmds_df$stress)[1]
  p2b <- ggplot2::ggplot(nmds_df, ggplot2::aes(NMDS1, NMDS2, color = macro_group)) +
    ggplot2::geom_point(size = 2.0, alpha = 0.85) +
    ggplot2::scale_color_manual(values = pal_macro) +
    ggplot2::labs(
      title = "Family-level scaffold composition (NMDS on Jaccard distance)",
      subtitle = paste0("Stress = ", round(stress_value, 3), " | Hulls by macro-group"),
      x = "NMDS1",
      y = "NMDS2",
      color = "Macro-group"
    ) +
    theme_pub(10)
  p2b <- add_group_hulls(p2b, nmds_df, "NMDS1", "NMDS2", "macro_group", alpha = 0.10)
}

pdf(out_fig2, width = 12.5, height = 7.5, useDingbats = FALSE)
print(p2)
if (!is.null(p2b)) print(p2b)
dev.off()

subtitle3 <- paste0(
  "Evidence = median pActivity with interquartile ranges. Effort filter for rho: n_records >= 20 and n_targets >= 3. ",
  "Spearman rho = ", ifelse(is.na(rho), "NA", round(rho, 2)),
  ", p = ", ifelse(is.na(p_rho), "NA", format.pval(p_rho, digits = 2)),
  " | residualized rho = ", ifelse(is.na(rho_partial), "NA", round(rho_partial, 2)),
  ", p = ", ifelse(is.na(p_partial), "NA", format.pval(p_partial, digits = 2))
)

land_plot <- as_base_df(
  land2 %>%
    dplyr::mutate(label = ifelse(dplyr::row_number() <= 20, family, NA_character_))
)

if (nrow(land_plot) > 0) {
  p3 <- ggplot2::ggplot(land_plot, ggplot2::aes(x = novelty_ratio, y = median_p_activity, color = macro_group)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = q25_p_activity, ymax = q75_p_activity),
      width = 0,
      alpha = 0.25,
      na.rm = TRUE
    ) +
    ggplot2::geom_point(ggplot2::aes(size = n_compounds), alpha = 0.85, na.rm = TRUE) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = label),
      size = 2.6,
      max.overlaps = 60,
      show.legend = FALSE,
      na.rm = TRUE
    ) +
    ggplot2::scale_color_manual(values = pal_macro) +
    ggplot2::scale_size_continuous(name = "Compounds per family") +
    ggplot2::labs(
      title = "Bioprospecting landscape (family level): novelty versus potency evidence",
      subtitle = subtitle3,
      x = "Structural novelty (scaffold/compound)",
      y = "Potency evidence (median pActivity; IQR)",
      color = "Macro-group"
    ) +
    theme_pub(10)
} else {
  p3 <- placeholder_plot(
    "Bioprospecting landscape (family level): novelty versus potency evidence",
    "No family-level bioactivity records were available for plotting."
  )
}

pdf(out_fig3, width = 12.5, height = 7.5, useDingbats = FALSE)
print(p3)
dev.off()

message("[INFO] Writing the statistical workbook.")
wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, "macro_groups_mapping")
openxlsx::writeData(wb, "macro_groups_mapping", export_df(groups))

openxlsx::addWorksheet(wb, "family_metrics_chem")
openxlsx::writeData(wb, "family_metrics_chem", export_df(family_metrics))

openxlsx::addWorksheet(wb, "macro_metrics_chem")
openxlsx::writeData(wb, "macro_metrics_chem", export_df(macro_metrics))

openxlsx::addWorksheet(wb, "univar_KW_eps2")
openxlsx::writeData(wb, "univar_KW_eps2", export_df(univar_kw_tbl))

openxlsx::addWorksheet(wb, "univar_Dunn_BH")
openxlsx::writeData(wb, "univar_Dunn_BH", export_df(univar_dunn_tbl))

openxlsx::addWorksheet(wb, "univar_CLD")
openxlsx::writeData(wb, "univar_CLD", export_df(univar_cld_tbl))

openxlsx::addWorksheet(wb, "PERMANOVA_all")
openxlsx::writeData(wb, "PERMANOVA_all", export_df(perma_all))

openxlsx::addWorksheet(wb, "BETADISPER_all")
openxlsx::writeData(wb, "BETADISPER_all", export_df(bd_all))

openxlsx::addWorksheet(wb, "PAIRWISE_PERMANOVA")
openxlsx::writeData(wb, "PAIRWISE_PERMANOVA", export_df(pw_all))

openxlsx::addWorksheet(wb, "indicator_scaffolds")
openxlsx::writeData(wb, "indicator_scaffolds", export_df(ind_scaf_tbl))

openxlsx::addWorksheet(wb, "indicator_classes")
openxlsx::writeData(wb, "indicator_classes", export_df(ind_class_tbl))

openxlsx::addWorksheet(wb, "bio_family")
openxlsx::writeData(wb, "bio_family", export_df(bio_family))

openxlsx::addWorksheet(wb, "landscape_family")
openxlsx::writeData(wb, "landscape_family", export_df(land2))

openxlsx::addWorksheet(wb, "top_priority_40")
openxlsx::writeData(wb, "top_priority_40", export_df(top_priority))

openxlsx::addWorksheet(wb, "FIG_paths")
openxlsx::writeData(
  wb,
  "FIG_paths",
  data.frame(
    Fig1 = out_fig1,
    Fig2 = out_fig2,
    Fig3 = out_fig3,
    Excel = out_xlsx,
    stringsAsFactors = FALSE
  )
)

openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)

message("[OK] Part III completed successfully.")
message("Excel: ", out_xlsx)
message("Fig1:  ", out_fig1)
message("Fig2:  ", out_fig2)
message("Fig3:  ", out_fig3)
