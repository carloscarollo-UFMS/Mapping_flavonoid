
`%||%` <- function(a, b) {
  if (is.null(a)) return(b)
  if (length(a) == 0) return(b)
  if (is.character(a) && length(a) == 1 && !nzchar(a)) return(b)
  a
}

as_base_df <- function(x) as.data.frame(x, stringsAsFactors = FALSE)

bind_rows_df <- function(...) as_base_df(dplyr::bind_rows(...))

bind_cols_df <- function(...) as_base_df(dplyr::bind_cols(...))

empty_df <- function(...) data.frame(..., stringsAsFactors = FALSE)

deframe_df <- function(x) {
  x <- as_base_df(x)
  if (ncol(x) < 2) stop("deframe_df requires at least two columns.")
  stats::setNames(x[[2]], x[[1]])
}

if (!exists("runtime")) stop("Object 'runtime' not found. This script must be called by the Main Pipeline.")
OUT_DIR <- runtime$OUT_DIR %||% NA_character_
if (is.na(OUT_DIR) || !nzchar(OUT_DIR) || !dir.exists(OUT_DIR)) stop("OUT_DIR does not exist: ", OUT_DIR)

cfg <- if (exists("cfg")) cfg else list()
cfg$auto_install_pkgs <- cfg$auto_install_pkgs %||% FALSE

tax_col <- cfg$analysis_tax_level %||% "family"
tax_col <- tolower(trimws(as.character(tax_col)))

out_dir <- file.path(OUT_DIR, "PartIII_ALL")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

fig_dir <- out_dir
out_fig_dir <- fig_dir

out_xlsx <- file.path(out_dir, "macro_groups_STATS_MASTER.xlsx")
out_fig1 <- file.path(fig_dir, "Fig1_ChemistryStats.pdf")
out_fig2 <- file.path(fig_dir, "Fig2_OrdinationStats.pdf")
out_fig3 <- file.path(fig_dir, "Fig3_BioprospectingLandscape.pdf")

set.seed(1)

pkgs <- c(
  "readr","readxl","openxlsx","janitor","dplyr","tidyr","stringr","tibble",
  "ggplot2","ggrepel","scales",
  "vegan","rstatix","multcompView","broom","rlang",
  "ComplexHeatmap","circlize","gridExtra","forcats","grid"
)

missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing)) {
  if (isTRUE(cfg$auto_install_pkgs)) {
    install.packages(missing, dependencies = TRUE)
  } else {
    stop(
      "Missing packages: ", paste(missing, collapse = ", "),
      "\nTip: set cfg$auto_install_pkgs <- TRUE in the controller if you want auto-install."
    )
  }
}

suppressPackageStartupMessages({
  library(readr); library(readxl); library(openxlsx); library(janitor)
  library(dplyr); library(tidyr); library(stringr); library(tibble)
  library(ggplot2); library(ggrepel); library(scales)
  library(vegan); library(rstatix); library(multcompView); library(broom)
  library(rlang)
  library(ComplexHeatmap); library(circlize); library(gridExtra); library(forcats); library(grid)
})

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

project_root <- getwd()

path_groups <- resolve_file(
  "groups.xlsx",
  c(
    cfg$path_groups %||% NA_character_,
    file.path(project_root, "inputs", "groups.xlsx"),
    file.path(project_root, "groups.xlsx"),
    file.path(OUT_DIR, "inputs", "groups.xlsx"),
    file.path(OUT_DIR, "groups.xlsx")
  ),
  required = TRUE
)

path_rdkit <- resolve_file(
  "lotus_flavonoids_rdkit_annotations.csv",
  c(
    cfg$path_rdkit %||% NA_character_,
    file.path(project_root, "inputs", "lotus_flavonoids_rdkit_annotations.csv"),
    file.path(project_root, "lotus_flavonoids_rdkit_annotations.csv"),
    file.path(OUT_DIR, "inputs", "lotus_flavonoids_rdkit_annotations.csv"),
    file.path(OUT_DIR, "lotus_flavonoids_rdkit_annotations.csv")
  ),
  required = TRUE
)

path_chembl <- resolve_file(
  "Lotus_Final_Database_v9_Fixed.xlsx",
  c(
    cfg$path_chembl %||% NA_character_,
    file.path(project_root, "inputs", "Lotus_Final_Database_v9_Fixed.xlsx"),
    file.path(project_root, "Lotus_Final_Database_v9_Fixed.xlsx"),
    file.path(OUT_DIR, "inputs", "Lotus_Final_Database_v9_Fixed.xlsx"),
    file.path(OUT_DIR, "Lotus_Final_Database_v9_Fixed.xlsx")
  ),
  required = TRUE
)

path_lotus <- NA_character_

if (!exists("lin_enriched")) {
  message("NOTE: lin_enriched not found in memory. Will try xlsx fallback (path_lotus).")
  
  path_lotus <- resolve_file(
    "LOTUS Excel (optional)",
    c(
      cfg$path_lotus %||% NA_character_,
      file.path(OUT_DIR, paste0(runtime$base_tag %||% "lotus_run", ".xlsx")),
      Sys.glob(file.path(OUT_DIR, "lotus_*.xlsx"))
    ),
    required = FALSE
  )
}

if (!exists("uni_enriched")) {
  message("NOTE: uni_enriched not found in memory. Sections that require uni_enriched will fail unless you load it.")
}


groups <- readxl::read_excel(path_groups, col_names = FALSE) %>%
  dplyr::transmute(
    family      = stringr::str_trim(as.character(.1)),
    color       = stringr::str_trim(as.character(.2)),
    macro_group = stringr::str_trim(as.character(.3))
  ) %>%
  dplyr::filter(!is.na(family), family != "", !is.na(macro_group), macro_group != "") %>%
  dplyr::filter(!tolower(family) %in% c("family", "familia")) %>%
  dplyr::mutate(
    color = ifelse(stringr::str_detect(color, "^#"), color, paste0("#", color)),
    color = stringr::str_replace_all(color, "\\s+", "")
  ) %>%
  dplyr::distinct(family, macro_group, .keep_all = TRUE)

pal_macro <- groups %>%
  dplyr::group_by(macro_group) %>%
  dplyr::summarise(color = dplyr::first(color), .groups = "drop") %>%
  deframe_df()

apg_order <- c(
  "Lycophytes","Ferns","Gymnosperms","ANA","Magnoliids",
  "Monocots","Monocots (Commelinids)","Basal Eudicots","Saxifragales",
  "Rosids (Fabids)","Rosids (Malvids)","Caryophy_Santalales",
  "Asterids (Basal)","Asterids (Lamiids)","Asterids (Campanulids)"
)
present <- unique(groups$macro_group)
macro_levels <- c(apg_order[apg_order %in% present], setdiff(sort(present), apg_order))


`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

get_lotus_lin <- function(sheet = "lin_enriched") {
  if (exists("lin_enriched", inherits = TRUE)) {
    x <- get("lin_enriched", inherits = TRUE)
    if (!inherits(x, "data.frame")) x <- as.data.frame(x)
    return(x)
  }
  
  p <- NA_character_
  if (exists("cfg", inherits = TRUE) && !is.null(cfg$path_lotus) && nzchar(as.character(cfg$path_lotus))) {
    p <- as.character(cfg$path_lotus)
  } else if (exists("path_lotus", inherits = TRUE) && !is.na(path_lotus) && nzchar(as.character(path_lotus))) {
    p <- as.character(path_lotus)
  }
  
  if (is.na(p) || !file.exists(p)) {
    stop(
      "lin_enriched not found in memory and no valid Excel path was provided.\n",
      "Fix: ensure Part I/Auto-Loader loaded 'lin_enriched' OR set cfg$path_lotus to the exported LOTUS Excel."
    )
  }
  
  readxl::read_excel(p, sheet = sheet, .name_repair = "minimal") |>
    janitor::clean_names()
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

wilson_ci <- function(k, n, z = 1.96) {
  if (is.na(k) || is.na(n) || n <= 0) return(c(NA_real_, NA_real_))
  phat <- k / n
  denom <- 1 + z^2 / n
  center <- (phat + z^2/(2*n)) / denom
  half <- (z * sqrt((phat*(1-phat) + z^2/(4*n))/n)) / denom
  c(max(0, center - half), min(1, center + half))
}

safe_num <- function(x) suppressWarnings(as.numeric(gsub("[^0-9.+-eE]", "", as.character(x))))

add_group_hulls <- function(p, df, x, y, grp, alpha = 0.10) {
  gvals <- unique(df[[grp]])
  hull_df <- list()
  for (g in gvals) {
    sub <- df[df[[grp]] == g & is.finite(df[[x]]) & is.finite(df[[y]]), , drop = FALSE]
    if (nrow(sub) < 3) next
    idx <- chull(sub[[x]], sub[[y]])
    hull_df[[as.character(g)]] <- sub[idx, , drop = FALSE]
  }
  hull <- dplyr::bind_rows(hull_df)
  if (nrow(hull) == 0) return(p)
  p + ggplot2::geom_polygon(
    data = hull,
    ggplot2::aes(x = .data[[x]], y = .data[[y]], fill = .data[[grp]]),
    alpha = alpha, color = NA, inherit.aes = FALSE, show.legend = FALSE
  )
}

cld_from_dunn <- function(dunn_df) {
  comps <- paste(dunn_df$group1, dunn_df$group2, sep = "-")
  pvals <- dunn_df$p.adj
  names(pvals) <- comps
  multcompView::multcompLetters(pvals)$Letters
}


apg_order <- c(
  "Lycophytes", "Ferns", "Gymnosperms",
  "ANA",
  "Magnoliids",
  "Monocots", "Monocots (Commelinids)",
  "Basal Eudicots",
  "Saxifragales",
  "Rosids (Fabids)", "Rosids (Malvids)",
  "Caryophy_Santalales",
  "Asterids (Basal)", "Asterids (Lamiids)", "Asterids (Campanulids)"
)

present <- unique(groups$macro_group)
macro_levels <- c(apg_order[apg_order %in% present], setdiff(sort(present), apg_order))


groups <- readxl::read_excel(path_groups, col_names = FALSE) %>%
  transmute(
    family      = str_trim(as.character(...1)),
    color       = str_trim(as.character(...2)),
    macro_group = str_trim(as.character(...3))
  ) %>%
  filter(!is.na(family), family != "", !is.na(macro_group), macro_group != "") %>%
  filter(!tolower(family) %in% c("family", "familia")) %>%
  mutate(
    color = ifelse(str_detect(color, "^#"), color, paste0("#", color)),
    color = str_replace_all(color, "\\s+", "")
  ) %>%
  distinct(family, macro_group, .keep_all = TRUE)

pal_macro <- groups %>%
  group_by(macro_group) %>%
  summarise(color = first(color), .groups = "drop") %>%
  deframe_df()

rdkit <- readr::read_csv(path_rdkit, show_col_types = FALSE, progress = FALSE) %>%
  clean_names() %>%
  mutate(
    family   = str_trim(as.character(family)),
    inchikey = str_trim(as.character(inchikey)),
    core14   = str_trim(as.character(core14))
  ) %>%
  left_join(groups %>% select(family, macro_group), by = "family")

fam_scaf <- rdkit %>%
  filter(!is.na(family), family != "",
         !is.na(macro_group), macro_group != "",
         !is.na(inchikey), inchikey != "",
         !is.na(core14), core14 != "") %>%
  distinct(family, macro_group, core14, inchikey)

family_metrics <- fam_scaf %>%
  group_by(family, macro_group) %>%
  summarise(
    n_compounds = n_distinct(inchikey),
    n_scaffolds = n_distinct(core14),
    novelty_ratio = n_scaffolds / pmax(n_compounds, 1),
    .groups = "drop"
  )

present <- unique(family_metrics$macro_group)
macro_levels <- c(apg_order[apg_order %in% present], setdiff(sort(present), apg_order))
family_metrics <- family_metrics %>%
  mutate(macro_group = factor(macro_group, levels = macro_levels))

theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = base_size - 1),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.title = element_text(face = "bold")
    )
}

plot_metric <- function(df, ycol, ylab, title, log_y = FALSE) {
  
  if (exists("macro_levels")) {
    df <- df %>% mutate(macro_group = factor(macro_group, levels = macro_levels))
  }
  
  df <- df %>%
    filter(!is.na(macro_group), macro_group != "") %>%
    filter(!is.na(.data[[ycol]]))
  
  n_df <- df %>%
    group_by(macro_group) %>%
    summarise(n = n(), .groups = "drop")
  
  lbl_map <- setNames(
    paste0(as.character(n_df$macro_group), "\n(n=", n_df$n, ")"),
    as.character(n_df$macro_group)
  )
  
  subtitle_text <- "Insufficient groups for statistical testing."
  cld_pos <- empty_df(macro_group = character(), y = numeric(), cld = character())
  
  n_groups <- n_distinct(df$macro_group)
  if (n_groups >= 2 && nrow(df) >= 10) {
    
    kw  <- rstatix::kruskal_test(df, reformulate("macro_group", ycol))
    eff <- rstatix::kruskal_effsize(df, reformulate("macro_group", ycol))
    
    subtitle_text <- paste0(
      "Kruskalâ€“Wallis p = ", format.pval(kw$p, digits = 2),
      " | effect size eps2 = ", round(eff$effsize, 3),
      " | Dunn posthoc: BH-FDR (letters shown when available)"
    )
    
    try({
      dunn <- rstatix::dunn_test(df, reformulate("macro_group", ycol), p.adjust.method = "BH")
      
      comps <- paste(dunn$group1, dunn$group2, sep = "-")
      pvals <- dunn$p.adj
      names(pvals) <- comps
      
      letters <- multcompView::multcompLetters(pvals)$Letters
      
      cld_tbl <- empty_df(
        macro_group = names(letters),
        cld = unname(letters)
      ) %>%
        mutate(macro_group = factor(macro_group, levels = levels(df$macro_group)))
      
      y_max <- df %>%
        group_by(macro_group) %>%
        summarise(y = max(.data[[ycol]], na.rm = TRUE), .groups = "drop") %>%
        mutate(macro_group = factor(macro_group, levels = levels(df$macro_group)))
      
      cld_pos <- y_max %>%
        left_join(cld_tbl, by = "macro_group") %>%
        mutate(y = y * 1.08)
    }, silent = TRUE)
  }
  
  p <- ggplot(df, aes(x = macro_group, y = .data[[ycol]], fill = macro_group)) +
    geom_violin(width = 0.9, alpha = 0.55, color = NA) +
    geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.9) +
    geom_point(position = position_jitter(width = 0.12), size = 0.75, alpha = 0.45) +
    scale_x_discrete(labels = lbl_map, drop = FALSE) +
    scale_fill_manual(values = pal_macro, guide = "none") +
    labs(title = title, subtitle = subtitle_text, x = NULL, y = ylab) +
    theme_pub(10)
  
  if (nrow(cld_pos) > 0 && all(c("macro_group","y","cld") %in% names(cld_pos))) {
    p <- p + geom_text(
      data = cld_pos,
      aes(x = macro_group, y = y, label = cld),
      inherit.aes = FALSE,
      size = 3.2,
      fontface = "bold"
    )
  }
  
  if (log_y) {
    p <- p + scale_y_continuous(trans = scales::pseudo_log_trans(base = 10))
  }
  
  return(p)
}

pA <- plot_metric(
  family_metrics, "n_compounds",
  "Unique compounds per family (pseudo-log10)",
  "A  Compound richness across macro-groups", log_y = TRUE
)

pB <- plot_metric(
  family_metrics, "n_scaffolds",
  "Unique scaffolds per family (pseudo-log10)",
  "B  Scaffold richness across macro-groups", log_y = TRUE
)

pC <- plot_metric(
  family_metrics, "novelty_ratio",
  "Novelty ratio (n_scaffolds / n_compounds)",
  "C  Structural novelty across macro-groups", log_y = FALSE
)

pdf(out_fig1, width = 13, height = 10, useDingbats = FALSE)
print(pA)
print(pB)
print(pC)
dev.off()

message("FIG 1 saved to: ", out_fig1)

fam_scaf_mat <- fam_scaf %>%
  mutate(val = 1L) %>%
  select(family, macro_group, core14, val) %>%
  distinct() %>%
  pivot_wider(names_from = core14, values_from = val, values_fill = 0L)

meta_fam <- fam_scaf_mat %>% select(family, macro_group) %>%
  mutate(macro_group = factor(macro_group, levels = macro_levels))
X <- fam_scaf_mat %>% select(-family, -macro_group)

keep_cols_X <- colSums(X, na.rm = TRUE) > 0
X <- X[, keep_cols_X, drop = FALSE]

keep_rows_X <- rowSums(X, na.rm = TRUE) > 0
X <- X[keep_rows_X, , drop = FALSE]
meta_fam <- meta_fam[keep_rows_X, , drop = FALSE]


dist_j <- vegdist(X, method = "jaccard", binary = TRUE)
perma <- adonis2(dist_j ~ macro_group, data = meta_fam, permutations = 999)
r2 <- perma$R2[1]; p_perma <- perma$`Pr(>F)`[1]

bd <- betadisper(dist_j, meta_fam$macro_group)
bd_test <- anova(bd)
p_disp <- bd_test$`Pr(>F)`[1]

pcoa <- cmdscale(dist_j, k = 2, eig = TRUE)
pcoa_df <- as.data.frame(pcoa$points)
colnames(pcoa_df) <- c("PCoA1","PCoA2")
pcoa_df <- bind_cols(meta_fam, pcoa_df)

hulls <- pcoa_df %>%
  group_by(macro_group) %>%
  filter(n() >= 3,
         n_distinct(paste0(round(PCoA1, 10), "_", round(PCoA2, 10))) >= 3) %>%
  slice(chull(PCoA1, PCoA2)) %>%
  ungroup()

ell_ok <- pcoa_df %>%
  group_by(macro_group) %>%
  summarise(n_pts = n(), .groups = "drop") %>%
  filter(n_pts >= 4) %>%
  pull(macro_group)
pcoa_ell <- pcoa_df %>% filter(macro_group %in% ell_ok)

pD <- ggplot(pcoa_df, aes(PCoA1, PCoA2, color = macro_group)) +
  geom_point(size = 2.0, alpha = 0.85) +
  scale_color_manual(values = pal_macro) +
  coord_equal() +
  labs(
    title = "Macro-groups in scaffold composition (Jaccard PCoA)",
    subtitle = paste0(
      "PERMANOVA: RÂ² = ", round(r2, 3), ", p = ", format.pval(p_perma, digits = 2),
      " | betadisper: p = ", format.pval(p_disp, digits = 2),
      " (interpret PERMANOVA with dispersion in mind)"
    ),
    x = "PCoA1", y = "PCoA2", color = "Macro-group"
  ) +
  theme_pub(10)

if (nrow(hulls) > 0) {
  pD <- pD +
    geom_polygon(
      data = hulls,
      aes(x = PCoA1, y = PCoA2, group = macro_group, fill = macro_group),
      inherit.aes = FALSE,
      alpha = 0.10, color = NA
    ) +
    scale_fill_manual(values = pal_macro, guide = "none")
}
if (nrow(pcoa_ell) > 0) {
  pD <- pD +
    stat_ellipse(
      data = pcoa_ell,
      aes(group = macro_group),
      type = "t", level = 0.68,
      linewidth = 0.6, alpha = 0.35
    )
}

pdf(out_fig2, width = 12, height = 6.5, useDingbats = FALSE)
print(pD)
dev.off()


if (!exists("macro_levels")) {
  macro_levels <- sort(unique(groups$macro_group))
}

if (!exists("pal_macro")) {
  pal_macro <- groups %>%
    dplyr::group_by(macro_group) %>%
    dplyr::summarise(color = dplyr::first(color), .groups = "drop") %>%
    deframe_df()
}

if (!exists("family_metrics")) {
  if (!exists("fam_scaf")) {
    stop("family_metrics not found and fam_scaf not found. Build fam_scaf first (family x core14 x inchikey).")
  }
  family_metrics <- fam_scaf %>%
    dplyr::group_by(family, macro_group) %>%
    dplyr::summarise(
      n_compounds = dplyr::n_distinct(inchikey),
      n_scaffolds = dplyr::n_distinct(core14),
      novelty_ratio = n_scaffolds / pmax(n_compounds, 1),
      .groups = "drop"
    ) %>%
    dplyr::mutate(macro_group = factor(macro_group, levels = macro_levels))
}

chembl <- readxl::read_excel(path_chembl, sheet = "MASTER_DATA") %>%
  janitor::clean_names() %>%
  mutate(
    inchikey = str_trim(as.character(inchikey)),
    standard_units = tolower(str_trim(as.character(standard_units))),
    standard_type  = str_trim(as.character(standard_type)),
    target_category_macro = str_trim(as.character(target_category_macro))
  )

if (!("family" %in% names(chembl))) stop("MASTER_DATA is missing required column: family")

if (is.list(chembl$family)) {
  chembl$family <- vapply(chembl$family, function(x) {
    if (length(x) == 0) return(NA_character_)
    paste(as.character(x), collapse = ";")
  }, character(1))
} else {
  chembl$family <- as.character(chembl$family)
}
chembl$family <- str_trim(chembl$family)

chembl2 <- chembl %>%
  filter(!is.na(inchikey), inchikey != "") %>%
  filter(!is.na(family), family != "") %>%
  tidyr::separate_rows(family, sep = ";") %>%
  mutate(family = str_trim(family)) %>%
  filter(!is.na(family), family != "") %>%
  left_join(groups %>% select(family, macro_group) %>% distinct(), by = "family") %>%
  filter(!is.na(macro_group), macro_group != "") %>%
  mutate(macro_group = factor(macro_group, levels = macro_levels))

chembl2 <- chembl2 %>%
  mutate(
    value_nm = case_when(
      !is.na(val_num) ~ as.numeric(val_num),
      standard_units == "nm" ~ as.numeric(standard_value),
      standard_units %in% c("um","Âµm") ~ as.numeric(standard_value) * 1000,
      TRUE ~ NA_real_
    ),
    p_activity = ifelse(!is.na(value_nm) & value_nm > 0, -log10(value_nm * 1e-9), NA_real_),
    is_active_10uM = ifelse(!is.na(value_nm) & value_nm <= 10000, TRUE, FALSE)
  )

bio_family <- chembl2 %>%
  group_by(family, macro_group) %>%
  summarise(
    n_records = n(),
    n_compounds_tested = n_distinct(inchikey),
    n_targets = n_distinct(target_category_macro),
    
    n_active = sum(is_active_10uM %in% TRUE, na.rm = TRUE),
    active_rate = mean(is_active_10uM, na.rm = TRUE),
    
    median_p_activity = suppressWarnings(median(p_activity, na.rm = TRUE)),
    q25_p_activity = suppressWarnings(quantile(p_activity, 0.25, na.rm = TRUE, names = FALSE)),
    q75_p_activity = suppressWarnings(quantile(p_activity, 0.75, na.rm = TRUE, names = FALSE)),
    best_p_activity = suppressWarnings(max(p_activity, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    median_p_activity = ifelse(is.infinite(median_p_activity), NA_real_, median_p_activity),
    best_p_activity   = ifelse(is.infinite(best_p_activity), NA_real_, best_p_activity)
  )

land <- family_metrics %>%
  mutate(macro_group = factor(macro_group, levels = macro_levels)) %>%
  left_join(bio_family, by = c("family", "macro_group"))

ci <- t(mapply(wilson_ci, land$n_active, land$n_records))
land$active_lo <- ci[, 1]
land$active_hi <- ci[, 2]

df_cor <- land %>%
  filter(!is.na(novelty_ratio), !is.na(median_p_activity),
         !is.na(n_records), !is.na(n_targets),
         n_records >= 20, n_targets >= 3)

if (nrow(df_cor) >= 10) {
  rho <- suppressWarnings(cor(df_cor$novelty_ratio, df_cor$median_p_activity, method = "spearman"))
  p_rho <- suppressWarnings(cor.test(df_cor$novelty_ratio, df_cor$median_p_activity, method = "spearman")$p.value)
} else {
  rho <- NA_real_
  p_rho <- NA_real_
}

subtitle_txt <- if (!is.na(rho)) {
  paste0(
    "Evidence = median pActivity per family (IQR bars). Effort filter: n_records>=20, n_targets>=3. ",
    "Spearman rho = ", round(rho, 2), ", p = ", format.pval(p_rho, digits = 2), "."
  )
} else {
  paste0(
    "Evidence = median pActivity per family (IQR bars). Effort filter: n_records>=20, n_targets>=3. ",
    "Insufficient families for Spearman (n<10). Consider relaxing to n_records>=10 & n_targets>=2."
  )
}

message("Families used for Spearman (df_cor): ", nrow(df_cor))

land <- land %>%
  mutate(
    score = scales::rescale(novelty_ratio, to = c(0, 1), na.rm = TRUE) +
      scales::rescale(median_p_activity, to = c(0, 1), na.rm = TRUE)
  ) %>%
  arrange(desc(score)) %>%
  mutate(label = ifelse(row_number() <= 20, family, NA_character_))

pBP <- ggplot(land, aes(x = novelty_ratio, y = median_p_activity, color = macro_group)) +
  geom_errorbar(aes(ymin = q25_p_activity, ymax = q75_p_activity),
                width = 0, alpha = 0.25, na.rm = TRUE) +
  geom_point(aes(size = n_compounds), alpha = 0.85, na.rm = TRUE) +
  ggrepel::geom_text_repel(aes(label = label),
                           size = 2.6, max.overlaps = 40,
                           show.legend = FALSE, na.rm = TRUE) +
  scale_color_manual(values = pal_macro) +
  scale_size_continuous(name = "Compounds per family") +
  labs(
    title = "Bioprospecting landscape (family level): novelty vs potency evidence",
    subtitle = subtitle_txt,
    x = "Structural novelty (scaffold/compound)",
    y = "Potency evidence (median pActivity; IQR)",
    color = "Macro-group"
  ) +
  theme_pub(10) +
  theme(legend.position = "right")

pdf(out_fig3, width = 12.5, height = 7.5, useDingbats = FALSE)
print(pBP)
dev.off()

message("FIG 3 saved to: ", out_fig3)


groups <- readxl::read_excel(path_groups, col_names = FALSE) %>%
  transmute(
    family      = str_trim(as.character(...1)),
    color       = str_trim(as.character(...2)),
    macro_group = str_trim(as.character(...3))
  ) %>%
  filter(!is.na(family), family != "", !is.na(macro_group), macro_group != "") %>%
  distinct(family, .keep_all = TRUE)

groups <- groups %>%
  mutate(color = ifelse(str_detect(color, "^#"), color, paste0("#", color))) %>%
  mutate(color = str_replace_all(color, "\\s+", ""))

macro_levels <- groups %>% distinct(macro_group) %>% pull(macro_group)
pal_macro <- groups %>% distinct(macro_group, color) %>% { setNames(.$color, .$macro_group) }

message("Loading RDKit annotations...")
rdkit <- readr::read_csv(path_rdkit, show_col_types = FALSE) %>%
  janitor::clean_names()

for (nm in c("family","genus","species","inchikey","smiles","class_np","core14")) {
  if (nm %in% names(rdkit)) rdkit[[nm]] <- str_trim(as.character(rdkit[[nm]]))
}

req_rdkit <- c("inchikey","smiles","family","class_np","core14")
miss_rdkit <- setdiff(req_rdkit, names(rdkit))
if (length(miss_rdkit)) stop("RDKit missing required columns: ", paste(miss_rdkit, collapse = ", "))

rdkit <- rdkit %>%
  filter(!is.na(family), family != "") %>%
  left_join(groups %>% select(family, macro_group) %>% distinct(), by = "family") %>%
  filter(!is.na(macro_group), macro_group != "") %>%
  mutate(macro_group = factor(macro_group, levels = macro_levels))

desc_candidates <- c(
  "mw","logp","tpsa","num_hba","num_hbd","num_rings","heavy_atom_number",
  "hba","hbd","rings"
)
bin_candidates <- c(
  "has_phenolic_oh","has_methoxy_aryl","has_prenyl_like","has_conj_carbonyl","has_probable_sugar",
  "contains_sugar"
)

desc_cols <- intersect(desc_candidates, names(rdkit))
bin_cols  <- intersect(bin_candidates, names(rdkit))

if ("hba" %in% names(rdkit) && !("num_hba" %in% names(rdkit))) rdkit$num_hba <- rdkit$hba
if ("hbd" %in% names(rdkit) && !("num_hbd" %in% names(rdkit))) rdkit$num_hbd <- rdkit$hbd
if ("rings" %in% names(rdkit) && !("num_rings" %in% names(rdkit))) rdkit$num_rings <- rdkit$rings

if (!("num_hba" %in% names(rdkit))) rdkit$num_hba <- NA_real_
if (!("num_hbd" %in% names(rdkit))) rdkit$num_hbd <- NA_real_
if (!("num_rings" %in% names(rdkit))) rdkit$num_rings <- NA_real_

for (nm in intersect(c("mw","logp","tpsa","num_hba","num_hbd","num_rings","heavy_atom_number"), names(rdkit))) {
  rdkit[[nm]] <- safe_num(rdkit[[nm]])
}

for (nm in bin_cols) {
  rdkit[[nm]] <- as.logical(rdkit[[nm]] %in% TRUE)
}

lotus_sheet <- "lin_enriched"

if (exists("lin_enriched")) {
  lotus <- janitor::clean_names(lin_enriched)
} else {
  lotus_sheet <- lotus_sheet %||% "lin_enriched"
  if (is.na(path_lotus) || !nzchar(path_lotus) || !file.exists(path_lotus)) {
    stop("lin_enriched not in memory and path_lotus is missing/invalid. ",
         "Load lin_enriched (parquet) or provide LOTUS xlsx.")
  }
  lotus <- readxl::read_excel(path_lotus, sheet = lotus_sheet) %>%
    janitor::clean_names()
}


for (nm in c("family","genus","species","inchikey","smiles","lotus_id")) {
  if (nm %in% names(lotus)) lotus[[nm]] <- str_trim(as.character(lotus[[nm]]))
}
if (!all(c("inchikey","family") %in% names(lotus))) {
  stop("LOTUS lin_enriched must include at least: inchikey, family")
}

lotus <- lotus %>%
  filter(!is.na(family), family != "") %>%
  left_join(groups %>% select(family, macro_group) %>% distinct(), by = "family") %>%
  filter(!is.na(macro_group), macro_group != "") %>%
  mutate(macro_group = factor(macro_group, levels = macro_levels))

message("Extracting physicochemical descriptors from LOTUS to enrich RDKit data...")
lotus_desc_lookup <- lotus %>%
  select(inchikey, 
         any_of(c("molecular_weight", "xlogp", "topo_psa", "heavy_atom_number"))) %>%
  mutate(
    mw_lotus   = if ("molecular_weight" %in% names(.)) safe_num(molecular_weight) else NA_real_,
    logp_lotus = if ("xlogp" %in% names(.)) safe_num(xlogp) else NA_real_,
    tpsa_lotus = if ("topo_psa" %in% names(.)) safe_num(topo_psa) else NA_real_,
    heavy_lotus= if ("heavy_atom_number" %in% names(.)) safe_num(heavy_atom_number) else NA_real_
  ) %>%
  group_by(inchikey) %>%
  summarise(
    mw_lotus    = median(mw_lotus, na.rm=TRUE),
    logp_lotus  = median(logp_lotus, na.rm=TRUE),
    tpsa_lotus  = median(tpsa_lotus, na.rm=TRUE),
    heavy_lotus = median(heavy_lotus, na.rm=TRUE),
    .groups = "drop"
  )

message("Collapsing RDKit to compound-level (unique InChIKey) ...")

rdkit_aug <- rdkit %>%
  left_join(lotus_desc_lookup, by = "inchikey")

chem_compounds <- rdkit_aug %>%
  group_by(inchikey) %>%
  summarise(
    smiles = first(na.omit(smiles)),
    class_np = first(na.omit(class_np)),
    core14   = first(na.omit(core14)),
    
    mw   = if ("mw" %in% names(rdkit)) median(mw, na.rm = TRUE) else median(mw_lotus, na.rm=TRUE),
    logp = if ("logp" %in% names(rdkit)) median(logp, na.rm = TRUE) else median(logp_lotus, na.rm=TRUE),
    tpsa = if ("tpsa" %in% names(rdkit)) median(tpsa, na.rm = TRUE) else median(tpsa_lotus, na.rm=TRUE),
    
    num_hba = median(num_hba, na.rm = TRUE),
    num_hbd = median(num_hbd, na.rm = TRUE),
    num_rings = median(num_rings, na.rm = TRUE),
    
    heavy_atom_number = if ("heavy_atom_number" %in% names(rdkit)) median(heavy_atom_number, na.rm = TRUE) else median(heavy_lotus, na.rm=TRUE),
    
    has_phenolic_oh   = if ("has_phenolic_oh" %in% names(rdkit)) any(has_phenolic_oh %in% TRUE) else NA,
    has_methoxy_aryl  = if ("has_methoxy_aryl" %in% names(rdkit)) any(has_methoxy_aryl %in% TRUE) else NA,
    has_prenyl_like   = if ("has_prenyl_like" %in% names(rdkit)) any(has_prenyl_like %in% TRUE) else NA,
    has_conj_carbonyl = if ("has_conj_carbonyl" %in% names(rdkit)) any(has_conj_carbonyl %in% TRUE) else NA,
    has_probable_sugar = if ("has_probable_sugar" %in% names(rdkit)) any(has_probable_sugar %in% TRUE) else NA,
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~replace(.x, is.infinite(.x), NA_real_)))


fam_comp <- rdkit %>%
  distinct(family, macro_group, inchikey)

fam_scaf <- rdkit %>%
  filter(!is.na(core14), core14 != "") %>%
  distinct(family, macro_group, core14, inchikey)

fam_class <- rdkit %>%
  filter(!is.na(class_np), class_np != "") %>%
  distinct(family, macro_group, class_np, inchikey)

message("Computing family-level chemistry metrics...")

family_metrics <- fam_scaf %>%
  group_by(family, macro_group) %>%
  summarise(
    n_compounds = n_distinct(inchikey),
    n_scaffolds = n_distinct(core14),
    novelty_ratio = n_scaffolds / pmax(n_compounds, 1),
    .groups = "drop"
  )

if (any(c("has_prenyl_like","has_probable_sugar","has_methoxy_aryl","has_phenolic_oh","has_conj_carbonyl") %in% names(chem_compounds))) {
  deco_by_family <- fam_comp %>%
    left_join(chem_compounds %>% select(inchikey, starts_with("has_")), by = "inchikey") %>%
    group_by(family) %>%
    summarise(
      prop_prenyl = mean(has_prenyl_like %in% TRUE, na.rm = TRUE),
      prop_sugar  = mean(has_probable_sugar %in% TRUE, na.rm = TRUE),
      prop_methoxy = mean(has_methoxy_aryl %in% TRUE, na.rm = TRUE),
      prop_phenolic_oh = mean(has_phenolic_oh %in% TRUE, na.rm = TRUE),
      prop_conj_carbonyl = mean(has_conj_carbonyl %in% TRUE, na.rm = TRUE),
      .groups = "drop"
    )
} else {
  deco_by_family <- empty_df(family = unique(family_metrics$family))
}

desc_present <- intersect(c("mw","logp","tpsa","num_hba","num_hbd","num_rings","heavy_atom_number"), names(chem_compounds))

desc_by_family <- fam_comp %>%
  left_join(chem_compounds %>% select(inchikey, all_of(desc_present)), by = "inchikey") %>%
  group_by(family) %>%
  summarise(across(all_of(desc_present), ~median(.x, na.rm = TRUE), .names = "med_{.col}"),
            .groups = "drop") %>%
  mutate(across(where(is.numeric), ~replace(.x, is.infinite(.x), NA_real_)))


family_metrics <- family_metrics %>%
  left_join(deco_by_family, by = "family") %>%
  left_join(desc_by_family, by = "family") %>%
  mutate(macro_group = factor(macro_group, levels = macro_levels))

macro_metrics <- family_metrics %>%
  group_by(macro_group) %>%
  summarise(
    n_families = n_distinct(family),
    n_compounds_median = median(n_compounds, na.rm = TRUE),
    n_compounds_IQR = IQR(n_compounds, na.rm = TRUE),
    n_scaffolds_median = median(n_scaffolds, na.rm = TRUE),
    novelty_ratio_median = median(novelty_ratio, na.rm = TRUE),
    .groups = "drop"
  )

message("Running univariate tests across macro-groups...")

univar_targets <- c(
  "n_compounds","n_scaffolds","novelty_ratio",
  "prop_prenyl","prop_sugar","prop_methoxy","prop_phenolic_oh","prop_conj_carbonyl",
  paste0("med_", desc_present)
)
univar_targets <- univar_targets[univar_targets %in% names(family_metrics)]

univar_kw <- list()
univar_dunn <- list()
univar_cld <- list()

for (m in univar_targets) {
  dfm <- family_metrics %>%
    filter(!is.na(macro_group), macro_group != "") %>%
    select(macro_group, value = all_of(m)) %>%
    filter(is.finite(value))
  
  if (nrow(dfm) < 12 || n_distinct(dfm$macro_group) < 2) next
  
  kw <- rstatix::kruskal_test(dfm, value ~ macro_group)
  eff <- rstatix::kruskal_effsize(dfm, value ~ macro_group)
  kw <- kw %>% mutate(metric = m, eps2 = eff$effsize)
  
  univar_kw[[m]] <- kw
  
  dunn <- rstatix::dunn_test(dfm, value ~ macro_group, p.adjust.method = "BH") %>%
    mutate(metric = m)
  univar_dunn[[m]] <- dunn
  
  letters <- tryCatch(cld_from_dunn(dunn), error = function(e) NULL)
  if (!is.null(letters)) {
    univar_cld[[m]] <- empty_df(metric = m, macro_group = names(letters), cld = unname(letters))
  }
}

univar_kw_tbl   <- bind_rows(univar_kw)
univar_dunn_tbl <- bind_rows(univar_dunn)
univar_cld_tbl  <- bind_rows(univar_cld)
message("Building incidence matrices...")

min_scaffold_family_support <- 3L      # keep scaffolds present in >= this many families
max_scaffold_features       <- 5000L   # hard cap to prevent huge matrices

scaf_freq <- fam_scaf %>%
  distinct(family, core14) %>%
  count(core14, name = "n_families") %>%
  arrange(desc(n_families))

keep_scaf <- scaf_freq %>%
  filter(n_families >= min_scaffold_family_support) %>%
  slice_head(n = max_scaffold_features) %>%
  pull(core14)

fam_scaf_reduced <- fam_scaf %>%
  filter(core14 %in% keep_scaf)

message("Scaffold reduction: kept ", length(keep_scaf), " scaffolds.")

mat_scaf <- fam_scaf_reduced %>%
  mutate(val = 1L) %>%
  select(family, macro_group, core14, val) %>%
  distinct() %>%
  pivot_wider(names_from = core14, values_from = val, values_fill = 0L)

meta_fam <- mat_scaf %>% select(family, macro_group)
X_scaf <- mat_scaf %>% select(-family, -macro_group)

mat_class <- fam_class %>%
  mutate(val = 1L) %>%
  select(family, macro_group, class_np, val) %>%
  distinct() %>%
  pivot_wider(names_from = class_np, values_from = val, values_fill = 0L)

meta_class <- mat_class %>% select(family, macro_group)
X_class <- mat_class %>% select(-family, -macro_group)

med_cols <- grep("^med_", names(family_metrics), value = TRUE)
X_desc <- family_metrics %>%
  select(family, macro_group, all_of(med_cols))

run_permanova_block <- function(X, meta, dist_method, binary = FALSE, tag = "X",
                                permutations = 999, run_nmds = TRUE) {
  
  keep <- complete.cases(meta$macro_group) & !is.na(meta$macro_group)
  X2 <- X[keep, , drop = FALSE]
  meta2 <- meta[keep, , drop = FALSE]
  
  meta2$macro_group <- droplevels(factor(meta2$macro_group))
  
  if (!binary) {
    X2 <- as.data.frame(X2)
    for (nm in names(X2)) X2[[nm]] <- safe_num(X2[[nm]])
    ok_cols <- vapply(X2, function(v) !all(is.na(v)), logical(1))
    X2 <- X2[, ok_cols, drop = FALSE]
    for (nm in names(X2)) {
      v <- X2[[nm]]
      if (any(is.na(v))) X2[[nm]][is.na(v)] <- median(v, na.rm = TRUE)
    }
    X2 <- as.matrix(X2)
  } else {
    X2 <- as.matrix(X2)
  }
  
  if (binary) {
    rs <- rowSums(X2 > 0, na.rm = TRUE)
    keep_rows <- rs > 0
    if (sum(!keep_rows) > 0) {
      message("Block ", tag, ": Dropping ", sum(!keep_rows), " empty rows (no features) to prevent NaN distances.")
      X2 <- X2[keep_rows, , drop = FALSE]
      meta2 <- meta2[keep_rows, , drop = FALSE]
      meta2$macro_group <- droplevels(factor(meta2$macro_group))
    }
  }
  
  if (nrow(X2) < 10 || ncol(X2) < 2 || nlevels(meta2$macro_group) < 2) {
    message("Block ", tag, ": Insufficient data/groups for PERMANOVA. Skipping.")
    return(list(tag=tag, permanova=empty_df(), betadisper=empty_df(), pairwise=empty_df(), ord_pcoa=empty_df(), ord_nmds=empty_df()))
  }
  
  dist_obj <- tryCatch({
    vegan::vegdist(X2, method = dist_method, binary = binary)
  }, error = function(e) {
    message("Error in vegdist for ", tag, ": ", e$message)
    return(NULL)
  })
  
  if (is.null(dist_obj) || any(is.na(dist_obj))) {
    message("Distance matrix contains NA/NaN or failed for ", tag, ". Skipping.")
    return(list(tag=tag, permanova=empty_df(), betadisper=empty_df(), pairwise=empty_df(), ord_pcoa=empty_df(), ord_nmds=empty_df()))
  }
  
  perma_tbl <- empty_df()
  tryCatch({
    perma <- vegan::adonis2(dist_obj ~ macro_group, data = meta2, permutations = permutations)
    perma_tbl <- broom::tidy(perma) %>% mutate(tag = tag)
  }, error = function(e) {
    message("PERMANOVA failed for ", tag, ": ", e$message)
  })
  
  bd_tbl <- empty_df()
  tryCatch({
    bd <- vegan::betadisper(dist_obj, meta2$macro_group)
    bd_tbl <- broom::tidy(anova(bd)) %>% mutate(tag = tag)
  }, error = function(e) {
    message("Betadisper failed for ", tag, ": ", e$message)
  })
  
  pw_tbl <- empty_df()
  tryCatch({
    levs <- levels(meta2$macro_group)
    if (length(levs) >= 2) {
      pw_list <- list()
      combos <- combn(levs, 2, simplify = FALSE)
      for (pair in combos) {
        g1 <- pair[1]; g2 <- pair[2]
        idx <- which(meta2$macro_group %in% c(g1, g2))
        if (length(idx) < 6) next
        
        d_ij <- as.dist(as.matrix(dist_obj)[idx, idx])
        m_ij <- meta2[idx, , drop = FALSE]
        m_ij$macro_group <- droplevels(factor(m_ij$macro_group))
        
        if (nlevels(m_ij$macro_group) < 2) next
        
        try({
          p_ij <- vegan::adonis2(d_ij ~ macro_group, data = m_ij, permutations = permutations)
          t_ij <- broom::tidy(p_ij) %>% filter(term == "macro_group") %>%
            mutate(tag = tag, group1 = g1, group2 = g2)
          pw_list[[paste(g1, g2)]] <- t_ij
        }, silent=TRUE)
      }
      pw_tbl <- bind_rows(pw_list)
      if (nrow(pw_tbl) > 0) pw_tbl <- pw_tbl %>% mutate(p_adj = p.adjust(p.value, method = "BH"))
    }
  }, error = function(e) NULL)
  
  ord_pcoa <- empty_df()
  tryCatch({
    pcoa <- cmdscale(dist_obj, k = 2, eig = TRUE)
    ord_pcoa <- as.data.frame(pcoa$points)
    colnames(ord_pcoa) <- c("Axis1", "Axis2")
    ord_pcoa <- bind_cols(meta2, ord_pcoa) %>% mutate(tag = tag)
  }, error = function(e) NULL)
  
  ord_nmds <- empty_df()
  if (isTRUE(run_nmds) && ncol(X2) < 3000) {
    tryCatch({
      nmds <- vegan::metaMDS(dist_obj, k=2, trymax=20, autotransform=FALSE, trace=FALSE)
      ord_nmds <- as.data.frame(nmds$points)
      colnames(ord_nmds) <- c("NMDS1", "NMDS2")
      ord_nmds <- bind_cols(meta2, ord_nmds) %>% mutate(tag = tag, stress = nmds$stress)
    }, error = function(e) NULL)
  }
  
  list(tag=tag, permanova=perma_tbl, betadisper=bd_tbl, pairwise=pw_tbl, ord_pcoa=ord_pcoa, ord_nmds=ord_nmds)
}


message("Running Scaffolds (Jaccard)...")
res_scaf <- run_permanova_block(X_scaf, meta_fam, dist_method = "jaccard", binary = TRUE, tag = "Scaffold_Jaccard")

message("Running Classes (Jaccard)...")
res_class <- run_permanova_block(X_class, meta_class, dist_method = "jaccard", binary = TRUE, tag = "Class_Jaccard")

message("Running Descriptors (Euclidean - Manual Fix)...")

meta_desc <- X_desc %>% select(family, macro_group)
Xd <- X_desc %>% select(-family, -macro_group)

Xd2 <- as.data.frame(Xd)
for (nm in names(Xd2)) Xd2[[nm]] <- safe_num(Xd2[[nm]])

ok_cols <- vapply(Xd2, function(v) !all(is.na(v)), logical(1))
Xd2 <- Xd2[, ok_cols, drop = FALSE]

for (nm in names(Xd2)) {
  v <- Xd2[[nm]]; if(any(is.na(v))) Xd2[[nm]][is.na(v)] <- median(v, na.rm=TRUE)
}

vars <- vapply(Xd2, function(x) var(x, na.rm=TRUE), numeric(1))
keep_var <- !is.na(vars) & vars > 1e-12 # Strict variance check
if(sum(!keep_var) > 0) message("Dropped ", sum(!keep_var), " constant descriptor columns.")
Xd2 <- Xd2[, keep_var, drop = FALSE]

if (ncol(Xd2) < 1) {
  message("WARNING: No valid descriptor columns remain (feature space empty). Skipping Euclidean analysis.")
  res_desc <- list(
    tag = "Descriptors_Euclidean",
    permanova = empty_df(term = "SKIPPED_NO_DATA", p.value = NA),
    betadisper = empty_df(),
    pairwise = empty_df(),
    ord_pcoa = empty_df(),
    ord_nmds = empty_df()
  )
} else {
  Xd2 <- scale(as.matrix(Xd2))
  Xd2[!is.finite(Xd2)] <- 0
  
  keep_rows <- !is.na(meta_desc$macro_group)
  Xd_final <- Xd2[keep_rows, , drop = FALSE]
  meta_final <- meta_desc[keep_rows, , drop = FALSE]
  
  meta_final$macro_group <- droplevels(factor(meta_final$macro_group))
  
  if (nrow(Xd_final) < 5 || nlevels(meta_final$macro_group) < 2) {
    message("WARNING: Not enough rows or groups for Descriptors analysis.")
    res_desc <- list(tag = "Descriptors_Euclidean", permanova=empty_df(), betadisper=empty_df(), pairwise=empty_df(), ord_pcoa=empty_df(), ord_nmds=empty_df())
  } else {
    dist_desc <- tryCatch({
      d <- vegan::vegdist(Xd_final, method = "euclidean")
      if (any(!is.finite(d))) stop("Distance matrix contains non-finite values")
      d
    }, error = function(e) {
      message("ERROR: vegdist failed on descriptors: ", e$message)
      NULL
    })
    
    if (is.null(dist_desc)) {
      message("Skipping Descriptors stats due to distance failure.")
      res_desc <- list(tag="Descriptors_Euclidean", permanova=empty_df(), betadisper=empty_df(), pairwise=empty_df(), ord_pcoa=empty_df(), ord_nmds=empty_df())
    } else {
      perma_tbl_d <- empty_df()
      tryCatch({
        perma_desc <- vegan::adonis2(dist_desc ~ macro_group, data = meta_final, permutations = 999)
        perma_tbl_d <- broom::tidy(perma_desc) %>% mutate(tag = "Descriptors_Euclidean")
      }, error = function(e) message("Descriptors PERMANOVA failed: ", e$message))
      
      bd_tbl_d <- empty_df()
      tryCatch({
        bd_desc <- vegan::betadisper(dist_desc, meta_final$macro_group)
        bd_tbl_d <- broom::tidy(anova(bd_desc)) %>% mutate(tag = "Descriptors_Euclidean")
      }, error = function(e) message("Descriptors Betadisper failed: ", e$message))
      
      ord_pcoa_d <- empty_df()
      tryCatch({
        pcoa_desc <- cmdscale(dist_desc, k = 2, eig = TRUE)
        pts <- as.data.frame(pcoa_desc$points)
        if (ncol(pts) >= 2) {
          colnames(pts) <- c("Axis1", "Axis2")
          ord_pcoa_d <- bind_cols(meta_final, pts) %>% mutate(tag = "Descriptors_Euclidean")
        }
      }, error = function(e) message("Descriptors PCoA failed: ", e$message))
      
      pw_tbl_d <- empty_df()
      tryCatch({
        levs <- levels(meta_final$macro_group)
        if(length(levs) >= 2) {
          pw_list <- list()
          combos <- combn(levs, 2, simplify = FALSE)
          for(pair in combos) {
            g1 <- pair[1]; g2 <- pair[2]
            idx <- which(meta_final$macro_group %in% c(g1, g2))
            if(length(idx) < 6) next
            d_sub <- as.dist(as.matrix(dist_desc)[idx, idx])
            m_sub <- meta_final[idx, , drop = FALSE]
            m_sub$macro_group <- droplevels(factor(m_sub$macro_group))
            if (nlevels(m_sub$macro_group) < 2) next
            
            try({
              p_sub <- vegan::adonis2(d_sub ~ macro_group, data = m_sub, permutations = 999)
              t_sub <- broom::tidy(p_sub) %>% filter(term == "macro_group") %>%
                mutate(tag = "Descriptors_Euclidean", group1 = g1, group2 = g2)
              pw_list[[paste(g1, g2)]] <- t_sub
            }, silent=TRUE)
          }
          pw_tbl_d <- bind_rows(pw_list)
          if(nrow(pw_tbl_d)>0) pw_tbl_d <- pw_tbl_d %>% mutate(p_adj = p.adjust(p.value, method = "BH"))
        }
      }, error = function(e) message("Descriptors Pairwise failed: ", e$message))
      
      res_desc <- list(
        tag = "Descriptors_Euclidean",
        permanova = perma_tbl_d, betadisper = bd_tbl_d, pairwise = pw_tbl_d,
        ord_pcoa = ord_pcoa_d, ord_nmds = empty_df()
      )
    }
  }
}

perma_all <- bind_rows(res_scaf$permanova, res_class$permanova, res_desc$permanova)
bd_all    <- bind_rows(res_scaf$betadisper, res_class$betadisper, res_desc$betadisper)
pw_all    <- bind_rows(res_scaf$pairwise,   res_class$pairwise,   res_desc$pairwise)

message("Block 8 completed successfully.")
ind_scaf_tbl <- empty_df()
ind_class_tbl <- empty_df()

if (requireNamespace("indicspecies", quietly = TRUE)) {
  suppressPackageStartupMessages(library(indicspecies))
  
  if (exists("X_scaf") && ncol(X_scaf) > 0) {
    scaf_freq <- colSums(X_scaf > 0)
    keep_scaf <- names(scaf_freq)[scaf_freq >= 5]
    
    if (length(keep_scaf) > 200) {
      keep_scaf <- names(sort(scaf_freq[keep_scaf], decreasing = TRUE))[1:200]
    }
    
    X_scaf_f  <- as.matrix(X_scaf[, keep_scaf, drop = FALSE])
    grp <- meta_fam$macro_group
    
    if (nrow(X_scaf_f) >= 10 && ncol(X_scaf_f) >= 5 && n_distinct(grp) >= 2) {
      message("Indicator analysis: scaffolds (top ", ncol(X_scaf_f), " features, 199 perms)...")
      ind_scaf <- tryCatch(
        indicspecies::multipatt(X_scaf_f, grp, func = "r.g", control = how(nperm = 199)),
        error = function(e) NULL
      )
      
      if (!is.null(ind_scaf)) {
        res <- as.data.frame(ind_scaf$sign)
        if (nrow(res) > 0) {
          res$feature <- rownames(res)
          ind_scaf_tbl <- res %>%
            as_base_df() %>%
            arrange(p.value) %>%
            mutate(tag = "indicator_scaffold")
        }
      }
    }
  }
  
  if (exists("X_class") && ncol(X_class) > 0) {
    class_freq <- colSums(X_class > 0)
    keep_class <- names(class_freq)[class_freq >= 3]
    X_class_f  <- as.matrix(X_class[, keep_class, drop = FALSE])
    grp2 <- meta_class$macro_group
    
    if (nrow(X_class_f) >= 10 && ncol(X_class_f) >= 3 && n_distinct(grp2) >= 2) {
      message("Indicator analysis: classes (", ncol(X_class_f), " features, 199 perms)...")
      ind_cls <- tryCatch(
        indicspecies::multipatt(X_class_f, grp2, func = "r.g", control = how(nperm = 199)),
        error = function(e) NULL
      )
      
      if (!is.null(ind_cls)) {
        res2 <- as.data.frame(ind_cls$sign)
        if (nrow(res2) > 0) {
          res2$feature <- rownames(res2)
          ind_class_tbl <- res2 %>%
            as_base_df() %>%
            arrange(p.value) %>%
            mutate(tag = "indicator_class")
        }
      }
    }
  }
  
} else {
  message("Package 'indicspecies' not installed; skipping indicator analyses.")
}

message("Loading ChEMBL MASTER_DATA...")
chembl <- readxl::read_excel(path_chembl, sheet = "MASTER_DATA") %>%
  janitor::clean_names() %>%
  mutate(
    inchikey = str_trim(as.character(inchikey)),
    standard_units = tolower(str_trim(as.character(standard_units))),
    standard_type  = str_trim(as.character(standard_type)),
    target_category_macro = str_trim(as.character(target_category_macro)),
    target_category_l3    = str_trim(as.character(target_category_l3))
  )

if (!("family" %in% names(chembl))) stop("MASTER_DATA is missing required column: family")

if (is.list(chembl$family)) {
  chembl$family <- vapply(chembl$family, function(x) {
    if (length(x) == 0) return(NA_character_)
    paste(as.character(x), collapse = ";")
  }, character(1))
} else {
  chembl$family <- as.character(chembl$family)
}
chembl$family <- str_trim(chembl$family)

chembl2 <- chembl %>%
  filter(!is.na(inchikey), inchikey != "") %>%
  filter(!is.na(family), family != "") %>%
  tidyr::separate_rows(family, sep = ";") %>%
  mutate(family = str_trim(family)) %>%
  filter(!is.na(family), family != "") %>%
  left_join(groups %>% select(family, macro_group) %>% distinct(), by = "family") %>%
  filter(!is.na(macro_group), macro_group != "") %>%
  mutate(macro_group = factor(macro_group, levels = macro_levels)) %>%
  mutate(
    val_num = safe_num(val_num),
    standard_value = safe_num(standard_value),
    value_nm = case_when(
      !is.na(val_num) ~ val_num,
      standard_units == "nm" ~ standard_value,
      standard_units %in% c("um","Âµm") ~ standard_value * 1000,
      TRUE ~ NA_real_
    ),
    p_activity = ifelse(!is.na(value_nm) & value_nm > 0, -log10(value_nm * 1e-9), NA_real_),
    is_active_10uM = ifelse(!is.na(value_nm) & value_nm <= 10000, TRUE, FALSE)
  )

bio_family <- chembl2 %>%
  group_by(family, macro_group) %>%
  summarise(
    n_records = n(),
    n_compounds_tested = n_distinct(inchikey),
    n_targets_macro = n_distinct(target_category_macro),
    n_targets_l3    = n_distinct(target_category_l3),
    n_targets       = n_distinct(target_category_macro),
    n_active = sum(is_active_10uM %in% TRUE, na.rm = TRUE),
    active_rate = mean(is_active_10uM, na.rm = TRUE),
    active_rate_cc = (n_active + 0.5) / (n_records + 1),
    median_p_activity = suppressWarnings(median(p_activity, na.rm = TRUE)),
    q25_p_activity = suppressWarnings(quantile(p_activity, 0.25, na.rm = TRUE, names = FALSE)),
    q75_p_activity = suppressWarnings(quantile(p_activity, 0.75, na.rm = TRUE, names = FALSE)),
    best_p_activity = suppressWarnings(max(p_activity, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    median_p_activity = ifelse(is.infinite(median_p_activity), NA_real_, median_p_activity),
    best_p_activity   = ifelse(is.infinite(best_p_activity), NA_real_, best_p_activity)
  )

land <- family_metrics %>%
  left_join(bio_family, by = c("family","macro_group"))

df_cor <- land %>%
  filter(!is.na(novelty_ratio), !is.na(median_p_activity),
         !is.na(n_records), !is.na(n_targets),
         n_records >= 20, n_targets >= 3)

rho <- NA_real_; p_rho <- NA_real_
if (nrow(df_cor) >= 10) {
  rho <- suppressWarnings(cor(df_cor$novelty_ratio, df_cor$median_p_activity, method = "spearman"))
  p_rho <- suppressWarnings(cor.test(df_cor$novelty_ratio, df_cor$median_p_activity, method = "spearman")$p.value)
}

df_pc <- land %>%
  filter(!is.na(novelty_ratio), !is.na(median_p_activity), !is.na(n_records)) %>%
  mutate(log_effort = log10(n_records + 1))

rho_partial <- NA_real_; p_partial <- NA_real_
if (nrow(df_pc) >= 20) {
  r1 <- residuals(lm(novelty_ratio ~ log_effort, data = df_pc))
  r2 <- residuals(lm(median_p_activity ~ log_effort, data = df_pc))
  rho_partial <- suppressWarnings(cor(r1, r2, method = "spearman"))
  p_partial <- suppressWarnings(cor.test(r1, r2, method = "spearman")$p.value)
}

land2 <- land %>%
  mutate(
    z_novelty = as.numeric(scale(novelty_ratio)),
    z_potency = as.numeric(scale(median_p_activity)),
    tier = case_when(
      z_novelty >= 1 & z_potency >= 1 ~ "Priority I: High novelty + High potency",
      z_novelty >= 1 & z_potency <  1 ~ "Hidden Gems: High novelty, lower evidence",
      z_novelty <  1 & z_potency >= 1 ~ "Stars: High evidence, lower novelty",
      TRUE ~ "Lower priority / baseline"
    ),
    score = scales::rescale(novelty_ratio, to = c(0,1), na.rm = TRUE) +
      scales::rescale(median_p_activity, to = c(0,1), na.rm = TRUE)
  ) %>%
  arrange(desc(score))

top_priority <- land2 %>%
  filter(!is.na(novelty_ratio), !is.na(median_p_activity)) %>%
  slice_head(n = 40)


plot_metric <- function(df, ycol, ylab, title, log_y = FALSE) {
  
  df <- df %>%
    filter(!is.na(macro_group), macro_group != "") %>%
    filter(is.finite(.data[[ycol]])) %>%
    mutate(macro_group = factor(macro_group, levels = macro_levels))
  
  n_df <- df %>% group_by(macro_group) %>% summarise(n = n(), .groups = "drop")
  lbl_map <- setNames(
    paste0(as.character(n_df$macro_group), "\n(n=", n_df$n, ")"),
    as.character(n_df$macro_group)
  )
  
  subtitle_text <- "Kruskalâ€“Wallis: not computed (insufficient data)."
  cld_pos <- empty_df(macro_group = character(), y = numeric(), cld = character())
  
  if (nrow(df) >= 12 && n_distinct(df$macro_group) >= 2) {
    kw  <- rstatix::kruskal_test(df, reformulate("macro_group", ycol))
    eff <- rstatix::kruskal_effsize(df, reformulate("macro_group", ycol))
    subtitle_text <- paste0(
      "Kruskalâ€“Wallis p = ", format.pval(kw$p, digits = 2),
      " | effect size eps2 = ", round(eff$effsize, 3),
      " | Dunn posthoc: BH-FDR (letters shown when available)"
    )
    
    dunn <- tryCatch(
      rstatix::dunn_test(df, reformulate("macro_group", ycol), p.adjust.method = "BH"),
      error = function(e) NULL
    )
    if (!is.null(dunn) && nrow(dunn) > 0) {
      letters <- tryCatch(cld_from_dunn(dunn), error = function(e) NULL)
      if (!is.null(letters)) {
        cld_tbl <- empty_df(macro_group = names(letters), cld = unname(letters)) %>%
          mutate(macro_group = factor(macro_group, levels = levels(df$macro_group)))
        y_max <- df %>% group_by(macro_group) %>%
          summarise(y = max(.data[[ycol]], na.rm = TRUE), .groups = "drop")
        cld_pos <- y_max %>%
          left_join(cld_tbl, by = "macro_group") %>%
          mutate(y = y * 1.08)
      }
    }
  }
  
  p <- ggplot(df, aes(x = macro_group, y = .data[[ycol]], fill = macro_group)) +
    geom_violin(width = 0.9, alpha = 0.55, color = NA) +
    geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.9) +
    geom_point(position = position_jitter(width = 0.12), size = 0.75, alpha = 0.45) +
    scale_x_discrete(labels = lbl_map, drop = FALSE) +
    scale_fill_manual(values = pal_macro, guide = "none") +
    labs(title = title, subtitle = subtitle_text, x = NULL, y = ylab) +
    theme_pub(10)
  
  if (nrow(cld_pos) > 0) {
    p <- p + geom_text(
      data = cld_pos,
      aes(x = macro_group, y = y, label = cld),
      inherit.aes = FALSE,
      size = 3.2,
      fontface = "bold"
    )
  }
  
  if (log_y) p <- p + scale_y_continuous(trans = scales::pseudo_log_trans(base = 10))
  p
}

pA <- plot_metric(family_metrics, "n_compounds",
                  "Unique compounds per family (pseudo-log10)",
                  "A  Compound richness across macro-groups", log_y = TRUE)

pB <- plot_metric(family_metrics, "n_scaffolds",
                  "Unique scaffolds per family (pseudo-log10)",
                  "B  Scaffold richness across macro-groups", log_y = TRUE)

pC <- plot_metric(family_metrics, "novelty_ratio",
                  "Novelty ratio (n_scaffolds / n_compounds)",
                  "C  Structural novelty across macro-groups", log_y = FALSE)

pdf(out_fig1, width = 13, height = 10, useDingbats = FALSE)
print(pA); print(pB); print(pC)
dev.off()

pcoa_df <- res_scaf$ord_pcoa %>%
  mutate(macro_group = factor(macro_group, levels = macro_levels))

perma_line <- perma_all %>%
  filter(tag == "Scaffold_Jaccard", term == "macro_group") %>%
  slice_head(n = 1)

bd_line <- bd_all %>%
  filter(tag == "Scaffold_Jaccard") %>%
  slice_head(n = 1)

subtitle2 <- paste0(
  "Distance: Jaccard (binary core14). PERMANOVA R2=",
  ifelse(nrow(perma_line)>0, round(perma_line$r.squared, 2), NA),
  ", p=", ifelse(nrow(perma_line)>0, format.pval(perma_line$p.value, digits = 2), NA),
  " | betadisper p=", ifelse(nrow(bd_line)>0, format.pval(bd_line$p.value, digits = 2), NA)
)

p2 <- ggplot(pcoa_df, aes(Axis1, Axis2, color = macro_group)) +
  geom_point(size = 2.0, alpha = 0.85) +
  scale_color_manual(values = pal_macro) +
  labs(
    title = "Family-level scaffold composition (PCoA on Jaccard distance)",
    subtitle = subtitle2,
    x = "PCoA1", y = "PCoA2", color = "Macro-group"
  ) +
  theme_pub(10)

p2 <- add_group_hulls(p2, pcoa_df, "Axis1", "Axis2", "macro_group", alpha = 0.10)

nmds_df <- res_scaf$ord_nmds
p2b <- NULL
if (nrow(nmds_df) > 0) {
  nmds_df <- nmds_df %>% mutate(macro_group = factor(macro_group, levels = macro_levels))
  stress <- unique(nmds_df$stress)[1]
  p2b <- ggplot(nmds_df, aes(NMDS1, NMDS2, color = macro_group)) +
    geom_point(size = 2.0, alpha = 0.85) +
    scale_color_manual(values = pal_macro) +
    labs(
      title = "Family-level scaffold composition (NMDS on Jaccard distance)",
      subtitle = paste0("Stress=", round(stress, 3), " | Hulls per macro-group"),
      x = "NMDS1", y = "NMDS2", color = "Macro-group"
    ) +
    theme_pub(10)
  p2b <- add_group_hulls(p2b, nmds_df, "NMDS1", "NMDS2", "macro_group", alpha = 0.10)
}

pdf(out_fig2, width = 12.5, height = 7.5, useDingbats = FALSE)
print(p2)
if (!is.null(p2b)) print(p2b)
dev.off()

subtitle3 <- paste0(
  "Evidence=median pActivity (IQR). Effort filter (rho): n_records>=20 & n_targets>=3. ",
  "Spearman rho=", ifelse(is.na(rho), "NA", round(rho, 2)),
  ", p=", ifelse(is.na(p_rho), "NA", format.pval(p_rho, digits = 2)),
  " | Effort-controlled residual rho=", ifelse(is.na(rho_partial), "NA", round(rho_partial, 2)),
  ", p=", ifelse(is.na(p_partial), "NA", format.pval(p_partial, digits = 2))
)

land_plot <- land2 %>%
  mutate(label = ifelse(row_number() <= 20, family, NA_character_))

p3 <- ggplot(land_plot, aes(x = novelty_ratio, y = median_p_activity, color = macro_group)) +
  geom_errorbar(aes(ymin = q25_p_activity, ymax = q75_p_activity),
                width = 0, alpha = 0.25, na.rm = TRUE) +
  geom_point(aes(size = n_compounds), alpha = 0.85, na.rm = TRUE) +
  ggrepel::geom_text_repel(aes(label = label),
                           size = 2.6, max.overlaps = 60,
                           show.legend = FALSE, na.rm = TRUE) +
  scale_color_manual(values = pal_macro) +
  scale_size_continuous(name = "Compounds per family") +
  labs(
    title = "Bioprospecting landscape (family level): novelty vs potency evidence",
    subtitle = subtitle3,
    x = "Structural novelty (scaffold/compound)",
    y = "Potency evidence (median pActivity; IQR)",
    color = "Macro-group"
  ) +
  theme_pub(10)

pdf(out_fig3, width = 12.5, height = 7.5, useDingbats = FALSE)
print(p3)
dev.off()

message("Writing Excel workbook...")
wb <- createWorkbook()

addWorksheet(wb, "macro_groups_mapping"); writeData(wb, "macro_groups_mapping", groups)
addWorksheet(wb, "family_metrics_chem");   writeData(wb, "family_metrics_chem", family_metrics)
addWorksheet(wb, "macro_metrics_chem");    writeData(wb, "macro_metrics_chem", macro_metrics)

addWorksheet(wb, "univar_KW_eps2");        writeData(wb, "univar_KW_eps2", univar_kw_tbl)
addWorksheet(wb, "univar_Dunn_BH");        writeData(wb, "univar_Dunn_BH", univar_dunn_tbl)
addWorksheet(wb, "univar_CLD");            writeData(wb, "univar_CLD", univar_cld_tbl)

addWorksheet(wb, "PERMANOVA_all");         writeData(wb, "PERMANOVA_all", perma_all)
addWorksheet(wb, "BETADISPER_all");        writeData(wb, "BETADISPER_all", bd_all)
addWorksheet(wb, "PAIRWISE_PERMANOVA");    writeData(wb, "PAIRWISE_PERMANOVA", pw_all)

addWorksheet(wb, "indicator_scaffolds");   writeData(wb, "indicator_scaffolds", ind_scaf_tbl)
addWorksheet(wb, "indicator_classes");     writeData(wb, "indicator_classes", ind_class_tbl)

addWorksheet(wb, "bio_family");            writeData(wb, "bio_family", bio_family)
addWorksheet(wb, "landscape_family");      writeData(wb, "landscape_family", land2)
addWorksheet(wb, "top_priority_40");       writeData(wb, "top_priority_40", top_priority)

addWorksheet(wb, "FIG_paths");             writeData(wb, "FIG_paths", empty_df(
  Fig1 = out_fig1, Fig2 = out_fig2, Fig3 = out_fig3, Excel = out_xlsx
))

saveWorkbook(wb, out_xlsx, overwrite = TRUE)

message("DONE.")
message("Excel: ", out_xlsx)
message("Fig1:  ", out_fig1)
message("Fig2:  ", out_fig2)
message("Fig3:  ", out_fig3)


out_dir <- file.path(OUT_DIR, "PartIII_ALL")
fig_dir <- out_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

xlsx_path <- file.path(out_dir, "macro_groups_STATS_MASTER.xlsx")


mapping <- read_excel(xlsx_path, sheet = "macro_groups_mapping")
pal_macro <- setNames(mapping$color, mapping$macro_group)

land <- read_excel(xlsx_path, sheet = "landscape_family")

ind_class <- read_excel(xlsx_path, sheet = "indicator_classes")

perma <- read_excel(xlsx_path, sheet = "PERMANOVA_all")

theme_nature <- function(base_size = 7) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(family = "sans", color = "black"),
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle = element_text(size = base_size, color = "gray30", margin = margin(b = 10)),
      axis.title = element_text(face = "bold", size = base_size),
      axis.text = element_text(color = "black", size = base_size),
      legend.position = "right",
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray92", size = 0.3),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      strip.text = element_text(face = "bold", size = base_size)
    )
}

message("Generating Figure 5 (radar plot)...")

z_cut <- 0.8

land <- land %>%
  mutate(
    is_priority = z_novelty > z_cut & z_potency > z_cut,
    label_text = ifelse(is_priority | n_compounds > 500, family, NA_character_)
  )

p_landscape <- ggplot(land, aes(x = z_novelty, y = z_potency)) +
  annotate("rect", xmin = z_cut, xmax = Inf, ymin = z_cut, ymax = Inf, 
           fill = "#E41A1C", alpha = 0.03) +
  
  geom_hline(yintercept = 0, linetype = "solid", color = "gray80", size = 0.3) +
  geom_vline(xintercept = 0, linetype = "solid", color = "gray80", size = 0.3) +
  geom_hline(yintercept = z_cut, linetype = "dashed", color = "gray50", size = 0.3) +
  geom_vline(xintercept = z_cut, linetype = "dashed", color = "gray50", size = 0.3) +
  
  geom_point(aes(color = macro_group, size = n_compounds), alpha = 0.75, stroke = 0.2) +
  
  annotate("text", x = max(land$z_novelty, na.rm=TRUE), y = max(land$z_potency, na.rm=TRUE), 
           label = "Priority Zone\n(Novelty + Potency)", 
           color = "#E41A1C", size = 2.5, hjust = 1, vjust = 1, fontface = "bold") +
  
  geom_text_repel(
    aes(label = label_text),
    size = 2.2,
    min.segment.length = 0,
    box.padding = 0.3,
    point.padding = 0.3,
    max.overlaps = 50,
    segment.color = "gray60",
    segment.size = 0.2,
    color = "black",
    bg.color = "white",
    bg.r = 0.15
  ) +
  
  scale_color_manual(values = pal_macro, name = "Lineage") +
  scale_size_continuous(range = c(1, 7), name = "Reported\nCompounds", breaks = c(10, 100, 1000)) +
  
  labs(
    title = "Bioprospecting Radar: Structural Novelty vs. Bioactive Potency",
    subtitle = "Family-level analysis. Top-right quadrant indicates high novelty and high historical potency.",
    x = "Structural Novelty (Z-score)",
    y = "Bioactivity Potency (Z-score of median pActivity)"
  ) +
  theme_nature() +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))

ggsave(file.path(out_fig_dir, "Fig5_Bioprospecting_Radar.pdf"), p_landscape, width = 180, height = 130, units = "mm")


message("Generating Figure 6 (fingerprint plot)...")

ind_long <- ind_class %>%
  filter(p.value < 0.1) %>% 
  pivot_longer(cols = starts_with("s."), names_to = "group_col", values_to = "is_associated") %>%
  filter(is_associated == 1) %>%
  mutate(macro_group = str_remove(group_col, "s\\."))

top_n_per_group <- 8 
top_features <- ind_long %>%
  group_by(macro_group) %>%
  arrange(desc(stat)) %>%
  slice_head(n = top_n_per_group) %>%
  ungroup()

if (nrow(top_features) < 5) {
  message("Warning: Few significant classes were detected. Relaxing filters...")
  ind_long <- ind_class %>%
    pivot_longer(cols = starts_with("s."), names_to = "group_col", values_to = "is_associated") %>%
    filter(is_associated == 1) %>%
    mutate(macro_group = str_remove(group_col, "s\\."))
  
  top_features <- ind_long %>%
    group_by(macro_group) %>%
    arrange(desc(stat)) %>%
    slice_head(n = top_n_per_group) %>%
    ungroup()
}

grp_levels <- sort(unique(top_features$macro_group))
top_features$macro_group <- factor(top_features$macro_group, levels = grp_levels)

feat_order <- top_features %>%
  mutate(grp_idx = as.numeric(macro_group)) %>%
  group_by(feature) %>%
  summarise(mean_grp = mean(grp_idx), max_stat = max(stat)) %>%
  arrange(mean_grp, desc(max_stat)) %>%
  pull(feature)

top_features$feature <- factor(top_features$feature, levels = feat_order)

p_fingerprint <- ggplot(top_features, aes(x = macro_group, y = feature)) +
  geom_tile(fill = NA, color = "gray95", size = 0.2) +
  
  geom_point(aes(size = stat, fill = p.value), shape = 21, color = "gray30", stroke = 0.2) +
  
  scale_size_continuous(range = c(2, 5), name = "Indicator Value\n(Association Strength)") +
  scale_fill_distiller(palette = "YlGnBu", direction = 1, name = "Significance\n(p-value)", limits = c(0, 0.1)) +
  
  labs(
    title = "Chemical Class Fingerprint by Plant Lineage",
    subtitle = "Top indicator chemical classes (IndVal) showing phylogenetic specialization.",
    x = NULL,
    y = NULL
  ) +
  theme_nature() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major = element_line(linetype = "dotted"),
    legend.position = "right"
  )

ggsave(file.path(out_fig_dir, "Fig6_Indicator_Fingerprint.pdf"), p_fingerprint, width = 180, height = 160, units = "mm")


message("Generating Figure 4 (variance paradox)...")

perma_clean <- perma %>%
  filter(!is.na(R2)) %>%
  mutate(
    Type_Label = case_when(
      tag == "Scaffold_Jaccard" ~ "Structural\n(Scaffolds)",
      tag == "Class_Jaccard" ~ "Structural\n(Classes)",
      tag == "Descriptors_Euclidean" ~ "Physicochemical\n(Properties)",
      TRUE ~ tag
    ),
    Type_Order = case_when(
      tag == "Class_Jaccard" ~ 1,
      tag == "Scaffold_Jaccard" ~ 2,
      tag == "Descriptors_Euclidean" ~ 3
    ),
    Significance = ifelse(p.value < 0.05, "Significant (p < 0.05)", "Not Significant")
  ) %>%
  arrange(Type_Order)

p_variance <- ggplot(perma_clean, aes(x = reorder(Type_Label, Type_Order), y = R2)) +
  geom_col(aes(fill = Significance), width = 0.6, alpha = 0.9) +
  
  geom_text(aes(label = scales::percent(R2, accuracy = 0.1)), 
            vjust = -0.5, size = 3, fontface = "bold") +
  
  annotate("segment", x = 1, xend = 2, y = 0.16, yend = 0.16, color = "black") +
  annotate("segment", x = 1, xend = 1, y = 0.155, yend = 0.16, color = "black") +
  annotate("segment", x = 2, xend = 2, y = 0.155, yend = 0.16, color = "black") +
  annotate("text", x = 1.5, y = 0.17, label = "Divergent Evolution\n(Distinct Architectures)", 
           size = 2.8, fontface = "italic", color = "black") +
  
  annotate("text", x = 3, y = 0.17, label = "Convergent Evolution\n(Overlapping Properties)", 
           size = 2.8, fontface = "italic", color = "gray40") +
  
  scale_fill_manual(values = c("Significant (p < 0.05)" = "#377EB8", "Not Significant" = "gray70")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.20), expand = c(0,0)) +
  
  labs(
    title = "The Phylogenetic Signal Paradox",
    subtitle = "Variance explained by lineage (PERMANOVA RÂ²)",
    y = "Phylogenetic Signal Strength (RÂ²)",
    x = NULL
  ) +
  theme_nature() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 8)
  )

ggsave(file.path(out_fig_dir, "Fig4_Variance_Paradox.pdf"), p_variance, width = 90, height = 90, units = "mm")

message("Process completed. Figures saved to: ", out_fig_dir)


out_fig_dir <- fig_dir

message("Loading groups...")

groups_raw <- readxl::read_excel(path_groups, col_names = FALSE, .name_repair = "unique")

groups <- groups_raw |>
  dplyr::select(1:3) |>
  rlang::set_names(c("family", "color", "macro_group")) |>
  dplyr::mutate(
    family      = stringr::str_trim(as.character(family)),
    color       = stringr::str_trim(as.character(color)),
    macro_group = stringr::str_trim(as.character(macro_group)),
    color       = dplyr::if_else(stringr::str_detect(color, "^#"), color, paste0("#", color))
  ) |>
  dplyr::filter(!is.na(family), nzchar(family),
                !is.na(macro_group), nzchar(macro_group)) |>
  dplyr::distinct(family, .keep_all = TRUE)

fam_keep <- unique(groups$family)
macro_levels <- sort(unique(groups$macro_group))
pal_macro <- stats::setNames(groups$color, groups$macro_group)

message("Loading LOTUS lin_enriched (prefer in-memory)...")

if (!exists("lin_enriched", inherits = TRUE)) {
  stop(
    "lin_enriched not found in memory.\n",
    "Fix: run Part I (or the controller auto-load parquet) BEFORE Part III."
  )
}

lotus_desc <- get("lin_enriched", inherits = TRUE) |>
  janitor::clean_names() |>
  dplyr::mutate(
    inchikey = stringr::str_trim(as.character(inchikey)),
    family   = stringr::str_trim(as.character(family))
  ) |>
  dplyr::filter(!is.na(inchikey), nzchar(inchikey),
                !is.na(family), nzchar(family),
                family %in% fam_keep) |>
  dplyr::left_join(groups |> dplyr::select(family, macro_group), by = "family") |>
  dplyr::filter(!is.na(macro_group), nzchar(macro_group)) |>
  dplyr::mutate(macro_group = factor(macro_group, levels = macro_levels))

need_lotus <- c("family","macro_group","inchikey","molecular_weight","xlogp","topo_psa","heavy_atom_number")
missing_lotus <- setdiff(need_lotus, names(lotus_desc))
if (length(missing_lotus)) stop("lin_enriched missing LOTUS columns: ", paste(missing_lotus, collapse = ", "))

lotus_key <- lotus_desc |>
  dplyr::select(dplyr::all_of(need_lotus)) |>
  dplyr::distinct(family, inchikey, .keep_all = TRUE)  # avoids exploding refs overweight

message("Loading RDKit annotations (recalculated)...")

rdkit_raw <- readr::read_csv(path_rdkit, show_col_types = FALSE) |>
  janitor::clean_names() |>
  dplyr::mutate(inchikey = stringr::str_trim(as.character(inchikey)))

need_rdkit <- c("inchikey","num_hba","num_hbd","num_rings")
missing_rdkit <- setdiff(need_rdkit, names(rdkit_raw))
if (length(missing_rdkit)) stop("RDKit annotations missing columns: ", paste(missing_rdkit, collapse = ", "))

rdkit_key <- rdkit_raw |>
  dplyr::group_by(inchikey) |>
  dplyr::summarise(
    num_hba   = suppressWarnings(stats::median(as.numeric(num_hba), na.rm = TRUE)),
    num_hbd   = suppressWarnings(stats::median(as.numeric(num_hbd), na.rm = TRUE)),
    num_rings = suppressWarnings(stats::median(as.numeric(num_rings), na.rm = TRUE)),
    has_phenolic_oh     = dplyr::if_any(dplyr::matches("^has_phenolic"), ~ any(.x == 1, na.rm = TRUE)),
    has_methoxy_aryl    = dplyr::if_any(dplyr::matches("^has_methoxy"), ~ any(.x == 1, na.rm = TRUE)),
    has_prenyl_like     = dplyr::if_any(dplyr::matches("^has_prenyl"), ~ any(.x == 1, na.rm = TRUE)),
    has_conj_carbonyl   = dplyr::if_any(dplyr::matches("^has_conj"), ~ any(.x == 1, na.rm = TRUE)),
    has_probable_sugar  = dplyr::if_any(dplyr::matches("^has_probable_sugar$|^has_sugar"), ~ any(.x == 1, na.rm = TRUE)),
    .groups = "drop"
  )

message("Joining LOTUS (trusted) + RDKit (recalculated) and aggregating by family...")

chem_fam <- lotus_key |>
  dplyr::left_join(rdkit_key, by = "inchikey") |>
  dplyr::group_by(family, macro_group) |>
  dplyr::summarise(
    MW         = stats::median(molecular_weight, na.rm = TRUE),
    LogP       = stats::median(xlogp, na.rm = TRUE),
    TPSA       = stats::median(topo_psa, na.rm = TRUE),
    HeavyAtoms = stats::median(heavy_atom_number, na.rm = TRUE),
    
    HBA   = stats::median(num_hba, na.rm = TRUE),
    HBD   = stats::median(num_hbd, na.rm = TRUE),
    Rings = stats::median(num_rings, na.rm = TRUE),
    
    N_InChIKey = dplyr::n_distinct(inchikey),
    .groups = "drop"
  ) |>
  dplyr::left_join(groups, by = c("family","macro_group")) |>
  dplyr::filter(stats::complete.cases(MW, LogP, TPSA, HBA, HBD, Rings))

plot_data <- chem_fam



message("Generating PCoA Biplot (Structure Only)...")

class_counts <- read_csv(path_rdkit, show_col_types = FALSE) %>% clean_names() %>%
  filter(!is.na(class_np), family %in% plot_data$family) %>%
  count(family, class_np) %>%
  pivot_wider(names_from = class_np, values_from = n, values_fill = 0) %>%
  column_to_rownames("family")

common_fams <- intersect(rownames(class_counts), plot_data$family)
class_mat <- class_counts[common_fams, ]
meta_pcoa <- plot_data %>% filter(family %in% common_fams) %>% arrange(match(family, rownames(class_mat)))

dist_struct <- vegdist(class_mat, method = "bray")
pcoa <- cmdscale(dist_struct, k = 2, eig = TRUE)
scores_pcoa <- as.data.frame(pcoa$points)
colnames(scores_pcoa) <- c("PCoA1", "PCoA2")
scores_pcoa <- bind_cols(meta_pcoa, scores_pcoa)

centroids_pcoa <- scores_pcoa %>%
  group_by(macro_group) %>%
  summarise(
    PCoA1 = mean(PCoA1, na.rm = TRUE),
    PCoA2 = mean(PCoA2, na.rm = TRUE),
    n_families = n()
  ) %>%
  ungroup()

fit_classes <- envfit(pcoa, class_mat, permutations = 999)
vec_classes <- as.data.frame(fit_classes$vectors$arrows * sqrt(fit_classes$vectors$r)) %>%
  rownames_to_column("Label") %>% 
  filter(fit_classes$vectors$pvals < 0.05) %>%
  mutate(Type = "Class")

if(nrow(vec_classes) > 0) colnames(vec_classes)[2:3] <- c("Dim1", "Dim2")

vectors_pcoa <- vec_classes

mult_pcoa <- 0.8 

p_pcoa <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray80") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray80") +
  
  geom_segment(data = vectors_pcoa, aes(x = 0, y = 0, xend = Dim1 * mult_pcoa, yend = Dim2 * mult_pcoa),
               arrow = arrow(length = unit(0.25, "cm")), color = "black", size = 0.7) +
  
  geom_text_repel(data = vectors_pcoa, aes(x = Dim1 * mult_pcoa, y = Dim2 * mult_pcoa, label = Label),
                  fontface = "italic", size = 3.5, color = "black", max.overlaps = 30, bg.color="white", bg.r=0.15) +
  
  geom_point(data = centroids_pcoa, aes(x = PCoA1, y = PCoA2, color = macro_group), 
             size = 6, alpha = 0.9, stroke = 1, shape = 16) +
  
  geom_text_repel(data = centroids_pcoa, aes(x = PCoA1, y = PCoA2, label = macro_group),
                  fontface = "bold", size = 4, color = "gray40", box.padding = 0.6) +
  
  scale_color_manual(values = pal_macro) +
  
  labs(
    title = "PCoA Biplot: Structural Divergence (Centroids)",
    subtitle = "Chemical class vectors projecting onto lineage space.",
    x = "PCoA 1", y = "PCoA 2"
  ) +
  theme_classic() +
  theme(legend.position = "none")

ggsave(file.path(out_fig_dir, "PCoA_Biplot_Structure_Chem_Only.pdf"), p_pcoa, width = 12, height = 9)

message("Chemical centroid biplots saved to: ", out_fig_dir)



out_fig_dir <- fig_dir  

groups <- read_excel(path_groups, col_names = FALSE) %>%
  transmute(
    family      = str_trim(as.character(...1)),
    color       = str_trim(as.character(...2)),
    macro_group = str_trim(as.character(...3))
  ) %>%
  filter(!is.na(family), !is.na(macro_group)) %>%
  distinct(family, .keep_all = TRUE) %>%
  mutate(color = ifelse(str_detect(color, "^#"), color, paste0("#", color)))

message("Loading LOTUS chemical descriptors (prefer in-memory lin_enriched)...")

lotus_desc <- get_lotus_lin("lin_enriched") %>%
  dplyr::mutate(
    inchikey = stringr::str_trim(as.character(.data$inchikey))
  ) %>%
  dplyr::rename(
    topo_psa = topoPSA
  ) %>%
  dplyr::select(inchikey, molecular_weight, xlogp, topo_psa, heavy_atom_number) %>%
  dplyr::mutate(
    dplyr::across(c(molecular_weight, xlogp, topo_psa, heavy_atom_number), suppressWarnings(as.numeric))
  ) %>%
  dplyr::filter(!is.na(inchikey), nzchar(inchikey)) %>%
  dplyr::group_by(inchikey) %>%
  dplyr::summarise(
    dplyr::across(dplyr::everything(), ~ stats::median(.x, na.rm = TRUE)),
    .groups = "drop"
  )

message("Loading RDKit annotations...")
rdkit <- read_csv(path_rdkit, show_col_types = FALSE) %>% 
  clean_names() %>%
  select(
    inchikey, family, num_hba, num_hbd, num_rings, class_np,
    has_prenyl_like, has_probable_sugar, has_methoxy_aryl, 
    has_phenolic_oh, has_conj_carbonyl
  )


message("Computing physicochemical property profiles...")
chem_data <- rdkit %>%
  left_join(lotus_desc, by = "inchikey") %>%
  filter(family %in% groups$family) %>%
  left_join(groups, by = "family") %>%
  filter(!is.na(macro_group)) %>%
  group_by(macro_group) %>%
  summarise(
    MW = median(molecular_weight, na.rm=TRUE),
    LogP = median(xlogp, na.rm=TRUE),
    TPSA = median(topo_psa, na.rm=TRUE),
    HBA = median(num_hba, na.rm=TRUE),
    HBD = median(num_hbd, na.rm=TRUE),
    Rings = median(num_rings, na.rm=TRUE),
    n_families = n_distinct(family), 
    .groups = "drop"
  ) %>%
  filter(complete.cases(.)) %>%
  column_to_rownames("macro_group")

mat_props <- as.matrix(chem_data %>% select(MW, LogP, TPSA, HBA, HBD, Rings))
mat_props_scaled <- scale(mat_props)

message("Computing chemical class profiles...")
class_profile <- rdkit %>%
  filter(!is.na(class_np)) %>%
  left_join(groups, by = "family") %>%
  filter(!is.na(macro_group)) %>%
  count(macro_group, class_np) %>%
  group_by(macro_group) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  select(macro_group, class_np, prop) %>%
  pivot_wider(names_from = class_np, values_from = prop, values_fill = 0) %>%
  column_to_rownames("macro_group")

common_groups <- intersect(rownames(mat_props_scaled), rownames(class_profile))
mat_classes <- as.matrix(class_profile[common_groups, ])
mat_classes_scaled <- scale(mat_classes)
mat_classes_scaled[is.na(mat_classes_scaled)] <- 0

message("Computing structural feature profiles...")
struct_profile <- rdkit %>%
  left_join(groups, by = "family") %>%
  filter(!is.na(macro_group)) %>%
  group_by(macro_group) %>%
  summarise(
    Prenyl = mean(has_prenyl_like, na.rm = TRUE),
    Sugar = mean(has_probable_sugar, na.rm = TRUE),
    Methoxy = mean(has_methoxy_aryl, na.rm = TRUE),
    Phenol_OH = mean(has_phenolic_oh, na.rm = TRUE),
    Conj_CO = mean(has_conj_carbonyl, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  column_to_rownames("macro_group")

mat_struct <- as.matrix(struct_profile[common_groups, ])
mat_struct_scaled <- scale(mat_struct)
mat_struct_scaled[is.na(mat_struct_scaled)] <- 0

mat_props_scaled <- mat_props_scaled[common_groups, ]
fam_counts <- chem_data[common_groups, "n_families"]

group_colors_vec <- groups %>% select(macro_group, color) %>% distinct() %>% deframe_df()
use_colors <- group_colors_vec[common_groups]
use_colors[is.na(use_colors)] <- "grey80"

col_props <- colorRamp2(c(-2, 0, 2), c("#2166AC", "#F7F7F7", "#B2182B")) 
col_struct <- colorRamp2(c(-1.5, 0, 2.5), c("white", "#E5F5F9", "#2CA25F"))
col_classes <- colorRamp2(c(-1, 0, 3), c("white", "#E0ECF4", "#8856A7"))

row_anot <- HeatmapAnnotation(
  "Families (n)" = anno_barplot(
    fam_counts, 
    width = unit(2.5, "cm"), 
    gp = gpar(fill = "gray40", col = NA),
    axis_param = list(side = "top", gp = gpar(fontsize = 8))
  ),
  "Lineage" = common_groups,
  col = list("Lineage" = use_colors),
  show_legend = TRUE,
  annotation_name_side = "top",
  which = "row"
)

message("Generating integrated heatmap...")

ht_props <- Heatmap(
  mat_props_scaled,
  name = "Props\n(Z-score)",
  col = col_props,
  
  cluster_rows = TRUE, 
  cluster_columns = TRUE,
  
  row_names_side = "left", 
  row_dend_side = "left",
  row_dend_width = unit(1.5, "cm"), 
  
  rect_gp = gpar(col = "white", lwd = 0.5),
  border = TRUE,
  row_names_gp = gpar(fontsize = 11, fontface = "bold"),
  column_names_gp = gpar(fontsize = 9),
  column_names_rot = 45,
  
  column_title = "Physicochemical\nProperties",
  column_title_gp = gpar(fontsize = 10, fontface = "bold")
)

ht_struct <- Heatmap(
  mat_struct_scaled,
  name = "Features\n(Z-score)",
  col = col_struct,
  
  cluster_rows = FALSE, 
  cluster_columns = TRUE,
  
  rect_gp = gpar(col = "white", lwd = 0.5),
  border = TRUE,
  show_row_names = FALSE, 
  column_names_gp = gpar(fontsize = 9),
  column_names_rot = 45,
  
  column_title = "Structural\nFeatures",
  column_title_gp = gpar(fontsize = 10, fontface = "bold")
)

ht_classes <- Heatmap(
  mat_classes_scaled,
  name = "Classes\n(Z-score)",
  col = col_classes,
  
  cluster_rows = FALSE, 
  cluster_columns = TRUE,
  
  rect_gp = gpar(col = "white", lwd = 0.5),
  border = TRUE,
  show_row_names = FALSE,
  column_names_gp = gpar(fontsize = 8, fontface = "italic"),
  column_names_rot = 45,
  
  right_annotation = row_anot,
  
  column_title = "Chemical\nClasses",
  column_title_gp = gpar(fontsize = 10, fontface = "bold")
)

ht_list <- ht_props + ht_struct + ht_classes


pdf_file <- file.path(out_fig_dir, "Heatmap_Combined_Analysis.pdf")

pdf(pdf_file, width = 18, height = 9)
ComplexHeatmap::draw(
  ht_list,
  column_title = "Integrated Chemical Analysis: Properties, Features & Classes",
  column_title_gp = grid::gpar(fontsize = 16, fontface = "bold"),
  gap = grid::unit(3, "mm"),
  padding = grid::unit(c(10, 10, 10, 25), "mm"),
  merge_legend = TRUE,
  heatmap_legend_side = "right"
)
dev.off()

message("Integrated heatmap saved to: ", pdf_file)


`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

stopifnot(exists("uni_enriched"))

tax_col <- if (exists("tax_col") && is.character(tax_col) && length(tax_col) == 1 && nzchar(tax_col)) {
  tolower(trimws(tax_col))
} else {
  "family"
}

if (!exists("map_tax_inchi")) {
  if (exists("lin_enriched")) {
    if (!all(c("inchikey", tax_col) %in% names(lin_enriched))) {
      stop("map_tax_inchi not found and lin_enriched lacks required columns: inchikey + ", tax_col)
    }
    map_tax_inchi <- lin_enriched %>%
      dplyr::select(inchikey, !!rlang::sym(tax_col)) %>%
      dplyr::distinct()
  } else {
    stop("map_tax_inchi not found. Provide it OR ensure lin_enriched is loaded to build it.")
  }
}

needed_uni <- c("inchikey", "murko_framework")
missing_uni <- setdiff(needed_uni, names(uni_enriched))
if (length(missing_uni)) {
  stop("uni_enriched is missing required columns: ", paste(missing_uni, collapse = ", "))
}

keep_cols <- intersect(
  names(uni_enriched),
  c("inchikey","murko_framework","hybrid_class","class_np","mw","tpsa","xlogp","hba","hbd","rings","fsp3")
)
if (!("inchikey" %in% keep_cols) || !("murko_framework" %in% keep_cols)) {
  stop("uni_enriched is missing inchikey/murko_framework after keep_cols. Available: ",
       paste(names(uni_enriched), collapse = ", "))
}

df_props <- uni_enriched %>%
  dplyr::select(dplyr::all_of(keep_cols)) %>%
  dplyr::distinct() %>%
  dplyr::inner_join(
    map_tax_inchi %>% dplyr::select(inchikey, !!rlang::sym(tax_col)),
    by = "inchikey"
  ) %>%
  dplyr::filter(!is.na(.data[[tax_col]]), nzchar(as.character(.data[[tax_col]])))

message("df_props columns: ", paste(names(df_props), collapse = ", "))

if (nrow(df_props) == 0) {
  message(">>> PART V SKIPPED: df_props is empty after join/filter.")
} else {
  
  message("\n>>> STARTING PART V: SCAFFOLD INNOVATION ANALYSIS...")
  
  if (!exists("taxa_keep")) {
    taxa_keep <- df_props %>%
      dplyr::distinct(.data[[tax_col]]) %>%
      dplyr::pull(1) %>%
      as.character()
  } else {
    taxa_keep <- as.character(taxa_keep)
    taxa_keep <- taxa_keep[!is.na(taxa_keep) & nzchar(taxa_keep)]
  }
  
  df_scaff <- df_props %>%
    dplyr::filter(!is.na(murko_framework), nzchar(murko_framework)) %>%
    dplyr::select(inchikey, murko_framework, Taxon = !!rlang::sym(tax_col))
  
  if (length(taxa_keep) > 0L) {
    df_scaff <- df_scaff %>% dplyr::filter(.data$Taxon %in% taxa_keep)
  }
  
  if (nrow(df_scaff) == 0) {
    message("      [Skipped] No Murcko frameworks found after filters.")
  } else {
    
    scaffold_stats <- df_scaff %>%
      dplyr::group_by(Taxon) %>%
      dplyr::summarise(
        N_Compounds = dplyr::n_distinct(inchikey),
        N_Scaffolds = dplyr::n_distinct(murko_framework),
        Innovation_Ratio = round(N_Scaffolds / pmax(N_Compounds, 1), 3),
        .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(Innovation_Ratio))
    
    scaff_distribution <- df_scaff %>%
      dplyr::distinct(murko_framework, Taxon) %>%
      dplyr::count(murko_framework, name = "Taxon_Count")
    
    unique_scaffolds <- scaff_distribution %>%
      dplyr::filter(Taxon_Count == 1) %>%
      dplyr::pull(murko_framework)
    
    exclusivity_stats <- df_scaff %>%
      dplyr::filter(murko_framework %in% unique_scaffolds) %>%
      dplyr::group_by(Taxon) %>%
      dplyr::summarise(
        N_Exclusive_Scaffolds = dplyr::n_distinct(murko_framework),
        Example_Exclusive_Scaffold = dplyr::first(murko_framework),
        .groups = "drop"
      )
    
    final_stats <- scaffold_stats %>%
      dplyr::left_join(exclusivity_stats, by = "Taxon") %>%
      dplyr::mutate(
        N_Exclusive_Scaffolds = tidyr::replace_na(N_Exclusive_Scaffolds, 0L),
        Pct_Exclusive = dplyr::if_else(
          N_Scaffolds > 0,
          round((N_Exclusive_Scaffolds / N_Scaffolds) * 100, 1),
          0
        )
      ) %>%
      dplyr::arrange(dplyr::desc(N_Exclusive_Scaffolds))
    
    tag <- runtime$base_tag %||% "lotus_run"
    out_scaff <- file.path(out_fig_dir, paste0(tag, "_STATS_Scaffold_Innovation.xlsx"))
    
    openxlsx::write.xlsx(
      x = list(Scaffold_Metrics = final_stats),
      file = out_scaff,
      overwrite = TRUE
    )
    
    message("      [Saved] Scaffold Innovation Stats: ", basename(out_scaff))
    message("      -> Top Exclusive Taxon: ", final_stats$Taxon[1], " (",
            final_stats$N_Exclusive_Scaffolds[1], " unique scaffolds)")
  }
  
  message(">>> PART V COMPLETED.")
}

if (!exists("final_stats") || is.null(final_stats) || nrow(final_stats) == 0) {
  message(">>> Innovation Panel SKIPPED: final_stats not available (or empty).")
} else {
  
  data_clean <- final_stats %>%
    dplyr::mutate(
      N_Compounds = as.numeric(N_Compounds),
      Pct_Exclusive = as.numeric(Pct_Exclusive),
      Taxon = stringr::str_trim(as.character(Taxon))
    )
  
  groups_raw <- readxl::read_excel(path_groups, col_names = FALSE, .name_repair = "minimal")
  
  groups_df <- empty_df(
    Taxon      = stringr::str_trim(as.character(groups_raw[[1]])),
    color_hex  = stringr::str_trim(as.character(groups_raw[[2]])),
    macro_group= stringr::str_trim(as.character(groups_raw[[3]]))
  ) %>%
    dplyr::filter(!is.na(Taxon), nzchar(Taxon), !is.na(macro_group), nzchar(macro_group)) %>%
    dplyr::distinct(Taxon, .keep_all = TRUE) %>%
    dplyr::mutate(color_hex = dplyr::if_else(stringr::str_detect(color_hex, "^#"),
                                             color_hex, paste0("#", color_hex)))
  
  official_palette <- groups_df %>%
    dplyr::distinct(macro_group, color_hex) %>%
    dplyr::filter(!is.na(macro_group), nzchar(macro_group), !is.na(color_hex), nzchar(color_hex)) %>%
    { setNames(.$color_hex, .$macro_group) }
  
  message("Palette loaded from groups.xlsx:")
  print(official_palette)
  
  data_plot <- data_clean %>%
    dplyr::left_join(groups_df, by = "Taxon") %>%
    dplyr::mutate(macro_group = dplyr::if_else(is.na(macro_group), "Other", macro_group)) %>%
    dplyr::filter(macro_group != "Other")
  
  if (nrow(data_plot) == 0) {
    message(">>> Innovation Panel SKIPPED: no taxa matched groups.xlsx mapping.")
  } else {
    
    label_subset <- data_plot %>%
      dplyr::filter(Pct_Exclusive > 40 & N_Compounds > 5)
    
    plot_b <- ggplot2::ggplot(data_plot, ggplot2::aes(x = N_Compounds, y = Pct_Exclusive)) +
      ggplot2::geom_point(ggplot2::aes(color = macro_group), size = 3.5, alpha = 0.8) +
      ggplot2::scale_color_manual(values = official_palette, name = "Lineage") +
      ggplot2::scale_x_log10(labels = scales::comma) +
      ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 1)) +
      ggrepel::geom_text_repel(
        data = label_subset,
        ggplot2::aes(label = Taxon, color = macro_group),
        size = 3, fontface = "bold",
        box.padding = 0.4, point.padding = 0.3,
        max.overlaps = 30, show.legend = FALSE
      ) +
      ggplot2::labs(
        title = "Exclusivity vs. Sampling Effort",
        subtitle = "Scaffold novelty across lineages (colors from groups.xlsx)",
        x = "Reported Compounds (Log Scale)",
        y = "Exclusive Scaffolds (%)"
      ) +
      ggplot2::theme_classic(base_size = 14) +
      ggplot2::theme(
        legend.position = "right",
        legend.title = ggplot2::element_text(face = "bold"),
        legend.text  = ggplot2::element_text(size = 9)
      )
    
    macro_stats <- data_plot %>%
      dplyr::group_by(macro_group) %>%
      dplyr::summarise(
        Mean_Exclusivity = mean(Pct_Exclusive, na.rm = TRUE),
        n_fams = dplyr::n(),
        .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(Mean_Exclusivity))
    
    plot_c <- ggplot2::ggplot(macro_stats,
                              ggplot2::aes(x = reorder(macro_group, Mean_Exclusivity),
                                           y = Mean_Exclusivity)) +
      ggplot2::geom_col(ggplot2::aes(fill = macro_group),
                        width = 0.7, color = "black", linewidth = 0.2) +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f%%", Mean_Exclusivity)),
                         hjust = -0.2, size = 3.5, fontface = "bold", color = "black") +
      ggplot2::scale_fill_manual(values = official_palette) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.4))) +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = "Lineage Innovation Ranking",
        subtitle = "Mean scaffold exclusivity per macro-group",
        x = NULL, y = "Mean Exclusive Scaffolds (%)"
      ) +
      ggplot2::theme_classic() +
      ggplot2::theme(legend.position = "none")
    
    ggplot2::ggsave(file.path(out_fig_dir, "Panel_B_Dynamic_Palette.pdf"),
                    plot = plot_b, width = 12, height = 8)
    ggplot2::ggsave(file.path(out_fig_dir, "Panel_C_Dynamic_Palette.pdf"),
                    plot = plot_c, width = 8, height = 7)
    
    message("Panels generated with the palette defined in groups.xlsx: ", out_fig_dir)
  }
}


library(readxl)
library(ggplot2)
library(dplyr)
library(ggrepel)

file_path <- out_xlsx

if(!file.exists(file_path)) stop("Required workbook not found: ", file_path)

df_landscape <- read_excel(file_path, sheet = "landscape_family")

df_calc <- df_landscape %>%
  filter(n_compounds >= 3) %>%
  mutate(
    n_exclusive_est = round(novelty_ratio * n_scaffolds),
    n_exclusive_est = if_else(n_exclusive_est > n_scaffolds, n_scaffolds, n_exclusive_est)
  )

total_scaffolds <- sum(df_calc$n_scaffolds, na.rm = TRUE)
total_exclusive <- sum(df_calc$n_exclusive_est, na.rm = TRUE)
global_rate <- total_exclusive / total_scaffolds

calculate_binom_p <- function(k, n, p) {
  if(is.na(k) || is.na(n) || n == 0) return(1.0)
  binom.test(x = k, n = n, p = p, alternative = "greater")$p.value
}

top_giants <- df_calc %>%
  arrange(desc(n_compounds)) %>%
  slice_head(n = 10) %>%
  pull(family)

plot_data <- df_calc %>%
  rowwise() %>%
  mutate(
    p_val = calculate_binom_p(n_exclusive_est, n_scaffolds, global_rate),
    log_p = -log10(p_val),
    
    tested_safe = ifelse(is.na(n_compounds_tested), 0, n_compounds_tested),
    gap_rate = 1 - (tested_safe / n_compounds),
    gap_rate = ifelse(gap_rate < 0, 0, gap_rate),
    
    is_significant = p_val < 0.05,
    is_giant = family %in% top_giants,
    
    is_highlight = is_significant | is_giant,
    
    full_label = paste0(family, "\n(", macro_group, ")"),
    
    label_text = if_else(is_highlight, full_label, NA_character_),
    
    font_face = dplyr::case_when(
      is_significant ~ "bold",
      is_giant ~ "italic",
      TRUE ~ "plain"
    )
  ) %>%
  ungroup() %>%
  arrange(is_highlight)

p_volcano <- ggplot(plot_data, aes(x = novelty_ratio, y = log_p)) +
  
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", size = 0.3) +
  geom_vline(xintercept = global_rate, linetype = "dotted", color = "gray50", size = 0.5) +
  
  annotate("text", x = 0.95, y = -log10(0.05) + 0.1, label = "Significance (p < 0.05)", 
           color = "red", size = 3, hjust = 1) +
  annotate("text", x = global_rate + 0.01, y = 0, label = "Global Avg", 
           color = "gray50", size = 3, angle = 90, hjust = 0, vjust = 0) +
  
  geom_point(data = subset(plot_data, !is_highlight),
             aes(size = n_compounds), 
             color = "gray90", fill = "gray95", shape = 21) +
  
  geom_point(data = subset(plot_data, is_highlight),
             aes(size = n_compounds, fill = gap_rate), 
             shape = 21, color = "black", stroke = 0.4) +
  
  geom_text_repel(aes(label = label_text, fontface = font_face), 
                  size = 2.5,
                  max.overlaps = 40,
                  box.padding = 0.6,
                  point.padding = 0.3,
                  force = 4,
                  min.segment.length = 0,
                  segment.color = "gray60",
                  segment.size = 0.3,
                  lineheight = 0.8) +
  
  scale_fill_gradientn(colors = c("#4575b4", "#ffffbf", "#d73027"),
                       name = "Untested Rate\n(Opportunity Gap)",
                       labels = scales::percent) +
  
  scale_size_continuous(trans = "log10", range = c(2, 9), 
                        name = "Compounds (n)",
                        breaks = c(10, 100, 1000, 5000)) +
  
  scale_x_continuous(labels = scales::percent, name = "Structural Innovation Rate") +
  scale_y_continuous(name = "Significance (-log10 p-value)") +
  
  theme_classic() +
  labs(
    title = "Discovery Volcano: Innovation vs. Opportunity",
    subtitle = "Bold = Statistically Novel (p < 0.05). Italic = Top 10 Largest Families. Color = % Untested."
  ) +
  theme(
    axis.title = element_text(face = "bold", size = 10),
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(face = "bold")
  )

print(p_volcano)
out_dir <- dirname(file_path)
ggsave(file.path(out_dir, "Fig3A_Volcano_Final_Clades.pdf"), 
       plot = p_volcano, width = 11, height = 9, device = cairo_pdf)

message("Figure saved to: ", file.path(out_dir, "Fig3A_Volcano_Final_Clades.pdf"))


library(readxl)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(ggpubr)

file_path <- out_xlsx

if(!file.exists(file_path)) {
  stop("Required workbook not found: ", file_path)
}

df_landscape <- read_excel(file_path, sheet = "landscape_family")

plot_data <- df_landscape %>%
  filter(n_compounds >= 3) %>%
  mutate(
    n_exclusive_est = round(novelty_ratio * n_scaffolds),
    p_val = purrr::map2_dbl(n_exclusive_est, n_scaffolds, 
                            ~binom.test(.x, .y, p=0.53, alternative="greater")$p.value),
    is_volcano_sig = p_val < 0.05,
    
    is_top10 = rank(-n_compounds) <= 10,
    label_fam = if_else(is_volcano_sig | is_top10, family, NA_character_)
  )

p_gradient <- ggplot(plot_data, aes(x = n_compounds, y = novelty_ratio)) +
  
  geom_smooth(method = "lm", formula = y ~ log10(x), 
              color = "gray30", fill = "gray85", alpha = 0.3, size = 0.8) +
  
  geom_point(aes(fill = novelty_ratio, size = n_compounds, color = is_volcano_sig), 
             shape = 21, stroke = 1.2) +
  
  geom_text_repel(aes(label = label_fam), 
                  size = 3.5,
                  fontface = ifelse(plot_data$is_volcano_sig, "bold", "italic"),
                  max.overlaps = 50,
                  box.padding = 0.6,
                  min.segment.length = 0,
                  segment.color = "gray50") +
  
  stat_cor(method = "spearman", label.x.npc = "center", label.y = 1.05, 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  
  scale_fill_gradientn(colors = c("#4575b4", "#ffffbf", "#d73027"), 
                       name = "Novelty Ratio", labels = scales::percent) +
  
  scale_color_manual(values = c("FALSE" = "gray40", "TRUE" = "black"), guide = "none") +
  
  scale_x_log10(name = "Sampling Effort (Number of Compounds, log scale)", 
                breaks = c(5, 10, 50, 100, 500, 1000, 5000)) +
  scale_y_continuous(name = "Structural Novelty Rate", labels = scales::percent, limits = c(0, 1.1)) +
  scale_size_continuous(range = c(3, 10), name = "Compounds (n)") +
  
  theme_bw() +
  labs(
    title = "Discovery Potential: Effort vs. Novelty Gradient",
    subtitle = "Color indicates novelty rate. Bold labels/Black borders indicate statistically significant innovators (Volcano p<0.05).",
    caption = "Top 10 largest families labeled in italics. Trend line represents expected diminishing returns."
  ) +
  theme(
    legend.position = "right",
    axis.title = element_text(face = "bold")
  )

print(p_gradient)
out_dir <- dirname(file_path)
ggsave(file.path(out_dir, "Fig3C_Bias_Gradient.pdf"), 
       plot = p_gradient, width = 10, height = 8, device = cairo_pdf)

message("Figure saved to: ", file.path(out_dir, "Fig3C_Bias_Gradient.pdf"))

library(readxl)
library(dplyr)
library(ggplot2)
library(ggrepel)

path_stats <- out_xlsx
path_raw <- path_chembl

message("Loading data...")
df_colors <- read_excel(path_stats, sheet = "macro_groups_mapping")
df_master <- read_excel(path_raw, sheet = "MASTER_DATA")

df_chem <- df_master %>%
  select(family, chemicalTaxonomyNPclassifierClass, fsp3, alogp, pActivity) %>%
  rename(chem_class = chemicalTaxonomyNPclassifierClass) %>%
  filter(!is.na(chem_class), !is.na(fsp3), !is.na(alogp))

class_stats <- df_chem %>%
  group_by(chem_class) %>%
  summarise(
    n_compounds = n(),
    med_fsp3 = median(fsp3, na.rm = TRUE),
    med_logp = median(alogp, na.rm = TRUE),
    med_potency = ifelse(all(is.na(pActivity)), NA, median(pActivity, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  filter(n_compounds >= 5)

cor_test <- cor.test(class_stats$med_fsp3, class_stats$med_potency, method = "spearman", use = "complete.obs")
cor_label <- paste0("Correlation (Fsp3 vs Potency):\nRho = ", round(cor_test$estimate, 2), " (p = ", round(cor_test$p.value, 2), ")\nNo significant link")

center_fsp3 <- median(class_stats$med_fsp3, na.rm = TRUE)
center_logp <- median(class_stats$med_logp, na.rm = TRUE)
potency_breaks <- quantile(class_stats$med_potency, probs = c(0.33, 0.66), na.rm = TRUE)

class_stats <- class_stats %>%
  mutate(
    potency_tier = case_when(
      is.na(med_potency) ~ "No Data",
      med_potency > potency_breaks[2] ~ "Elite Potency",
      med_potency > potency_breaks[1] ~ "Moderate",
      TRUE ~ "Baseline"
    ),
    label_txt = chem_class, 
    font_face = ifelse(potency_tier == "Elite Potency", "bold", "plain")
  )


p_rel <- ggplot(class_stats, aes(x = med_fsp3, y = med_logp)) +
  
  geom_vline(xintercept = center_fsp3, linetype = "dashed", color = "gray80") +
  geom_hline(yintercept = center_logp, linetype = "dashed", color = "gray80") +
  
  annotate("label", x = 0, y = max(class_stats$med_logp), label = cor_label, 
           hjust = 0, vjust = 1, size = 3, fill = "white", alpha = 0.8, color = "black") +
  
  geom_point(aes(fill = potency_tier, size = n_compounds), 
             shape = 21, color = "black", stroke = 0.5, alpha = 0.9) +
  
  geom_text_repel(aes(label = label_txt, fontface = font_face), 
                  size = 3.5, 
                  box.padding = 0.5, 
                  force = 30,
                  max.overlaps = Inf,
                  min.segment.length = 0) +
  
  scale_fill_manual(values = c(
    "Elite Potency" = "#D7191C",   
    "Moderate" = "#FDAE61",        
    "Baseline" = "#ABD9E9",        
    "No Data" = "gray95"
  ), name = "Relative Potency") +
  
  scale_size_continuous(trans = "log10", range = c(4, 12), name = "Compounds (n)") +
  scale_x_continuous(name = paste0("Complexity (Median Fsp3) - Center: ", round(center_fsp3, 2))) +
  scale_y_continuous(name = paste0("Lipophilicity (Median LogP) - Center: ", round(center_logp, 2))) +
  
  theme_bw() +
  labs(
    title = "Chemotaxonomic Constellation",
    subtitle = "Mapping all flavonoid classes by complexity and lipophilicity.",
    caption = "All classes with n >= 5 compounds shown."
  ) +
  theme(axis.title = element_text(face = "bold"), legend.position = "right")

out_dir <- dirname(path_stats)
ggsave(file.path(out_dir, "Fig3C_Constellation_AllLabels.pdf"), 
       plot = p_rel, width = 12, height = 10, device = cairo_pdf)

message("Figure saved with all labels: ", file.path(out_dir, "Fig3C_Constellation_AllLabels.pdf"))


try(detach("package:ComplexHeatmap", unload=TRUE), silent=TRUE)

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(viridis)
library(circlize)
library(ComplexHeatmap)
library(grid)

path_stats <- out_xlsx
path_raw <- path_chembl

message("Loading data...")
df_colors <- read_excel(path_stats, sheet = "macro_groups_mapping")
macro_palette <- setNames(df_colors$color, df_colors$macro_group)
fam_map <- df_colors %>% select(family, macro_group) %>% distinct()
df_master <- read_excel(path_raw, sheet = "MASTER_DATA")

df_bio <- df_master %>%
  filter(!is.na(pActivity), !is.na(Target_Category_Macro)) %>%
  left_join(fam_map, by = "family") %>%
  filter(!is.na(macro_group)) %>%
  filter(!Target_Category_Macro %in% c("Unclassified", "Miscellaneous", "ADMET"))

message("Preparing Panel A...")

heatmap_data <- df_bio %>%
  group_by(macro_group, Target_Category_Macro) %>%
  summarise(med_potency = median(pActivity), .groups = "drop") %>%
  pivot_wider(names_from = Target_Category_Macro, values_from = med_potency) %>%
  as.data.frame()

rownames(heatmap_data) <- heatmap_data$macro_group
mat <- as.matrix(heatmap_data[, -1]) 

clade_counts <- df_bio %>%
  group_by(macro_group) %>%
  summarise(n = n(), .groups="drop")
row_counts <- clade_counts$n[match(rownames(mat), clade_counts$macro_group)]

row_ha <- HeatmapAnnotation(
  "Hits (n)" = anno_barplot(row_counts, 
                            gp = gpar(fill = "gray40", col = NA), 
                            width = unit(1.5, "cm"),
                            border = FALSE),
  which = "row",
  annotation_name_side = "top",
  annotation_name_gp = gpar(fontsize = 8)
)

mat_for_clustering <- mat
mat_for_clustering[is.na(mat_for_clustering)] <- min(mat, na.rm = TRUE)
dend_row <- hclust(dist(mat_for_clustering, method = "euclidean"), method = "ward.D2")
dend_col <- hclust(dist(t(mat_for_clustering), method = "euclidean"), method = "ward.D2")

col_fun <- colorRamp2(c(5, 6, 7), c("black", "#B63679", "#FCFDBF")) 

ht <- Heatmap(mat, 
              name = "pActivity", 
              col = col_fun,
              na_col = "gray95",
              cluster_rows = dend_row,
              cluster_columns = dend_col,
              rect_gp = gpar(col = "white", lwd = 0.5),
              row_names_gp = gpar(fontsize = 9, fontface = "bold"),
              column_names_gp = gpar(fontsize = 9, fontface = "bold"),
              column_names_rot = 45,
              right_annotation = row_ha,
              column_title = "A. Phylogenetic Bioactivity Atlas",
              column_title_gp = gpar(fontsize = 12, fontface = "bold"),
              heatmap_legend_param = list(direction = "horizontal", title_position = "topcenter")
)

message("Preparing Panel B...")

class_owner <- df_bio %>%
  filter(!is.na(chemicalTaxonomyNPclassifierClass)) %>%
  group_by(chemicalTaxonomyNPclassifierClass, macro_group) %>%
  summarise(n = n(), .groups="drop_last") %>%
  top_n(1, n) %>%
  rename(dominant_clade = macro_group) %>%
  select(chemicalTaxonomyNPclassifierClass, dominant_clade)

bubble_data <- df_bio %>%
  filter(!is.na(chemicalTaxonomyNPclassifierClass)) %>%
  inner_join(class_owner, by = "chemicalTaxonomyNPclassifierClass") %>%
  group_by(chemicalTaxonomyNPclassifierClass, dominant_clade, Target_Category_Macro) %>%
  summarise(med_potency = median(pActivity), count = n(), .groups = "drop") %>%
  filter(count >= 3) %>%
  arrange(desc(med_potency)) %>%
  group_by(Target_Category_Macro) %>%
  slice_head(n = 6)

p_bubble <- ggplot(bubble_data, aes(x = Target_Category_Macro, y = chemicalTaxonomyNPclassifierClass)) +
  geom_point(aes(size = count, fill = dominant_clade), 
             shape = 21, color = "black", stroke = 0.3, alpha = 0.8) +
  scale_fill_manual(values = macro_palette, name = "Lineage") +
  scale_size_continuous(range = c(2, 6), name = "Hits") +
  theme_bw() +
  labs(title = "B. Chemical Drivers of Potency", x = "", y = "") +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 12)
  )

message("Preparing Panel C...")

ridge_data <- df_bio %>%
  group_by(macro_group) %>%
  filter(n() >= 10) %>%
  ungroup() %>%
  mutate(macro_group = reorder(macro_group, pActivity, median))

p_ridge <- ggplot(ridge_data, aes(x = pActivity, y = macro_group, fill = macro_group)) +
  geom_density_ridges(scale = 2, rel_min_height = 0.01, alpha = 0.8, color = "white", size=0.2) +
  geom_vline(xintercept = 6, linetype = "dashed", alpha = 0.5) +
  scale_fill_manual(values = macro_palette) +
  theme_ridges() +
  labs(title = "C. Evolutionary Potency Distribution", x = "Bioactivity (pActivity)", y = "") +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 12)
  )

message("Drawing final A4 Composite...")

out_file <- file.path(dirname(path_stats), "Figure4_Full_Composite_A4_Fixed.pdf")

cairo_pdf(out_file, width = 8.27, height = 11.69)

grid.newpage()
pushViewport(viewport(layout = grid.layout(3, 1, heights = c(1.2, 1, 0.8))))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(ht, newpage = FALSE, heatmap_legend_side = "bottom", 
     padding = unit(c(10, 2, 2, 2), "mm"))
popViewport()

print(p_bubble, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))

print(p_ridge, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))

dev.off()

message("Success! File saved at: ", out_file)


library(readxl)
library(dplyr)
library(tidyr)
library(ggraph)
library(tidygraph)
library(igraph)
library(stringr)

path_master <- path_chembl
path_stats <- out_xlsx

message("Loading and processing data...")
df_rdkit <- read.csv(path_rdkit)
df_master <- read_excel(path_master, sheet = "MASTER_DATA")

df_merged <- df_master %>%
  filter(!is.na(pActivity), !is.na(Target_Category_Macro)) %>%
  filter(!Target_Category_Macro %in% c("Unclassified", "Miscellaneous", "ADMET")) %>%
  mutate(inchikey_clean = str_trim(inchikey)) %>%
  inner_join(df_rdkit, by = c("inchikey" = "inchikey")) %>%
  mutate(
    Motif_Type = case_when(
      has_prenyl_like %in% c("True", TRUE) ~ "Prenylated",
      has_methoxy_aryl %in% c("True", TRUE) ~ "Methoxylated",
      has_probable_sugar %in% c("True", TRUE) ~ "Glycosylated",
      TRUE ~ "Simple"
    ),
    Node_Label = str_squish(paste(Motif_Type, class_np)) 
  )

top_compounds_per_node <- df_merged %>%
  group_by(Node_Label) %>%
  arrange(desc(pActivity)) %>% 
  slice(1) %>% 
  ungroup()

top_15_nodes <- top_compounds_per_node %>%
  arrange(desc(pActivity)) %>%
  head(15)

edges <- df_merged %>%
  group_by(Node_Label, Target_Category_Macro) %>%
  summarise(weight = max(pActivity), n_hits = n(), .groups = "drop") %>% 
  filter(n_hits >= 3) %>% 
  rename(from = Node_Label, to = Target_Category_Macro)

chem_nodes <- df_merged %>%
  group_by(Node_Label) %>%
  summarise(
    lead_potency = max(pActivity),
    motif = first(Motif_Type), 
    count = n(), 
    .groups = "drop"
  ) %>%
  mutate(name = str_squish(Node_Label)) %>%
  filter(name %in% edges$from) %>%
  left_join(top_15_nodes %>% select(Node_Label, best_inchikey = inchikey), by = c("name" = "Node_Label")) %>%
  mutate(type = "Chemical")

target_nodes <- edges %>% select(name = to) %>% distinct() %>% 
  mutate(lead_potency = NA, motif = "Target", count = NA, best_inchikey = NA, type = "Target")

nodes <- bind_rows(chem_nodes, target_nodes) %>% distinct(name, .keep_all = TRUE)

graph <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE) %>%
  activate(nodes) %>%
  mutate(degree = centrality_degree()) %>%
  mutate(
    label_category = case_when(
      type == "Target" ~ "Target",
      !is.na(best_inchikey) ~ "Top15_Red",
      degree > quantile(degree, 0.90, na.rm=TRUE) ~ "Hub_Simple",
      TRUE ~ "None"
    ),
    final_label = case_when(
      label_category == "Target" ~ name,
      label_category == "Top15_Red" ~ paste0(name, "\n", substr(best_inchikey, 1, 9), "..."),
      label_category == "Hub_Simple" ~ name,
      TRUE ~ NA_character_
    )
  )

set.seed(42) 
p_net <- ggraph(graph, layout = 'fr') + 
  
  geom_edge_link(aes(width = weight, alpha = n_hits, color = weight)) + 
  scale_edge_color_gradient(low = "gray92", high = "gray10", guide = "none") +
  scale_edge_width(range = c(0.6, 2.5), guide="none") + 
  scale_edge_alpha(range = c(0.4, 0.9), guide="none") +
  
  geom_node_point(aes(filter = type == "Target", size = degree), 
                  shape = 22, fill = "black", color = "white", stroke = 1.5) +
  
  geom_node_point(aes(filter = type == "Chemical", 
                      fill = lead_potency, 
                      shape = motif, 
                      size = degree), 
                  color = "black", stroke = 0.5, alpha = 0.9) +
  
  geom_node_text(aes(label = final_label, filter = label_category == "Target"), 
                 repel = TRUE, fontface = "bold", size = 5, bg.color = "white") +
  
  geom_node_text(aes(label = final_label, filter = label_category == "Top15_Red"), 
                 repel = TRUE, size = 3.5, fontface = "bold", color = "red4",
                 bg.color = "white", bg.r = 0.1, min.segment.length = 0, force = 30) +
  
  geom_node_text(aes(label = final_label, filter = label_category == "Hub_Simple"), 
                 repel = TRUE, size = 3, fontface = "plain", color = "gray20") +
  
  scale_fill_gradientn(
    colors = c("blue3", "#FAA000", "red4"), 
    name = "Example Potency", 
    na.value = "gray95"
  ) +
  
  scale_shape_manual(values = c(
    "Simple" = 21, "Prenylated" = 24, "Methoxylated" = 23, "Glycosylated" = 25, "Target" = 22
  )) +
  
  scale_size_continuous(range = c(4, 13), guide="none") +
  theme_void() +
  theme(legend.position = "right", plot.title = element_text(face="bold", size=16)) +
  labs(title = "Figure 5: Structure-Activity Relationship Network")

out_file <- file.path(dirname(path_stats), "Fig5_SAR_Network_Final_Colors.pdf")
ggsave(out_file, plot = p_net, width = 15, height = 12, device = cairo_pdf)

message("Network figure saved with lead-potency coloring and the custom palette: ", out_file)


library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(stringr)

path_master <- path_chembl
out_dir <- file.path(OUT_DIR, "PartIII_ALL")

message("Loading input tables...")
df_rdkit <- read.csv(path_rdkit)
df_master <- read_excel(path_master, sheet = "MASTER_DATA")

df_stat <- df_master %>%
  filter(!is.na(pActivity), !is.na(Target_Category_Macro)) %>%
  filter(!Target_Category_Macro %in% c("Unclassified", "Miscellaneous", "ADMET")) %>%
  mutate(inchikey_clean = str_trim(inchikey)) %>%
  inner_join(df_rdkit, by = c("inchikey" = "inchikey")) %>%
  mutate(
    Motif_Type = case_when(
      has_prenyl_like == "True" | has_prenyl_like == TRUE ~ "Prenylated",
      has_methoxy_aryl == "True" | has_methoxy_aryl == TRUE ~ "Methoxylated",
      has_probable_sugar == "True" | has_probable_sugar == TRUE ~ "Glycosylated",
      TRUE ~ "Simple"
    )
  )

df_stat$Motif_Type <- factor(df_stat$Motif_Type, levels = c("Simple", "Glycosylated", "Methoxylated", "Prenylated"))

summary_stats <- df_stat %>%
  group_by(Target_Category_Macro, Motif_Type) %>%
  summarise(
    N_Compounds = n(),
    Median_Potency = median(pActivity, na.rm = TRUE),
    Mean_Potency = mean(pActivity, na.rm = TRUE),
    SD_Potency = sd(pActivity, na.rm = TRUE),
    .groups = "drop"
  )

baseline_stats <- summary_stats %>%
  filter(Motif_Type == "Simple") %>%
  select(Target_Category_Macro, Base_Median = Median_Potency)

final_stats_table <- summary_stats %>%
  left_join(baseline_stats, by = "Target_Category_Macro") %>%
  mutate(
    Potency_Gain_vs_Simple = Median_Potency - Base_Median,
    Fold_Change_Linear = 10^(Potency_Gain_vs_Simple)
  ) %>%
  arrange(Target_Category_Macro, desc(Median_Potency))

write.csv(final_stats_table, file.path(out_dir, "Table_Statistics_Potency_Gain.csv"), row.names = FALSE)

stat_test <- df_stat %>%
  group_by(Target_Category_Macro) %>%
  wilcox_test(pActivity ~ Motif_Type, ref.group = "Simple") %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>%
  mutate(xy.position = 9)

write.csv(stat_test, file.path(out_dir, "Table_Significance_Tests.csv"), row.names = FALSE)


p_stats <- ggboxplot(df_stat, x = "Motif_Type", y = "pActivity",
                     fill = "Motif_Type", 
                     facet.by = "Target_Category_Macro",
                     outlier.shape = NA,
                     alpha = 0.6) +
  
  geom_jitter(aes(color = Motif_Type), alpha = 0.2, width = 0.2, size = 0.5) +
  
  stat_compare_means(method = "wilcox.test", ref.group = "Simple", 
                     label = "p.signif", label.y = 10, hide.ns = TRUE) +
  
  scale_fill_manual(values = c("Simple"="gray80", "Glycosylated"="#440154", "Methoxylated"="#21908C", "Prenylated"="#FDE725")) +
  scale_color_manual(values = c("Simple"="gray60", "Glycosylated"="#440154", "Methoxylated"="#21908C", "Prenylated"="#FDE725")) +
  
  labs(title = "Statistical Impact of Structural Motifs on Bioactivity",
       subtitle = "Wilcoxon test vs. Simple scaffold (* p<0.05, **** p<0.0001)",
       y = "Potency (pActivity)", x = "Structural Class") +
  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none",
        strip.background = element_rect(fill = "white", color="black"),
        strip.text = element_text(face="bold"))

ggsave(file.path(out_dir, "Fig_Statistical_Boxplot.pdf"), plot = p_stats, width = 12, height = 8)

message("Analysis completed.")
message("1. Potency gain table: Table_Statistics_Potency_Gain.csv")
message("2. Significance tests: Table_Significance_Tests.csv")
message("3. Figure: Fig_Statistical_Boxplot.pdf")

