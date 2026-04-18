# Tree_APG_ITOL: phylogeny reconstruction and APG clade export for iTOL.

`%||%` <- function(a, b) {
  if (is.null(a)) return(b)
  if (length(a) == 0) return(b)
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
      "install.packages(c(", paste(sprintf("\"%s\"", missing), collapse = ", "), "), dependencies = TRUE)\n"
    )
  }
  invisible(TRUE)
}

normalize_path <- function(x) normalizePath(x, winslash = "/", mustWork = FALSE)

find_project_root <- function(start = getwd()) {
  candidates <- unique(normalize_path(c(start, file.path(start, ".."), file.path(start, "../.."))))
  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "scripts", "APGIV_family_order_clades_WorldFlora.csv"))) {
      return(candidate)
    }
    if (dir.exists(file.path(candidate, "scripts")) && dir.exists(file.path(candidate, "outputs"))) {
      return(candidate)
    }
  }
  normalize_path(start)
}

find_first_existing <- function(paths) {
  paths <- as.character(paths)
  paths <- paths[!is.na(paths) & nzchar(paths)]
  hit <- paths[file.exists(paths)]
  if (length(hit)) normalize_path(hit[1]) else NA_character_
}

resolve_file <- function(label, candidates, required = TRUE) {
  hit <- find_first_existing(candidates)
  if (is.na(hit) && isTRUE(required)) {
    stop("Required file not found [", label, "]. Tried:\n- ", paste(candidates, collapse = "\n- "))
  }
  hit
}

rename_first_match <- function(df, pattern, new_name) {
  if (new_name %in% names(df)) return(df)
  hits <- grep(pattern, names(df), ignore.case = TRUE)
  if (!length(hits)) {
    stop("Required column not found for pattern: ", pattern)
  }
  names(df)[hits[1]] <- new_name
  df
}

trim_chr <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x
}

normalize_taxon <- function(x) {
  x <- trim_chr(x)
  keep <- !is.na(x)
  x[keep] <- paste0(toupper(substr(x[keep], 1, 1)), tolower(substr(x[keep], 2, nchar(x[keep]))))
  x
}

extract_genus <- function(labels) {
  vapply(
    strsplit(as.character(labels), "_", fixed = TRUE),
    function(parts) parts[1],
    FUN.VALUE = character(1)
  )
}

kv <- function(key, value) paste(c(key, value), collapse = "\t")

translate_clade <- function(raw_clade, translation_map) {
  raw_clade <- trim_chr(raw_clade)
  translated <- rep("Other", length(raw_clade))

  keep_raw <- !is.na(raw_clade) & raw_clade %in% unname(translation_map)
  translated[keep_raw] <- raw_clade[keep_raw]

  needs_map <- !is.na(raw_clade) & raw_clade %in% names(translation_map)
  translated[needs_map] <- unname(translation_map[raw_clade[needs_map]])

  translated
}

get_final_clade <- function(family_name, node1, node2, node3, node4, order_name) {
  if (family_name == "Namaceae") return("Asterids (Lamiids)")
  if (family_name == "Vitaceae") return("Rosids (Basal)")
  if (family_name == "Cynomoriaceae") return("Saxifragales")

  if (family_name %in% c(
    "Ophioglossaceae", "Helminthostachyaceae", "Botrychiaceae", "Psilotaceae",
    "Pteridaceae", "Equisetaceae", "Marattiaceae", "Osmundaceae",
    "Hymenophyllaceae", "Gleicheniaceae", "Cyatheaceae", "Dryopteridaceae",
    "Polypodiaceae", "Aspleniaceae", "Salviniaceae", "Marsileaceae",
    "Dicksoniaceae"
  )) {
    return("Ferns")
  }

  if (family_name %in% c(
    "Pinaceae", "Araucariaceae", "Podocarpaceae", "Sciadopityaceae", "Taxaceae",
    "Cephalotaxaceae", "Cupressaceae", "Cycadaceae", "Zamiaceae",
    "Ginkgoaceae", "Gnetaceae", "Ephedraceae", "Welwitschiaceae"
  )) {
    return("Gymnosperms")
  }

  if (family_name %in% c("Lycopodiaceae", "Selaginellaceae", "Isoetaceae")) {
    return("Lycophytes")
  }

  if (family_name %in% c(
    "Polytrichaceae", "Aulacomniaceae", "Hypnaceae", "Mniaceae", "Dicranaceae",
    "Bryaceae", "Leucobryaceae", "Marchantiaceae", "Sphagnaceae"
  )) {
    return("Bryophytes")
  }

  for (node in list(node4, node3, node2, node1)) {
    node <- trim_chr(node)
    if (!is.na(node)) return(node)
  }

  order_name <- trim_chr(order_name)
  if (identical(order_name, "Boraginales")) return("Asterids (Lamiids)")
  if (identical(order_name, "Vitales")) return("Rosids (Basal)")
  if (identical(order_name, "Saxifragales")) return("Saxifragales")

  "Other"
}

write_tree_colors <- function(df, path_out) {
  lines <- c(
    "TREE_COLORS",
    "SEPARATOR TAB",
    "DATA",
    paste(df$ID, "branch", df$Color, "normal", sep = "\t"),
    paste(df$ID, "label", df$Color, "bold", "1.5", sep = "\t")
  )
  writeLines(lines, path_out)
}

write_color_strip <- function(df, palette, path_out) {
  used_colors <- palette[names(palette) %in% unique(df$FinalClade)]
  header <- c(
    "DATASET_COLORSTRIP",
    "SEPARATOR TAB",
    kv("DATASET_LABEL", "APG IV Clades"),
    kv("COLOR", "#000000"),
    kv("STRIP_WIDTH", "100"),
    kv("LEGEND_TITLE", "Major Clades"),
    kv("LEGEND_SHAPES", rep("1", length(used_colors))),
    kv("LEGEND_COLORS", as.character(used_colors)),
    kv("LEGEND_LABELS", names(used_colors)),
    "DATA"
  )
  body <- paste(df$ID, df$Color, df$FinalClade, sep = "\t")
  writeLines(c(header, body), path_out)
}

cfg <- if (exists("cfg")) cfg else list()

pkgs <- c("ape", "V.PhyloMaker2", "readxl", "readr")
require_pkgs(pkgs)

suppressPackageStartupMessages({
  library(ape)
  library(V.PhyloMaker2)
  library(readxl)
  library(readr)
})

project_root <- normalize_path(cfg$project_root %||% find_project_root())
outputs_dir <- file.path(project_root, "outputs")

lotus_candidates <- c(
  cfg$path_lotus %||% NA_character_,
  if (dir.exists(outputs_dir)) {
    list.files(
      outputs_dir,
      pattern = "lotus_kingdom_Plantae.*\\.xlsx$",
      recursive = TRUE,
      full.names = TRUE
    )
  } else {
    character(0)
  }
)

path_lotus <- resolve_file("LOTUS workbook", lotus_candidates, required = TRUE)
path_apg <- resolve_file(
  "APG IV reference table",
  c(
    cfg$path_apg %||% NA_character_,
    file.path(project_root, "scripts", "APGIV_family_order_clades_WorldFlora.csv"),
    file.path(project_root, "APGIV_family_order_clades_WorldFlora.csv"),
    file.path(getwd(), "APGIV_family_order_clades_WorldFlora.csv")
  ),
  required = TRUE
)

out_dir <- normalize_path(cfg$out_dir %||% file.path(project_root, "phylo_outputs_FINAL"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("----------------------------------------------------------------")
message(">>> [1/5] Loading source tables")
message("----------------------------------------------------------------")

lin_raw <- as_base_df(
  readxl::read_excel(path_lotus, sheet = "lin_enriched", .name_repair = "minimal")
)
lin_raw <- rename_first_match(lin_raw, "^family$", "family_name")
lin_raw <- rename_first_match(lin_raw, "^genus$", "genus_name")
lin_raw$family_name <- normalize_taxon(lin_raw$family_name)
lin_raw$genus_name <- normalize_taxon(lin_raw$genus_name)

keep_rows <- !is.na(lin_raw$family_name) & nzchar(lin_raw$family_name) &
  !is.na(lin_raw$genus_name) & nzchar(lin_raw$genus_name)

sp_list <- lin_raw[keep_rows, c("family_name", "genus_name"), drop = FALSE]
sp_list <- sp_list[!duplicated(sp_list$family_name), , drop = FALSE]
sp_list$species_dummy <- paste0(sp_list$genus_name, "_sp")
sp_list <- as_base_df(sp_list)

if (nrow(sp_list) == 0) {
  stop("No valid family/genus records were found in sheet 'lin_enriched'.")
}

apg_db <- as_base_df(readr::read_csv(path_apg, show_col_types = FALSE))
if (!"Family" %in% names(apg_db)) {
  stop("Column 'Family' was not found in the APG reference table.")
}

for (col_name in c("Order", "Node.1", "Node.2", "Node.3", "Node.4")) {
  if (!col_name %in% names(apg_db)) apg_db[[col_name]] <- NA_character_
}

apg_db$Family <- normalize_taxon(apg_db$Family)

message("----------------------------------------------------------------")
message(">>> [2/5] Assigning APG clades")
message("----------------------------------------------------------------")

match_idx <- match(sp_list$family_name, apg_db$Family)
classified_df <- cbind(
  sp_list,
  apg_db[match_idx, c("Order", "Node.1", "Node.2", "Node.3", "Node.4"), drop = FALSE]
)
classified_df <- as_base_df(classified_df)

classified_df$raw_clade <- vapply(
  seq_len(nrow(classified_df)),
  function(i) {
    get_final_clade(
      family_name = classified_df$family_name[i],
      node1 = classified_df$Node.1[i],
      node2 = classified_df$Node.2[i],
      node3 = classified_df$Node.3[i],
      node4 = classified_df$Node.4[i],
      order_name = classified_df$Order[i]
    )
  },
  FUN.VALUE = character(1)
)

clade_translation <- c(
  "Fabids" = "Rosids (Fabids)",
  "Malvids" = "Rosids (Malvids)",
  "Lamiids" = "Asterids (Lamiids)",
  "Campanulids" = "Asterids (Campanulids)",
  "Superrosids" = "Saxifragales",
  "Superasterids" = "Caryophyllales",
  "Asterids" = "Asterids (Basal)",
  "Rosids" = "Rosids (Basal)",
  "Commelinids" = "Monocots (Commelinids)",
  "Monocots" = "Monocots",
  "Eudicots" = "Basal Eudicots",
  "Magnoliids" = "Magnoliids",
  "Basal Angiosperms" = "Basal Angiosperms",
  "Gymnosperms" = "Gymnosperms",
  "Ferns" = "Ferns",
  "Lycophytes" = "Lycophytes",
  "Bryophytes" = "Bryophytes",
  "Saxifragales" = "Saxifragales",
  "Santalales" = "Santalales",
  "Caryophyllales" = "Caryophyllales"
)

classified_df$FinalClade <- translate_clade(classified_df$raw_clade, clade_translation)
classified_df <- as_base_df(classified_df)

message("----------------------------------------------------------------")
message(">>> [3/5] Building phylogeny")
message("----------------------------------------------------------------")

sp_input <- data.frame(
  species = sp_list$species_dummy,
  genus = sp_list$genus_name,
  family = sp_list$family_name,
  stringsAsFactors = FALSE
)

sp_input$family[sp_input$genus %in% c("Helminthostachys", "Botrychium", "Ophioglossum", "Mankyua")] <- "Ophioglossaceae"
sp_input$family[sp_input$family %in% c("Helminthostachyaceae", "Botrychiaceae")] <- "Ophioglossaceae"
sp_input$family[sp_input$family %in% c("Namaceae", "Heliotropiaceae", "Cordiaceae", "Ehretiaceae")] <- "Boraginaceae"
sp_input$family[sp_input$genus == "Cephalotaxus"] <- "Taxaceae"
sp_input$family[sp_input$genus %in% c("Sambucus", "Viburnum")] <- "Adoxaceae"
sp_input$family[sp_input$genus == "Cleome"] <- "Cleomaceae"
sp_input <- as_base_df(sp_input)

data("GBOTB.extended.TPL", package = "V.PhyloMaker2", envir = environment())
data("nodes.info.1.TPL", package = "V.PhyloMaker2", envir = environment())

phylo_result <- V.PhyloMaker2::phylo.maker(
  sp_input,
  tree = GBOTB.extended.TPL,
  nodes = nodes.info.1.TPL,
  scenarios = "S3"
)

tree_temp <- phylo_result$scenario.3
if (is.null(tree_temp)) {
  stop("PhyloMaker did not return a scenario.3 tree.")
}

message("----------------------------------------------------------------")
message(">>> [4/5] Restoring original family labels")
message("----------------------------------------------------------------")

tip_genera <- extract_genus(tree_temp$tip.label)
genus_to_family <- stats::setNames(sp_list$family_name, sp_list$genus_name)
restored_labels <- unname(genus_to_family[tip_genera])
missing_labels <- is.na(restored_labels) | !nzchar(restored_labels)
restored_labels[missing_labels] <- tree_temp$tip.label[missing_labels]

tree_temp$tip.label <- restored_labels
tree_final <- ape::ladderize(tree_temp)

path_newick <- file.path(out_dir, "newick_end.nwk")
ape::write.tree(tree_final, file = path_newick)

message("----------------------------------------------------------------")
message(">>> [5/5] Writing iTOL outputs")
message("----------------------------------------------------------------")

clade_palette <- c(
  "Bryophytes" = "#7F7F7F",
  "Lycophytes" = "#7F7F7F",
  "Ferns" = "#505050",
  "Gymnosperms" = "#000000",
  "Basal Angiosperms" = "#E5C494",
  "Magnoliids" = "#FFD92F",
  "Monocots" = "#7570B3",
  "Monocots (Commelinids)" = "#7570B3",
  "Basal Eudicots" = "#FC8D62",
  "Saxifragales" = "#1B9E77",
  "Rosids (Basal)" = "#8FBC8F",
  "Rosids (Fabids)" = "#66A61E",
  "Rosids (Malvids)" = "#E7298A",
  "Caryophyllales" = "#E31A1C",
  "Santalales" = "#FF7F00",
  "Asterids (Basal)" = "#FF7F00",
  "Asterids (Lamiids)" = "#377EB8",
  "Asterids (Campanulids)" = "#984EA3",
  "Other" = "#000000"
)

export_df <- data.frame(ID = tree_final$tip.label, stringsAsFactors = FALSE)
match_export <- match(export_df$ID, classified_df$family_name)
export_df$FinalClade <- classified_df$FinalClade[match_export]
export_df$FinalClade[is.na(export_df$FinalClade)] <- "Other"
export_df$Color <- unname(clade_palette[export_df$FinalClade])
export_df$Color[is.na(export_df$Color)] <- "#000000"
export_df <- as_base_df(export_df)

utils::write.csv(export_df, file.path(out_dir, "APG_clade_assignments.csv"), row.names = FALSE)
write_tree_colors(export_df, file.path(out_dir, "iTOL_tree_colors.txt"))
write_color_strip(export_df, clade_palette, file.path(out_dir, "iTOL_clade_strip.txt"))

check_families <- c("Namaceae", "Vitaceae", "Cynomoriaceae", "Helminthostachyaceae", "Ophioglossaceae")
check_df <- export_df[export_df$ID %in% check_families, c("ID", "FinalClade"), drop = FALSE]
if (nrow(check_df) > 0) {
  message("[INFO] Selected family assignments:")
  print(check_df, row.names = FALSE)
}

message("Completed. Files written to: ", out_dir)
