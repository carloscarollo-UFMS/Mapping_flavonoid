# Taxonomic distribution, structural diversity, and functional evidence of selected flavonoid subclasses across vascular plant families

This repository contains the analytical workflow, processed datasets, audit files, statistical outputs, and figure source files associated with the manuscript:

> **Taxonomic distribution, structural diversity, and functional evidence of selected flavonoid subclasses across vascular plant families**

The project integrates taxonomically normalized occurrence records from LOTUS, RDKit-based structural annotation, a branch-length-free family-level visualization tree, and curated ChEMBL bioactivity records. The analyses were designed to distinguish documented chemical patterns from unequal sampling and selective functional testing.

## Study scope

The analysis focuses on 17 selected and harmonized flavonoid subclasses:

- 2-Arylbenzofurans
- Anthocyanidins
- Aurones
- Chalcones
- Coumestans
- Dihydroflavonols
- Flavan-3-ols
- Flavandiols
- Flavanones
- Flavans
- Flavonolignans
- Isoflavanones
- Isoflavones
- Neoflavonoids
- Proanthocyanidins
- Pterocarpans
- Rotenoids

Flavones and flavonols were excluded a priori because together they represented 53.37% of the source-level LOTUS flavonoid universe and could obscure variation among the selected subclasses. This analytical exclusion does not imply lower biological, ecological, or chemotaxonomic importance.

## Validated analytical dataset

The manuscript-level workflow produced:

- 10,863 structures assigned to the 17 selected and harmonized subclasses
- 10,679 structures with retrievable plant occurrence records
- 10,496 taxonomically resolved compounds
- 10,469 compounds eligible for the primary structural analyses
- 33,109 deduplicated documentary occurrence records
- 25,900 distinct compound–taxon associations
- 213 vascular plant families in the family-level visualization tree
- 131 families with at least 10 eligible compounds in the primary family-level inferential analyses
- 2,064 compounds mapped to ChEMBL identifiers
- 1,003 compounds with eligible primary quantitative bioactivity records
- 1,000 compounds contributing to the final functional aggregation

The primary family-level threshold was **n ≥ 10 compounds**, with sensitivity analyses at **n ≥ 5** and **n ≥ 20**.

## Analytical workflow

The workflow is organized into five connected stages:

1. **LOTUS extraction and taxonomic normalization**  
   Extraction from MongoDB, World Flora Online-assisted taxonomic resolution, occurrence deduplication, and audit-table generation.

2. **Bioactivity integration and structural-input preparation**  
   ChEMBL identifier mapping, activity retrieval and filtering, functional-domain assignment, binary-evidence classification, and preparation of the RDKit input table.

3. **RDKit structural annotation**  
   Structure standardization, quality control, physicochemical descriptors, Bemis–Murcko scaffolds, structural motifs, and compound-level audit files.

4. **Family-level visualization tree**  
   Deterministic representative-species selection, V.PhyloMaker2 placement, APG-based lineage annotation, topology auditing, and iTOL exports. The tree is used as a family-level visualization framework rather than as a dated phylogeny.

5. **Statistical analyses and figures**  
   Family-level multivariate analyses, sampling-adjusted scaffold richness, functional evidence summaries, Fsp3 and motif–potency models, sensitivity analyses, and the descriptive structure–bioactivity evidence network.

## Repository structure

```text
.
├── README.md
├── CITATION.cff
├── Main_Pipeline_end.R
├── scripts/
│   ├── Part_I_Extraction_v2_1.R
│   ├── Part_II_Bio_iTOL_Prep_v2_1.R
│   ├── 3_annotate_flavonoids_rdkit_auto.py
│   ├── Tree_APG_ITOL_v2_1.R
│   ├── Part_III_Figures_Stats_v2_1.R
│   └── APGIV_family_order_clades_WorldFlora.csv
├── outputs/
│   └── lotus_kingdom_Plantae_20260622/
│       ├── PartI_ALL/
│       ├── PartII_ALL/
│       ├── RDKit_ALL/
│       ├── Tree_APG_iTOL/
│       ├── PartIII_ALL/
│       ├── lotus_kingdom_Plantae_20260622_lin_enriched.parquet
│       ├── lotus_kingdom_Plantae_20260622_uni_enriched.parquet
│       └── lotus_kingdom_Plantae_20260622_lin_compound_species.parquet
└── phylo_outputs_FINAL/
    ├── family_tree_s3.nwk
    ├── newick_end.nwk
    ├── family_representatives.csv
    ├── APG_clade_assignments.csv
    └── tree and taxonomic audit files
```

## Principal scripts

### `Main_Pipeline_end.R`

Creates the `cfg` and `runtime` objects, defines the run directory, and controls execution of the R modules. Review the configuration block before running the workflow, particularly:

- database and collection names;
- taxonomic scope;
- input and output paths;
- module switches;
- the primary threshold of `analysis_min_compounds_per_taxon = 10L`;
- the sensitivity thresholds `c(5L, 10L, 20L)`;
- the validated random seed `20260622L`.

### `scripts/Part_I_Extraction_v2_1.R`

Performs LOTUS extraction, hierarchical taxonomic resolution, deduplication, and audit exports.

### `scripts/Part_II_Bio_iTOL_Prep_v2_1.R`

Maps compounds to ChEMBL, retrieves and filters bioactivity records, assigns functional domains, generates binary evidence summaries, and exports the input used by RDKit.

### `scripts/3_annotate_flavonoids_rdkit_auto.py`

Standardizes molecular structures and generates compound- and occurrence-level structural annotations. The main compact compound-level output is:

```text
RDKit_ALL/lotus_flavonoids_rdkit_compounds.csv
```

### `scripts/Tree_APG_ITOL_v2_1.R`

Constructs the branch-length-free family-level visualization tree and exports Newick, iTOL, APG-assignment, representative-species, substitution, exclusion, and topology-audit files.

### `scripts/Part_III_Figures_Stats_v2_1.R`

Runs the validated manuscript analyses. The script requires a primary family threshold of **n ≥ 10**, uses sensitivity thresholds of 5, 10, and 20 compounds, applies 999 permutations in multivariate tests, and uses the random seed `20260622`.

## Software requirements

The validated analytical run used:

- R 4.5.3
- Python 3.14.5
- RDKit 2026.03.3
- pandas 3.0.3
- MongoDB Community Server and MongoDB Database Tools for reconstruction from the raw LOTUS collections

Important R packages include:

```text
arrow, here, mongolite, jsonlite, dplyr, tidyr, readr, readxl,
writexl, openxlsx, stringr, stringi, progress, httr, purrr,
tibble, ggplot2, ggrepel, patchwork, scales, vegan, mgcv,
rstatix, multcompView, broom, sandwich, lmtest, ComplexHeatmap,
circlize, ggridges, igraph, tidygraph, ggraph, gridExtra, ape,
V.PhyloMaker2
```

Install the Python dependencies with:

```bash
python -m pip install pandas rdkit
```

Install CRAN R packages as needed with `install.packages()`. `V.PhyloMaker2` should be installed according to its official distribution instructions when it is not available from the configured CRAN repository.

## Reproducing the analyses

### 1. Clone the repository

```bash
git clone https://github.com/carloscarollo-UFMS/Mapping_flavonoid.git
cd Mapping_flavonoid
```

### 2. Configure the R controller

Open `Main_Pipeline_end.R` and verify database settings, external-data paths, module switches, the run tag, and the output directory.

The R modules depend on the `cfg` and `runtime` objects created by the controller. Running the modules in isolated R sessions without reproducing these objects is not recommended.

### 3. Run Parts I and II when reconstructing from source data

Execute the controller with the relevant module switches enabled:

```bash
Rscript Main_Pipeline_end.R
```

Part I requires a local MongoDB instance containing the configured LOTUS collection and the World Flora Online classification backbone. Part II additionally requires internet access to query ChEMBL unless the validated cached or processed outputs are used.

### 4. Run the RDKit annotation step

```bash
python scripts/3_annotate_flavonoids_rdkit_auto.py \
  --input outputs/lotus_kingdom_Plantae_20260622/PartII_ALL/lotus_kingdom_Plantae_20260622__flavonoids_for_rdkit.csv \
  --output-dir outputs/lotus_kingdom_Plantae_20260622/RDKit_ALL
```

The script exports compact compound-level annotations, occurrence-level annotations, invalid-structure and conflict audits, structural-QC records, and run metadata.

### 5. Generate or audit the family-level tree

Run `Tree_APG_ITOL_v2_1.R` through the project configuration used by the controller. The validated tree and all corresponding audit files are also provided in `phylo_outputs_FINAL/`.

### 6. Run the final statistical analyses

Set the Part III module to active in `Main_Pipeline_end.R` and execute:

```bash
Rscript Main_Pipeline_end.R
```

The final outputs are written to:

```text
outputs/lotus_kingdom_Plantae_20260622/PartIII_ALL/
```

## Reproduction from the committed processed data

Users who do not need to repeat the raw MongoDB extraction can begin from the committed processed Parquet, RDKit, taxonomic, tree, and ChEMBL-derived tables. These files support inspection and reproduction of the downstream manuscript statistics and figures without repeating every database query.

Public-database content can change over time. Exact reconstruction of the validated run should therefore use the committed processed tables and metadata rather than newly retrieved ChEMBL or taxonomic records.

## External and large data files

Raw MongoDB dumps and the complete World Flora Online classification backbone are not stored directly in this repository because of their size. Their download location and required local paths are documented in this README and in `Main_Pipeline_end.R`; no separate `external_data/` directory is required.

The occurrence-level file `lotus_flavonoids_rdkit_annotations.csv` is also omitted from the normal GitHub repository because it exceeds the standard single-file limit. The smaller compound-level file and the corresponding metadata and audit tables are provided for downstream reproduction.

When external files are used, record their version, retrieval date, and cryptographic hash. Do not commit local credentials, `.Renviron`, MongoDB connection strings, or private access tokens.

## Main outputs

### Part I

- taxonomically normalized occurrence tables;
- WFO resolution summaries and conflict audits;
- deduplication summaries;
- family-coverage summaries;
- processed Parquet datasets.

### Part II

- ChEMBL mapping and filtering flow;
- curated quantitative bioactivity tables;
- target and assay functional-domain classifications;
- binary-evidence summaries;
- compound–family, compound–species, and compound–subclass maps;
- RDKit input file and run metadata.

### RDKit

- compound-level structural descriptors and motifs;
- Bemis–Murcko scaffolds;
- primary structural-eligibility flag;
- review and exclusion audits;
- software and run metadata.

### Family-level tree

- validated Newick files;
- deterministic representative-species table;
- APG lineage assignments;
- iTOL annotation files;
- backbone substitutions, exclusions, and topology corrections;
- run metadata.

### Part III

- manuscript figure components;
- family-level chemical and structural matrices;
- global and pairwise PERMANOVA results;
- multivariate-dispersion tests;
- sampling-adjusted scaffold-richness models;
- univariate tests and multiple-testing correction;
- functional-domain and subclass summaries;
- Fsp3 and motif–potency models;
- descriptive network edges and nodes;
- sensitivity analyses, model objects, and audit files.

## Journal supplementary material

The formal Supplementary Information, supplementary tables, and supplementary figures are distributed through the journal article webpage and are not duplicated in this repository. The repository retains the processed tables, audit files, model outputs, and figure source files required to trace and reproduce the reported analyses.

After publication, the article DOI and direct link to the journal-hosted supplementary material should be added to this section.

## Interpretation limits

This repository documents the available evidence in LOTUS and ChEMBL; it does not represent a complete inventory of plant metabolomes or biological activities. Important limitations include:

- unequal phytochemical documentation among families and subclasses;
- selective testing and reporting of bioactivity;
- heterogeneous assay, target, and measurement contexts;
- incomplete taxonomic or database coverage;
- descriptive rather than causal interpretation of the structure–bioactivity evidence network.

The absence of indexed ChEMBL evidence should be interpreted as a functional annotation gap, not as evidence of inactivity.

## Citation

A machine-readable citation file is provided as `CITATION.cff`. Until the associated article receives its final bibliographic information, cite the repository using its GitHub URL and release version.

Repository:

```text
https://github.com/carloscarollo-UFMS/Mapping_flavonoid
```

When the manuscript is published, its DOI and complete citation should be added to this section and to `CITATION.cff`.

## Contact

For questions about the workflow or data organization, open an issue in this repository or contact the corresponding author through the institutional information provided in the associated manuscript.
