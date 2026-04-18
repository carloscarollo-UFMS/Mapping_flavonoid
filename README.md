# Comparative flavonoid chemotaxonomy across plant lineages reveals structural specialization and target-dependent functional patterns

This repository contains the code and supporting data resources associated with the manuscript:

**Comparative flavonoid chemotaxonomy across plant lineages reveals structural specialization and target-dependent functional patterns**

The workflow integrates:

- MongoDB-based extraction from a curated LOTUS subset
- taxonomic curation and normalization with World Flora Online
- bioactivity integration
- RDKit-based structural annotation
- macro-group statistical analyses
- phylogeny/iTOL-oriented formatting

The analytical workflow is implemented in **R** and **Python**.

---

## Repository contents

### Scripts

- `scripts/Main_Pipeline_end.R`
- `scripts/Part I_Extraction.R`
- `scripts/Part II - Bio_iTOL_Prep.R`
- `scripts/3_annotate_flavonoids_rdkit.py`
- `scripts/Part III - Figures_Stats.R`
- `scripts/Tree_APG_ITOL.R`

### Required files

- `inputs/classification.tsv`
- `mongo/subset_minor_flavonoids_V.bson`
- `mongo/subset_minor_flavonoids_V.metadata.json`

### Additional files for full reproduction of Part II global context

- `mongo/lotusUniqueNaturalProduct.bson`
- `mongo/lotusUniqueNaturalProduct.metadata.json`

---

## Recommended repository structure

```text
.
├─ README.md
├─ scripts/
│  ├─ Main_Pipeline_end.R
│  ├─ Part I_Extraction.R
│  ├─ Part II - Bio_iTOL_Prep.R
│  ├─ 3_annotate_flavonoids_rdkit.py
│  ├─ Part III - Figures_Stats.R
│  └─ Tree_APG_ITOL.R
├─ inputs/
│  └─ classification.tsv
├─ mongo/
│  ├─ subset_minor_flavonoids_V.bson
│  ├─ subset_minor_flavonoids_V.metadata.json
│  ├─ lotusUniqueNaturalProduct.bson
│  └─ lotusUniqueNaturalProduct.metadata.json
├─ phylo_outputs_FINAL/
│  └─ APG_clade_assignments.csv
└─ outputs/
```

---

## Software requirements

- **R**
- **Python 3**
- **MongoDB Community Server**
- **MongoDB Database Tools** (`mongorestore`)

---

## R dependencies

The workflow uses packages including:

- `here`
- `arrow`
- `dplyr`
- `mongolite`
- `jsonlite`
- `progress`
- `stringr`
- `stringi`
- `writexl`
- `readr`
- `readxl`
- `httr`
- `ggplot2`
- `grid`
- `gridExtra`
- `scales`
- `openxlsx`
- `janitor`
- `ggrepel`
- `vegan`
- `rstatix`
- `multcompView`
- `broom`
- `ape`
- `V.PhyloMaker2`

Example installation block:

```r
install.packages(c(
  "here", "arrow", "dplyr", "mongolite", "jsonlite", "progress",
  "stringr", "stringi", "writexl", "readr", "readxl", "httr",
  "ggplot2", "gridExtra", "scales", "openxlsx", "janitor",
  "ggrepel", "vegan", "rstatix", "multcompView", "broom", "ape"
))
```

If `V.PhyloMaker2` is not available from your default CRAN mirror, install it separately according to its official instructions.

---

## Python dependencies

The RDKit annotation step requires:

- `pandas`
- `rdkit`

Example installation:

```bash
pip install pandas rdkit
```

---

## MongoDB configuration

The workflow is configured around:

- **database:** `lotusdb`
- **primary collection:** `subset_minor_flavonoids_V`

Example restore command:

```bash
mongorestore --db lotusdb --collection subset_minor_flavonoids_V ./mongo
```

The scripts use the environment variable `LOTUS_MONGO_URL`. A typical local value is:

```text
mongodb://127.0.0.1:27017
```

---

## Important note about Part II

`subset_minor_flavonoids_V` is sufficient for the subset-based workflow and for Part I.

However, **Part II is not fully reproduced with this subset alone**. One section of Part II also queries the additional MongoDB collection `lotusUniqueNaturalProduct` in the same database (`lotusdb`) to recover global occurrence context across plant families.

Therefore:

- the subset collection is sufficient for the core extraction and downstream subset-based analyses;
- full reproduction of the global-context step in Part II requires `lotusUniqueNaturalProduct`;
- if `lotusUniqueNaturalProduct` is not installed locally, the global occurrence layer used in Part II will be incomplete.

---

## World Flora Online input

Taxonomic normalization in Part I uses the World Flora Online backbone through `classification.tsv`.

Place the file at:

```text
inputs/classification.tsv
```

Alternatively, provide an explicit path through `cfg$wfo_csv_path` in `Main_Pipeline_end.R`.

Example:

```r
wfo_csv_path = "C:/Users/Carollo/Documents/Maira_Bio/data/wfo/classification.tsv"
```

---

## Part III macro-group mapping

Part III requires a family-level mapping table with three columns:

1. `family`
2. `color`
3. `macro_group`

In the current version of the script, this file is handled automatically:

- if `groups.xlsx` already exists, it is used directly;
- if it does not exist, Part III searches for `APG_clade_assignments.csv`;
- if found, Part III automatically generates `inputs/groups.xlsx` and continues normally.

The expected source for this automatic generation is:

```text
phylo_outputs_FINAL/APG_clade_assignments.csv
```

This file is produced by `Tree_APG_ITOL.R` and provides the family identity, clade assignment, and color mapping needed by Part III.

---

## Workflow logic

The current workflow is organized as follows:

1. **Part I** extracts, curates, and normalizes records from MongoDB
2. **Part II** prepares the flavonoid subset for structural annotation and integrates bioactivity metadata
3. **Python RDKit annotation** computes structural descriptors and motif flags
4. **Tree_APG_ITOL** generates APG/clade assignments
5. **Part III** generates macro-group statistics and figures

Recommended execution order:

1. `Part I_Extraction.R`
2. `Part II - Bio_iTOL_Prep.R`
3. `3_annotate_flavonoids_rdkit.py`
4. `Tree_APG_ITOL.R`
5. `Part III - Figures_Stats.R`

`Tree_APG_ITOL.R` should be run before Part III when `groups.xlsx` is not already available, because it generates `APG_clade_assignments.csv`, which Part III can convert automatically into the required macro-group mapping table.

---

## Main pipeline script

`Main_Pipeline_end.R` centralizes:

- taxonomic scope
- database settings
- module execution switches
- output directory logic

Users intending to run the workflow through the main driver script should review the configuration block and activate the desired modules before execution.

Part II and Part III expect `cfg` and `runtime` objects created by the main pipeline. Running the full workflow through a single R session is therefore recommended.

---

## Expected outputs

The workflow generates run-specific outputs under `outputs/`, including:

- curated analytical tables
- Excel summaries
- RDKit input and annotation tables
- statistical result tables
- PDF figures
- phylogeny/iTOL-compatible files

Typical outputs include:

### Part I
- run-specific LOTUS Excel workbook
- `lin_enriched`
- `uni_enriched`
- exported baseline datasets

### Part II
- `PartII_ALL/<base_tag>__flavonoids_for_rdkit.csv`
- `PartII_ALL/<base_tag>_BIO_A_Summary.xlsx`
- `PartII_ALL/<base_tag>_BIO_C_Global_Context.xlsx`
- `PartII_ALL/<base_tag>_Lotus_Final_Database_v9_Fixed.xlsx`

### Python step
- `lotus_flavonoids_rdkit_annotations.csv`

### Tree/APG
- `phylo_outputs_FINAL/APG_clade_assignments.csv`

### Part III
- `PartIII_ALL/macro_groups_STATS_MASTER.xlsx`
- `PartIII_ALL/Fig1_ChemistryStats.pdf`
- `PartIII_ALL/Fig2_OrdinationStats.pdf`
- `PartIII_ALL/Fig3_BioprospectingLandscape.pdf`

---

## Notes on reproducibility

This repository provides the resources needed to reproduce the subset-based comparative flavonoid chemotaxonomy workflow associated with the manuscript.

Complete reproduction of the **global occurrence context** used in Part II additionally depends on the MongoDB collection `lotusUniqueNaturalProduct`.

The macro-group mapping used by Part III is automatically reproducible from the APG assignment export, removing the need for a manually supplied `groups.xlsx` file when `APG_clade_assignments.csv` is available.

---

## Citation

If you use this repository, please cite the associated manuscript and any linked data resources.
