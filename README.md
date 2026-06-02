# T2D Single-Cell Multiomics

This repository provides the cleaned analysis code for a human pancreatic islet single-cell multiome study across control (CON), prediabetes (PRE), and type 2 diabetes (T2D) states.

The codebase contains publication-oriented workflows for scRNA-seq and scATAC-seq preprocessing, WNN atlas construction, cell type annotation, machine-learning-based feature selection, DE method benchmarking, functional enrichment, ATAC regulatory analysis, cell communication analysis, pseudotime analysis, and supplementary table generation.

## Overview

The analysis workflow is organized into ordered modules:

1. `00_atlas_construction`
   Load multiome data, perform RNA/ATAC preprocessing, construct the WNN atlas, run dimensionality reduction and clustering, annotate cell types, and generate atlas-level visualizations.

2. `01_feature_matrix_preparation`
   Export RNA expression matrices and ATAC peak matrices for downstream machine-learning analyses.

3. `02_ml_feature_selection`
   Run XGBoost/BorutaShap models and identify stage-transition associated genes and peaks, referred to as STAGE genes and STAGE peaks.

4. `03_model_evaluation`
   Summarize model performance, feature stability, and feature overlap across cell types and disease-stage comparisons.

5. `04_de_method_benchmark`
   Compare single-cell-level Wilcoxon testing with donor-level pseudobulk DESeq2 to evaluate conventional differential expression strategies.

6. `05_functional_enrichment`
   Perform RNA-based GO/KEGG enrichment and ATAC-associated GREAT enrichment analyses.

7. `06_atac_motif_peak_gene`
   Analyze chromVAR motif activity, motif-to-STAGE peak mapping, peak-to-gene links, and TF–peak–gene regulatory networks.

8. `07_cell_communication`
   Run MultiNicheNet ligand–receptor and ligand–target analyses.

9. `08_pseudotime`
   Perform scRiskCell, monocle3, and tradeSeq analyses to characterize trajectory-associated molecular dynamics.

10. `09_figures_and_tables`
    Generate final summary plots and supplementary tables.

## Repository Structure

```text
configs/                  Path and analysis configuration files
envs/                     Python and R dependency lists
scripts/                  Ordered analysis modules
src/                      Shared helper functions and utilities
data/                     Placeholder directory; input data are not included
results/                  Placeholder directory; outputs are regenerated locally
docs/                     Workflow documentation, file manifest, and cleanup report
```

## Data Availability

TODO: Add GEO/accession and data download instructions.

Raw sequencing data, processed matrices, Seurat objects, model files, generated figures, and large result tables are not included in this repository.

Users should prepare the required input files locally and configure their paths in:

```text
configs/config.yaml
```

The default paths in `configs/config.yaml` are repository-relative placeholders, such as:

```text
data/processed
results
```

## Installation

### Python environment

Install Python dependencies with:

```bash
pip install -r envs/requirements.txt
```

Alternatively, create the conda environment:

```bash
conda env create -f envs/environment.yml
conda activate t2d_single_cell_multiomics
```

### R environment

R package dependencies are listed in:

```text
envs/R_packages.txt
```

Bioconductor packages should be installed with Bioconductor-aware tools, such as `BiocManager`.

## Configuration

Before running the workflow, edit:

```text
configs/config.yaml
```

to specify local paths for input data, processed objects, external resources, and output directories.

The configuration file controls paths for major inputs and outputs, including:

```text
raw data
processed data
Seurat objects
RNA expression matrices
ATAC peak matrices
model outputs
DE benchmark results
functional enrichment results
ATAC regulatory analysis outputs
cell communication results
pseudotime results
supplementary tables
```

## Step-by-Step Usage

Run scripts from the repository root in the following order.

### 1. Atlas construction

```bash
Rscript scripts/00_atlas_construction/00_run_multiome_atlas_construction.R
```

### 2. Feature matrix preparation

```bash
Rscript scripts/01_feature_matrix_preparation/00_export_rna_protein_coding_matrix.R
Rscript scripts/01_feature_matrix_preparation/01_filter_atac_peak_matrix.R
python3 scripts/01_feature_matrix_preparation/02_convert_atac_csv_to_parquet.py
```

### 3. Machine-learning feature selection

```bash
python3 scripts/02_ml_feature_selection/00_train_rna_xgboost_borutashap.py
python3 scripts/02_ml_feature_selection/01_train_atac_xgboost_borutashap.py
```

### 4. DE method benchmark

```bash
Rscript scripts/04_de_method_benchmark/00_run_de_method_benchmark.R
```

Additional downstream modules are documented in:

```text
docs/workflow.md
```

## Expected Outputs

Outputs are written to the configured `results/` subdirectories. These may include:

```text
model performance summaries
STAGE gene and peak tables
DE benchmark outputs
functional enrichment results
motif activity summaries
peak-to-gene link tables
TF–peak–gene regulatory networks
MultiNicheNet communication results
pseudotime-associated dynamic gene tables
figures
supplementary tables
```

Generated outputs are not tracked by Git and should be regenerated locally by running the scripts.

## Large Files Not Included

This repository intentionally excludes large files and generated outputs, including:

```text
raw sequencing data
processed count matrices
Seurat objects
AnnData-like objects
parquet matrices
trained model files
intermediate analysis objects
generated figures
generated supplementary tables
large CSV/TSV/XLSX result files
logs and cache files
```

This design keeps the repository lightweight and focused on reproducible analysis code.

## Documentation

Additional documentation is available in:

```text
docs/workflow.md
docs/code_overview.md
docs/file_manifest.md
docs/cleanup_report.md
```

These files describe the full workflow, script-level inputs and outputs, original file sources, and cleanup decisions.

## Citation

TODO: Add manuscript citation after publication.

## License

TODO: Add license information.

## Contact

TODO: Add contact information.
