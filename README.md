# T2D Single-Cell Multiomics

Cleaned analysis code release for a human pancreatic islet single-cell multiome study across CON, PRE, and T2D. The repository contains publication-oriented scripts for scRNA-seq/scATAC-seq preprocessing, WNN atlas construction, cell type annotation, machine-learning feature selection, DE method benchmarking, functional enrichment, ATAC regulatory analyses, cell communication, pseudotime analysis, and supplementary table generation.

## Workflow Overview

1. `00_atlas_construction`: load multiome data, run RNA/ATAC preprocessing, WNN integration, clustering, annotation, and atlas-level plots.
2. `01_feature_matrix_preparation`: export RNA and ATAC feature matrices for machine learning.
3. `02_ml_feature_selection`: run XGBoost/BorutaShap and export STAGE genes/peaks.
4. `03_model_evaluation`: summarize model performance and feature overlap.
5. `04_de_method_benchmark`: compare single-cell Wilcoxon and donor-level pseudobulk DESeq2.
6. `05_functional_enrichment`: run RNA GO/KEGG and ATAC GREAT enrichment.
7. `06_atac_motif_peak_gene`: run chromVAR, motif/STAGE peak mapping, peak-to-gene links, and TF-peak-gene networks.
8. `07_cell_communication`: run MultiNicheNet ligand-receptor analyses.
9. `08_pseudotime`: run scRiskCell, monocle3, and tradeSeq analyses.
10. `09_figures_and_tables`: generate final summary plots and supplementary tables.

## Repository Structure

```text
configs/                  Path and analysis configuration
envs/                     Python/R dependency lists
scripts/                  Ordered analysis modules
src/                      Shared utilities
examples/                 Small non-sensitive format examples
data/                     Placeholder only; data are not included
results/                  Placeholder only; outputs are regenerated
docs/                     Workflow, manifest, and cleanup documentation
```

## Data Availability

TODO: Add GEO/accession and data download instructions.

Raw data, processed matrices, Seurat objects, model files, generated figures, and large result tables are not included. Configure local paths in `configs/config.yaml` before running the workflow.

## Installation

Python dependencies:

```bash
pip install -r envs/requirements.txt
```

Conda environment:

```bash
conda env create -f envs/environment.yml
conda activate t2d_single_cell_multiomics
```

R dependencies are listed in `envs/R_packages.txt`. Install Bioconductor packages with Bioconductor-aware tooling such as `BiocManager`.

## Configuration

Edit `configs/config.yaml` to point to prepared local data and output directories. The default paths are repository-relative placeholders such as `data/processed` and `results`.

## Step-By-Step Usage

Run scripts from the repository root in this order:

```bash
Rscript scripts/00_atlas_construction/00_run_multiome_atlas_construction.R
Rscript scripts/01_feature_matrix_preparation/00_export_rna_protein_coding_matrix.R
Rscript scripts/01_feature_matrix_preparation/01_filter_atac_peak_matrix.R
python3 scripts/01_feature_matrix_preparation/02_convert_atac_csv_to_parquet.py
python3 scripts/02_ml_feature_selection/00_train_rna_xgboost_borutashap.py
python3 scripts/02_ml_feature_selection/01_train_atac_xgboost_borutashap.py
Rscript scripts/04_de_method_benchmark/00_run_de_method_benchmark.R
```

Continue with the downstream modules documented in `docs/workflow.md`.

## Expected Outputs

Outputs are written under configured `results` subdirectories, including model summaries, DE benchmark outputs, enrichment tables, motif activity summaries, regulatory network tables, communication results, pseudotime tables, figures, and supplementary tables.

## Large Files Not Included

The repository intentionally excludes raw data, processed data, Seurat/AnnData-like objects, parquet matrices, trained models, generated figures, and generated supplementary tables.

## Citation

TODO: Add manuscript citation after publication.

## License

TODO: Add license information.

## Contact

TODO: Add contact information.

