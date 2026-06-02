# Cleanup Report

## Summary

Created a cleaned GitHub-ready repository with 37 main analysis scripts and 58 files total. The cleaned repository keeps code, configuration, documentation, dependency lists, and two small example CSV files only. Large data, generated results, object files, figures, model artifacts, logs, and caches were excluded.

## Created Directories

- `configs`
- `envs`
- `scripts/00_atlas_construction`
- `scripts/01_feature_matrix_preparation`
- `scripts/02_ml_feature_selection`
- `scripts/03_model_evaluation`
- `scripts/04_de_method_benchmark`
- `scripts/05_functional_enrichment`
- `scripts/06_atac_motif_peak_gene`
- `scripts/07_cell_communication`
- `scripts/08_pseudotime`
- `scripts/09_figures_and_tables`
- `src/utils`
- `src/preprocessing`
- `src/modeling`
- `src/enrichment`
- `src/atac_analysis`
- `src/communication`
- `src/visualization`
- `examples`
- `data`
- `results`
- `docs`

No historical-version directory was created.

## Included Core Scripts

- Atlas construction: `scripts/00_atlas_construction/00_run_multiome_atlas_construction.R`
- Feature matrix preparation: three RNA/ATAC export and conversion scripts.
- ML feature selection: RNA and ATAC XGBoost/BorutaShap scripts.
- Model evaluation: RNA/ATAC UpSet and model performance scripts.
- DE method benchmark: `scripts/04_de_method_benchmark/00_run_de_method_benchmark.R`
- Functional enrichment: RNA DE, GO/KEGG, GREAT, and enrichment plotting scripts.
- ATAC regulatory analysis: chromVAR, motif/STAGE peak mapping, peak-to-gene linking, TF-peak-gene network construction and plotting scripts.
- Cell communication: MultiNicheNet run and plotting scripts.
- Pseudotime: scRiskCell, monocle3, tradeSeq, and pseudotime plotting scripts.
- Figures and tables: supplementary table generation scripts.

## Parameterized Scripts

Repeated downstream script families were consolidated into representative parameterized public scripts:

- GO/KEGG enrichment
- GREAT enrichment
- chromVAR motif activity
- peak-to-gene links
- motif-to-STAGE peak mapping
- TF-peak-gene network construction
- MultiNicheNet communication
- scRiskCell scoring
- monocle3/tradeSeq pseudotime analysis
- supplementary table generation

The cleanup preserved cell type labels, comparison names, thresholds, model parameters, statistical methods, input/output fields, and analysis order from the original scripts. Path handling was changed to use `configs/config.yaml` and repository-relative locations.

## Excluded Scripts And Files

Excluded categories:

- Raw data
- Processed data
- Seurat-like objects and other serialized R objects
- Parquet matrices
- Trained model files
- Generated figures
- Generated result tables
- Generated supplementary tables
- Logs and cache files
- Temporary no-extension files
- Historical or duplicate scripts

Historical versions marked in `docs/file_manifest.md`:

- `multinichenet copy.R`
- `multinichenet_v2.R`
- `plot_v2/plot.R`
- Reference WNN script not selected as the public atlas source

## Example Files

Kept two small example files:

- `examples/stage_genes_example.csv`
- `examples/stage_peaks_example.csv`

These examples demonstrate file formats only. Full results are excluded and should be regenerated from scripts.

## Path And Config Cleanup

Added:

- `configs/config.yaml`
- `src/utils/config.R`
- `src/utils/config.py`

Replaced old local absolute paths in public scripts with config-derived variables and repository-relative paths for:

- raw data
- processed data
- Seurat object
- RNA expression parquet
- ATAC peak parquet
- model outputs
- DE benchmark outputs
- functional enrichment outputs
- GREAT outputs
- STAGE gene/peak outputs
- ATAC regulatory outputs
- communication outputs
- pseudotime outputs
- supplementary table outputs
- external NicheNet resources

Chinese comments and messages were removed from public scripts. Greek cell type symbols were preserved where they are part of biological labels.

## Generated Documentation

- `README.md`
- `docs/code_overview.md`
- `docs/workflow.md`
- `docs/file_manifest.md`
- `docs/cleanup_report.md`
- `data/README.md`
- `results/README.md`
- `examples/README.md`
- Module README files for functional enrichment, ATAC regulatory analysis, cell communication, and pseudotime.

## Generated Environment Files

- `envs/requirements.txt`
- `envs/environment.yml`
- `envs/R_packages.txt`

Dependencies were inferred from actual `import`, `library()`, and `require()` statements. No reliable package versions were found in the scanned scripts, so package names were listed without pinned versions.

## Quality Checks

Commands and results:

```bash
du -sh .
# 940K .
```

```bash
find . -type f -size +5M -print
# No files found.
```

Local absolute path scan:

```bash
# Passed: matches only in docs/file_manifest.md original path entries.
```

```bash
find . \( -name "*.qs" -o -name "*.rds" -o -name "*.RDS" -o -name "*.parquet" -o -name "*.pkl" -o -name "*.joblib" -o -name "*.xlsx" -o -name "*.log" -o -name "*.pdf" -o -name "*.png" \) -print
# No files found.
```

```bash
Rscript -e 'files <- list.files("scripts", pattern="[.]R$", recursive=TRUE, full.names=TRUE); for (f in files) parse(f)'
# Passed.
```

```bash
python3 -m py_compile $(find scripts src -type f -name "*.py")
# Passed.
```

```bash
grep -R "00_preprocessing_wNN/\|01_cell_annotation/\|04_atac_motif_peak_gene/02_run_great_enrichment.R\|07_figures_and_tables/02_run_go_kegg_enrichment.R\|archive/" README.md docs/workflow.md docs/code_overview.md
# No matches.
```

## Remaining Manual Confirmation

- Add license information.
- Add manuscript citation after publication.
- Add contact information.
- Add GEO/accession and data download instructions.
- Confirm whether the two example CSV files are acceptable for public release.
- Confirm external NicheNet resource download instructions.
- Confirm any manuscript-specific figure/table naming conventions.

## GitHub Pre-Upload Checklist

- Review `configs/config.yaml` placeholders.
- Review `docs/file_manifest.md` for provenance and excluded-file status.
- Review example CSV files for public-release suitability.
- Confirm no generated large outputs are staged.
- Run `git status` before committing.

## Recommended Git Commands

```bash
cd clean_github_repo

git init
git status
git add .
git commit -m "Initial cleaned code release"
git branch -M main
git remote add origin <YOUR_GITHUB_REPOSITORY_URL>
git push -u origin main
```
