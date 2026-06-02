# Workflow

Run all commands from the repository root after editing `configs/config.yaml`.

## 00 Atlas Construction

```bash
Rscript scripts/00_atlas_construction/00_run_multiome_atlas_construction.R
```

Builds the multiome atlas, runs WNN integration, clustering, annotation, and atlas-level visualization.

## 01 Feature Matrix Preparation

```bash
Rscript scripts/01_feature_matrix_preparation/00_export_rna_protein_coding_matrix.R
Rscript scripts/01_feature_matrix_preparation/01_filter_atac_peak_matrix.R
python3 scripts/01_feature_matrix_preparation/02_convert_atac_csv_to_parquet.py
```

Exports RNA and ATAC feature matrices for ML.

## 02 ML Feature Selection

```bash
python3 scripts/02_ml_feature_selection/00_train_rna_xgboost_borutashap.py
python3 scripts/02_ml_feature_selection/01_train_atac_xgboost_borutashap.py
```

Runs XGBoost/BorutaShap and exports STAGE genes/peaks.

## 03 Model Evaluation

```bash
Rscript scripts/03_model_evaluation/00_plot_rna_stage_feature_upset.R
Rscript scripts/03_model_evaluation/01_plot_atac_stage_feature_upset.R
python3 scripts/03_model_evaluation/02_summarize_rna_model_performance.py
python3 scripts/03_model_evaluation/03_summarize_atac_model_performance.py
```

Summarizes model performance and STAGE feature overlap.

## 04 DE Method Benchmark

```bash
Rscript scripts/04_de_method_benchmark/00_run_de_method_benchmark.R
```

Compares single-cell Wilcoxon and donor-level pseudobulk DESeq2 before downstream interpretation. This module supports DE benchmark results and related supplementary tables; it is not the main ML training workflow.

## 05 Functional Enrichment

```bash
Rscript scripts/05_functional_enrichment/00_run_rna_de_for_enrichment.R
Rscript scripts/05_functional_enrichment/01_run_rna_go_kegg_enrichment.R
Rscript scripts/05_functional_enrichment/02_run_atac_great_enrichment.R
Rscript scripts/05_functional_enrichment/03_plot_functional_enrichment.R
Rscript scripts/05_functional_enrichment/04_plot_kegg_matrix_bubble.R
```

Runs RNA GO/KEGG and ATAC GREAT enrichment.

## 06 ATAC Motif Peak Gene

```bash
Rscript scripts/06_atac_motif_peak_gene/00_run_chromvar_motif_activity.R
Rscript scripts/06_atac_motif_peak_gene/01_map_motifs_to_stage_peaks.R
Rscript scripts/06_atac_motif_peak_gene/02_run_peak_to_gene_links.R
Rscript scripts/06_atac_motif_peak_gene/03_build_tf_peak_gene_network.R
Rscript scripts/06_atac_motif_peak_gene/04_plot_tf_peak_gene_network.R
Rscript scripts/06_atac_motif_peak_gene/05_plot_tf_peak_gene_sankey.R
Rscript scripts/06_atac_motif_peak_gene/06_plot_tf_gene_circle_network.R
```

Runs chromVAR, motif/STAGE peak mapping, peak-to-gene links, and TF-peak-gene regulatory network analysis.

## 07 Cell Communication

```bash
Rscript scripts/07_cell_communication/00_run_multinichenet.R
Rscript scripts/07_cell_communication/01_plot_multinichenet_results.R
```

Runs MultiNicheNet ligand-receptor prioritization and plotting.

## 08 Pseudotime

```bash
python3 scripts/08_pseudotime/00_compute_scriskcell_scores.py
Rscript scripts/08_pseudotime/01_run_monocle3_trajectory.R
Rscript scripts/08_pseudotime/02_run_tradeseq_dynamic_genes.R
Rscript scripts/08_pseudotime/03_plot_pseudotime_results.R
```

Runs scRiskCell, monocle3, tradeSeq, and pseudotime plotting.

## 09 Figures And Tables

```bash
python3 scripts/09_figures_and_tables/00_make_supp_table5_rna_stage_genes.py
python3 scripts/09_figures_and_tables/01_make_supp_table5_atac_stage_peaks.py
python3 scripts/09_figures_and_tables/02_make_differential_motif_activity_table.py
python3 scripts/09_figures_and_tables/03_make_table8_peak_to_gene_links.py
python3 scripts/09_figures_and_tables/04_make_table9_candidate_motif_peak_gene_triplets.py
python3 scripts/09_figures_and_tables/05_make_supp_table10_multinichenet_prioritized_lr.py
python3 scripts/09_figures_and_tables/06_make_supp_table11_multinichenet_ligand_target_links.py
python3 scripts/09_figures_and_tables/07_make_supp_table12_tradeseq_pseudotime_dynamic_genes.py
```

Generates final supplementary tables from regenerated outputs.

