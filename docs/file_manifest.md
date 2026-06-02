# File Manifest

Status values:

- `included`: copied into the cleaned repository.
- `example`: copied as a small non-sensitive example file.
- `excluded`: not copied.
- `excluded_historical_version`: historical, duplicate, or previous-attempt script; not copied.

| Original path | New path | File type | Status | Module | Manual confirmation | Notes |
|---|---|---|---|---|---|---|
| `/home/jywei/wnn_v3/plot/wnn.R` | `scripts/00_atlas_construction/00_run_multiome_atlas_construction.R` | R script | included | 00_atlas_construction | no | Confirmed by user as public WNN/atlas source script. |
| `/home/jywei/wnn_v3/wnn.R` |  | R script | excluded_historical_version | 00_atlas_construction | no | Reference version only; not used as main public workflow. |
| `/home/jywei/ml_code/RNA/export_expr_protein_coding_filtered.R` | `scripts/01_feature_matrix_preparation/00_export_rna_protein_coding_matrix.R` | R script | included | 01_feature_matrix_preparation | no | RNA protein-coding matrix export. |
| `/home/jywei/ml_code/ATAC/atac_filter.R` | `scripts/01_feature_matrix_preparation/01_filter_atac_peak_matrix.R` | R script | included | 01_feature_matrix_preparation | no | ATAC feature matrix filtering. |
| `/home/jywei/ml_code/ATAC/csv2par.py` | `scripts/01_feature_matrix_preparation/02_convert_atac_csv_to_parquet.py` | Python script | included | 01_feature_matrix_preparation | no | ATAC CSV-to-parquet conversion. |
| `/home/jywei/ml_code/RNA/RNA_ml.py` | `scripts/02_ml_feature_selection/00_train_rna_xgboost_borutashap.py` | Python script | included | 02_ml_feature_selection | no | RNA XGBoost/BorutaShap workflow. |
| `/home/jywei/ml_code/ATAC/ATAC_ml.py` | `scripts/02_ml_feature_selection/01_train_atac_xgboost_borutashap.py` | Python script | included | 02_ml_feature_selection | no | ATAC XGBoost/BorutaShap workflow. |
| `/home/jywei/ml_code/upset/RNA_upset.R` | `scripts/03_model_evaluation/00_plot_rna_stage_feature_upset.R` | R script | included | 03_model_evaluation | no | RNA STAGE feature overlap plotting. |
| `/home/jywei/ml_code/upset/ATAC_upset.R` | `scripts/03_model_evaluation/01_plot_atac_stage_feature_upset.R` | R script | included | 03_model_evaluation | no | ATAC STAGE feature overlap plotting. |
| `/home/jywei/supplementary_tables/make_supp_table4_ml_performance.py` | `scripts/03_model_evaluation/02_summarize_rna_model_performance.py` | Python script | included | 03_model_evaluation | no | RNA model performance summary. |
| `/home/jywei/supplementary_tables/make_supp_table4_ATAC_performance.py` | `scripts/03_model_evaluation/03_summarize_atac_model_performance.py` | Python script | included | 03_model_evaluation | no | ATAC model performance summary. |
| `/home/jywei/wnn_v3/de_method_benchmark/de_method_benchmark.R` | `scripts/04_de_method_benchmark/00_run_de_method_benchmark.R` | R script | included | 04_de_method_benchmark | no | Independent DE method benchmark module. |
| `/home/jywei/downstream_v3/*/GO_RNA/DE.R` | `scripts/05_functional_enrichment/00_run_rna_de_for_enrichment.R` | R script family | included | 05_functional_enrichment | no | Parameterized from representative GO_RNA DE scripts. |
| `/home/jywei/downstream_v3/*/GO_RNA/GO.R` | `scripts/05_functional_enrichment/01_run_rna_go_kegg_enrichment.R` | R script family | included | 05_functional_enrichment | no | Parameterized from representative GO/KEGG scripts. |
| `/home/jywei/downstream_v3/*/great/great.R` | `scripts/05_functional_enrichment/02_run_atac_great_enrichment.R` | R script family | included | 05_functional_enrichment | no | Parameterized GREAT/localGREAT workflow. |
| `/home/jywei/downstream_v3/*/GO_RNA/plot.R` | `scripts/05_functional_enrichment/03_plot_functional_enrichment.R` | R script family | included | 05_functional_enrichment | no | Functional enrichment plotting. |
| `/home/jywei/downstream_v3/acinar_CON_PRE/GO_RNA/acinar_KEGG_matrix_bubble_plot.R` | `scripts/05_functional_enrichment/04_plot_kegg_matrix_bubble.R` | R script | included | 05_functional_enrichment | no | Selected KEGG matrix bubble plot. |
| `/home/jywei/downstream_v3/*/chromvar/chromvar.R` | `scripts/06_atac_motif_peak_gene/00_run_chromvar_motif_activity.R` | R script family | included | 06_atac_motif_peak_gene | no | Parameterized chromVAR motif activity workflow. |
| `/home/jywei/downstream_v3/*/mlpeaks_chromvar_mapping/map_chromvar_to_MLpeaks.R` | `scripts/06_atac_motif_peak_gene/01_map_motifs_to_stage_peaks.R` | R script family | included | 06_atac_motif_peak_gene | no | Motif-to-STAGE peak mapping. |
| `/home/jywei/downstream_v3/*/linkpeak/*.R` | `scripts/06_atac_motif_peak_gene/02_run_peak_to_gene_links.R` | R script family | included | 06_atac_motif_peak_gene | no | Peak-to-gene linking. |
| `/home/jywei/downstream_v3/*/netlink/netlink.R` | `scripts/06_atac_motif_peak_gene/03_build_tf_peak_gene_network.R` | R script family | included | 06_atac_motif_peak_gene | no | TF-peak-gene triplet construction. |
| `/home/jywei/downstream_v3/*/netlink/plot.R` | `scripts/06_atac_motif_peak_gene/04_plot_tf_peak_gene_network.R` | R script family | included | 06_atac_motif_peak_gene | no | Regulatory network plotting. |
| `/home/jywei/downstream_v3/*/netlink/sankey.R` | `scripts/06_atac_motif_peak_gene/05_plot_tf_peak_gene_sankey.R` | R script family | included | 06_atac_motif_peak_gene | no | Sankey plotting. |
| `/home/jywei/downstream_v3/*/netlink/cirlce.R` | `scripts/06_atac_motif_peak_gene/06_plot_tf_gene_circle_network.R` | R script family | included | 06_atac_motif_peak_gene | no | Circle network plotting; original filename typo retained only in source reference. |
| `/home/jywei/downstream_v3/*/multinichenet/multinichenet.R` | `scripts/07_cell_communication/00_run_multinichenet.R` | R script family | included | 07_cell_communication | no | Parameterized MultiNicheNet workflow. |
| `/home/jywei/downstream_v3/*/multinichenet/plot/plot.R` | `scripts/07_cell_communication/01_plot_multinichenet_results.R` | R script family | included | 07_cell_communication | no | MultiNicheNet plotting. |
| `/home/jywei/downstream_v3/riskcell/*/1.py` | `scripts/08_pseudotime/00_compute_scriskcell_scores.py` | Python script family | included | 08_pseudotime | no | Parameterized scRiskCell workflow. |
| `/home/jywei/downstream_v3/riskcell_monocle/*/1.R` | `scripts/08_pseudotime/01_run_monocle3_trajectory.R` | R script family | included | 08_pseudotime | no | Parameterized monocle3 trajectory workflow. |
| `/home/jywei/downstream_v3/riskcell_monocle/beta/union_native/1.R` | `scripts/08_pseudotime/02_run_tradeseq_dynamic_genes.R` | R script | included | 08_pseudotime | no | tradeSeq dynamic gene analysis. |
| `/home/jywei/downstream_v3/riskcell_monocle/acinar/plot_beautified_from_qs.R` | `scripts/08_pseudotime/03_plot_pseudotime_results.R` | R script | included | 08_pseudotime | no | Pseudotime plotting. |
| `/home/jywei/supplementary_tables/make_supp_table5_RNA_stage_genes.py` | `scripts/09_figures_and_tables/00_make_supp_table5_rna_stage_genes.py` | Python script | included | 09_figures_and_tables | no | Supplementary STAGE gene table generation. |
| `/home/jywei/supplementary_tables/make_supp_table5_ATAC_stage_peaks.py` | `scripts/09_figures_and_tables/01_make_supp_table5_atac_stage_peaks.py` | Python script | included | 09_figures_and_tables | no | Supplementary STAGE peak table generation. |
| `/home/jywei/supplementary_tables/make_supp_table_differential_motif_activity.py` | `scripts/09_figures_and_tables/02_make_differential_motif_activity_table.py` | Python script | included | 09_figures_and_tables | no | Differential motif activity table. |
| `/home/jywei/supplementary_tables/make_table8_peak_to_gene_links.py` | `scripts/09_figures_and_tables/03_make_table8_peak_to_gene_links.py` | Python script | included | 09_figures_and_tables | no | Peak-to-gene supplementary table. |
| `/home/jywei/supplementary_tables/make_table9_candidate_motif_peak_gene_triplets.py` | `scripts/09_figures_and_tables/04_make_table9_candidate_motif_peak_gene_triplets.py` | Python script | included | 09_figures_and_tables | no | Regulatory triplet supplementary table. |
| `/home/jywei/supplementary_tables/make_supp_table10_multinichenet_prioritized_lr.py` | `scripts/09_figures_and_tables/05_make_supp_table10_multinichenet_prioritized_lr.py` | Python script | included | 09_figures_and_tables | no | MultiNicheNet LR supplementary table. |
| `/home/jywei/supplementary_tables/make_supp_table11_multinichenet_ligand_target_links.py` | `scripts/09_figures_and_tables/06_make_supp_table11_multinichenet_ligand_target_links.py` | Python script | included | 09_figures_and_tables | no | Ligand-target supplementary table. |
| `/home/jywei/supplementary_tables/make_supp_table12_tradeseq_pseudotime_dynamic_genes.py` | `scripts/09_figures_and_tables/07_make_supp_table12_tradeseq_pseudotime_dynamic_genes.py` | Python script | included | 09_figures_and_tables | no | tradeSeq supplementary table. |
| `/home/jywei/ml_code/RNA/STAGE_genes/Ductal_cell__PRE_vs_T2D__test_shap_geq0p6.csv` | `examples/stage_genes_example.csv` | CSV | example | examples | yes | Small format example only. |
| `/home/jywei/ml_code/ATAC/STAGE_peaks/Ductal_cell/CON_vs_T2D_shap_marker_genes.csv` | `examples/stage_peaks_example.csv` | CSV | example | examples | yes | Small format example only. |
| `/home/jywei/downstream_v3/alpha_CON_PRE/multinichenet/multinichenet copy.R` |  | R script | excluded_historical_version | 07_cell_communication | no | Explicitly excluded by user. |
| `/home/jywei/downstream_v3/alpha_PRE_T2D/multinichenet/multinichenet_v2.R` |  | R script | excluded_historical_version | 07_cell_communication | no | Explicitly excluded by user. |
| `/home/jywei/downstream_v3/beta_CON_PRE/multinichenet/plot_v2/plot.R` |  | R script | excluded_historical_version | 07_cell_communication | no | Explicitly excluded by user. |
| `/home/jywei/wnn_v3/atac_peaks_autosomes.parquet` |  | parquet | excluded | data | no | Large processed ATAC matrix. |
| `/home/jywei/wnn_v3/expr_protein_coding_filtered.parquet` |  | parquet | excluded | data | no | Large processed RNA matrix. |
| `/home/jywei/downstream_v3/**/*.qs` |  | QS/R object | excluded | results | no | Large Seurat/chromVAR/MultiNicheNet/pseudotime objects. |
| `/home/jywei/downstream_v3/**/*.rds` |  | RDS object | excluded | results | no | Large or generated intermediate objects. |
| `/home/jywei/downstream_v3/**/*.parquet` |  | parquet | excluded | results | no | Intermediate matrices. |
| `/home/jywei/downstream_v3/**/*.png` |  | figure | excluded | results | no | Generated figures. |
| `/home/jywei/downstream_v3/**/*.pdf` |  | figure | excluded | results | no | Generated figures. |
| `/home/jywei/downstream_v3/**/*.log` |  | log | excluded | logs | no | Runtime logs. |
| `/home/jywei/supplementary_tables/tables/*.csv` |  | generated table | excluded | 09_figures_and_tables | yes | Generated supplementary outputs; regenerate from scripts. |
| `/home/jywei/supplementary_tables/tables/*.xlsx` |  | generated table | excluded | 09_figures_and_tables | yes | Generated supplementary outputs; regenerate from scripts. |
| `/home/jywei/supplementary_tables/__pycache__/*.pyc` |  | cache | excluded | cache | no | Python bytecode cache. |

