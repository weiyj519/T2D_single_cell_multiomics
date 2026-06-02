#!/usr/bin/env Rscript

# Script name: 00_run_multinichenet.R
# Purpose: Cleaned public workflow script for the T2D single-cell multiome analysis.
# Input: Paths and analysis settings are read from configs/config.yaml.
# Output: Module-specific outputs are written under the configured results directory.
# Main steps: Load configuration, run the original analysis logic, and export reproducible outputs.
# Example command: Rscript scripts/07_cell_communication/00_run_multinichenet.R

source(file.path(getwd(), "src", "utils", "config.R"))
config <- load_config()
raw_data_dir <- resolve_config_path(config, "paths.raw_data_dir")
processed_data_dir <- resolve_config_path(config, "paths.processed_data_dir")
seurat_object_path <- resolve_config_path(config, "paths.seurat_object_path")
rna_expression_parquet <- resolve_config_path(config, "paths.rna_expression_parquet")
atac_peak_parquet <- resolve_config_path(config, "paths.atac_peak_parquet")
result_dir <- ensure_dir(resolve_config_path(config, "paths.result_dir"))
figure_dir <- ensure_dir(resolve_config_path(config, "paths.figure_dir"))
supplementary_table_dir <- ensure_dir(resolve_config_path(config, "paths.supplementary_table_dir"))
model_result_dir <- ensure_dir(resolve_config_path(config, "paths.model_result_dir"))
rna_model_result_dir <- ensure_dir(resolve_config_path(config, "paths.rna_model_result_dir"))
atac_model_result_dir <- ensure_dir(resolve_config_path(config, "paths.atac_model_result_dir"))
de_benchmark_result_dir <- ensure_dir(resolve_config_path(config, "paths.de_benchmark_result_dir"))
enrichment_result_dir <- ensure_dir(resolve_config_path(config, "paths.enrichment_result_dir"))
go_result_dir <- ensure_dir(resolve_config_path(config, "paths.go_result_dir"))
great_result_dir <- ensure_dir(resolve_config_path(config, "paths.great_result_dir"))
stage_gene_dir <- ensure_dir(resolve_config_path(config, "paths.stage_gene_dir"))
stage_peak_dir <- ensure_dir(resolve_config_path(config, "paths.stage_peak_dir"))
downstream_result_dir <- ensure_dir(resolve_config_path(config, "paths.downstream_result_dir"))
atac_regulatory_result_dir <- ensure_dir(resolve_config_path(config, "paths.atac_regulatory_result_dir"))
communication_result_dir <- ensure_dir(resolve_config_path(config, "paths.communication_result_dir"))
pseudotime_result_dir <- ensure_dir(resolve_config_path(config, "paths.pseudotime_result_dir"))

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(SingleCellExperiment)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(multinichenetr)
  library(nichenetr)
  library(muscat)
})

# =========================
# 0)  / 
# =========================
qs_path <- seurat_object_path
ml_path <- "results/stage_features/stage_genes/Acinar_cell__CON_vs_PRE__test_shap_geq0p6.csv"

outdir <- "results/downstream/acinar_CON_PRE/multinichenet"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# receiver（Seurat meta  celltype ）
# acinar  celltype /，，
receiver_raw_candidates <- c(
  "Acinar cell", "acinar cell", "Acinar_cell", "acinar_cell",
  "Acinar", "acinar"
)

# 
groups_keep <- c("CON", "PRE")

# MultiNicheNet  contrast ：，
# （CON-PRE  PRE-CON）
contrasts_oi <- c("'CON-PRE','PRE-CON'")

# （DE）（ DE ）
logFC_threshold <- 0.25
p_val_threshold <- 0.05
p_val_adj <- FALSE

# 
fraction_cutoff <- 0.05
min_cells <- 10
min_sample_prop <- 0.50  #  celltype ： >=50% 
seed <- 1
set.seed(seed)
covariates = NA
batches = NA
# ligand activity 
ligand_activity_down <- TRUE      #  ligand activity
top_n_ligands <- 250
top_n_target <- 250

#  ML ， mean_abs_SHAP ， geneset 
#（MultiNicheNet ：DE geneset  1/200  1/10 ）
trim_ml_by_abs_shap <- FALSE
top_n_ml <- 300

# NicheNet v2 （ vignette  Zenodo  RDS ）
# ， RDS 
#lr_network_url <- "https://zenodo.org/record/7074291/files/lr_network_human_allInfo_30112033.rds"
#ligand_target_matrix_url <- "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_human_allInfo_30112033.rds"
# ：
lr_network_path <- "data/external/lr_network_human_allInfo_30112033.rds"  # 
ligand_target_matrix_path <- "data/external/ligand_target_matrix_nsga2r_final.rds"  # 

# =========================
# 1)  Seurat； CON/PRE； SCE
# =========================
cat("[1] Loading Seurat object...\n")
obj <- qread(qs_path)
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)

stopifnot(all(c("orig.ident", "group", "celltype") %in% colnames(obj@meta.data)))

obj <- subset(obj, subset = group %in% groups_keep)
obj$group <- factor(obj$group, levels = groups_keep)

#  "sample" （ vignette  sample_id="sample"）
obj$sample <- obj$orig.ident

#  SCE
sce <- as.SingleCellExperiment(obj, assay = "RNA")

# 
colData(sce)$group <- obj$group
colData(sce)$sample <- obj$sample
colData(sce)$celltype <- obj$celltype

# ：//，
colData(sce)$group <- make.names(as.character(colData(sce)$group))
colData(sce)$celltype <- make.names(as.character(colData(sce)$celltype))
colData(sce)$sample <- make.names(as.character(colData(sce)$sample))

# （ ML  SCE ）
sce <- makenames_SCE(sce)

group_id <- "group"
sample_id <- "sample"
celltype_id <- "celltype"

receiver_oi <- make.names(receiver_raw_candidates[1])

# =========================
# 2)  sender  receiver
# =========================
all_celltypes <- colData(sce) %>% as.data.frame() %>% pull(.data[[celltype_id]]) %>% unique() %>% sort()

# receiver （）
receiver_hits <- intersect(all_celltypes, make.names(receiver_raw_candidates))
if (length(receiver_hits) == 0) {
  hint <- sort(unique(all_celltypes[grepl("acinar", all_celltypes, ignore.case = TRUE)]))
  stop(" acinar receiver celltype； celltype ： ",
       if (length(hint) == 0) paste(utils::head(all_celltypes, 10), collapse = ", ") else paste(hint, collapse = ", "))
}
receivers_oi <- receiver_hits[1]
receiver_oi <- receivers_oi

# senders： receiver （/， all_celltypes）
senders_oi <- all_celltypes

cat("[2] Receivers:", paste(receivers_oi, collapse=", "), "\n")
cat("[2] Senders:", paste(senders_oi, collapse=", "), "\n")

# Contrast ：group  contrast “”
contrast_tbl <- tibble(
  contrast = c("CON-PRE", "PRE-CON"),
  group    = c(make.names("CON"), make.names("PRE"))
)

# =========================
# 3)  
# =========================
cat("[3] Abundance info...\n")
abundance_info <- get_abundance_info(
  sce = sce,
  sample_id = sample_id,
  group_id = group_id,
  celltype_id = celltype_id,
  min_cells = min_cells,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi
)




ggsave(file.path(outdir, "abund_plot_sample.png"),
       plot = abundance_info$abund_plot_sample,
       width = 10, height = 6, dpi = 300)

ggsave(file.path(outdir, "abund_plot_group.png"),
       plot = abundance_info$abund_plot_group,
       width = 10, height = 6, dpi = 300)

ggsave(file.path(outdir, "abund_barplot.png"),
       plot = abundance_info$abund_barplot,
       width = 10, height = 6, dpi = 300)


#1)  group × celltype  min_cells
library(dplyr)
library(SummarizedExperiment)

abundance_df_summarized <- abundance_info$abundance_data %>% 
  mutate(keep = as.logical(keep)) %>% 
  group_by(group_id, celltype_id) %>% 
  summarise(samples_present = sum(keep), .groups = "drop")

#2)  condition-specific celltypes  absent celltypes（）
celltypes_absent_one_condition <- abundance_df_summarized %>% 
  filter(samples_present == 0) %>% 
  pull(celltype_id) %>% 
  unique()

celltypes_present_one_condition <- abundance_df_summarized %>% 
  filter(samples_present >= 2) %>% 
  pull(celltype_id) %>% 
  unique()

condition_specific_celltypes <- intersect(
  celltypes_absent_one_condition, 
  celltypes_present_one_condition
)

total_nr_conditions <- SummarizedExperiment::colData(sce)[, group_id] %>% 
  unique() %>% 
  length() 


absent_celltypes <- abundance_df_summarized %>% 
  filter(samples_present < 2) %>% 
  group_by(celltype_id) %>% 
  summarise(n = dplyr::n(), .groups = "drop") %>% 
  filter(n == total_nr_conditions) %>% 
  pull(celltype_id)
  
cat("condition-specific celltypes:\n")
print(condition_specific_celltypes)

cat("absent celltypes:\n")
print(absent_celltypes)

#3)  senders/receivers， SCE（）
analyse_condition_specific_celltypes <- FALSE

if (analyse_condition_specific_celltypes == TRUE) {
  senders_oi   <- setdiff(senders_oi, absent_celltypes)
  receivers_oi <- setdiff(receivers_oi, absent_celltypes)
} else {
  senders_oi   <- setdiff(senders_oi, union(absent_celltypes, condition_specific_celltypes))
  receivers_oi <- setdiff(receivers_oi, union(absent_celltypes, condition_specific_celltypes))
}

# ： receiver  senders 
senders_oi <- unique(c(senders_oi, receivers_oi))

#  SCE  celltypes
sce <- sce[, SummarizedExperiment::colData(sce)[, celltype_id] %in% c(senders_oi, receivers_oi)]

cat("Kept senders:\n"); print(senders_oi)
cat("Kept receivers:\n"); print(receivers_oi)
cat("Remaining celltypes in SCE:\n")
print(sort(unique(SummarizedExperiment::colData(sce)[, celltype_id])))

#4)（），
readr::write_csv(abundance_df_summarized,
                 file.path(outdir, "abundance_df_summarized.csv"))


# =========================
# 4) ：（fraction）
# =========================
cat("[4] Gene expression fractions...\n")
frq_list <- get_frac_exprs(
  sce = sce,
  sample_id = sample_id,
  celltype_id = celltype_id,
  group_id = group_id,
  fraction_cutoff = fraction_cutoff,
  min_cells = min_cells,
  min_sample_prop = min_sample_prop
)

# receiver （ background）
expressed_receiver <- frq_list$expressed_df %>%
  filter(celltype == receiver_oi, expressed == TRUE) %>%
  pull(gene) %>% unique()

cat("[4] #background expressed genes (receiver): ", length(expressed_receiver), "\n")

# =========================
# 5)  NicheNet v2 （）
# =========================

# =========================
# 5)  NicheNet v2 （，）
# =========================
cat("[5] Loading NicheNet v2 networks...\n")
lr_network <- readRDS(lr_network_path)
ligand_target_matrix <- readRDS(ligand_target_matrix_path)


# lr_network ；ligand_target_matrix 

# 1)  gene （ makenames_SCE ）
lr_network_all <- lr_network %>%
  mutate(
    ligand   = make.names(ligand),
    receptor = make.names(receptor)
  )

colnames(ligand_target_matrix) <- make.names(colnames(ligand_target_matrix))
rownames(ligand_target_matrix) <- make.names(rownames(ligand_target_matrix))

# 2)  (ligand, receptor) 
lr_network_use <- lr_network_all %>% distinct(ligand, receptor)

# 3) ： ligand_target_matrix  ligands（），
lr_network <- lr_network_use %>% filter(ligand %in% colnames(ligand_target_matrix))
ligand_target_matrix_use <- ligand_target_matrix[, unique(lr_network$ligand), drop = FALSE]

# 4) 
ligands_all   <- unique(lr_network$ligand)
receptors_all <- unique(lr_network$receptor)

cat("[5] lr_network edges:", nrow(lr_network), "\n")
cat("[5] ligands:", length(ligands_all), " receptors:", length(receptors_all), "\n")
cat("[5] ligand_target_matrix dim:", paste(dim(ligand_target_matrix_use), collapse=" x "), "\n")



# =========================
# 6)  + （）
# =========================
cat("[6] Processing abundance/expression info...\n")

abundance_expression_info <- process_abundance_expression_info(
  sce = sce,
  sample_id = sample_id,
  group_id = group_id,
  celltype_id = celltype_id,
  min_cells = min_cells,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network_use,   # 5 distinct(ligand,receptor)
  frq_list = frq_list,
  abundance_info = abundance_info
)

abundance_expression_info$celltype_info$pb_df %>% head()
abundance_expression_info$celltype_info$pb_df_group %>% head()
abundance_expression_info$sender_receiver_info$pb_df %>% head()
abundance_expression_info$sender_receiver_info$pb_df_group %>% head()

# =========================
# 7)  celltype  pseudobulk DE（）
# =========================
cat("[7] Pseudobulk DE via muscat...\n")
DE_info <- get_DE_info(
  sce = sce,
  sample_id = sample_id,
  group_id = group_id,
  celltype_id = celltype_id,
  contrasts_oi = contrasts_oi,
  min_cells = min_cells,
  expressed_df = frq_list$expressed_df,
  batches = batches, covariates = covariates
)

DE_info$celltype_de$de_output_tidy %>% head()
celltype_de <- DE_info$celltype_de$de_output_tidy


ggsave(
  filename = file.path(outdir, "DE_hist_pvals.png"),
  plot     = DE_info$hist_pvals,
  width    = 18, height = 5, dpi = 500
)




# sender/receiver  DE （）

sender_receiver_de = multinichenetr::combine_sender_receiver_de(
  sender_de = celltype_de,
  receiver_de = celltype_de,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network
)
sender_receiver_de %>% head(20)

# =========================
# 8)  ML genes  geneset_of_interest（）
# =========================
sig_receiver <- DE_info$celltype_de$de_output_tidy %>%
  filter(cluster_id == receiver_oi, contrast == "CON-PRE") %>%
  filter(p_val <= p_val_threshold, abs(logFC) >= logFC_threshold) %>%
  arrange(p_val)

write_csv(sig_receiver, file.path(outdir, "DE_acinar_CON_vs_PRE_sig.csv"))
 

# =========================
# 8) Ligand activity ( wrapper) + ML genes  geneset_oi
#    ： DE  background， ML “” geneset
# =========================
cat("[8] Ligand activity via official get_ligand_activities_targets_DEgenes() with ML-masked DE...\n")

# 8.1  ML genes  receiver 
ml_df <- readr::read_csv(ml_path, show_col_types = FALSE)
stopifnot("gene" %in% colnames(ml_df))

ml_df <- ml_df %>% mutate(gene = make.names(gene))

if (trim_ml_by_abs_shap) {
  stopifnot("mean_abs_SHAP" %in% colnames(ml_df))
  ml_df <- ml_df %>% arrange(desc(mean_abs_SHAP)) %>% slice_head(n = top_n_ml)
}

# receiver （4 expressed_receiver）
ml_genes_all <- intersect(unique(ml_df$gene), expressed_receiver)
cat("[8] ML genes in receiver expressed background:", length(ml_genes_all), "\n")
readr::write_csv(tibble(gene = ml_genes_all),
                 file.path(outdir, "ML_genes_in_receiver_background.csv"))

# 8.2 “”： ML  geneset_oi（）
# ：get_ligand_activities_targets_DEgenes()  logFC_threshold + p_val_threshold  receiver_de  geneset
#  ML  p_val/p_adj 1，logFC0 -> 
celltype_de_mlmasked <- celltype_de %>%
  dplyr::mutate(
    is_ml = .data$gene %in% ml_genes_all,
    p_val = ifelse(.data$is_ml, .data$p_val, 1),
    p_adj = ifelse(.data$is_ml, .data$p_adj, 1),
    logFC = ifelse(.data$is_ml, .data$logFC, 0)
  ) %>%
  dplyr::select(-is_ml)

# 8.3 ： geneset/background （）
geneset_assessment_ml <- contrast_tbl$contrast %>%
  lapply(
    multinichenetr::process_geneset_data,
    receiver_de = celltype_de_mlmasked,
    logFC_threshold = logFC_threshold,
    p_val_adj = p_val_adj,
    p_val_threshold = p_val_threshold
  ) %>% bind_rows()
readr::write_csv(geneset_assessment_ml, file.path(outdir, "geneset_assessment_MLmasked.csv"))
print(geneset_assessment_ml)

# 8.4  wrapper： ligand_activities + targets（ activity_scaled / direction）
verbose <- TRUE
cores_system <- 8
n.cores <- min(cores_system, length(unique(celltype_de_mlmasked$cluster_id)))

ligand_activities_targets_DEgenes <- suppressMessages(suppressWarnings(
  multinichenetr::get_ligand_activities_targets_DEgenes(
    receiver_de = celltype_de_mlmasked,
    receivers_oi = intersect(receivers_oi, celltype_de$cluster_id %>% unique()),
    ligand_target_matrix = ligand_target_matrix_use,  # 5
    logFC_threshold = logFC_threshold,
    p_val_threshold = p_val_threshold,
    p_val_adj = p_val_adj,
    top_n_target = top_n_target,
    verbose = verbose,
    n.cores = n.cores
  )
))
ligand_activities_targets_DEgenes$ligand_activities %>% head(20)
# 8.5 
ligand_activities_tbl <- ligand_activities_targets_DEgenes$ligand_activities
readr::write_csv(ligand_activities_tbl,
                 file.path(outdir, "ligand_activities_targets_MLmasked.csv"))

cat("[8] ligand_activities rows:", nrow(ligand_activities_tbl), "\n")
print(ligand_activities_tbl %>% dplyr::arrange(desc(activity_scaled)) %>% head(20))



# =========================
# 9) Prioritization:  sender–ligand–receiver–receptor
# =========================
cat("[9] Prioritization via generate_prioritization_tables()...\n")

# “ targets  ligand activity”（）， FALSE
# “ targets  ligand activity”， TRUE（）
ligand_activity_down <- FALSE

# 1) sender_receiver_tbl： sender-receiver 
sender_receiver_tbl <- sender_receiver_de %>%
  dplyr::distinct(sender, receiver)

# 2) grouping_tbl：（ sample  group ）
metadata_combined <- SummarizedExperiment::colData(sce) %>%
  tibble::as_tibble()

if(!is.na(batches)[1]) {
  grouping_tbl <- metadata_combined %>%
    dplyr::select(all_of(c(sample_id, group_id, batches))) %>%
    dplyr::distinct()
  colnames(grouping_tbl) <- c("sample","group",batches)
} else {
  grouping_tbl <- metadata_combined %>%
    dplyr::select(all_of(c(sample_id, group_id))) %>%
    dplyr::distinct()
  colnames(grouping_tbl) <- c("sample","group")
}

# 3)  prioritization
prioritization_tables <- suppressMessages(
  multinichenetr::generate_prioritization_tables(
    sender_receiver_info = abundance_expression_info$sender_receiver_info,
    sender_receiver_de = sender_receiver_de,
    ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
    contrast_tbl = contrast_tbl,
    sender_receiver_tbl = sender_receiver_tbl,
    grouping_tbl = grouping_tbl,
    scenario = "regular",               # （）
    fraction_cutoff = fraction_cutoff,  # 
    abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
    abundance_data_sender   = abundance_expression_info$abundance_data_sender,
    ligand_activity_down = ligand_activity_down
  )
)

# 4) 
names(prioritization_tables)
prioritization_tables$group_prioritization_tbl %>% head(20)
# 5)  priority_table（）
#    “”： data.frame/tibble  csv
for(nm in names(prioritization_tables)) {
  obj <- prioritization_tables[[nm]]
  if(is.data.frame(obj)) {
    readr::write_csv(obj, file.path(outdir, paste0("prioritization_", nm, ".csv")))
  }
}

# 6) ， 20 （，）
if("prioritization_tbl" %in% names(prioritization_tables)) {
  prioritization_tables$prioritization_tbl %>% dplyr::arrange(dplyr::desc(prioritization_score)) %>% head(20)
}
if("prioritization_table" %in% names(prioritization_tables)) {
  prioritization_tables$prioritization_table %>% dplyr::arrange(dplyr::desc(prioritization_score)) %>% head(20)
}
if("prioritization_scores" %in% names(prioritization_tables)) {
  prioritization_tables$prioritization_scores %>% head(20)
}

cat("[9] Prioritization done. CSVs saved to: ", outdir, "\n", sep="")




# =========================
# 10) Optional: LR-target correlation across samples
# =========================
cat("[10] LR-target correlation inference...\n")

# receivers_oi： receiver（ receiver）
receivers_cor <- prioritization_tables$group_prioritization_tbl$receiver %>% unique()

lr_target_prior_cor <- multinichenetr::lr_target_prior_cor_inference(
  receivers_oi = receivers_cor,
  abundance_expression_info = abundance_expression_info,
  celltype_de = celltype_de,
  grouping_tbl = grouping_tbl,
  prioritization_tables = prioritization_tables,
  ligand_target_matrix = ligand_target_matrix_use,  # 5
  logFC_threshold = logFC_threshold,
  p_val_threshold = p_val_threshold,
  p_val_adj = p_val_adj
)

# 
names(lr_target_prior_cor)

# =========================
# 11) Save all key MultiNicheNet outputs (lite) as .qs
# =========================
cat("[11] Building and saving multinichenet_output...\n")

multinichenet_output <- list(
  celltype_info = abundance_expression_info$celltype_info,
  celltype_de = celltype_de,
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  prioritization_tables = prioritization_tables,
  grouping_tbl = grouping_tbl,
  lr_target_prior_cor = lr_target_prior_cor
)

# ：
multinichenet_output_lite <- multinichenetr::make_lite_output(multinichenet_output)

#  qs（）
qs_out <- file.path(outdir, "multinichenet_output_lite.qs")
qs::qsave(multinichenet_output_lite, qs_out, preset = "high")

cat("[11] Saved: ", qs_out, "\n", sep="")

#
library(dplyr)
library(RColorBrewer)
# 1️⃣  top50
prioritized_tbl_oi_all <- get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables,
  top_n = 50,
  rank_per_group = FALSE
)

# 2️⃣ join 
prioritized_tbl_oi <- multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
  dplyr::filter(.data$id %in% prioritized_tbl_oi_all$id) %>%
  dplyr::distinct(.data$id, .data$sender, .data$receiver, .data$ligand, .data$receptor, .data$group) %>%
  dplyr::left_join(
    prioritized_tbl_oi_all %>%
      dplyr::select(.data$id, .data$group, .data$prioritization_score, .data$prioritization_rank),
    by = c("id", "group")
  )

# 3️⃣ NA → 0
prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] <- 0

# 4️⃣  receiver（α cell -> receiver_oi）
prioritized_tbl_oi_receiver <- prioritized_tbl_oi %>%
  filter(receiver == receiver_oi)

# 
stopifnot(nrow(prioritized_tbl_oi_receiver) > 0)
prioritized_tbl_oi_receiver <- prioritized_tbl_oi_receiver %>%
  mutate(
    sender = recode(sender,
      "α.cell" = "Alpha.cell",
      "β.cell" = "Beta.cell",
      "δ.cell" = "Delta.cell"
    ),
    receiver = recode(receiver,
      "α.cell" = "Alpha.cell",
      "β.cell" = "Beta.cell",
      "δ.cell" = "Delta.cell"
    )
  )
senders_receivers <- union(
  prioritized_tbl_oi_receiver$sender,
  prioritized_tbl_oi_receiver$receiver
) %>% sort()
n_sr <- length(senders_receivers)
base_pal <- RColorBrewer::brewer.pal(min(11, max(3, n_sr)), "Spectral")
cols <- setNames(colorRampPalette(base_pal)(n_sr), senders_receivers)

pdf_file <- file.path(outdir, paste0("circos_top50_global_", receiver_oi, ".pdf"))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

grDevices::pdf(pdf_file, width = 9, height = 9, useDingbats = FALSE)
multinichenetr::make_circos_group_comparison(
  prioritized_tbl_oi_receiver,
  colors_sender = cols,
  colors_receiver = cols
)
grDevices::dev.off()





# n_sr <- length(senders_receivers)
# base_pal <- RColorBrewer::brewer.pal(min(11, max(3, n_sr)), "Spectral")
# cols <- setNames(colorRampPalette(base_pal)(n_sr), senders_receivers)

#   pdf_file <- file.path(outdir, paste0("circos_top50_global_", receiver_oi, ".pdf"))
#   dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

#   grDevices::pdf(pdf_file, width = 9, height = 9, useDingbats = FALSE)
#   multinichenetr::make_circos_group_comparison(
#     prioritized_tbl_oi_global,
#     colors_sender = cols,
#     colors_receiver = cols
#   )
#   grDevices::dev.off()




# =========================
# 1)  group  Top N
# =========================
library(dplyr)
library(RColorBrewer)
library(multinichenetr)

# 1)  top50（ receiver ）

prioritized_tbl_oi_all <- get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables,
  top_n = 50,
  rank_per_group = TRUE
)
print(table(prioritized_tbl_oi_all$group, useNA = "ifany"))

prioritized_tbl_oi = 
  multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% prioritized_tbl_oi_all$id) %>%
  distinct(id, sender, receiver, ligand, receptor, group) %>% 
  left_join(prioritized_tbl_oi_all)
prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0

senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), prioritized_tbl_oi$receiver %>% unique()) %>% sort()

colors_sender = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
colors_receiver = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
pdf_file <- file.path(outdir, "circos_grouptop50.pdf")

grDevices::pdf(pdf_file, width = 9, height = 9, useDingbats = FALSE)

circos_list <- multinichenetr::make_circos_group_comparison(
  prioritized_tbl_oi,
  colors_sender,
  colors_receiver
)

grDevices::dev.off()







groups_use <- c("CON", "PRE")

for (group_oi in groups_use) {

  # 1)  group  top50
  prioritized_tbl_group <- get_top_n_lr_pairs(
    multinichenet_output$prioritization_tables,
    top_n = 50,
    groups_oi = group_oi
  )

  # 2) ， group  sender/receiver/ligand/receptor 
  prioritized_tbl_group_full <-
    multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
    dplyr::filter(group == group_oi, id %in% prioritized_tbl_group$id) %>%
    dplyr::distinct(id, sender, receiver, ligand, receptor, group) %>%
    dplyr::left_join(
      prioritized_tbl_group %>%
        dplyr::select(id, group, prioritization_score, prioritization_rank),
      by = c("id", "group")
    )

  # 2.1) ， legend  ...
  prioritized_tbl_group_full <- prioritized_tbl_group_full %>%
    dplyr::mutate(
      sender = dplyr::recode(
        sender,
        "α.cell" = "Alpha.cell",
        "β.cell" = "Beta.cell",
        "δ.cell" = "Delta.cell"
      ),
      receiver = dplyr::recode(
        receiver,
        "α.cell" = "Alpha.cell",
        "β.cell" = "Beta.cell",
        "δ.cell" = "Delta.cell"
      )
    )

  # 3)  group  sender / receiver
  senders_receivers <- union(
    unique(prioritized_tbl_group_full$sender),
    unique(prioritized_tbl_group_full$receiver)
  ) %>% sort()

  n_sr <- length(senders_receivers)
  base_pal <- RColorBrewer::brewer.pal(min(11, max(3, n_sr)), "Spectral")
  cols <- setNames(colorRampPalette(base_pal)(n_sr), senders_receivers)

  # 4) 
  pdf_file <- file.path(outdir, paste0("circos_top50_", group_oi, ".pdf"))
  grDevices::pdf(pdf_file, width = 9, height = 9, useDingbats = FALSE)

  multinichenetr::make_circos_group_comparison(
    prioritized_tbl_group_full,
    colors_sender = cols,
    colors_receiver = cols
  )

  grDevices::dev.off()
}


# prioritized_tbl_oi_all <- get_top_n_lr_pairs(
#   multinichenet_output$prioritization_tables,
#   top_n = 50,
#   rank_per_group = TRUE
# )
# print(table(prioritized_tbl_oi_all$group, useNA = "ifany"))
# # 2)  id + group join （）
# prioritized_tbl_oi <-
#   multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
#   filter(id %in% prioritized_tbl_oi_all$id) %>%
#   distinct(id, sender, receiver, ligand, receptor, group) %>%
#   inner_join(   #  inner_join，； left_join 
#     prioritized_tbl_oi_all %>% dplyr::select(id, group, prioritization_score, prioritization_rank),
#     by = c("id", "group")
#   )

# # 3)  receiver（： 50 ）
# prioritized_tbl_oi_receiver <- prioritized_tbl_oi %>%
#   dplyr::filter(.data$receiver == .env$receiver_oi)
# print(table(prioritized_tbl_oi_receiver$group, useNA = "ifany"))
# stopifnot(nrow(prioritized_tbl_oi_receiver) > 0)

# # 4) 
# senders_receivers <- sort(unique(c(prioritized_tbl_oi_receiver$sender, prioritized_tbl_oi_receiver$receiver)))
# n_sr <- length(senders_receivers)
# cols <- setNames(colorRampPalette(brewer.pal(min(11, max(3, n_sr)), "Spectral"))(n_sr), senders_receivers)

# # 5) 
# pdf_file <- file.path(outdir, "circos_grouptop50.pdf")
# grDevices::pdf(pdf_file, width = 9, height = 9, useDingbats = FALSE)
# make_circos_group_comparison(
#   prioritized_tbl_oi,
#   colors_sender = cols,
#   colors_receiver = cols
# )
# grDevices::dev.off()
# file.info(pdf_file)[, c("size", "mtime")]


# #legend
# pdf_file <- file.path(outdir, "circos_top50_unicode.pdf")
# dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# if (capabilities("cairo")) {
#   grDevices::cairo_pdf(pdf_file, width = 9, height = 9, family = "DejaVu Sans")
# } else {
#   grDevices::pdf(pdf_file, width = 9, height = 9, useDingbats = FALSE)
# }

# ： circos_list  replayPlot
# ： device 
circos_list <- multinichenetr::make_circos_group_comparison(
  prioritized_tbl_oi_receiver,
  colors_sender = cols,
  colors_receiver = cols
)

grDevices::dev.off()
file.info(pdf_file)[, c("size", "mtime")]




##########################
#####bubble plot #########
##########################
#top50
library(dplyr)
library(multinichenetr)

group_oi <- "PRE"

# 1)  overall top50
prioritized_tbl_overall50 <- multinichenetr::get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables,
  top_n = 50,
  rank_per_group = FALSE
)

print(table(prioritized_tbl_overall50$group, useNA = "ifany"))

# 2) join （sender/receiver/ligand/receptor）
prioritized_tbl_overall50_full <-
  multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% prioritized_tbl_overall50$id) %>%
  distinct(id, sender, receiver, ligand, receptor, group) %>%
  left_join(
    prioritized_tbl_overall50 %>% dplyr::select(id, group, prioritization_score, prioritization_rank),
    by = c("id", "group")
  ) %>%
  mutate(prioritization_score = ifelse(is.na(prioritization_score), 0, prioritization_score))

# 3)  group + receiver
prioritized_tbl_overall50_receiver <- prioritized_tbl_overall50_full %>%
  filter(group == group_oi, receiver == receiver_oi)

cat("n(overall50 ∩", group_oi, "∩ receiver=", receiver_oi, ") = ", nrow(prioritized_tbl_overall50_receiver), "\n")

# 4) bubble plot
plot_overall50 <- multinichenetr::make_sample_lr_prod_activity_plots_Omnipath(
  multinichenet_output$prioritization_tables,
  prioritized_tbl_overall50_receiver %>%
    inner_join(lr_network_all, by = c("ligand", "receptor"))
)

ggplot2::ggsave(
  filename = file.path(outdir, paste0("bubble_overallTop50_in_", group_oi, "_", receiver_oi, "_receiver.png")),
  plot = plot_overall50,
  width = 16, height = 12, dpi = 500
)
ggplot2::ggsave(
  filename = file.path(outdir, paste0(
    "bubble_overallTop50_in_", group_oi, "_", receiver_oi, "_receiver.pdf"
  )),
  plot = plot_overall50,
  width = 16, height = 12,
  device = cairo_pdf   # ：/
)
##########################
#top50
library(dplyr)
library(multinichenetr)

group_oi <- "PRE"

# 1)  overall top50
prioritized_tbl_overall50 <- multinichenetr::get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables,
  top_n = 50,
  rank_per_group = FALSE
)

print(table(prioritized_tbl_overall50$group, useNA = "ifany"))

# 2) join （sender/receiver/ligand/receptor）
prioritized_tbl_overall50_full <-
  multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% prioritized_tbl_overall50$id) %>%
  distinct(id, sender, receiver, ligand, receptor, group) %>%
  left_join(
    prioritized_tbl_overall50 %>% dplyr::select(id, group, prioritization_score, prioritization_rank),
    by = c("id", "group")
  ) %>%
  mutate(prioritization_score = ifelse(is.na(prioritization_score), 0, prioritization_score))

# 3)  group + receiver
prioritized_tbl_overall50_receiver <- prioritized_tbl_overall50_full %>%
  filter(group == group_oi, receiver == receiver_oi)

cat("n(overall50 ∩", group_oi, "∩ receiver=", receiver_oi, ") = ", nrow(prioritized_tbl_overall50_receiver), "\n")

# 4) bubble plot
plot_overall50 <- multinichenetr::make_sample_lr_prod_activity_plots_Omnipath(
  multinichenet_output$prioritization_tables,
  prioritized_tbl_overall50_receiver %>%
    inner_join(lr_network_all, by = c("ligand", "receptor"))
)

ggplot2::ggsave(
  filename = file.path(outdir, paste0("bubble_overallTop50_in_", group_oi, "_", receiver_oi, "_receiver.png")),
  plot = plot_overall50,
  width = 16, height = 12, dpi = 500
)
ggplot2::ggsave(
  filename = file.path(outdir, paste0(
    "bubble_overallTop50_in_", group_oi, "_", receiver_oi, "_receiver.pdf"
  )),
  plot = plot_overall50,
  width = 16, height = 12,
  device = cairo_pdf   # ：/
)

#######################
#top50


group_oi <- "CON"

prioritized_tbl_group_top50 <- multinichenetr::get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables,
  50,
  groups_oi = group_oi
)

prioritized_tbl_group_top50_full <-
  multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
  dplyr::filter(id %in% prioritized_tbl_group_top50$id, group == group_oi) %>%
  dplyr::distinct(id, sender, receiver, ligand, receptor, group) %>%
  dplyr::left_join(
    prioritized_tbl_group_top50 %>%
      dplyr::select(id, group, prioritization_score, prioritization_rank),
    by = c("id", "group")
  ) %>%
  dplyr::filter(receiver == receiver_oi)
plot_group_top50 <- multinichenetr::make_sample_lr_prod_activity_plots_Omnipath(
  multinichenet_output$prioritization_tables,
  prioritized_tbl_group_top50_full %>% inner_join(lr_network_all, by = c("ligand","receptor"))
)

ggplot2::ggsave(
  file.path(outdir, paste0("bubble_groupTop50_", group_oi, "_", receiver_oi, "_receiver.png")),
  plot_group_top50,
  width = 16, height = 12, dpi = 500
)
ggplot2::ggsave(
  file.path(outdir, paste0(
    "bubble_groupTop50_", group_oi, "_", receiver_oi, "_receiver.pdf"
  )),
  plot_group_top50,
  width = 16, height = 12,
  device = cairo_pdf   # ：/
)

group_oi <- "PRE"

prioritized_tbl_group_top50 <- multinichenetr::get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables,
  50,
  groups_oi = group_oi
)
prioritized_tbl_group_top50_full <-
  multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% prioritized_tbl_group_top50$id, group == group_oi) %>%
  distinct(id, sender, receiver, ligand, receptor, group) %>%
  left_join(
    prioritized_tbl_group_top50 %>% dplyr::select(id, group, prioritization_score, prioritization_rank),
    by = c("id","group")
  ) %>%
  filter(receiver == receiver_oi)  # 

plot_group_top50 <- multinichenetr::make_sample_lr_prod_activity_plots_Omnipath(
  multinichenet_output$prioritization_tables,
  prioritized_tbl_group_top50_full %>% inner_join(lr_network_all, by = c("ligand","receptor"))
)

ggplot2::ggsave(
  file.path(outdir, paste0("bubble_groupTop50_", group_oi, "_", receiver_oi, "_receiver.png")),
  plot_group_top50,
  width = 16, height = 12, dpi = 500
)
ggplot2::ggsave(
  file.path(outdir, paste0(
    "bubble_groupTop50_", group_oi, "_", receiver_oi, "_receiver.pdf"
  )),
  plot_group_top50,
  width = 16, height = 12,
  device = cairo_pdf   # ：/
)



###
#0）： prioritized LR （： top50， receiver=receiver_oi）


receiver_focus <- receiver_oi
top_n_use <- 50

#  top50（CON 50 + PRE 50）
prioritized_tbl_oi_all <- multinichenetr::get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables,
  top_n = top_n_use,
  rank_per_group = TRUE,
  receivers_oi = receiver_focus
)

# 
print(table(prioritized_tbl_oi_all$group, useNA = "ifany"))
stopifnot(nrow(prioritized_tbl_oi_all) > 0)

# 1： NicheNet  target（）

lr_target_prior = prioritized_tbl_oi_all %>% inner_join(
        multinichenet_output$ligand_activities_targets_DEgenes$ligand_activities %>%
          distinct(ligand, target, direction_regulation, contrast) %>% inner_join(contrast_tbl) %>% ungroup() 
        ) 
lr_target_df = lr_target_prior %>% distinct(group, sender, receiver, ligand, receptor, id, target, direction_regulation) 

lr_target_df %>% filter(target %in% union(lr_network$ligand, lr_network$receptor))

network_prior = infer_intercellular_regulatory_network(lr_target_df, prioritized_tbl_oi_all)
network_prior$links %>% head()
network_prior$nodes %>% head()
# 5) （，）
readr::write_csv(network_prior$nodes, file.path(outdir, "network_nodes.csv"))
readr::write_csv(network_prior$links, file.path(outdir, "network_links.csv"))



# 2： + LR–target （）
#： lr_target_prior_cor = lr_target_prior_cor_inference(...)  multinichenet_output$lr_target_prior_cor
lr_target_prior_cor_filtered = 
  multinichenet_output$prioritization_tables$group_prioritization_tbl$group %>% unique() %>% 
  lapply(function(group_oi){
    lr_target_prior_cor_filtered = multinichenet_output$lr_target_prior_cor %>%
      inner_join(
        multinichenet_output$ligand_activities_targets_DEgenes$ligand_activities %>%
          distinct(ligand, target, direction_regulation, contrast)
        ) %>% 
      inner_join(contrast_tbl) %>% filter(group == group_oi)
    
    lr_target_prior_cor_filtered_up = lr_target_prior_cor_filtered %>% 
      filter(direction_regulation == "up") %>% 
      filter( (rank_of_target < top_n_target) & (pearson > 0.33))
    
    lr_target_prior_cor_filtered_down = lr_target_prior_cor_filtered %>% 
      filter(direction_regulation == "down") %>% 
      filter( (rank_of_target < top_n_target) & (pearson < -0.33))
    lr_target_prior_cor_filtered = bind_rows(
      lr_target_prior_cor_filtered_up, 
      lr_target_prior_cor_filtered_down
      )
}) %>% bind_rows()

lr_target_df = lr_target_prior_cor_filtered %>% 
  distinct(group, sender, receiver, ligand, receptor, id, target, direction_regulation) 
network_cor = infer_intercellular_regulatory_network(lr_target_df, prioritized_tbl_oi_all)
network_cor$links %>% head()
network_cor$nodes %>% head()
# 5) （，）
readr::write_csv(network_cor$nodes, file.path(outdir, "network_cor_nodes.csv"))
readr::write_csv(network_cor$links, file.path(outdir, "network_cor_links.csv"))



colors_sender[receiver_oi] = "pink" # the  original yellow background with white font is not very readable
network_graph = visualize_network(network_prior, colors_sender)
network_graph$plot





suppressPackageStartupMessages({
  library(dplyr)
  library(RColorBrewer)
  library(multinichenetr)
})

sanitize_celltype <- function(x){
  x %>%
    gsub("β\\.cell", "beta.cell",  .) %>%
    gsub("α\\.cell", "alpha.cell", .) %>%
    gsub("δ\\.cell", "delta.cell", .)
}

make_celltype_palette <- function(celltypes){
  celltypes <- sort(unique(celltypes))
  n <- length(celltypes)

  # ； celltype ，
  base_cols <- c(
    "#6E8FB2", "#7DA494", "#E5A79A", "#C16E71",
    "#ABC8E5", "#D8A0C1", "#9F8DB8", "#D0D08A"
  )

  pal <- if (n <= length(base_cols)) base_cols[seq_len(n)] else grDevices::colorRampPalette(base_cols)(n)
  setNames(pal, celltypes)
}

save_plot_any <- function(p, pdf_file, png_file, w = 12, h = 10, dpi = 300){
  # PDF： cairo_pdf（）； pdf()
  if ("cairo_pdf" %in% ls(getNamespace("grDevices"))) {
    grDevices::cairo_pdf(pdf_file, width = w, height = h)
  } else {
    grDevices::pdf(pdf_file, width = w, height = h, useDingbats = FALSE)
  }
  if (inherits(p, "ggplot")) print(p) else grDevices::replayPlot(p)
  grDevices::dev.off()

  # PNG
  grDevices::png(png_file, width = w, height = h, units = "in", res = dpi)
  if (inherits(p, "ggplot")) print(p) else grDevices::replayPlot(p)
  grDevices::dev.off()
}

plot_network_by_group <- function(network_obj, outdir, prefix){
  stopifnot(all(c("nodes","links") %in% names(network_obj)))

  # ---- ：nodes  node （ row.names ）----
  # （）
  net_nodes <- network_obj$nodes %>%
  dplyr::mutate(
    node = sanitize_celltype(node),
    celltype = sanitize_celltype(celltype)
  ) %>%
  dplyr::group_by(node) %>%
  dplyr::summarise(
    celltype  = dplyr::first(.data$celltype),
    gene      = dplyr::first(.data$gene),
    type_gene = paste(sort(unique(.data$type_gene)), collapse = "|"),
    .groups = "drop"
  )

  # links （ nodes$node ）
  net_links <- network_obj$links %>%
    mutate(
      sender_ligand   = sanitize_celltype(sender_ligand),
      receiver_target = sanitize_celltype(receiver_target)
    )

  groups <- sort(unique(net_links$group))
  if (length(groups) == 0) {
    message("[WARN] ", prefix, ": no groups found in links.")
    return(invisible(NULL))
  }

  # （ group ）
  cols_sender <- make_celltype_palette(net_nodes$celltype)

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  for (g in groups) {
    links_g <- net_links %>% filter(group == g)
    if (nrow(links_g) == 0) {
      message("[WARN] ", prefix, " group=", g, ": 0 links, skip plotting.")
      next
    }

    network_g <- list(nodes = net_nodes, links = links_g)

    graph_g <- multinichenetr::visualize_network(network_g, cols_sender)
    p <- graph_g$plot

    pdf_file <- file.path(outdir, paste0(prefix, "_", g, ".pdf"))
    png_file <- file.path(outdir, paste0(prefix, "_", g, ".png"))

    save_plot_any(p, pdf_file, png_file, w =12.5, h = 8, dpi = 500)
    message("[OK] saved: ", pdf_file)
    message("[OK] saved: ", png_file)
  }

  invisible(NULL)
}
plot_network_by_group(network_prior, outdir, prefix = "intercellular_network_prior")
plot_network_by_group(network_cor,   outdir, prefix = "intercellular_network_cor")



#  LR （A→B）
#  target  LR  ligand/receptor ，
# （A  B，B  A）
#  LR “”
network <- network_prior
network$prioritized_lr_interactions
prioritized_tbl_oi_network = prioritized_tbl_oi_all %>% inner_join(
  network$prioritized_lr_interactions)
prioritized_tbl_oi_network
group_oi = "CON"
prioritized_tbl_oi_CON = prioritized_tbl_oi_network %>% filter(group == group_oi)

plot_oi = make_sample_lr_prod_activity_plots_Omnipath(
  multinichenet_output$prioritization_tables, 
  prioritized_tbl_oi_CON %>% inner_join(lr_network_all)
  )
ggplot2::ggsave(
  filename = file.path(outdir, paste0("bubble_networkTop_in_", group_oi, "_", receiver_oi, "_receiver.png")),
  plot = plot_oi,
  width = 16, height = 12, dpi = 500
)
ggplot2::ggsave(
  filename = file.path(outdir, paste0(
    "bubble_networkTop_in_", group_oi, "_", receiver_oi, "_receiver.pdf"
  )),
  plot = plot_oi,
  width = 16, height = 12,
  device = cairo_pdf   # ：/
)
group_oi = "PRE"
prioritized_tbl_oi_PRE = prioritized_tbl_oi_network %>% filter(group == group_oi)

plot_oi = make_sample_lr_prod_activity_plots_Omnipath(
  multinichenet_output$prioritization_tables, 
  prioritized_tbl_oi_PRE %>% inner_join(lr_network_all)
  )
ggplot2::ggsave(
  filename = file.path(outdir, paste0("bubble_networkTop_in_", group_oi, "_", receiver_oi, "_receiver.png")),
  plot = plot_oi,
  width = 16, height = 12, dpi = 500
)
ggplot2::ggsave(
  filename = file.path(outdir, paste0(
    "bubble_networkTop_in_", group_oi, "_", receiver_oi, "_receiver.pdf"
  )),
  plot = plot_oi,
  width = 16, height = 12,
  device = cairo_pdf   # ：/
)






# prioritized LR（sender–ligand→receiver–receptor） “”：
#  ligand activity ，
#   target genes  activity，
#  target “”
group_oi    <- "CON"      #  "PRE"
receiver_oi <- receiver_focus   # （make.names ）

prioritized_tbl_oi_group_50 <- multinichenetr::get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables,
  top_n      = 50,
  groups_oi  = group_oi,
  receivers_oi = receiver_oi
)

# 
print(prioritized_tbl_oi_group_50 %>% dplyr::count(.data$group))

stopifnot(nrow(prioritized_tbl_oi_group_50) > 0)

combined_plot <- multinichenetr::make_ligand_activity_target_plot(
  group_oi,
  receiver_oi,
  prioritized_tbl_oi_group_50,
  multinichenet_output$prioritization_tables,
  multinichenet_output$ligand_activities_targets_DEgenes,
  contrast_tbl,
  multinichenet_output$grouping_tbl,
  multinichenet_output$celltype_info,     # <- 8：
  ligand_target_matrix,                   # <- 9
  plot_legend = FALSE
)
ggplot2::ggsave(
  filename = file.path(outdir, paste0("ligand_activity_target_top50_", group_oi, "_", receiver_oi, ".png")),
  plot = combined_plot$combined_plot,
  width = 13, height = 10, dpi = 500
)
ggplot2::ggsave(
  filename = file.path(outdir, paste0(
    "ligand_activity_target_top50_", group_oi, "_", receiver_oi, ".pdf"
  )),
  plot = combined_plot$combined_plot,
  width = 13, height = 10,
  device = cairo_pdf   # ：/
)
group_oi    <- "PRE"      #  "CON"
receiver_oi <- receiver_focus   # （make.names ）

prioritized_tbl_oi_group_50 <- multinichenetr::get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables,
  top_n      = 50,
  groups_oi  = group_oi,
  receivers_oi = receiver_oi
)

# 
print(prioritized_tbl_oi_group_50 %>%  dplyr::count(group))
stopifnot(nrow(prioritized_tbl_oi_group_50) > 0)

combined_plot <- multinichenetr::make_ligand_activity_target_plot(
  group_oi,
  receiver_oi,
  prioritized_tbl_oi_group_50,
  multinichenet_output$prioritization_tables,
  multinichenet_output$ligand_activities_targets_DEgenes,
  contrast_tbl,
  multinichenet_output$grouping_tbl,
  multinichenet_output$celltype_info,     # <- 8：
  ligand_target_matrix,                   # <- 9
  plot_legend = FALSE
)
ggplot2::ggsave(
  filename = file.path(outdir, paste0("ligand_activity_target_top50_", group_oi, "_", receiver_oi, ".png")),
  plot = combined_plot$combined_plot,
  width = 14, height = 10, dpi = 500
)
ggplot2::ggsave(
  filename = file.path(outdir, paste0(
    "ligand_activity_target_top50_", group_oi, "_", receiver_oi, ".pdf"
  )),
  plot = combined_plot$combined_plot,
  width = 14, height = 10,
  device = cairo_pdf   # ：/
)





