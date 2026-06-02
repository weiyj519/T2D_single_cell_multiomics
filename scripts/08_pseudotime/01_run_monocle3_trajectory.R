#!/usr/bin/env Rscript

# Script name: 01_run_monocle3_trajectory.R
# Purpose: Cleaned public workflow script for the T2D single-cell multiome analysis.
# Input: Paths and analysis settings are read from configs/config.yaml.
# Output: Module-specific outputs are written under the configured results directory.
# Main steps: Load configuration, run the original analysis logic, and export reproducible outputs.
# Example command: Rscript scripts/08_pseudotime/01_run_monocle3_trajectory.R

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

# =========================================================
# Monocle3 native trajectory on stage-selected acinar cells
# trajectory genes = ML-union genes
# merged from union.R + monocle3.R (Strategy B only)
# =========================================================

suppressPackageStartupMessages({
  library(qs)
  library(Matrix)
  library(Seurat)
  library(monocle3)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(Signac)
  library(data.table)
})

# =========================================================
# 0. 
# =========================================================
out_dir <- "results/downstream/riskcell_monocle/acinar"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# 1. 
# =========================================================
selected_path <- "results/downstream/riskcell/acinar/acinar_selected_cells_for_pseudotime.tsv"

ml_path_con_pre <- "results/model_outputs/rna/Acinar_cell/CON_vs_PRE_shap_marker_genes.csv"
ml_path_pre_t2d <- "results/model_outputs/rna/Acinar_cell/PRE_vs_T2D_shap_marker_genes.csv"
ml_path_con_t2d <- "results/model_outputs/rna/Acinar_cell/CON_vs_T2D_shap_marker_genes.csv"

seurat_path <- seurat_object_path

#  align_cds（ donor ）
DO_ALIGN <- TRUE

# =========================================================
# 2. helper 
# =========================================================
get_layer_mat <- function(obj, assay = "RNA", layer = c("counts", "data")) {
  layer <- match.arg(layer)
  tryCatch(
    LayerData(obj, assay = assay, layer = layer),
    error = function(e) GetAssayData(obj, assay = assay, slot = layer)
  )
}

save_plot <- function(p, filename, width = 6, height = 5, dpi = 300) {
  ggsave(filename, p, width = width, height = height, dpi = dpi)
}

rasterize_points <- function(p, dpi = 300) {
  if (!requireNamespace("ggrastr", quietly = TRUE)) {
    message("[WARN] ggrastr not installed, save as vector points.")
    return(p)
  }
  ggrastr::rasterise(p, layers = "Point", dpi = dpi)
}

read_ml_table <- function(path, set_name) {
  if (!file.exists(path)) stop("File not found: ", path)

  df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)

  gene_col <- intersect(c("gene", "Gene", "feature", "Feature"), colnames(df))[1]
  if (is.na(gene_col)) stop("No gene column found in: ", path)

  weight_col <- intersect(
    c("mean_abs_SHAP", "mean_abs_shap", "abs_SHAP", "abs_shap",
      "mean_SHAP", "SHAP", "importance", "gain"),
    colnames(df)
  )[1]

  out <- data.frame(
    gene = as.character(df[[gene_col]]),
    raw_weight = if (is.na(weight_col)) 1 else suppressWarnings(as.numeric(df[[weight_col]])),
    stringsAsFactors = FALSE
  )

  out <- out %>%
    filter(!is.na(gene), gene != "") %>%
    group_by(gene) %>%
    summarise(raw_weight = max(raw_weight, na.rm = TRUE), .groups = "drop")

  bad_weight <- !is.finite(out$raw_weight)
  out$raw_weight[bad_weight] <- 1

  if (nrow(out) == 0) stop("No valid genes found in: ", path)

  if (max(out$raw_weight, na.rm = TRUE) == 0) {
    out$raw_weight <- 1
  }

  out$weight <- out$raw_weight / max(out$raw_weight, na.rm = TRUE)
  out$gene_set <- set_name
  out
}

choose_root_nodes_by_group <- function(cds, group_col = "group", early_group = "CON") {
  cv <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  cv <- as.matrix(cv[colnames(cds), , drop = FALSE])
  if (ncol(cv) > 1) cv <- cv[, 1, drop = FALSE]

  part_vec <- as.character(partitions(cds))
  grp_vec  <- as.character(colData(cds)[, group_col])
  pg_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name

  root_nodes <- c()

  for (pt in unique(part_vec)) {
    cells_part <- colnames(cds)[part_vec == pt]
    grp_part   <- grp_vec[part_vec == pt]

    early_cells <- cells_part[grp_part == early_group]
    cells_use <- if (length(early_cells) > 0) early_cells else cells_part

    tab <- table(cv[cells_use, 1])
    if (length(tab) == 0) next

    winner <- names(which.max(tab))[1]
    winner_num <- suppressWarnings(as.numeric(winner))
    if (!is.na(winner_num) && winner_num <= length(pg_nodes)) {
      root_nodes <- c(root_nodes, pg_nodes[winner_num])
    } else {
      root_nodes <- c(root_nodes, winner)
    }
  }

  unique(root_nodes)
}

# =========================================================
# 3.  stage-selected acinar cells
# =========================================================
cat("\n========== 1. Loading stage-selected acinar cells ==========\n")

sel <- fread(selected_path)

combined <- qread(seurat_path)
Idents(combined) <- combined$celltype
acinar_obj <- subset(combined, idents = "Acinar cell")
rm(combined); gc()

cells_keep <- intersect(colnames(acinar_obj), sel$Cell_id)
acinar_obj <- subset(acinar_obj, cells = cells_keep)

acinar_obj$Stage_selected <- sel$Stage_selected[match(colnames(acinar_obj), sel$Cell_id)]
acinar_obj$Stage_selected <- factor(acinar_obj$Stage_selected, levels = c("CON", "PRE", "T2D"))

cat("Selected acinar cells:", ncol(acinar_obj), "\n")
cat("\n[Stage_selected distribution]\n")
print(table(acinar_obj$Stage_selected, useNA = "ifany"))
cat("\n[Original group distribution]\n")
print(table(acinar_obj$group, useNA = "ifany"))
cat("\n[Sample distribution]\n")
print(table(acinar_obj$orig.ident, useNA = "ifany"))

DefaultAssay(acinar_obj) <- "RNA"
acinar_obj <- tryCatch({
  JoinLayers(acinar_obj)
}, error = function(e) {
  message("[INFO] JoinLayers skipped: ", conditionMessage(e))
  acinar_obj
})

# =========================================================
# 4.  ML  union
# =========================================================
cat("\n========== 2. Reading ML-union genes ==========\n")

tbl_con_pre <- read_ml_table(ml_path_con_pre, "CON_vs_PRE")
tbl_pre_t2d <- read_ml_table(ml_path_pre_t2d, "PRE_vs_T2D")
tbl_con_t2d <- read_ml_table(ml_path_con_t2d, "CON_vs_T2D")

ml_all <- bind_rows(tbl_con_pre, tbl_pre_t2d, tbl_con_t2d)

ml_union_df <- ml_all %>%
  group_by(gene) %>%
  summarise(
    n_sets = n_distinct(gene_set),
    sets = paste(sort(unique(gene_set)), collapse = ";"),
    max_raw_weight = max(raw_weight, na.rm = TRUE),
    max_weight = max(weight, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sets), desc(max_raw_weight), gene)

ml_union <- unique(ml_union_df$gene)

cat("\n[ML set sizes]\n")
print(table(ml_all$gene_set))
cat("\n[ML union size before expression filtering]:", length(ml_union), "\n")

write.csv(
  ml_union_df,
  file.path(out_dir, "acinar_ml_union_gene_membership.csv"),
  row.names = FALSE
)

# =========================================================
# 5.  ML-union genes 
# =========================================================
cat("\n========== 3. Restricting matrix to ML-union genes ==========\n")

counts_mat_full <- get_layer_mat(acinar_obj, assay = "RNA", layer = "counts")
counts_mat_full <- as(counts_mat_full, "dgCMatrix")

genes_use <- intersect(ml_union, rownames(counts_mat_full))
cat("[ML union genes found in acinar RNA matrix]:", length(genes_use), "\n")

if (length(genes_use) < 10) {
  stop("Too few ML-union genes found in expression matrix: ", length(genes_use))
}

counts_mat <- counts_mat_full[genes_use, , drop = FALSE]
counts_mat <- counts_mat[Matrix::rowSums(counts_mat) > 0, , drop = FALSE]

cat("[Genes retained after removing all-zero rows]:", nrow(counts_mat), "\n")

if (nrow(counts_mat) < 10) {
  stop("Too few non-zero genes left after filtering: ", nrow(counts_mat))
}

cell_meta <- acinar_obj@meta.data[colnames(acinar_obj), , drop = FALSE]
cell_meta$group <- factor(cell_meta$group, levels = c("CON", "PRE", "T2D"))
cell_meta$Stage_selected <- factor(cell_meta$Stage_selected, levels = c("CON", "PRE", "T2D"))

gene_meta <- data.frame(
  gene_short_name = rownames(counts_mat),
  row.names = rownames(counts_mat),
  stringsAsFactors = FALSE
)

write.csv(
  data.frame(gene = rownames(counts_mat)),
  file.path(out_dir, "acinar_ml_union_genes_used_for_trajectory.csv"),
  row.names = FALSE
)

# =========================================================
# 6.  Monocle3 CDS
# =========================================================
cat("\n========== 4. Building CDS ==========\n")

cds <- new_cell_data_set(
  expression_data = counts_mat,
  cell_metadata   = cell_meta,
  gene_metadata   = gene_meta
)

num_dim_use <- min(50, nrow(counts_mat) - 1, ncol(counts_mat) - 1)
if (!is.finite(num_dim_use) || num_dim_use < 2) {
  stop("num_dim_use is too small: ", num_dim_use)
}
cat("[num_dim used in preprocess_cds]:", num_dim_use, "\n")

# =========================================================
# 7. Monocle3 native： +  + UMAP +  + 
# =========================================================
cat("\n========== 5. Monocle3 trajectory on ML-union genes ==========\n")

set.seed(1)
cds <- preprocess_cds(cds, num_dim = num_dim_use)

if (DO_ALIGN &&
    "orig.ident" %in% colnames(colData(cds)) &&
    length(unique(colData(cds)$orig.ident)) > 1) {
  cat("[INFO] Running align_cds on orig.ident\n")
  cds <- align_cds(cds, alignment_group = "orig.ident")
}

umap_neighbors <- min(50, max(10, floor(ncol(cds) / 100)))
cat("[UMAP n_neighbors]:", umap_neighbors, "\n")

cds <- reduce_dimension(
  cds,
  reduction_method = "UMAP",
  preprocess_method = "PCA",
  umap.metric = "cosine",
  umap.min_dist = 0.1,
  umap.n_neighbors = 30,
  verbose = TRUE
)

cds <- cluster_cells(
  cds,
  reduction_method = "UMAP",
  cluster_method = "leiden",
  k = 20,
  resolution = 0.00001,
  verbose = TRUE
)

cat("\n[Cluster distribution]\n")
print(table(clusters(cds)))
cat("\n[Partition distribution]\n")
print(table(partitions(cds)))

cds <- learn_graph(
  cds,
  use_partition = TRUE,
  verbose = TRUE
)

# =========================================================
# 8.  root： Stage_selected = CON
# =========================================================
cat("\n========== 6. Choosing root ==========\n")

manual_root_nodes <- NULL
if (is.null(manual_root_nodes)) {
  root_nodes <- choose_root_nodes_by_group(
    cds,
    group_col = "group",
    early_group = "CON"
  )
} else {
  root_nodes <- manual_root_nodes
}

cat("[Chosen root nodes]\n")

print(root_nodes)

cds <- order_cells(cds, root_pr_nodes = root_nodes)

# =========================================================
# 9. 
# =========================================================
cat("\n========== 7. Plotting ==========\n")

stage_colors <- c(CON = "#D4D4D4", PRE = "#F4B7AD", T2D = "#90A4C4")

colData(cds)$Stage_selected <- factor(colData(cds)$Stage_selected,
                                      levels = c("CON", "PRE", "T2D"))
colData(cds)$group <- factor(colData(cds)$group,
                             levels = c("CON", "PRE", "T2D"))

p_stage <- plot_cells(
  cds,
  color_cells_by = "Stage_selected",
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = FALSE,
  graph_label_size = 3,
  cell_size = 0.5
) +
  scale_color_manual(values = stage_colors, name = "Stage_selected") +
  ggtitle("Stage-selected Acinar cells (ML-union genes + Monocle3 native)") +
  theme(legend.position = "right")
p_stage <- rasterize_points(p_stage, dpi = 400)
save_plot(p_stage, file.path(out_dir, "umap_by_stage.pdf"), 7, 5.5)

p_group <- plot_cells(
  cds,
  color_cells_by = "group",
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = FALSE,
  graph_label_size = 3,
  cell_size = 0.5
) +
  scale_color_manual(values = stage_colors, name = "group") +
  ggtitle("") +
  theme(legend.position = "right")
p_group <- rasterize_points(p_group, dpi = 400)
save_plot(p_group, file.path(out_dir, "umap_by_group.pdf"), 6, 5.5)

p_cluster <- plot_cells(
  cds,
  color_cells_by = "cluster",
  label_groups_by_cluster = TRUE,
  label_leaves = TRUE,
  label_branch_points = FALSE,
  cell_size = 0.5
) + ggtitle("")
p_cluster <- rasterize_points(p_cluster, dpi = 400)
save_plot(p_cluster, file.path(out_dir, "umap_by_cluster.pdf"), 6, 5.5)

p_pt <- plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = FALSE,
  cell_size = 0.5
) +
  ggtitle("") +
  scale_color_viridis_c()
p_pt <- rasterize_points(p_pt, dpi = 400)
save_plot(p_pt, file.path(out_dir, "umap_by_pseudotime.pdf"), 6, 5.5)

# =========================================================
# 10. pseudotime 
# =========================================================
cat("\n========== 8. Pseudotime summary ==========\n")

pt_vec <- pseudotime(cds)

pt_df <- data.frame(
  cell = colnames(cds),
  pseudotime = as.numeric(pt_vec),
  Stage_selected = colData(cds)$Stage_selected,
  group = colData(cds)$group,
  orig.ident = colData(cds)$orig.ident,
  partition = as.character(partitions(cds)),
  cluster = as.character(clusters(cds)),
  stringsAsFactors = FALSE
)

pt_finite <- pt_df[is.finite(pt_df$pseudotime), ]

cat("[Finite pseudotime counts]\n")
print(table(is.finite(pt_df$pseudotime)))

cat("\n[Pseudotime summary by Stage_selected]\n")
print(tapply(pt_finite$pseudotime, pt_finite$Stage_selected, summary))

cat("\n[Pseudotime median by Stage_selected]\n")
medians <- tapply(pt_finite$pseudotime, pt_finite$Stage_selected, median)
print(medians)

if (medians["CON"] < medians["PRE"] && medians["PRE"] < medians["T2D"]) {
  cat("\n✅ SUCCESS: pseudotime median order is CON < PRE < T2D\n")
} else {
  cat("\n⚠️ WARNING: pseudotime median order is NOT CON < PRE < T2D\n")
  cat("   Order:", names(sort(medians)), "\n")
}

cat("\n[Kruskal-Wallis test]\n")
print(kruskal.test(pseudotime ~ Stage_selected, data = pt_finite))

cat("\n[Pairwise Wilcoxon test]\n")
print(pairwise.wilcox.test(pt_finite$pseudotime, pt_finite$Stage_selected, p.adjust.method = "BH"))

p_violin <- ggplot(pt_finite, aes(x = Stage_selected, y = pseudotime, fill = Stage_selected)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.12, outlier.size = 0.3, fill = "white", alpha = 0.7) +
  scale_fill_manual(values = stage_colors) +
  labs(title = NULL, x = "Stage", y = "Pseudotime") +
  theme_classic() +
  theme(legend.position = "none")
save_plot(p_violin, file.path(out_dir, "pseudotime_violin_by_stage.pdf"), 4, 4)

p_donor <- ggplot(pt_finite, aes(x = reorder(orig.ident, pseudotime, median),
                                 y = pseudotime, fill = Stage_selected)) +
  geom_boxplot(outlier.size = 0.3) +
  scale_fill_manual(values = stage_colors) +
  labs(title = NULL, x = "Donor", y = "Pseudotime") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_donor, file.path(out_dir, "pseudotime_boxplot_by_donor.pdf"), 8, 5)

write.csv(pt_df, file.path(out_dir, "pseudotime_table.csv"), row.names = FALSE)


# =========================================================
# 11. tradeSeq dynamic genes on Monocle3 pseudotime
#  pseudotime_table.csv 
# =========================================================
suppressPackageStartupMessages({
  library(tradeSeq)
  library(SingleCellExperiment)
  library(BiocParallel)
})

cat("\n========== 9. tradeSeq dynamic gene analysis ==========\n")

out_dir_ts <- file.path(out_dir, "tradeSeq")
dir.create(out_dir_ts, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------
# 11.1 
# ---------------------------------------------------------
# TRUE  = 
# FALSE =  ML-union genes 
USE_FULL_GENESET_FOR_TRADESEQ <- FALSE

#  partition， "1"，
MANUAL_PARTITION <- NULL

#  evaluateK；， FALSE
RUN_EVALUATEK <- TRUE

#  evaluateK，
NKNOTS_DEFAULT <- 6

# tradeSeq 
QVALUE_CUTOFF <- 0.05
DETECT_FRAC_CUTOFF <- 0.05

# ---------------------------------------------------------
# 11.2  partition
# ---------------------------------------------------------
if (!is.null(MANUAL_PARTITION)) {
  target_partition <- as.character(MANUAL_PARTITION)
} else {
  target_partition <- names(sort(table(pt_finite$partition), decreasing = TRUE))[1]
}

cells_use <- pt_finite$cell[pt_finite$partition == target_partition]

cat("[target partition]:", target_partition, "\n")
cat("[cells used]:", length(cells_use), "\n")

if (length(cells_use) < 50) {
  warning("Cells in target partition are quite few: ", length(cells_use))
}

meta_use <- pt_finite[match(cells_use, pt_finite$cell), , drop = FALSE]
stopifnot(all(meta_use$cell == cells_use))

# ---------------------------------------------------------
# 11.3 
# ---------------------------------------------------------
if (USE_FULL_GENESET_FOR_TRADESEQ) {
  cat("[tradeSeq gene set] using FULL genes\n")
  counts_source <- counts_mat_full
} else {
  cat("[tradeSeq gene set] using ML-union genes only\n")
  counts_source <- counts_mat
}

counts_ts <- counts_source[, cells_use, drop = FALSE]

#  pseudotime 
ord <- order(meta_use$pseudotime)
meta_use  <- meta_use[ord, , drop = FALSE]
counts_ts <- counts_ts[, meta_use$cell, drop = FALSE]

# ---------------------------------------------------------
# 11.4 
# ---------------------------------------------------------
min_cells_detected <- max(10, ceiling(0.05 * ncol(counts_ts)))
keep_gene <- Matrix::rowSums(counts_ts > 0) >= min_cells_detected
counts_ts <- counts_ts[keep_gene, , drop = FALSE]

cat("[genes kept after filtering]:", nrow(counts_ts), "\n")
cat("[cells kept]:", ncol(counts_ts), "\n")

if (nrow(counts_ts) < 20) {
  stop("Too few genes left for tradeSeq: ", nrow(counts_ts))
}

# ---------------------------------------------------------
# 11.5  tradeSeq 
#  lineage
# ---------------------------------------------------------
pseudotime_mat <- matrix(meta_use$pseudotime, ncol = 1)
rownames(pseudotime_mat) <- meta_use$cell
colnames(pseudotime_mat) <- "Lineage1"

cellWeights <- matrix(1, nrow = nrow(meta_use), ncol = 1)
rownames(cellWeights) <- meta_use$cell
colnames(cellWeights) <- "Lineage1"

# ---------------------------------------------------------
# 11.6 evaluateK（）
# ： nGenes， 500 
# ---------------------------------------------------------
nknots_use <- NKNOTS_DEFAULT

if (RUN_EVALUATEK) {
  cat("\n========== 9.1 evaluateK ==========\n")

  nGenes_eval <- min(200, nrow(counts_ts))

  set.seed(1)
  icMat <- evaluateK(
    counts = counts_ts,
    pseudotime = pseudotime_mat,
    cellWeights = cellWeights,
    k = 3:7,
    nGenes = nGenes_eval,
    plot = FALSE,
    verbose = TRUE,
    parallel = FALSE
  )

  write.csv(
    icMat,
    file.path(out_dir_ts, "tradeSeq_evaluateK.csv"),
    row.names = TRUE
  )

  cat("[evaluateK finished] nGenes =", nGenes_eval, "\n")
  cat("[Using nknots]:", nknots_use, "\n")
} else {
  cat("[Skip evaluateK] using nknots =", nknots_use, "\n")
}

# ---------------------------------------------------------
# 11.7 fitGAM
# ---------------------------------------------------------
cat("\n========== 9.2 fitGAM ==========\n")

set.seed(1)
sce_ts <- fitGAM(
  counts = counts_ts,
  pseudotime = pseudotime_mat,
  cellWeights = cellWeights,
  nknots = nknots_use,
  verbose = TRUE,
  parallel = FALSE
)


# ---------------------------------------------------------
# 11.8  sce_ts
#  plotSmoothers(pointCol = "Stage_selected") 
# ---------------------------------------------------------
colData(sce_ts)$Stage_selected <- factor(
  meta_use$Stage_selected,
  levels = c("CON", "PRE", "T2D")
)
colData(sce_ts)$group <- factor(
  meta_use$group,
  levels = c("CON", "PRE", "T2D")
)
colData(sce_ts)$orig.ident <- meta_use$orig.ident
colData(sce_ts)$partition  <- meta_use$partition
colData(sce_ts)$cluster    <- meta_use$cluster

cat("\n[colData(sce_ts)]\n")
print(colnames(colData(sce_ts)))
cat("\n[Stage_selected distribution in sce_ts]\n")
print(table(colData(sce_ts)$Stage_selected, useNA = "ifany"))

# ---------------------------------------------------------
# 11.9 associationTest：
# ---------------------------------------------------------
cat("\n========== 9.3 associationTest ==========\n")

assoc_res <- associationTest(sce_ts)
assoc_res <- as.data.frame(assoc_res)
assoc_res$gene <- rownames(assoc_res)
assoc_res$qvalue <- p.adjust(assoc_res$pvalue, method = "BH")

detect_frac_vec <- Matrix::rowMeans(counts_ts > 0)
mean_count_vec  <- Matrix::rowMeans(counts_ts)

assoc_res$detect_frac <- detect_frac_vec[assoc_res$gene]
assoc_res$mean_count  <- mean_count_vec[assoc_res$gene]

assoc_res <- assoc_res %>%
  dplyr::arrange(qvalue, dplyr::desc(waldStat))

if (!exists("QVALUE_CUTOFF")) QVALUE_CUTOFF <- 0.05
if (!exists("DETECT_FRAC_CUTOFF")) DETECT_FRAC_CUTOFF <- 0.05

write.csv(
  assoc_res,
  file.path(out_dir_ts, "tradeSeq_associationTest_all.csv"),
  row.names = FALSE
)

assoc_sig <- assoc_res %>%
  dplyr::filter(
    !is.na(qvalue),
    qvalue < QVALUE_CUTOFF,
    detect_frac > DETECT_FRAC_CUTOFF
  ) %>%
  dplyr::arrange(qvalue, dplyr::desc(waldStat))

write.csv(
  assoc_sig,
  file.path(out_dir_ts, "tradeSeq_associationTest_sig.csv"),
  row.names = FALSE
)

cat("[significant dynamic genes]:", nrow(assoc_sig), "\n")
print(head(assoc_sig[, c("gene", "waldStat", "pvalue", "qvalue", "detect_frac")], 20))

write.table(
  assoc_sig$gene,
  file.path(out_dir_ts, "tradeSeq_dynamic_gene_list.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# ---------------------------------------------------------
# 11.10 ：/
# ---------------------------------------------------------
cat("\n========== 9.4 startVsEndTest ==========\n")

sev_res <- startVsEndTest(sce_ts)
sev_res <- as.data.frame(sev_res)
sev_res$gene <- rownames(sev_res)
sev_res$qvalue <- p.adjust(sev_res$pvalue, method = "BH")
sev_res$detect_frac <- detect_frac_vec[sev_res$gene]
sev_res$mean_count  <- mean_count_vec[sev_res$gene]

sev_res <- sev_res %>%
  dplyr::arrange(qvalue, dplyr::desc(waldStat))

if (!exists("QVALUE_CUTOFF")) QVALUE_CUTOFF <- 0.05
if (!exists("DETECT_FRAC_CUTOFF")) DETECT_FRAC_CUTOFF <- 0.05

write.csv(
  sev_res,
  file.path(out_dir_ts, "tradeSeq_startVsEndTest_all.csv"),
  row.names = FALSE
)

sev_sig <- sev_res %>%
  dplyr::filter(
    !is.na(qvalue),
    qvalue < QVALUE_CUTOFF,
    detect_frac > DETECT_FRAC_CUTOFF
  ) %>%
  dplyr::arrange(qvalue, dplyr::desc(waldStat))

write.csv(
  sev_sig,
  file.path(out_dir_ts, "tradeSeq_startVsEndTest_sig.csv"),
  row.names = FALSE
)

cat("[significant start-vs-end genes]:", nrow(sev_sig), "\n")

# ---------------------------------------------------------
# 11.11  smoother （ Stage_selected ）
# ---------------------------------------------------------
cat("\n========== 9.5 plotting tradeSeq smoothers ==========\n")

plot_tradeseq_gene <- function(gene, sce_obj, counts_mat, stage_colors) {
  p <- tradeSeq::plotSmoothers(
    models = sce_obj,
    counts = counts_mat,
    gene = gene,
    pointCol = "Stage_selected",
    curvesCols = "black"
  ) +
    scale_color_manual(
      values = stage_colors,
      breaks = c("CON", "PRE", "T2D"),
      drop = FALSE,
      name = "Stage_selected"
    ) +
    labs(
      title = gene,
      x = "Pseudotime",
      y = "Log(expression + 1)"
    ) +
    theme_bw(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 13),
      legend.position = "right"
    )

  p <- rasterize_points(p, dpi = 400)
  return(p)
}

#  associationTest 
candidate_genes <- assoc_sig$gene

# 
manual_genes <- c("NR5A2",
                  "RB1CC1",
                  "MECOM",
                
                  "GMDS",
                 
                  "MT1E",
                  "IGF1R")

manual_genes_present <- manual_genes[manual_genes %in% rownames(counts_ts)]
manual_genes_missing <- setdiff(manual_genes, manual_genes_present)

if (length(manual_genes_missing) > 0) {
  cat("[Missing requested genes]:", paste(manual_genes_missing, collapse = ", "), "\n")
}

if (length(manual_genes_present) > 0) {
  genes_to_plot <- manual_genes_present
} else {
  genes_to_plot <- head(candidate_genes, 12)
}

cat("[Genes selected for smoother plots]\n")
print(genes_to_plot)

if (length(genes_to_plot) > 0) {
  plot_list_ts <- lapply(
    genes_to_plot,
    function(g) plot_tradeseq_gene(g, sce_ts, counts_ts, stage_colors)
  )

  p_ts_panel <- patchwork::wrap_plots(plot_list_ts, ncol = 2, guides = "collect") &
    theme(
      legend.position = "right",
      strip.text = element_text(size = 14),
      plot.title = element_text(size = 18),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 13)
    )

  save_plot(
    p_ts_panel,
    file.path(out_dir_ts, "tradeSeq_stage_colored_smoothers.pdf"),
    width = 10,
    height = max(8, ceiling(length(genes_to_plot) / 2) * 2.8),
    dpi = 300
  )
}

# =========================================================
# 12.  qs 
# =========================================================
cat("\n========== 10. Saving all objects into one .qs ==========" ,"\n")

all_object_names <- ls(envir = .GlobalEnv, all.names = TRUE)
all_results <- mget(all_object_names, envir = .GlobalEnv, inherits = FALSE)
qs_file <- file.path(out_dir, "acinar_all_objects.qs")

qsave(
  all_results,
  qs_file
)

cat("Saved:", qs_file, "\n")




