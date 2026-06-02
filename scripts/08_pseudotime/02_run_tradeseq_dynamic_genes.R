#!/usr/bin/env Rscript

# Script name: 02_run_tradeseq_dynamic_genes.R
# Purpose: Cleaned public workflow script for the T2D single-cell multiome analysis.
# Input: Paths and analysis settings are read from configs/config.yaml.
# Output: Module-specific outputs are written under the configured results directory.
# Main steps: Load configuration, run the original analysis logic, and export reproducible outputs.
# Example command: Rscript scripts/08_pseudotime/02_run_tradeseq_dynamic_genes.R

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
# Monocle3 native trajectory on stage-selected beta cells
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
out_dir <- "results/downstream/riskcell_monocle/beta/union_native"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# 1. 
# =========================================================
selected_path <- "results/downstream/riskcell/beta_selected_cells_for_pseudotime.tsv"

ml_path_con_pre <- "results/model_outputs/rna/β_cell/CON_vs_PRE_shap_marker_genes.csv"
ml_path_pre_t2d <- "results/model_outputs/rna/β_cell/PRE_vs_T2D_shap_marker_genes.csv"
ml_path_con_t2d <- "results/model_outputs/rna/β_cell/CON_vs_T2D_shap_marker_genes.csv"

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

umap_box_theme <- function(base_size = 18) {
  theme_classic(base_size = base_size) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.3),
      aspect.ratio = 1,
      axis.line = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 1.1),
      axis.ticks.length = grid::unit(0.22, "cm"),
      axis.title = element_text(size = 22, face = "bold", color = "black"),
      axis.text = element_text(size = 18, face = "bold", color = "black"),
      legend.title = element_text(size = 18, face = "bold", color = "black"),
      legend.text = element_text(size = 16, face = "bold", color = "black"),
      plot.title = element_text(size = 22, face = "bold", hjust = 0.5, color = "black"),
      strip.text = element_text(size = 18, face = "bold", color = "black"),
      text = element_text(face = "bold", color = "black")
    )
}

style_umap_plot <- function(p, text_size = 5) {
  p <- p + umap_box_theme()

  for (i in seq_along(p$layers)) {
    geom_class <- class(p$layers[[i]]$geom)[1]
    if (geom_class %in% c("GeomText", "GeomLabel", "GeomTextRepel", "GeomLabelRepel")) {
      p$layers[[i]]$aes_params$fontface <- "bold"
      layer_size <- p$layers[[i]]$aes_params$size
      if (is.null(layer_size) || (is.numeric(layer_size) && layer_size < text_size)) {
        p$layers[[i]]$aes_params$size <- text_size
      }
    }
  }

  p
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
# 3.  stage-selected beta cells
# =========================================================
cat("\n========== 1. Loading stage-selected beta cells ==========\n")

sel <- fread(selected_path)

combined <- qread(seurat_path)
Idents(combined) <- combined$celltype
beta_obj <- subset(combined, idents = "β cell")
rm(combined); gc()

cells_keep <- intersect(colnames(beta_obj), sel$Cell_id)
beta_obj <- subset(beta_obj, cells = cells_keep)

beta_obj$Stage_selected <- sel$Stage_selected[match(colnames(beta_obj), sel$Cell_id)]
beta_obj$Stage_selected <- factor(beta_obj$Stage_selected, levels = c("CON", "PRE", "T2D"))

cat("Selected beta cells:", ncol(beta_obj), "\n")
cat("\n[Stage_selected distribution]\n")
print(table(beta_obj$Stage_selected, useNA = "ifany"))
cat("\n[Original group distribution]\n")
print(table(beta_obj$group, useNA = "ifany"))
cat("\n[Sample distribution]\n")
print(table(beta_obj$orig.ident, useNA = "ifany"))

DefaultAssay(beta_obj) <- "RNA"
beta_obj <- tryCatch({
  JoinLayers(beta_obj)
}, error = function(e) {
  message("[INFO] JoinLayers skipped: ", conditionMessage(e))
  beta_obj
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
  file.path(out_dir, "beta_ml_union_gene_membership.csv"),
  row.names = FALSE
)

# =========================================================
# 5.  ML-union genes 
# =========================================================
cat("\n========== 3. Restricting matrix to ML-union genes ==========\n")

counts_mat_full <- get_layer_mat(beta_obj, assay = "RNA", layer = "counts")
counts_mat_full <- as(counts_mat_full, "dgCMatrix")

genes_use <- intersect(ml_union, rownames(counts_mat_full))
cat("[ML union genes found in beta RNA matrix]:", length(genes_use), "\n")

if (length(genes_use) < 10) {
  stop("Too few ML-union genes found in expression matrix: ", length(genes_use))
}

counts_mat <- counts_mat_full[genes_use, , drop = FALSE]
counts_mat <- counts_mat[Matrix::rowSums(counts_mat) > 0, , drop = FALSE]

cat("[Genes retained after removing all-zero rows]:", nrow(counts_mat), "\n")

if (nrow(counts_mat) < 10) {
  stop("Too few non-zero genes left after filtering: ", nrow(counts_mat))
}

cell_meta <- beta_obj@meta.data[colnames(beta_obj), , drop = FALSE]
cell_meta$group <- factor(cell_meta$group, levels = c("CON", "PRE", "T2D"))
cell_meta$Stage_selected <- factor(cell_meta$Stage_selected, levels = c("CON", "PRE", "T2D"))

gene_meta <- data.frame(
  gene_short_name = rownames(counts_mat),
  row.names = rownames(counts_mat),
  stringsAsFactors = FALSE
)

write.csv(
  data.frame(gene = rownames(counts_mat)),
  file.path(out_dir, "beta_ml_union_genes_used_for_trajectory.csv"),
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
  umap.min_dist = 0.3,
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

stage_colors <- c(CON = "#B7CCDD", PRE = "#E8D7A9", T2D = "#C48986")

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
  graph_label_size = 5,
  cell_size = 0.5
) +
  scale_color_manual(values = stage_colors, name = "Stage_selected")
p_stage <- style_umap_plot(p_stage) + theme(legend.position = "right")
p_stage <- rasterize_points(p_stage, dpi = 400)
save_plot(p_stage, file.path(out_dir, "umap_by_stage.pdf"), 7, 5.5)

p_group <- plot_cells(
  cds,
  color_cells_by = "group",
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = FALSE,
  graph_label_size = 5,
  cell_size = 0.5
) +
  scale_color_manual(values = stage_colors, name = "group")
p_group <- style_umap_plot(p_group) + theme(legend.position = "right")
p_group <- rasterize_points(p_group, dpi = 400)
save_plot(p_group, file.path(out_dir, "umap_by_group.pdf"), 6, 5.5)

p_cluster <- plot_cells(
  cds,
  color_cells_by = "cluster",
  label_groups_by_cluster = TRUE,
  label_leaves = TRUE,
  label_branch_points = FALSE,
  group_label_size = 6,
  graph_label_size = 5,
  cell_size = 0.5
)
p_cluster <- style_umap_plot(p_cluster)
p_cluster <- rasterize_points(p_cluster, dpi = 400)
save_plot(p_cluster, file.path(out_dir, "umap_by_cluster.pdf"), 6, 5.5)

p_pt <- plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = FALSE,
  graph_label_size = 5,
  cell_size = 0.5
) +
  scale_color_viridis_c()
p_pt <- style_umap_plot(p_pt)
p_pt <- rasterize_points(p_pt, dpi = 400)
save_plot(p_pt, file.path(out_dir, "umap_by_pseudotime.pdf"), 6, 5.5)

p_points <- plot_cells(
  cds,
  color_cells_by = "Stage_selected",
  label_groups_by_cluster = FALSE,
  label_branch_points = FALSE,
  label_roots = FALSE,
  label_leaves = FALSE,
  label_principal_points = TRUE,
  graph_label_size = 5,
  cell_size = 0.5
) +
  scale_color_manual(values = stage_colors, name = "Stage_selected")
p_points <- style_umap_plot(p_points) + theme(legend.position = "right")
p_points <- rasterize_points(p_points, dpi = 400)
save_plot(p_points, file.path(out_dir, "umap_principal_points.pdf"), 7, 5.5)

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
# 13. ONLY partition 1: dynamic analysis on the main trajectory
# =========================================================
cat("\n========== 10. Partition 1 only: dynamic analysis ==========\n")

out_dir_p1 <- file.path(out_dir, "partition1_only")
dir.create(out_dir_p1, recursive = TRUE, showWarnings = FALSE)

cat("\n[Full partition distribution]\n")
print(table(partitions(cds)))

cat("\n[Partition x original group]\n")
print(table(partitions(cds), colData(cds)$group))

cat("\n[Partition x Stage_selected]\n")
print(table(partitions(cds), colData(cds)$Stage_selected))

main_partition <- "1"
cells_p1 <- colnames(cds)[as.character(partitions(cds)) == main_partition]

if (length(cells_p1) < 50) {
  stop("Too few cells in partition ", main_partition, ": ", length(cells_p1))
}

cds_p1 <- cds[, cells_p1]

cat("\n[Retained cells in partition 1]: ", ncol(cds_p1), "\n", sep = "")
cat("[Retained genes in partition 1]: ", nrow(cds_p1), "\n", sep = "")

cat("\n[Partition 1 group distribution]\n")
print(table(colData(cds_p1)$group, useNA = "ifany"))

cat("\n[Partition 1 Stage_selected distribution]\n")
print(table(colData(cds_p1)$Stage_selected, useNA = "ifany"))


# ---------------------------------------------------------
# 13.2  partition 1  learn_graph
#  partition ，
# ---------------------------------------------------------
cat("\n========== 10.1 Re-learning graph on partition 1 ==========\n")

cds_p1 <- learn_graph(
  cds_p1,
  use_partition = FALSE,
  verbose = TRUE
)

# ---------------------------------------------------------
# 13.3  root， partition 1  CON 
#  Stage_selected  group
# ---------------------------------------------------------
cat("\n========== 10.2 Choosing root for partition 1 ==========\n")

manual_root_nodes_p1 <- NULL
#  root，：
# manual_root_nodes_p1 <- c("Y_21")

if (is.null(manual_root_nodes_p1)) {
  root_nodes_p1 <- choose_root_nodes_by_group(
    cds_p1,
    group_col = "Stage_selected",
    early_group = "CON"
  )
} else {
  root_nodes_p1 <- manual_root_nodes_p1
}

cat("[Chosen root nodes for partition 1]\n")
print(root_nodes_p1)

if (length(root_nodes_p1) == 0) {
  stop("No root nodes found for partition 1. Please inspect the graph and set manual_root_nodes_p1.")
}

cds_p1 <- order_cells(cds_p1, root_pr_nodes = root_nodes_p1)

# ---------------------------------------------------------
# 13.4 1（partition 1）
# ---------------------------------------------------------
cat("\n========== 10.3 Partition 1 pseudotime violin by stage ==========" ,"\n")

pt_p1_plot <- data.frame(
  cell = colnames(cds_p1),
  pseudotime = as.numeric(pseudotime(cds_p1)),
  Stage_selected = colData(cds_p1)$Stage_selected,
  stringsAsFactors = FALSE
)

pt_p1_plot <- pt_p1_plot %>%
  filter(is.finite(pseudotime)) %>%
  mutate(Stage_selected = factor(Stage_selected, levels = c("CON", "PRE", "T2D")))

if (nrow(pt_p1_plot) == 0) {
  warning("No finite pseudotime values in partition 1 for violin plot.")
} else {
  p_violin_p1 <- ggplot(pt_p1_plot, aes(x = Stage_selected, y = pseudotime, fill = Stage_selected)) +
    geom_violin(trim = FALSE, alpha = 0.7, color = NA) +
    geom_boxplot(width = 0.12, outlier.size = 0.3, fill = "white", alpha = 0.8) +
    scale_fill_manual(values = stage_colors) +
    labs(x = "Stage", y = "Pseudotime") +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_blank()
    )

  save_plot(
    p_violin_p1,
    file.path(out_dir, "pseudotime_violin_by_stage_partition1.pdf"),
    width = 4.5,
    height = 4.5,
    dpi = 300
  )

  write.csv(
    pt_p1_plot,
    file.path(out_dir, "pseudotime_partition1_table.csv"),
    row.names = FALSE
  )

  cat("[Partition 1 pseudotime summary by stage]\n")
  print(tapply(pt_p1_plot$pseudotime, pt_p1_plot$Stage_selected, summary))

  cat("\n[Partition 1 Kruskal-Wallis test]\n")
  print(kruskal.test(pseudotime ~ Stage_selected, data = pt_p1_plot))

  cat("\n[Partition 1 Pairwise Wilcoxon test]\n")
  print(pairwise.wilcox.test(pt_p1_plot$pseudotime, pt_p1_plot$Stage_selected, p.adjust.method = "BH"))
}

# =========================================================
# 14. tradeSeq on partition 1 main trajectory
#    append after cds_p1 has been ordered
# =========================================================
cat("\n========== 11. tradeSeq on partition 1 ==========\n")

suppressPackageStartupMessages({
  library(tradeSeq)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(edgeR)
  library(ggplot2)
  library(patchwork)
})

out_dir_ts <- file.path(out_dir_p1, "tradeSeq")
dir.create(out_dir_ts, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------
# 14.1 helper
# ---------------------------------------------------------
get_p_col <- function(df) {
  pcols <- grep("^pvalue", colnames(df), value = TRUE)
  if (length(pcols) == 0) stop("No pvalue column found.")
  pcols[1]
}

make_result_df <- function(res_mat) {
  df <- as.data.frame(res_mat)
  df$gene <- rownames(df)
  df <- df[, c("gene", setdiff(colnames(df), "gene")), drop = FALSE]
  pcol <- get_p_col(df)
  df$qvalue <- p.adjust(df[[pcol]], method = "BH")
  df <- df[order(df$qvalue, df[[pcol]]), , drop = FALSE]
  rownames(df) <- NULL
  df
}

# ---------------------------------------------------------
# 14.2 extract count matrix + pseudotime for finite cells
# ---------------------------------------------------------
cat("\n========== 11.1 Preparing counts and pseudotime ==========\n")

pt_p1 <- pseudotime(cds_p1)
keep_cells_ts <- is.finite(pt_p1)

cat("[Total cells in cds_p1]:", ncol(cds_p1), "\n")
cat("[Finite pseudotime cells kept for tradeSeq]:", sum(keep_cells_ts), "\n")

if (sum(keep_cells_ts) < 100) {
  stop("Too few finite-pseudotime cells for tradeSeq: ", sum(keep_cells_ts))
}

cds_ts <- cds_p1[, keep_cells_ts]

counts_ts <- SummarizedExperiment::assay(cds_ts, "counts")
counts_ts <- as(counts_ts, "dgCMatrix")

pt_ts <- as.numeric(pseudotime(cds_ts))

#  lineage： pseudotime +  1  cellWeights
pseudotime_mat <- matrix(
  pt_ts,
  ncol = 1,
  dimnames = list(colnames(cds_ts), "Lineage1")
)

cellWeights_mat <- matrix(
  1,
  nrow = length(pt_ts),
  ncol = 1,
  dimnames = list(colnames(cds_ts), "Lineage1")
)

stopifnot(
  ncol(counts_ts) == nrow(pseudotime_mat),
  nrow(pseudotime_mat) == nrow(cellWeights_mat)
)

cat("[Genes for tradeSeq]:", nrow(counts_ts), "\n")
cat("[Cells for tradeSeq]:", ncol(counts_ts), "\n")

# ---------------------------------------------------------
# 14.3 optional: evaluateK
#  evaluateK  nknots
# ---------------------------------------------------------
cat("\n========== 11.2 evaluateK ==========\n")

RUN_EVALUATE_K <- TRUE
k_grid <- 4:7
nknots_use <- 6   # ； evaluateK 

if (RUN_EVALUATE_K) {
  set.seed(1)

  pdf(file.path(out_dir_ts, "evaluateK_k4to7.pdf"), width = 6, height = 5)
  icMat <- tradeSeq::evaluateK(
    counts = counts_ts,
    pseudotime = pseudotime_mat,
    cellWeights = cellWeights_mat,
    k = k_grid,
    nGenes = min(200, nrow(counts_ts)),
    plot = TRUE,
    verbose = TRUE
  )
  dev.off()

  write.csv(
    icMat,
    file.path(out_dir_ts, "evaluateK_k4to7_matrix.csv"),
    row.names = TRUE
  )

  cat("[evaluateK done] Current nknots_use =", nknots_use, "\n")
}

# ---------------------------------------------------------
# 14.4 fitGAM
# ---------------------------------------------------------
cat("\n========== 11.3 fitGAM ==========\n")

set.seed(1)
sce_ts <- tradeSeq::fitGAM(
  counts = counts_ts,
  pseudotime = pseudotime_mat,
  cellWeights = cellWeights_mat,
  nknots = nknots_use,
  verbose = TRUE,
  parallel = FALSE,
  sce = TRUE
)

# ，
colData(sce_ts)$Stage_selected <- colData(cds_ts)$Stage_selected
colData(sce_ts)$group <- colData(cds_ts)$group
colData(sce_ts)$orig.ident <- colData(cds_ts)$orig.ident
colData(sce_ts)$monocle_pseudotime <- pt_ts

saveRDS(sce_ts, file.path(out_dir_ts, "tradeSeq_fitGAM_sce.rds"))

# ---------------------------------------------------------
# 14.5 associationTest:  pseudotime 
# ---------------------------------------------------------
cat("\n========== 11.4 associationTest ==========\n")

assoc_res <- tradeSeq::associationTest(sce_ts)
assoc_df <- make_result_df(assoc_res)

write.csv(
  assoc_df,
  file.path(out_dir_ts, "associationTest_results.csv"),
  row.names = FALSE
)

assoc_sig <- assoc_df[!is.na(assoc_df$qvalue) & assoc_df$qvalue < 0.05, , drop = FALSE]

write.csv(
  assoc_sig,
  file.path(out_dir_ts, "associationTest_sig_q005.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(gene = assoc_sig$gene),
  file.path(out_dir_ts, "associationTest_dynamic_genes_q005.csv"),
  row.names = FALSE
)

cat("[associationTest significant genes, q < 0.05]:", nrow(assoc_sig), "\n")

# ---------------------------------------------------------
# 14.6 startVsEndTest:  vs 
# ---------------------------------------------------------
cat("\n========== 11.5 startVsEndTest ==========\n")

startend_res <- tradeSeq::startVsEndTest(sce_ts)
startend_df <- make_result_df(startend_res)

write.csv(
  startend_df,
  file.path(out_dir_ts, "startVsEndTest_results.csv"),
  row.names = FALSE
)

startend_sig <- startend_df[!is.na(startend_df$qvalue) & startend_df$qvalue < 0.05, , drop = FALSE]

write.csv(
  startend_sig,
  file.path(out_dir_ts, "startVsEndTest_sig_q005.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(gene = startend_sig$gene),
  file.path(out_dir_ts, "startVsEndTest_dynamic_genes_q005.csv"),
  row.names = FALSE
)

cat("[startVsEndTest significant genes, q < 0.05]:", nrow(startend_sig), "\n")

# ---------------------------------------------------------
# 14.7  Monocle3 graph_test 
# ---------------------------------------------------------
cat("\n========== 11.6 overlap with Monocle3 graph_test ==========\n")

if (exists("dynamic_genes_p1")) {
  overlap_assoc_graph <- intersect(assoc_sig$gene, dynamic_genes_p1)
  overlap_startend_graph <- intersect(startend_sig$gene, dynamic_genes_p1)

  write.csv(
    data.frame(gene = overlap_assoc_graph),
    file.path(out_dir_ts, "overlap_associationTest_vs_graphTest.csv"),
    row.names = FALSE
  )

  write.csv(
    data.frame(gene = overlap_startend_graph),
    file.path(out_dir_ts, "overlap_startVsEnd_vs_graphTest.csv"),
    row.names = FALSE
  )

  cat("[Overlap: associationTest ∩ graph_test]:", length(overlap_assoc_graph), "\n")
  cat("[Overlap: startVsEndTest ∩ graph_test]:", length(overlap_startend_graph), "\n")
}

# ---------------------------------------------------------
# 14.8  tradeSeq smoother 
#  associationTest  graph_test 
# ---------------------------------------------------------
cat("\n========== 11.7 plotting tradeSeq smoothers ==========\n")

candidate_genes <- assoc_sig$gene
if (exists("dynamic_genes_p1")) {
  overlap_genes <- intersect(candidate_genes, dynamic_genes_p1)
  candidate_genes <- unique(c(overlap_genes, assoc_sig$gene))
}

candidate_genes <- candidate_genes[candidate_genes %in% rownames(counts_ts)]

# 
manual_genes <- c("INS",  "TTR",   "LIMCH1", "PDE4D",
                  "PRUNE2",
                  "CPNE4")
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



# 
colData(sce_ts)$Stage_selected <- factor(
  colData(sce_ts)$Stage_selected,
  levels = c("CON", "PRE", "T2D")
)
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

  p <- rasterize_points(p, dpi = 400)  # 
  p
}

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
    file.path(out_dir_ts, "tradeSeq__genes2.pdf"),
    width = 10,
    height = max(8, ceiling(length(genes_to_plot) / 2) * 3.2),
    dpi = 300
  )

  
}

# =========================================================
# 15.  qs 
# =========================================================
cat("\n========== 12. Saving all objects into one .qs ==========" ,"\n")

all_object_names <- ls(envir = .GlobalEnv, all.names = TRUE)
all_results <- mget(all_object_names, envir = .GlobalEnv, inherits = FALSE)
qs_file <- file.path(out_dir, "union_native_all_objects.qs")

qsave(
  all_results,
  qs_file
)

cat("Saved:", qs_file, "\n")

#  qs  .GlobalEnv
loaded_results <- qread(qs_file)
if (!is.list(loaded_results)) {
  stop("Loaded object is not a list: ", qs_file)
}
list2env(loaded_results, envir = .GlobalEnv)
cat("Reloaded objects into .GlobalEnv:", length(loaded_results), "\n")
