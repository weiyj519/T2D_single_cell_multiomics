#!/usr/bin/env Rscript

# Script name: 03_plot_pseudotime_results.R
# Purpose: Cleaned public workflow script for the T2D single-cell multiome analysis.
# Input: Paths and analysis settings are read from configs/config.yaml.
# Output: Module-specific outputs are written under the configured results directory.
# Main steps: Load configuration, run the original analysis logic, and export reproducible outputs.
# Example command: Rscript scripts/08_pseudotime/03_plot_pseudotime_results.R

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
# Standalone beautified plotting from saved Lineage1 objects
# - Does not modify or rerun 1.R
# - Input: acinar_all_objects.qs
# - Output: boxed UMAPs + binned pseudotime heatmap for manual genes
# =========================================================

candidate_libs <- c(
  "python",
  "python"
)
candidate_libs <- candidate_libs[dir.exists(candidate_libs)]
.libPaths(unique(c(candidate_libs, .libPaths())))

suppressPackageStartupMessages({
  library(qs)
  library(Matrix)
  library(monocle3)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(pheatmap)
  library(SummarizedExperiment)
})

out_dir <- "results/downstream/riskcell_monocle/acinar"
qs_file <- file.path(out_dir, "acinar_all_objects.qs")
plot_dir <- file.path(out_dir, "beautified_from_qs")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(qs_file)) {
  stop("Saved qs object not found: ", qs_file)
}

cat("[Loading saved objects]\n", qs_file, "\n", sep = "")
loaded_results <- qs::qread(qs_file)
if (!is.list(loaded_results)) {
  stop("Loaded qs object is not a list.")
}
list2env(loaded_results, envir = .GlobalEnv)

manual_genes <- c("ASAP1", "PTPRG", "PDE5A", "CHGA", "PDE3B", "PDE4D")

if (!exists("stage_colors")) {
  stage_colors <- c(CON = "#4683B4", PRE = "#67B3AD", T2D = "#F3AF44")
}

save_plot_both <- function(p, basename, width = 6, height = 5, dpi = 300) {
  ggsave(
    filename = file.path(plot_dir, paste0(basename, ".pdf")),
    plot = p,
    width = width,
    height = height,
    dpi = dpi
  )
  ggsave(
    filename = file.path(plot_dir, paste0(basename, ".png")),
    plot = p,
    width = width,
    height = height,
    dpi = dpi
  )
}

rasterize_points2 <- function(p, dpi = 400) {
  if (!requireNamespace("ggrastr", quietly = TRUE)) {
    return(p)
  }
  ggrastr::rasterise(p, layers = "Point", dpi = dpi)
}

boxed_umap_theme <- theme(
  aspect.ratio = 1,
  panel.border = element_rect(color = "black", fill = NA, linewidth = 1.25),
  axis.line = element_blank(),
  axis.ticks = element_line(color = "black", linewidth = 0.45),
  axis.text = element_text(color = "black", size = 10),
  axis.title = element_text(color = "black", size = 12),
  panel.grid = element_blank(),
  panel.background = element_rect(fill = "white", color = NA),
  plot.background = element_rect(fill = "white", color = NA),
  legend.title = element_text(size = 11),
  legend.text = element_text(size = 10)
)

# =========================================================
# 1. Re-draw UMAPs with a thick full box
# =========================================================

if (!exists("cds")) {
  stop("Object `cds` was not found in the saved qs file.")
}

cat("[Drawing boxed UMAPs]\n")

colData(cds)$Stage_selected <- factor(
  colData(cds)$Stage_selected,
  levels = c("CON", "PRE", "T2D")
)
colData(cds)$group <- factor(
  colData(cds)$group,
  levels = c("CON", "PRE", "T2D")
)

p_stage <- plot_cells(
  cds,
  color_cells_by = "Stage_selected",
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = FALSE,
  graph_label_size = 3,
  cell_size = 0.5
) +
  scale_color_manual(values = stage_colors, name = "Stage_selected", drop = FALSE) +
  ggtitle("Stage-selected acinar cells") +
  boxed_umap_theme +
  theme(legend.position = "right")
p_stage <- rasterize_points2(p_stage)
save_plot_both(p_stage, "umap_by_stage_boxed", width = 7, height = 5.5)

p_group <- plot_cells(
  cds,
  color_cells_by = "group",
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = FALSE,
  graph_label_size = 3,
  cell_size = 0.5
) +
  scale_color_manual(values = stage_colors, name = "group", drop = FALSE) +
  ggtitle(NULL) +
  boxed_umap_theme +
  theme(legend.position = "right")
p_group <- rasterize_points2(p_group)
save_plot_both(p_group, "umap_by_group_boxed", width = 6, height = 5.5)

p_cluster <- plot_cells(
  cds,
  color_cells_by = "cluster",
  label_groups_by_cluster = TRUE,
  label_leaves = TRUE,
  label_branch_points = FALSE,
  cell_size = 0.5
) +
  ggtitle(NULL) +
  boxed_umap_theme +
  theme(legend.position = "right")
p_cluster <- rasterize_points2(p_cluster)
save_plot_both(p_cluster, "umap_by_cluster_boxed", width = 6, height = 5.5)

p_pt <- plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = FALSE,
  cell_size = 0.5
) +
  ggtitle(NULL) +
  scale_color_viridis_c(name = "Pseudotime") +
  boxed_umap_theme +
  theme(legend.position = "right")
p_pt <- rasterize_points2(p_pt)
save_plot_both(p_pt, "umap_by_pseudotime_boxed", width = 6, height = 5.5)
