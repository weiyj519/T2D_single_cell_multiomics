#!/usr/bin/env Rscript

# Script name: 01_plot_atac_stage_feature_upset.R
# Purpose: Cleaned public workflow script for the T2D single-cell multiome analysis.
# Input: Paths and analysis settings are read from configs/config.yaml.
# Output: Module-specific outputs are written under the configured results directory.
# Main steps: Load configuration, run the original analysis logic, and export reproducible outputs.
# Example command: Rscript scripts/03_model_evaluation/01_plot_atac_stage_feature_upset.R

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
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(UpSetR)
})

# ========== 1) ：ATAC （ *_shap_marker_genes.csv） ==========
base_dir <- "results/model_outputs/atac"
files <- list.files(
  base_dir,
  pattern = "_shap_marker_genes\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)

#  celltype （ 3 ）
cell_order <- c("α_cell","β_cell","δ_cell","Acinar_cell","Ductal_cell","Immune_cell","Stellate_cell","PP_cell")

# Set size（）： celltype 
cell_colors <- c(
  "α_cell" = "#cdb4d7",
  "β_cell" = "#f68b6a",
  "δ_cell" = "#50b688",
  "Acinar_cell" = "#5b2d90",
  "Ductal_cell" = "#ffda67",
  "Immune_cell" = "#b7c47d",
  "Stellate_cell" = "#ce6091",
  "PP_cell" = "#99bae2"
)

# ========== 2) ： +  cell/task ==========
read_gene_list <- function(f) {
  df <- suppressMessages(readr::read_csv(f, show_col_types = FALSE))
  #  gene ；（）
  gene_col <- intersect(names(df), c("gene","Gene","GENE","feature","Feature","peak","Peak","PEAK"))[1]
  if (is.na(gene_col)) gene_col <- names(df)[1]
  genes <- df[[gene_col]] %>%
    as.character() %>%
    str_trim() %>%
    .[. != "" & !is.na(.)] %>%
    unique()
  genes
}

parse_meta <- function(f) {
  b <- basename(f)
  # : /.../β_cell/CON_vs_T2D_shap_marker_genes.csv
  cell <- basename(dirname(f))
  task <- str_replace(b, "_shap_marker_genes\\.csv$", "")
  task <- str_replace(task, "\\.csv$", "")
  tibble(
    file = f,
    cell = cell,
    task = task
  )
}

allowed_tasks <- c("CON_vs_PRE", "PRE_vs_T2D", "CON_vs_T2D")
meta <- map_dfr(files, parse_meta) %>%
  filter(cell %in% cell_order, task %in% allowed_tasks) %>%
  mutate(cell = factor(cell, levels = cell_order))

if (nrow(meta) == 0) {
  stop("No input files found. Check base_dir and filename pattern: *_shap_marker_genes.csv")
}

# ========== 3)  task  UpSet ==========
plot_one_task <- function(task_name, out_prefix,
                          nintersects = 30,
                          point_size = 3.2,
                          line_size  = 0.8,
                          text_scale = 2.1,
                          title_size = 18) {

  plot_title <- str_replace(task_name, "_vs_", " vs ")

  sub <- meta %>% filter(task == task_name) %>% arrange(cell)

  gene_sets <- setNames(
    map(sub$file, read_gene_list),
    as.character(sub$cell)
  )

  set_sizes <- lengths(gene_sets)
  sets_order <- names(sort(set_sizes, decreasing = TRUE))
  sets_plot <- rev(sets_order)
  sets_bar_colors <- unname(cell_colors[sets_plot])
  sets_bar_colors[is.na(sets_bar_colors)] <- "grey70"

  # UpSetR  data.frame（），fromList 
  upset_df <- UpSetR::fromList(gene_sets)

  #  set size（）
  set_size_tbl <- tibble(
    cell = names(set_sizes),
    set_size = as.integer(set_sizes)
  ) %>%
    arrange(desc(set_size))
  readr::write_csv(set_size_tbl, paste0(out_prefix, "_set_size.csv"))

  # （ + ）
  membership_tbl <- upset_df %>%
    tibble::rownames_to_column("feature")

  intersection_tbl <- membership_tbl %>%
    group_by(across(all_of(sets_plot))) %>%
    summarise(
      intersection_size = dplyr::n(),
      features = paste(sort(feature), collapse = ";"),
      .groups = "drop"
    ) %>%
    mutate(n_sets = rowSums(across(all_of(sets_plot)))) %>%
    filter(n_sets > 0) %>%
    arrange(desc(intersection_size), desc(n_sets))

  readr::write_csv(intersection_tbl, paste0(out_prefix, "_intersection_all.csv"))
  readr::write_csv(
    dplyr::slice_head(intersection_tbl, n = nintersects),
    paste0(out_prefix, "_intersection_top", nintersects, ".csv")
  )

  cairo_pdf(
    paste0(out_prefix, ".pdf"),
    width = 9.5,
    height = 7,
    family = "Arial"
  )
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(family = "Arial", font = 2, cex = 1.2)
  print(
    UpSetR::upset(
      upset_df,
      sets = sets_plot,
      keep.order = TRUE,
      order.by = "freq",
      nintersects = nintersects,
      point.size = point_size,
      line.size = line_size,
      text.scale = text_scale,
      sets.bar.color = sets_bar_colors,
      mainbar.y.label = "Intersection size",
      sets.x.label = "Set size"
    )
  )
  grid::grid.text(
    plot_title,
    x = 0.5,
    y = grid::unit(0.98, "npc"),
    gp = grid::gpar(fontface = "bold", fontsize = title_size, fontfamily = "Arial")
  )
  dev.off()

  png(
    paste0(out_prefix, ".png"),
    width = 2200,
    height = 1600,
    res = 300,
    type = "cairo",
    family = "Arial"
  )
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(family = "Arial", font = 2, cex = 1.2)
  print(
    UpSetR::upset(
      upset_df,
      sets = sets_plot,
      keep.order = TRUE,
      order.by = "freq",
      nintersects = nintersects,
      point.size = point_size,
      line.size = line_size,
      text.scale = text_scale,
      sets.bar.color = sets_bar_colors,
      mainbar.y.label = "Intersection size",
      sets.x.label = "Set size"
    )
  )
  grid::grid.text(
    plot_title,
    x = 0.5,
    y = grid::unit(0.98, "npc"),
    gp = grid::gpar(fontface = "bold", fontsize = title_size, fontfamily = "Arial")
  )
  dev.off()

  message(
    "[OK] ", task_name,
    " -> ", out_prefix,
    ".pdf/.png + _set_size.csv + _intersection_all.csv + _intersection_top",
    nintersects,
    ".csv"
  )
}

dir.create("results/downstream/upset/ATAC_upset_out", showWarnings = FALSE)

plot_one_task("CON_vs_PRE",  "results/downstream/upset/ATAC_upset_out/UpSet_CON_vs_PRE_annotated")
plot_one_task("PRE_vs_T2D",  "results/downstream/upset/ATAC_upset_out/UpSet_PRE_vs_T2D_annotated")
plot_one_task("CON_vs_T2D",  "results/downstream/upset/ATAC_upset_out/UpSet_CON_vs_T2D_annotated")




