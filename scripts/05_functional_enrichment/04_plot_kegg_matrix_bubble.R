#!/usr/bin/env Rscript

# Script name: 04_plot_kegg_matrix_bubble.R
# Purpose: Cleaned public workflow script for the T2D single-cell multiome analysis.
# Input: Paths and analysis settings are read from configs/config.yaml.
# Output: Module-specific outputs are written under the configured results directory.
# Main steps: Load configuration, run the original analysis logic, and export reproducible outputs.
# Example command: Rscript scripts/05_functional_enrichment/04_plot_kegg_matrix_bubble.R

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

# ==============================================================================
# Acinar KEGG matrix bubble plot
# Compare: CON vs PRE and PRE vs T2D
# Run with:
#   python run -n Renv Rscript acinar_KEGG_matrix_bubble_plot.R
# ==============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))

required_pkgs <- c("ggplot2", "dplyr", "forcats", "viridis", "tibble", "scales")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(forcats)
  library(viridis)
})

base_dir <- "results/downstream/acinar_CON_PRE"

kegg_files <- tibble::tibble(
  Comparison = factor(
    c("Acinar\nCON vs PRE", "Acinar\nPRE vs T2D"),
    levels = c("Acinar\nCON vs PRE", "Acinar\nPRE vs T2D")
  ),
  File = c(
    "results/downstream/acinar_CON_PRE/GO_RNA/KEGG.csv",
    "results/downstream/acinar_PRE_T2D/GO_RNA/KEGG.csv"
  )
)

# （ Comparison ）
target_pathways <- list(
  "Acinar\nCON vs PRE" = c(
    "Longevity regulating pathway",
    "Thyroid hormone synthesis",
    "mTOR signaling pathway",
    "FoxO signaling pathway",
    "Insulin secretion"
  ),
  "Acinar\nPRE vs T2D" = c(
    "mTOR signaling pathway",
    "FoxO signaling pathway",
    "IL-17 signaling pathway",
    "Th17 cell differentiation",
    "Relaxin signaling pathway"
  )
)

parse_gene_ratio <- function(x) {
  vapply(strsplit(as.character(x), "/", fixed = TRUE), function(z) {
    if (length(z) != 2 || any(is.na(suppressWarnings(as.numeric(z))))) {
      return(NA_real_)
    }
    as.numeric(z[1]) / as.numeric(z[2])
  }, numeric(1))
}

read_kegg <- function(file, comparison) {
  read.csv(file, stringsAsFactors = FALSE, check.names = FALSE) %>%
    mutate(
      Comparison = comparison,
      GeneRatioNumeric = parse_gene_ratio(GeneRatio),
      pvalue = as.numeric(pvalue),
      log10_pvalue = -log10(pmax(pvalue, .Machine$double.xmin))
    )
}

parsed_kegg <- bind_rows(
  Map(read_kegg, kegg_files$File, as.character(kegg_files$Comparison))
)

# “” Comparison x Pathway ； CSV ，（）
requested_pairs <- bind_rows(lapply(as.character(kegg_files$Comparison), function(comp) {
  tibble::tibble(
    Comparison = comp,
    Description = target_pathways[[comp]]
  )
}))

plot_data <- requested_pairs %>%
  left_join(parsed_kegg, by = c("Comparison", "Description")) %>%
  mutate(
    Comparison = factor(
      Comparison,
      levels = c("Acinar\nCON vs PRE", "Acinar\nPRE vs T2D")
    )
  )

# （ KEGG.csv ）
missing_tbl <- plot_data %>%
  filter(is.na(ID)) %>%
  distinct(Comparison, Description)
if (nrow(missing_tbl) > 0) {
  message("[WARN]  KEGG.csv ，：")
  apply(missing_tbl, 1, function(x) message("  - ", x[[1]], ": ", x[[2]]))
}

pathway_order <- unique(c(
  target_pathways[["Acinar\nCON vs PRE"]],
  target_pathways[["Acinar\nPRE vs T2D"]]
))

plot_data <- plot_data %>%
  mutate(Description = factor(Description, levels = rev(pathway_order)))

wrap_labels <- function(x, width = 36) {
  vapply(x, function(label) paste(strwrap(label, width = width), collapse = "\n"), character(1))
}

write.csv(
  plot_data %>%
    arrange(Comparison, pvalue) %>%
    select(
      Comparison, category, subcategory, ID, Description, GeneRatio,
      GeneRatioNumeric, pvalue, p.adjust, qvalue, Count, geneID
    ),
  file = file.path(base_dir, "acinar_KEGG_selected_for_matrix_bubble.csv"),
  row.names = FALSE
)

p <- ggplot(plot_data, aes(x = Comparison, y = Description)) +
  geom_point(aes(size = GeneRatioNumeric, color = log10_pvalue), alpha = 0.92, na.rm = TRUE) +
  scale_color_gradientn(
    colors = c("#377EB8", "#ABDDA4", "#F46D43", "#A50026"),
    name = "-log10(P.value)"
  ) +
  scale_size_continuous(
    range = c(3.2, 9),
    name = "Gene Ratio",
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_y_discrete(labels = wrap_labels) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    text = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold", size = 14, color = "black"),
    axis.text.y = element_text(face = "bold", size = 16, color = "black", lineheight = 0.92),
    axis.ticks = element_line(linewidth = 1.1, color = "black"),
    axis.ticks.length = unit(0.18, "cm"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.4),
    panel.grid.major = element_line(color = "grey90", linetype = "dashed", linewidth = 0.7),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(face = "bold", size = 12),
    legend.key.height = unit(0.75, "cm"),
    plot.margin = margin(10, 16, 10, 10)
  )

pdf_file <- file.path(base_dir, "acinar_KEGG_top10_matrix_bubble_plot.pdf")
png_file <- file.path(base_dir, "acinar_KEGG_top10_matrix_bubble_plot.png")

ggsave(pdf_file, p, width = 8.2, height = 6, device = cairo_pdf)
ggsave(png_file, p, width = 9.2, height = 8.4, dpi = 300)

print(p)
message("Saved: ", pdf_file)
message("Saved: ", png_file)
