#!/usr/bin/env Rscript

# Script name: 03_plot_functional_enrichment.R
# Purpose: Cleaned public workflow script for the T2D single-cell multiome analysis.
# Input: Paths and analysis settings are read from configs/config.yaml.
# Output: Module-specific outputs are written under the configured results directory.
# Main steps: Load configuration, run the original analysis logic, and export reproducible outputs.
# Example command: Rscript scripts/05_functional_enrichment/03_plot_functional_enrichment.R

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
  library(ggplot2)
  library(tidytext)   # reorder_within / scale_y_reordered
  library(grid)       # unit
})

#====================#
# 1. 
#====================#
indir <- "results/downstream/beta_CON_PRE/GO_RNA"
outdir <- "results/downstream/beta_CON_PRE/GO_RNA/plot"
top_n        <- 5      # 10
p_cutoff     <- 0.05     #  pvalue 
wrap_width   <- 55       # term 
include_kegg <- TRUE    # FALSE = GO（2）
                         # TRUE  = GO + KEGG 
#====================#
# 2. 
#====================#
read_res <- function(path, onto) {
  if (!file.exists(path)) {
    message("[WARN] : ", path)
    return(NULL)
  }
  df <- read_csv(path, show_col_types = FALSE)
  if (nrow(df) == 0) {
    message("[WARN] : ", path)
    return(NULL)
  }
  df$ONTOLOGY <- onto
  df
}

GO_BP <- read_res(file.path(indir, "GO_BP.csv"), "BP")
GO_CC <- read_res(file.path(indir, "GO_CC.csv"), "CC")
GO_MF <- read_res(file.path(indir, "GO_MF.csv"), "MF")
KEGG  <- read_res(file.path(indir, "KEGG.csv"),  "KEGG")

combined_df <- bind_rows(GO_BP, GO_CC, GO_MF, KEGG)

if (nrow(combined_df) == 0) {
  stop("")
}

#====================#
# 3. 
#====================#
plot_df <- combined_df %>%
  filter(!is.na(pvalue), pvalue > 0, pvalue < p_cutoff) %>%
  mutate(
    Category = recode(
      ONTOLOGY,
      BP   = "Biological Process",
      CC   = "Cellular Component",
      MF   = "Molecular Function",
      KEGG = "KEGG Pathway"
    ),
    Category = factor(
      Category,
      levels = c(
        "Biological Process",
        "Cellular Component",
        "Molecular Function",
        "KEGG Pathway"
      )
    )
  ) %>%
  group_by(Category) %>%
  arrange(pvalue, .by_group = TRUE) %>%
  slice_head(n = top_n) %>%
  ungroup() %>%
  mutate(
    logP = -log10(pvalue),
    Description_wrap = str_wrap(Description, width = wrap_width)
  )

if (nrow(plot_df) == 0) {
  stop(" pvalue ， p_cutoff")
}

#====================#
# 4. 
#====================#
pal <- c(
  "Biological Process" = "#E77C73",
  "Cellular Component" = "#67B84B",
  "Molecular Function" = "#7693E6",
  "KEGG Pathway"       = "#A97AE8"
)

#====================#
# 5. （）
#====================#
p_bar <- ggplot(
  plot_df,
  aes(
    x = logP,
    y = reorder_within(Description_wrap, logP, Category),
    fill = Category
  )
) +
  geom_col(width = 0.72, color = NA) +
  facet_grid(
    Category ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_y_reordered() +
  scale_fill_manual(values = pal) +
  labs(
    x = "-log10(p-value)",
    y = NULL,
    fill = "Category"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),

    strip.background = element_rect(fill = "#F0F0F0", color = "black", linewidth = 0.7),
    strip.text.y = element_text(size = 11, face = "bold", angle = 270),

    axis.text.y = element_text(color = "black", size = 12, lineheight = 0.92),
    axis.text.x = element_text(color = "black", size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),

    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),

    panel.spacing.y = unit(0.25, "lines"),
    plot.margin = margin(8, 14, 8, 8)
  )

p_bar

ggsave(
  file.path(outdir, "GO_KEGG_barplot_facet_pvalue_top5.pdf"),
  p_bar, width = 7.2, height = 7.5, dpi = 600
)

#====================#
# 6. 
#====================#

#  GeneRatio  "a/b" 
parse_ratio <- function(x) {
  sapply(strsplit(as.character(x), "/"), function(xx) {
    if (length(xx) == 2) {
      as.numeric(xx[1]) / as.numeric(xx[2])
    } else {
      NA_real_
    }
  })
}

plot_df_dot <- plot_df %>%
  mutate(
    GeneRatio_num = parse_ratio(GeneRatio),
    Count = as.numeric(Count)
  ) %>%
  filter(!is.na(GeneRatio_num), !is.na(Count))

if (nrow(plot_df_dot) == 0) {
  stop("， GeneRatio  Count ")
}

p_dot <- ggplot(
  plot_df_dot,
  aes(
    x = GeneRatio_num,
    y = reorder_within(Description_wrap, logP, Category),
    size = Count,
    color = logP
  )
) +
  geom_point(alpha = 0.9) +
  facet_grid(
    Category ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_y_reordered() +
  scale_color_gradientn(
    colors = c("#FFB9B9", "#FF5A5A", "#D60000"),
    name = "-log10(p-value)"
  ) +
  scale_size_continuous(
    name = "Gene Count",
    range = c(3.5, 8)
  ) +
  labs(
    x = "Gene Ratio",
    y = NULL,
    color = "-log10(p-value)",
    size = "Gene Count"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey88", linewidth = 0.35),

    strip.background = element_rect(fill = "#F0F0F0", color = "black", linewidth = 0.7),
    strip.text.y = element_text(size = 11, face = "bold", angle = 270),

    axis.text.y = element_text(color = "black", size = 12, lineheight = 0.92),
    axis.text.x = element_text(color = "black", size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),

    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),

    panel.spacing.y = unit(0.25, "lines"),
    plot.margin = margin(8, 14, 8, 8)
  )

p_dot

ggsave(
  file.path(outdir, "GO_KEGG_dotplot_facet_pvalue_top5.pdf"),
  p_dot, width = 7.2, height = 7.5, dpi = 600
)

