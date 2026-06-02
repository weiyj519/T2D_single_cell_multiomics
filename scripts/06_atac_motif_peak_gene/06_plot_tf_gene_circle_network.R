#!/usr/bin/env Rscript

# Script name: 06_plot_tf_gene_circle_network.R
# Purpose: Cleaned public workflow script for the T2D single-cell multiome analysis.
# Input: Paths and analysis settings are read from configs/config.yaml.
# Output: Module-specific outputs are written under the configured results directory.
# Main steps: Load configuration, run the original analysis logic, and export reproducible outputs.
# Example command: Rscript scripts/06_atac_motif_peak_gene/06_plot_tf_gene_circle_network.R

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
  library(tidyverse)
  library(scales)
  library(ggrepel)
})

# ============================================================
# 0. 
# ============================================================

in_csv <- "results/downstream/beta_CON_PRE/netlink/TF_peak_gene_triplets_RNA_consistent_abs_rnaLFC0.25_padj0.05.csv"

outdir <- "results/downstream/beta_CON_PRE/netlink/circle_network"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# 
top_tf_n <- 10              #  TF
top_gene_per_tf <- 5      #  TF  target genes
min_strength_rna <- 0.02      # ， 0.01 / 0.02，
use_only_pass_sig <- TRUE  #  pass_sig == TRUE 
node_size_factor <- 1.5   # 
balance_tf_directions <- TRUE #  CON_up/PRE_up  TF
direction_radius <- 0.1      #  CON_up/PRE_up 
tf_radius <- 0.66            #  TF 
gene_radius <- 1.02           #  gene （ TF）
gene_label_radius <- 1.08     #  gene （ gene_radius）
direction_levels <- c("CON_up", "PRE_up", "Unknown")

# ============================================================
# 1. 
# ============================================================

raw <- read_csv(in_csv, show_col_types = FALSE)

need_cols <- c(
  "TF", "TF_gene", "peak", "gene",
  "TF_direction", "peak_direction", "gene_direction",
  "strength", "strength_rna",
  "rna_lfc", "rna_padj",
  "mean_abs_SHAP", "link_score", "link_fdr", "pass_sig"
)

miss <- setdiff(need_cols, colnames(raw))
if (length(miss) > 0) {
  warning("：", paste(miss, collapse = ", "))
}

dominant <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return("Unknown")
  names(sort(table(x), decreasing = TRUE))[1]
}

df <- raw %>%
  mutate(
    TF_node   = if ("TF_gene" %in% colnames(raw)) as.character(TF_gene) else as.character(TF),
    gene_node = as.character(gene),

    peak = as.character(peak),

    TF_direction   = as.character(TF_direction),
    peak_direction = as.character(peak_direction),
    gene_direction = as.character(gene_direction),

    strength = as.numeric(strength),
    strength_rna = as.numeric(strength_rna),
    mean_abs_SHAP = as.numeric(mean_abs_SHAP),
    link_score = as.numeric(link_score),
    link_fdr = as.numeric(link_fdr),
    rna_lfc = as.numeric(rna_lfc),
    rna_padj = as.numeric(rna_padj),

    pass_sig = as.logical(pass_sig)
  ) %>%
  filter(
    !is.na(TF_node), TF_node != "",
    !is.na(gene_node), gene_node != "",
    !is.na(peak), peak != ""
  )

if (use_only_pass_sig && "pass_sig" %in% colnames(df)) {
  df <- df %>% filter(pass_sig %in% TRUE)
}

df <- df %>%
  filter(!is.na(strength_rna)) %>%
  filter(strength_rna >= min_strength_rna)

message("[INFO] Triplets after filtering: ", nrow(df))

# ============================================================
# 2.  TF
#    ： gene peak  strength_rna
# ============================================================

tf_rank_all <- df %>%
  group_by(TF_node) %>%
  summarise(
    n_genes = n_distinct(gene_node),
    n_peaks = n_distinct(peak),
    total_strength = sum(abs(strength_rna), na.rm = TRUE),
    TF_direction = dominant(TF_direction),
    chromvar_lfc = suppressWarnings(first(na.omit(chromvar_lfc))),
    .groups = "drop"
  ) %>%
  arrange(desc(n_genes), desc(total_strength), desc(n_peaks))

if (balance_tf_directions) {
  n_tf_direction <- max(n_distinct(tf_rank_all$TF_direction), 1)
  top_tf_n_per_direction <- ceiling(top_tf_n / n_tf_direction)

  tf_rank <- tf_rank_all %>%
    group_by(TF_direction) %>%
    slice_head(n = top_tf_n_per_direction) %>%
    ungroup() %>%
    mutate(
      TF_direction = factor(TF_direction, levels = direction_levels)
    ) %>%
    arrange(TF_direction, desc(n_genes), desc(total_strength), desc(n_peaks)) %>%
    slice_head(n = top_tf_n) %>%
    mutate(TF_direction = as.character(TF_direction))
} else {
  tf_rank <- tf_rank_all %>%
    slice_head(n = top_tf_n)
}

tf_keep <- tf_rank$TF_node

message("[INFO] Selected TFs: ", paste(tf_keep, collapse = ", "))

# ============================================================
# 2.1 TF （ TF ）
# ============================================================

tf_colors <- setNames(scales::hue_pal()(length(tf_keep)), tf_keep)
direction_colors <- c(CON_up = "#0072B2", PRE_up = "#D55E00", Unknown = "grey70")
node_colors <- c(tf_colors, direction_colors)

# ============================================================
# 3.  TF-peak-gene  TF-gene 
#     peak  TF-gene 
# ============================================================

edges_all <- df %>%
  filter(TF_node %in% tf_keep) %>%
  group_by(TF_node, gene_node) %>%
  summarise(
    n_peaks = n_distinct(peak),
    weight = sum(abs(strength_rna), na.rm = TRUE),
    max_strength = max(abs(strength_rna), na.rm = TRUE),
    mean_SHAP = mean(mean_abs_SHAP, na.rm = TRUE),

    TF_direction = dominant(TF_direction),
    gene_direction = dominant(gene_direction),
    peak_direction = dominant(peak_direction),

    rna_lfc = suppressWarnings(first(na.omit(rna_lfc))),
    rna_padj = suppressWarnings(first(na.omit(rna_padj))),
    .groups = "drop"
  ) %>%
  arrange(TF_node, desc(weight), desc(n_peaks))

#  TF  top target genes
edges <- edges_all %>%
  group_by(TF_node) %>%
  slice_max(order_by = weight, n = top_gene_per_tf, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    #  gene  TF 
    gene_id = paste(TF_node, gene_node, sep = "::"),
    direction_id = TF_direction
  )

message("[INFO] Edges used for plot: ", nrow(edges))
message("[INFO] Gene nodes used for plot (TF-gene pairs): ", nrow(edges))
message("[INFO] Unique genes (for reference): ", n_distinct(edges$gene_node))

write_csv(edges, file.path(outdir, "circle_network_edges_TF_gene.csv"))
write_csv(tf_rank, file.path(outdir, "circle_network_selected_TFs.csv"))

# ============================================================
# 4. 
# ============================================================

tf_nodes <- edges %>%
  group_by(TF_node) %>%
  summarise(
    name = first(TF_node),
    type = "TF",
    direction = dominant(TF_direction),
    n_edges = n(),
    n_genes = n_distinct(gene_node),
    n_peaks = sum(n_peaks, na.rm = TRUE),
    node_weight = sum(weight, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(group = name) %>%
  select(name, type, group, direction, n_edges, n_genes, n_peaks, node_weight)

gene_nodes <- edges %>%
  transmute(
    #  TF-gene  gene （ TF ）
    gene_id = gene_id,
    name = gene_node,
    type = "gene",
    group = TF_node,
    direction = gene_direction,
    direction_id = direction_id,
    n_edges = 1,
    n_genes = 1,
    n_peaks = n_peaks,
    node_weight = weight,
    parent_TF = TF_node
  )

# ============================================================
# 5.  circle layout
#    genes ，TF ，CON_up/PRE_up 
# ============================================================

# -----  gene  -----
gene_nodes <- gene_nodes %>%
  mutate(
    parent_TF = factor(parent_TF, levels = tf_keep)
  ) %>%
  arrange(parent_TF, factor(direction, levels = direction_levels), desc(node_weight))

gene_nodes <- gene_nodes %>%
  mutate(group = factor(group, levels = tf_keep))

n_gene <- nrow(gene_nodes)

theta_gene <- seq(pi / 2, pi / 2 - 2 * pi, length.out = n_gene + 1)[-(n_gene + 1)]

gene_nodes <- gene_nodes %>%
  mutate(
    theta = theta_gene,
    x = gene_radius * cos(theta),
    y = gene_radius * sin(theta),

    label_x = gene_label_radius * cos(theta),
    label_y = gene_label_radius * sin(theta),

    angle_raw = theta * 180 / pi,
    angle = if_else(label_x < 0, angle_raw + 180, angle_raw),
    hjust = if_else(label_x < 0, 1, 0),

    plot_size = rescale(log1p(node_weight), to = c(2.2, 5.2) * node_size_factor)
  )

# -----  TF  -----
tf_nodes <- tf_nodes %>%
  mutate(
    name = factor(name, levels = tf_keep)
  ) %>%
  arrange(name)

tf_nodes <- tf_nodes %>%
  mutate(group = factor(as.character(group), levels = tf_keep))

n_tf <- nrow(tf_nodes)

theta_tf <- seq(pi / 2, pi / 2 - 2 * pi, length.out = n_tf + 1)[-(n_tf + 1)]

tf_nodes <- tf_nodes %>%
  mutate(
    theta = theta_tf,
    x = tf_radius * cos(theta),
    y = tf_radius * sin(theta),

    label_x = x,
    label_y = y,

    angle = 0,
    hjust = 0.5,

    direction_id = direction,
    plot_size = rescale(log1p(node_weight), to = c(10, 28) * node_size_factor)
  )

# -----  CON_up/PRE_up  -----
direction_nodes_base <- tf_nodes %>%
  group_by(direction) %>%
  summarise(
    name = first(direction),
    type = "direction",
    group = first(direction),
    direction = first(direction),
    direction_id = first(direction),
    n_edges = sum(n_edges, na.rm = TRUE),
    n_genes = sum(n_genes, na.rm = TRUE),
    n_peaks = sum(n_peaks, na.rm = TRUE),
    node_weight = sum(node_weight, na.rm = TRUE),
    .groups = "drop"
  )

direction_nodes <- direction_nodes_base %>%
  left_join(
    tf_nodes %>%
      group_by(direction) %>%
      summarise(
        theta = atan2(mean(sin(theta)), mean(cos(theta))),
        .groups = "drop"
      ),
    by = "direction"
  ) %>%
  arrange(factor(direction, levels = direction_levels)) %>%
  mutate(
    x = direction_radius * cos(theta),
    y = direction_radius * sin(theta),

    label_x = x,
    label_y = y,

    angle = 0,
    hjust = 0.5,

    plot_size = rescale(log1p(node_weight), to = c(7.0, 11.0) * node_size_factor)
  )

nodes <- bind_rows(
  gene_nodes %>% select(gene_id, direction_id, name, type, group, direction, x, y, label_x, label_y, angle, hjust, plot_size, node_weight, parent_TF),
  direction_nodes %>% select(direction_id, name, type, group, direction, x, y, label_x, label_y, angle, hjust, plot_size, node_weight),
  tf_nodes %>% select(direction_id, name, type, group, direction, x, y, label_x, label_y, angle, hjust, plot_size, node_weight)
)

write_csv(nodes, file.path(outdir, "circle_network_nodes.csv"))

# ============================================================
# 6. 
# ============================================================

edge_plot <- edges %>%
  left_join(
    tf_nodes %>% select(TF_node = name, x_tf = x, y_tf = y),
    by = "TF_node"
  ) %>%
  left_join(
    direction_nodes %>% select(direction_id, x_direction = x, y_direction = y),
    by = "direction_id"
  ) %>%
  left_join(
    gene_nodes %>% select(gene_id, x_gene = x, y_gene = y),
    by = "gene_id"
  ) %>%
  mutate(
    edge_width = rescale(log1p(weight), to = c(0.15, 0.85)),
    edge_alpha = rescale(log1p(weight), to = c(0.20, 0.55))
  )

# ============================================================
# 8. 
# ============================================================

p <- ggplot() +

  # CON_up/PRE_up -> TF
  geom_curve(
    data = edge_plot,
    aes(
      x = x_direction, y = y_direction,
      xend = x_tf, yend = y_tf,
      linewidth = edge_width,
      alpha = edge_alpha
    ),
    curvature = 0.08,
    color = "grey50",
    lineend = "round"
  ) +

  # TF -> gene
  geom_curve(
    data = edge_plot,
    aes(
      x = x_tf, y = y_tf,
      xend = x_gene, yend = y_gene,
      linewidth = edge_width,
      alpha = edge_alpha
    ),
    curvature = 0.08,
    color = "grey50",
    lineend = "round"
  ) +

  scale_linewidth_identity() +
  scale_alpha_identity() +

  #  gene 
  geom_point(
    data = gene_nodes,
    aes(x = x, y = y, size = plot_size, fill = group),
    shape = 21,
    color = "white",
    stroke = 0.25,
    alpha = 0.95
  ) +

  #  TF 
  geom_point(
    data = tf_nodes,
    aes(x = x, y = y, size = plot_size, fill = group),
    shape = 21,
    color = "white",
    stroke = 0.45,
    alpha = 0.75
  ) +

  #  CON_up/PRE_up 
  geom_point(
    data = direction_nodes,
    aes(x = x, y = y, size = plot_size, fill = direction),
    shape = 21,
    color = "white",
    stroke = 0.35,
    alpha = 0.95
  ) +

  scale_size_identity() +

  # TF 
  geom_text(
    data = tf_nodes,
    aes(x = x, y = y, label = name),
    size = 4.3,
    fontface = "bold",
    color = "black"
  ) +

  # CON_up/PRE_up 
  geom_text(
    data = direction_nodes,
    aes(x = label_x, y = label_y, label = name),
    size = 3.5,
    fontface = "bold",
    color = "black"
  ) +

  # gene 
  geom_text(
    data = gene_nodes,
    aes(
      x = label_x,
      y = label_y,
      label = name,
      angle = angle,
      hjust = hjust
    ),
    size = 4,
    color = "black",
    fontface = "bold"
  ) +

  scale_fill_manual(values = node_colors, name = "TF", breaks = tf_keep, drop = FALSE, guide = "none") +

  coord_equal(
    xlim = c(-1.35, 1.35),
    ylim = c(-1.35, 1.35),
    clip = "off"
  ) +

  theme_void(base_size = 12) +
  theme(
    legend.position = "none",
    plot.margin = margin(20, 40, 20, 40)
  )

print(p)

ggsave(
  file.path(outdir, "beta_CON_PRE_TF_gene_circle_network.png"),
  p,
  width = 8.5,
  height = 8.5,
  dpi = 500,
  bg = "white"
)

ggsave(
  file.path(outdir, "beta_CON_PRE_TF_gene_circle_network.pdf"),
  p,
  width = 8.5,
  height = 8.5,
  bg = "white"
)

message("[INFO] Done.")
message("[INFO] Output dir: ", outdir)
