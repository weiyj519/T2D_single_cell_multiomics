#!/usr/bin/env Rscript

# Script name: 04_plot_tf_peak_gene_network.R
# Purpose: Cleaned public workflow script for the T2D single-cell multiome analysis.
# Input: Paths and analysis settings are read from configs/config.yaml.
# Output: Module-specific outputs are written under the configured results directory.
# Main steps: Load configuration, run the original analysis logic, and export reproducible outputs.
# Example command: Rscript scripts/06_atac_motif_peak_gene/04_plot_tf_peak_gene_network.R

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
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(scales)
})

# ：（ netlink ：Rscript plot.R）
outdir <- "results/downstream/acinar_CON_PRE/netlink_mlpeak"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ：
# -  TF->gene （ TF_peak_gene_triplets_RNA_consistent.csv ）
# - “ TF +  gene”：
#   - TF ： TF_direction （PRE_up / CON_up）
#   - gene ： gene （）
#   - ： gene （edge_color_by_gene=TRUE）


top_n_sankey_each_dir <- 20

# =========================
# “”（/，）
# =========================
# 1) （CON_up / PRE_up / T2D_up） triplets（ strength_rna  N）
#    -  Inf 
max_triplets_per_dir <- 30

# 2)  TF  gene（ TF-gene importance  N）
#    -  Inf 
max_genes_per_tf_in_plot <- 50

circ_mean <- function(theta) atan2(mean(sin(theta)), mean(cos(theta)))

plot_tf_gene_circle_unique <- function(df,
                                       tf_col = "TF",
                                       gene_col = "target",
                                       weight_col = "importance",
                                       tf_dir_df = NULL,
                                       tf_dir_cols = c(PRE_up = "#7EE2FA", CON_up = "#FF9999", T2D_up = "#FF9999"),
                                       top_n_per_tf = 50,
                                       weight_transform = function(x) log1p(x),
                                       tf_cols = NULL,
                                       r_tf = 0.45,
                                       r_gene = 1.1,
                                       label_expand = 1.08,
                                       edge_alpha = 0.35,
                                       edge_color = "grey50",
                                       edge_color_by_tf = FALSE,
                                       edge_color_by_gene = TRUE,
                                       edge_width_range = c(0.8, 1.6),
                                       tf_size_range = c(8, 26),
                                       gene_size_range = c(5, 9),
                                       gene_label_size = 5.5,
                                       tf_label_size = 8,
                                       max_gene_labels = Inf,
                                       x_expand = 1.25,
                                       spread_tf_angles = TRUE) {

  dat <- df %>%
    transmute(
      TF   = as.character(.data[[tf_col]]),
      gene = as.character(.data[[gene_col]]),
      w    = as.numeric(.data[[weight_col]])
    ) %>%
    filter(!is.na(TF), !is.na(gene), is.finite(w), w > 0) %>%
    group_by(TF, gene) %>%
    summarise(w = sum(w), .groups = "drop") %>%
    mutate(wp = weight_transform(w)) %>%
    filter(is.finite(wp), wp > 0)

  if (is.finite(top_n_per_tf) && top_n_per_tf > 0) {
    dat <- dat %>%
      group_by(TF) %>%
      slice_max(order_by = wp, n = top_n_per_tf, with_ties = FALSE) %>%
      ungroup()
  }

  tf_order <- dat %>% group_by(TF) %>% summarise(tf_w = sum(wp), .groups="drop") %>%
    arrange(desc(tf_w)) %>% pull(TF)

  dat <- dat %>% mutate(TF = factor(TF, levels = tf_order))

  gene_primary <- dat %>%
    group_by(gene) %>%
    slice_max(order_by = wp, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(gene, TF_primary = as.character(TF))

  gene_order <- dat %>%
    left_join(gene_primary, by="gene") %>%
    arrange(factor(TF_primary, levels=tf_order), desc(wp), gene) %>%
    pull(gene) %>%
    unique()

  nG <- length(gene_order)
  gene_theta <- tibble(
    gene = gene_order,
    theta = seq(0, 2*pi, length.out = nG + 1)[1:nG]
  ) %>%
    mutate(
      x = r_gene * cos(theta),
      y = r_gene * sin(theta)
    ) %>%
    left_join(gene_primary, by="gene")

  gene_cols <- setNames(scales::hue_pal(l = 65, c = 80)(nG), gene_order)

  tf_theta <- dat %>%
    left_join(gene_theta %>% select(gene, theta), by="gene") %>%
    group_by(TF) %>%
    summarise(theta = circ_mean(theta), tf_w = sum(wp), .groups="drop") %>%
    mutate(theta_raw = theta)

  if (spread_tf_angles) {
    nTF <- nrow(tf_theta)
    tf_theta <- tf_theta %>%
      arrange(theta_raw) %>%
      mutate(theta = (row_number() - 1) / nTF * 2 * pi)
  } else {
    tf_theta <- tf_theta %>% mutate(theta = theta_raw)
  }

  tf_theta <- tf_theta %>%
    mutate(
      x = r_tf * cos(theta),
      y = r_tf * sin(theta)
    )

  if (!is.null(tf_dir_df)) {
    tf_dir_df2 <- tf_dir_df %>%
      transmute(TF = as.character(.data[[tf_col]]), TF_direction = as.character(TF_direction)) %>%
      distinct(TF, TF_direction)
    tf_theta <- tf_theta %>%
      left_join(tf_dir_df2, by = "TF")
  } else {
    tf_theta <- tf_theta %>% mutate(TF_direction = NA_character_)
  }

  tf_theta <- tf_theta %>%
    mutate(col_dir = dplyr::case_when(
      !is.na(TF_direction) & TF_direction %in% names(tf_dir_cols) ~ unname(tf_dir_cols[TF_direction]),
      TRUE ~ "#BDBDBD"
    ))

  tfs <- tf_order
  if (is.null(tf_cols)) {
    n_base <- min(8, max(3, length(tfs)))
    base_pal <- if (requireNamespace("RColorBrewer", quietly = TRUE)) {
      RColorBrewer::brewer.pal(n_base, "Set2")
    } else {
      scales::hue_pal(l = 60, c = 70)(n_base)
    }
    tf_cols <- if (length(tfs) <= length(base_pal)) {
      setNames(base_pal[seq_along(tfs)], tfs)
    } else {
      setNames(hue_pal(l = 60, c = 70)(length(tfs)), tfs)
    }
  } else {
    miss <- setdiff(tfs, names(tf_cols))
    if (length(miss) > 0) {
      tf_cols <- c(tf_cols, setNames(hue_pal(l = 60, c = 70)(length(miss)), miss))
    }
    tf_cols <- tf_cols[tfs]
  }

  tf_theta <- tf_theta %>%
    mutate(col = tf_cols[as.character(TF)],
           size = rescale(tf_w, to=tf_size_range))

  gene_nodes <- dat %>%
    group_by(gene) %>%
    summarise(gene_w = sum(wp), .groups="drop") %>%
    right_join(gene_theta, by="gene") %>%
    mutate(
      col = gene_cols[gene],
      size = rescale(gene_w, to=gene_size_range)
    )

  n_labels <- if (is.infinite(max_gene_labels)) {
    nrow(gene_nodes)
  } else {
    max(0, min(as.integer(max_gene_labels), nrow(gene_nodes)))
  }

  gene_labels <- gene_nodes %>%
    arrange(desc(gene_w)) %>%
    slice_head(n = n_labels) %>%
    mutate(
      ang = theta * 180/pi,
      ang2 = ifelse(ang > 90 & ang < 270, ang + 180, ang),
      hjust = ifelse(ang > 90 & ang < 270, 1, 0),
      xlab = label_expand * x,
      ylab = label_expand * y
    )

  edges <- dat %>%
    left_join(tf_theta %>% transmute(TF, x1=x, y1=y, tf_col=col), by="TF") %>%
    left_join(gene_nodes %>% transmute(gene, x2=x, y2=y, gene_col=col), by="gene") %>%
    mutate(
      ew = rescale(wp, to=edge_width_range),
      ecol = dplyr::case_when(
        edge_color_by_tf ~ tf_col,
        edge_color_by_gene ~ gene_col,
        TRUE ~ edge_color
      )
    )

  ggplot() +
    geom_curve(
      data = edges,
      aes(x=x1, y=y1, xend=x2, yend=y2, linewidth=ew),
      curvature = 0.12,
      colour = edges$ecol,
      alpha = edge_alpha,
      lineend = "round"
    ) +
    scale_linewidth_identity() +
    geom_point(
      data = gene_nodes,
      aes(x=x, y=y, size=size),
      shape=21, fill=gene_nodes$col,
      colour="white", stroke=0.35, alpha=0.90
    ) +
    geom_point(
      data = tf_theta,
      aes(x=x, y=y, size=size),
      shape=21, fill=tf_theta$col_dir,
      colour="white", stroke=0.45, alpha=0.90
    ) +
    scale_size_identity() +
    geom_text(
      data = tf_theta,
      aes(x=x, y=y, label=as.character(TF)),
      fontface="bold", size=tf_label_size, colour="black"
    ) +
    geom_text(
      data = gene_labels,
      aes(x=xlab, y=ylab, label=gene, angle=ang2, hjust=hjust),
      fontface="bold", size=gene_label_size,
      colour = gene_labels$col
    ) +
    coord_equal(xlim=c(-x_expand, x_expand), ylim=c(-x_expand, x_expand), clip="off") +
    theme_void(base_size = 12)
}



# =========================
# ： triplets ->  TF-gene ->  TF_direction -> 
# =========================

# tri_all：“”TF-peak-gene 
# （）， CSV 
tri_all <- NULL
if (exists("TF_peak_gene_consistent")) {
  tri_all <- TF_peak_gene_consistent
} else {
  tri_csv <- file.path(outdir, "TF_peak_gene_triplets_RNA_consistent.csv")
  if (!file.exists(tri_csv)) {
    stop(" triplets ：", tri_csv, "\n netlink.R  TF_peak_gene_triplets_RNA_consistent.csv")
  }
  tri_all <- readr::read_csv(
    tri_csv,
    show_col_types = FALSE
  )
}

#  CSV （ TopN）
tri_sel <- tri_all %>%
  filter(!is.na(TF_direction)) %>%
  filter(!is.na(strength_rna), is.finite(strength_rna)) %>%
  distinct(TF_gene, peak, gene, .keep_all = TRUE) %>%
  arrange(desc(strength_rna))

readr::write_csv(
  tri_sel,
  file.path(outdir, "TF_peak_gene_triplets_RNA_consistent_ALL_for_circle.csv")
)

dirs_present <- tri_sel %>%
  dplyr::pull(TF_direction) %>%
  as.character() %>%
  unique() %>%
  na.omit() %>%
  as.character()

preferred_order <- c("CON_up", "PRE_up", "T2D_up")
dirs_to_plot <- c(intersect(preferred_order, dirs_present), setdiff(dirs_present, preferred_order))

if (length(dirs_to_plot) == 0) {
  stop("tri_sel  TF_direction（ NA ）")
}

plot_one_dir <- function(dir_label) {
  tri_dir <- tri_sel %>%
    filter(TF_direction == dir_label) %>%
    arrange(desc(strength_rna))

  if (is.finite(max_triplets_per_dir) && max_triplets_per_dir > 0) {
    tri_dir <- tri_dir %>% slice_head(n = as.integer(max_triplets_per_dir))
  }

  if (nrow(tri_dir) == 0) return(invisible(NULL))

  #  TF->gene （ TF-gene  peak ）
  tf_gene <- tri_dir %>%
    transmute(TF = TF_gene, target = gene, importance = strength_rna) %>%
    group_by(TF, target) %>%
    summarise(importance = sum(importance, na.rm = TRUE), .groups = "drop")

  #  TF_direction ， tf_dir_df 
  tf_dir <- tri_dir %>%
    filter(!is.na(TF_gene), !is.na(TF_direction)) %>%
    group_by(TF = TF_gene, TF_direction) %>%
    summarise(w = sum(strength_rna, na.rm = TRUE), .groups = "drop") %>%
    group_by(TF) %>%
    slice_max(order_by = w, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    dplyr::select(TF, TF_direction)

  p <- plot_tf_gene_circle_unique(
    df = tf_gene,
    tf_col = "TF",
    gene_col = "target",
    weight_col = "importance",
    tf_dir_df = tf_dir,
    top_n_per_tf = max_genes_per_tf_in_plot,
    weight_transform = function(x) log1p(x),
    r_tf = 0.8,
    tf_label_size = 6,
    edge_color_by_tf = FALSE,
    max_gene_labels = Inf,
    x_expand = 1.28
  )

  safe_label <- gsub("[^A-Za-z0-9_]+", "_", dir_label)
  out_png <- file.path(outdir, paste0("TF_gene_circle_unique_", safe_label, ".png"))
  ggsave(out_png, p, width =14 , height = 14, dpi = 600, bg = "white")
  cat("[INFO] Saved circle plot: ", out_png, " (n_triplets=", nrow(tri_dir), ")\n", sep = "")
  p
}

plots <- lapply(dirs_to_plot, plot_one_dir)

# （Rscript ）
if (length(plots) > 0 && !is.null(plots[[1]])) {
  print(plots[[1]])
}
