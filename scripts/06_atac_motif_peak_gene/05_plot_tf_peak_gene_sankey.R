#!/usr/bin/env Rscript

# Script name: 05_plot_tf_peak_gene_sankey.R
# Purpose: Cleaned public workflow script for the T2D single-cell multiome analysis.
# Input: Paths and analysis settings are read from configs/config.yaml.
# Output: Module-specific outputs are written under the configured results directory.
# Main steps: Load configuration, run the original analysis logic, and export reproducible outputs.
# Example command: Rscript scripts/06_atac_motif_peak_gene/05_plot_tf_peak_gene_sankey.R

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

# install.packages("ggsankeyfier")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggsankeyfier)
  library(scales)
})

# =========================================================
# 1)  / 
# =========================================================
infile <- "results/downstream/acinar_CON_PRE/netlink/TF_peak_gene_triplets_RNA_consistent_abs_rnaLFC0.25_padj0.05.csv"
outdir <- "results/downstream/acinar_CON_PRE/netlink"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# =========================================================
# 2) 
# =========================================================
to_logical_safe <- function(x) {
  if (is.logical(x)) return(x)
  x <- tolower(trimws(as.character(x)))
  x %in% c("true", "t", "1", "yes", "y")
}

choose_col <- function(dat, preferred, candidates, label) {
  if (preferred %in% colnames(dat)) return(preferred)

  hit <- candidates[candidates %in% colnames(dat)]
  if (length(hit) > 0) {
    message("[INFO] : ", hit[1], "  ", label)
    return(hit[1])
  }

  stop(
    " ", label, " \n",
    ": ", preferred, "\n",
    ": ", paste(candidates, collapse = ", "), "\n",
    ":\n", paste(colnames(dat), collapse = ", ")
  )
}

# =========================================================
# 3) 
# =========================================================
plot_sankey_tf_peak_gene_combined_dirs <- function(
    triplet_csv,
    out_file,
    title = NULL,
    top_n_each_dir = Inf,
    max_tf_each_dir = 20,
    weight_col = "strength_rna",
    direction_col = "gene_direction",
    tf_col = "TF_gene",
    keep_dirs = c("CON_up", "PRE_up"),
    filter_pass_sig = TRUE,
    use_prop = TRUE,
    normalize_by_dir = FALSE,
    show_peak_labels = FALSE,
    min_gene_total = 0.03,
    width = 12,
    height = 8,
    dpi = 500
) {

  # -------------------------
  # 3.1 
  # -------------------------
  dat <- read_csv(triplet_csv, show_col_types = FALSE)

  tf_col <- choose_col(
    dat, tf_col,
    candidates = c("TF_gene", "TF", "motif"),
    label = "TF"
  )

  direction_col <- choose_col(
    dat, direction_col,
    candidates = c("gene_direction", "TF_direction", "peak_direction", "direction"),
    label = "direction"
  )

  peak_col <- choose_col(
    dat, "peak",
    candidates = c("peak", "Peak", "peak_id"),
    label = "peak"
  )

  gene_col <- choose_col(
    dat, "gene",
    candidates = c("gene", "gene_name", "symbol", "SYMBOL"),
    label = "gene"
  )

  if (!weight_col %in% colnames(dat)) {
    stop(
      ": ", weight_col,
      "\n:\n", paste(colnames(dat), collapse = ", ")
    )
  }

  # -------------------------
  # 3.2 
  # -------------------------
  dat <- dat %>%
    mutate(
      TF_plot   = as.character(.data[[tf_col]]),
      peak_plot = as.character(.data[[peak_col]]),
      gene_plot = as.character(.data[[gene_col]]),
      direction = as.character(.data[[direction_col]]),
      weight_raw = abs(as.numeric(.data[[weight_col]]))
    ) %>%
    filter(
      !is.na(TF_plot), TF_plot != "",
      !is.na(peak_plot), peak_plot != "",
      !is.na(gene_plot), gene_plot != "",
      !is.na(direction), direction != "",
      is.finite(weight_raw), !is.na(weight_raw), weight_raw > 0
    )

  if (filter_pass_sig && "pass_sig" %in% colnames(dat)) {
    dat <- dat %>%
      mutate(pass_sig_logical = to_logical_safe(pass_sig)) %>%
      filter(pass_sig_logical)
  }

  dat <- dat %>%
    filter(direction %in% keep_dirs)

  if (nrow(dat) == 0) {
    stop("， pass_sig")
  }

  # -------------------------
  # 3.3  top N
  #  top_n_each_dir = Inf，，
  # “”，
  # -------------------------
  tri <- dat %>%
    group_by(direction, TF_plot, peak_plot, gene_plot) %>%
    summarise(flow_raw = sum(weight_raw, na.rm = TRUE), .groups = "drop")

  if (!is.null(max_tf_each_dir) && is.finite(max_tf_each_dir)) {
    keep_tfs <- tri %>%
      group_by(direction, TF_plot) %>%
      summarise(tf_total = sum(flow_raw, na.rm = TRUE), .groups = "drop") %>%
      group_by(direction) %>%
      arrange(desc(tf_total), .by_group = TRUE) %>%
      slice_head(n = max_tf_each_dir) %>%
      ungroup()

    tri <- tri %>%
      semi_join(keep_tfs, by = c("direction", "TF_plot"))
  }

  if (is.finite(top_n_each_dir)) {
    tri <- tri %>%
      group_by(direction) %>%
      arrange(desc(flow_raw), .by_group = TRUE) %>%
      slice_head(n = top_n_each_dir) %>%
      ungroup()
  } else {
    tri <- tri %>%
      arrange(direction, desc(flow_raw))
  }

  if (nrow(tri) == 0) {
    stop("TopN ")
  }

  if (!is.null(min_gene_total) && is.finite(min_gene_total) && min_gene_total > 0) {
    keep_genes <- tri %>%
      group_by(gene_plot) %>%
      summarise(gene_total = sum(flow_raw, na.rm = TRUE), .groups = "drop") %>%
      filter(gene_total >= min_gene_total) %>%
      pull(gene_plot)

    tri <- tri %>%
      filter(gene_plot %in% keep_genes)
  }

  if (nrow(tri) == 0) {
    stop("Gene ， min_gene_total")
  }

  # -------------------------
  # 3.4 
  # -------------------------
  if (use_prop) {
    if (normalize_by_dir) {
      tri <- tri %>%
        group_by(direction) %>%
        mutate(flow = flow_raw / sum(flow_raw)) %>%
        ungroup()
    } else {
      tri <- tri %>%
        mutate(flow = flow_raw / sum(flow_raw))
    }
  } else {
    tri <- tri %>%
      mutate(flow = flow_raw)
  }

  # 
  selected_csv <- file.path(
    dirname(out_file),
    paste0(tools::file_path_sans_ext(basename(out_file)), "_selected_triplets.csv")
  )
  write_csv(tri %>% arrange(direction, desc(flow_raw)), selected_csv)
  message("[SAVED selected triplets] ", selected_csv)

  # -------------------------
  # 3.5 
  # Direction -> TF -> Peak -> Gene
  # -------------------------
  tri_plot <- tri %>%
    mutate(
      Direction = direction,
      TF   = TF_plot,
      Peak = peak_plot,
      Gene = gene_plot
    ) %>%
    select(direction, Direction, TF, Peak, Gene, flow)

  sank <- pivot_stages_longer(
    tri_plot,
    stages_from = c("Direction", "TF", "Peak", "Gene"),
    values_from = "flow",
    additional_aes_from = "direction"
  )

  # -------------------------
  # 3.6 
  # -------------------------
  if (!"direction" %in% colnames(sank)) {
    stop("pivot_stages_longer  direction ， CON_up / PRE_up ")
  }

  sank <- sank %>%
    mutate(
      stage = as.character(stage),
      node = as.character(node),
      direction = as.character(direction),
      stage = factor(stage, levels = c("Direction", "TF", "Peak", "Gene")),
      direction = factor(direction, levels = keep_dirs),
      node_label = node
    )

  if (!show_peak_labels) {
    sank <- sank %>%
      mutate(node_label = if_else(stage == "Peak", "", node_label))
  }

  if (any(is.na(sank$direction))) {
    stop(" Sankey edge  direction， CON_up / PRE_up ")
  }

  # -------------------------
  # 3.7 
  # -------------------------
  dir_cols <- c(
    "CON_up" = "#62B197",
    "PRE_up" = "#E18E6D"
  )

  missing_dirs <- setdiff(keep_dirs, names(dir_cols))
  if (length(missing_dirs) > 0) {
    extra_cols <- scales::hue_pal()(length(missing_dirs))
    names(extra_cols) <- missing_dirs
    dir_cols <- c(dir_cols, extra_cols)
  }
  dir_cols <- dir_cols[keep_dirs]

  pos <- position_sankey(
    v_space = "auto",
    order = "ascending",
    align = "justify"
  )

  pos_text <- position_sankey(
    v_space = "auto",
    order = "ascending",
    align = "justify",
    nudge_x = 0.10
  )

  # -------------------------
  # 3.8 
  # -------------------------
  p <- ggplot(
    sank,
    aes(
      x = stage,
      y = flow,
      group = node,
      connector = connector,
      edge_id = edge_id
    )
  ) +
    geom_sankeyedge(
      aes(fill = direction),
      alpha = 0.62,
      color = "transparent",
      position = pos
    ) +
    geom_sankeynode(
      fill = "white",
      color = "black",
      linewidth = 0.45,
      position = pos
    ) +
    geom_text(
      aes(label = node_label),
      stat = "sankeynode",
      position = pos_text,
      hjust = 0,
      vjust = 0,
      size = 4.8,
      fontface = "bold",
      color = "gray20",
      check_overlap = FALSE
    ) +
    scale_fill_manual(
      values = dir_cols,
      breaks = keep_dirs,
      drop = FALSE,
      na.translate = FALSE,
      name = "Direction"
    ) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 22) +
    theme(
      text = element_text(face = "bold"),
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 18, face = "bold", color = "gray20"),
      axis.title = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 15, face = "bold"),
      plot.title = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )

  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }

  ggsave(
    filename = out_file,
    plot = p,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )

  message("[SAVED plot] ", out_file)

  invisible(p)
}

# =========================================================
# 4) 
# =========================================================

# 1：
# ， CON_up / PRE_up 
plot_sankey_tf_peak_gene_combined_dirs(
  triplet_csv = infile,
  out_file = file.path(outdir, "TF_peak_gene_sankey_acinar_CON_PRE_combined_global_all_sig.pdf"),
  title = NULL,
  top_n_each_dir = Inf,
  max_tf_each_dir = 10,
  weight_col = "strength_rna",
  direction_col = "gene_direction",
  tf_col = "TF_gene",
  keep_dirs = c("CON_up", "PRE_up"),
  filter_pass_sig = TRUE,
  use_prop = TRUE,
  normalize_by_dir = FALSE,
  show_peak_labels = FALSE,
  min_gene_total = 0.03,
  width = 12,
  height = 8,
  dpi = 500
)



# 2：
# ，
plot_sankey_tf_peak_gene_combined_dirs(
  triplet_csv = infile,
  out_file = file.path(outdir, "TF_peak_gene_sankey_acinar_CON_PRE_combined_withinDir_all_sig.pdf"),
  title = NULL,
  top_n_each_dir = Inf,
  max_tf_each_dir = 10,
  weight_col = "strength_rna",
  direction_col = "gene_direction",
  tf_col = "TF_gene",
  keep_dirs = c("CON_up", "PRE_up"),
  filter_pass_sig = TRUE,
  use_prop = TRUE,
  normalize_by_dir = TRUE,
  show_peak_labels = FALSE,
  min_gene_total = 0.03,
  width = 12,
  height = 8,
  dpi = 500
)



message("\n[DONE] ：\n")
