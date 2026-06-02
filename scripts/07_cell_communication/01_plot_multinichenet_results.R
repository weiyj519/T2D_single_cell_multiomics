#!/usr/bin/env Rscript

# Script name: 01_plot_multinichenet_results.R
# Purpose: Cleaned public workflow script for the T2D single-cell multiome analysis.
# Input: Paths and analysis settings are read from configs/config.yaml.
# Output: Module-specific outputs are written under the configured results directory.
# Main steps: Load configuration, run the original analysis logic, and export reproducible outputs.
# Example command: Rscript scripts/07_cell_communication/01_plot_multinichenet_results.R

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
  library(tidyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(ggrepel)
  library(grid)
})

# =========================================================
# 0. 
# =========================================================

outdir <- "results/downstream/acinar_CON_PRE/multinichenet"

plot_dir <- file.path(outdir, "plot")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# 
prior_csv <- file.path(outdir, "prioritization_group_prioritization_tbl.csv")
ligtar_csv <- file.path(outdir, "ligand_activities_targets_MLmasked.csv")

# receiver ：multinichenet  make.names （ Acinar cell -> Acinar.cell）
receiver_focus <- "Acinar.cell"

groups_keep <- c("CON", "PRE")

# 
top_n_heatmap_per_group <- 15   # CON/PRE  top N， 2N
top_n_dumbbell <- 25
top_n_ligands_per_group <- 8
top_targets_per_ligand <- 8
plot_font_family <- "Arial"
plot_base_size <- 18

# =========================================================
# 1. 
# =========================================================

clean_label <- function(x) {
  x %>%
    as.character() %>%
    # （ α.cell， alpha cell）
    str_replace_all("α\\.cell", "α cell") %>%
    str_replace_all("β\\.cell", "β cell") %>%
    str_replace_all("δ\\.cell", "δ cell") %>%
    str_replace_all(regex("\\balpha\\s*cell\\b", ignore_case = TRUE), "α cell") %>%
    str_replace_all(regex("\\bbeta\\s*cell\\b", ignore_case = TRUE), "β cell") %>%
    str_replace_all(regex("\\bdelta\\s*cell\\b", ignore_case = TRUE), "δ cell") %>%
    str_replace_all("\\.", " ")
}

detect_col <- function(df, candidates, required = TRUE, label = "column") {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) {
    if (required) {
      stop("Cannot find ", label, ". Tried: ", paste(candidates, collapse = ", "))
    } else {
      return(NA_character_)
    }
  }
  hit[1]
}

make_lr_id <- function(df) {
  if ("id" %in% colnames(df)) {
    as.character(df$id)
  } else {
    paste(df$sender, df$receiver, df$ligand, df$receptor, sep = "|")
  }
}

# =========================================================
# 2.  prioritization 
# =========================================================

stopifnot(file.exists(prior_csv))

prior0 <- read_csv(prior_csv, show_col_types = FALSE)

need_basic <- c("sender", "receiver", "ligand", "receptor", "group")
miss_basic <- setdiff(need_basic, colnames(prior0))

if (length(miss_basic) > 0) {
  stop("prioritization ：", paste(miss_basic, collapse = ", "))
}

score_col <- detect_col(
  prior0,
  candidates = c(
    "prioritization_score",
    "prioritization_score_group",
    "score",
    "priority_score",
    "final_score",
    "aggregate_rank_score"
  ),
  required = TRUE,
  label = "prioritization score column"
)

prior <- prior0 %>%
  mutate(
    sender = as.character(sender),
    receiver = as.character(receiver),
    ligand = as.character(ligand),
    receptor = as.character(receptor),
    group = as.character(group),
    score = as.numeric(.data[[score_col]]),
    lr_id = make_lr_id(.),
    sender_label = clean_label(sender),
    receiver_label = clean_label(receiver),
    ligand_label = clean_label(ligand),
    receptor_label = clean_label(receptor),
    lr_label = paste0(ligand_label, " - ", receptor_label),
    sender_lr_label = paste0(sender_label, " | ", ligand_label, " - ", receptor_label)
  ) %>%
  filter(group %in% groups_keep) %>%
  filter(receiver == receiver_focus) %>%
  filter(!is.na(score))

cat("[INFO] prioritization rows after receiver/group filter:", nrow(prior), "\n")
cat("[INFO] score column used:", score_col, "\n")

if (nrow(prior) == 0) {
  stop(
    "prior  receiver_focus ：",
    paste(unique(prior0$receiver), collapse = ", ")
  )
}

# =========================================================
# Figure 1. LR 
# ：CON/PRE
# ：L-R pair
# ：sender 
# =========================================================

# CON  PRE  top N，
lr_scores_by_group <- prior %>%
  group_by(group, lr_id) %>%
  summarise(max_score = max(score, na.rm = TRUE), .groups = "drop") %>%
  filter(group %in% groups_keep) %>%
  arrange(group, desc(max_score), lr_id) %>%
  group_by(group) %>%
  mutate(rank_in_group = row_number()) %>%
  ungroup()

target_total_lr <- length(groups_keep) * top_n_heatmap_per_group
selected_lr <- character(0)

for (g in groups_keep) {
  selected_lr <- c(
    selected_lr,
    lr_scores_by_group %>%
      filter(group == g, rank_in_group <= top_n_heatmap_per_group) %>%
      pull(lr_id)
  )
}

selected_lr <- unique(selected_lr)

rank_cursor <- top_n_heatmap_per_group + 1

max_rank <- lr_scores_by_group %>%
  group_by(group) %>%
  summarise(max_rank = max(rank_in_group), .groups = "drop") %>%
  summarise(max(max_rank), .groups = "drop") %>%
  pull()

while (length(selected_lr) < target_total_lr && rank_cursor <= max_rank) {
  for (g in groups_keep) {
    if (length(selected_lr) >= target_total_lr) break

    cand <- lr_scores_by_group %>%
      filter(group == g, rank_in_group == rank_cursor) %>%
      pull(lr_id)

    cand <- cand[!cand %in% selected_lr]

    if (length(cand) > 0) {
      selected_lr <- c(selected_lr, cand[1])
    }
  }

  selected_lr <- unique(selected_lr)
  rank_cursor <- rank_cursor + 1
}

top_lr_ids <- selected_lr

hm_df <- prior %>%
  filter(lr_id %in% top_lr_ids) %>%
  group_by(lr_id, group) %>%
  summarise(
    score = max(score, na.rm = TRUE),
    sender_label = first(sender_label),
    lr_label = first(lr_label),
    .groups = "drop"
  ) %>%
  complete(lr_id, group = groups_keep, fill = list(score = 0)) %>%
  group_by(lr_id) %>%
  mutate(
    sender_label = sender_label[!is.na(sender_label)][1],
    lr_label = lr_label[!is.na(lr_label)][1],
    max_score = max(score, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  filter(!is.na(sender_label), !is.na(lr_label)) %>%
  arrange(sender_label, desc(max_score), lr_label) %>%
  mutate(
    sender_facet_label = factor(
      sender_label,
      levels = unique(sender_label)
    ),

    #  sender  L-R 
    lr_y_id = paste(sender_label, lr_label, sep = "___"),
    lr_y_id = factor(lr_y_id, levels = rev(unique(lr_y_id))),

    group = factor(group, levels = groups_keep)
  )

cat("[INFO] LR heatmap rows:", nrow(hm_df), "\n")
cat("[INFO] LR pairs shown:", length(unique(hm_df$lr_y_id)), "\n")
cat("[INFO] Senders shown:", length(unique(hm_df$sender_label)), "\n")

p_lr_heatmap <- ggplot(hm_df, aes(x = group, y = lr_y_id, fill = score)) +
  geom_tile(
    color = "white",
    linewidth = 0.45,
    width = 0.92,
    height = 0.92
  ) +
  facet_grid(
    rows = vars(sender_facet_label),
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_y_discrete(
    labels = function(x) str_replace(x, "^.*___", "")
  ) +
  scale_fill_gradient(
    low = "grey95",
    high = "#2878B5",
    name = "Prioritization\nscore"
  ) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = NULL,
    y = NULL
  ) +
  theme_classic(
    base_size = plot_base_size,
    base_family = plot_font_family
  ) +
  theme(
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.7
    ),
    # sender （）
    panel.spacing.y = unit(0.22, "lines"),

    # sender 
    strip.placement = "outside",
    strip.background.y = element_rect(
      fill = "grey90",
      color = "black",
      linewidth = 0.7
    ),
    strip.text.y.left = element_text(
      face = "bold",
      size = 14,
      angle = 0,
      margin = margin(r = 5, l = 5)
    ),

    axis.text.x = element_text(
      face = "bold",
      size = 18
    ),
    axis.text.y = element_text(
      face = "bold",
      size = 12
    ),
    axis.ticks.y = element_blank(),

    plot.title = element_blank(),
    plot.subtitle = element_blank(),

    legend.title = element_text(
      face = "bold",
      size = 14
    ),
    legend.text = element_text(size = 13),
    legend.position = "bottom",
    legend.location = "plot",
    legend.justification = "center",
    legend.box.just = "center"
  )

ggsave(
  file.path(plot_dir, "Fig_LR_priority_heatmap_senderBox_CON_PRE_acinarReceiver.pdf"),
  p_lr_heatmap,
  device = cairo_pdf,
  width =4.5,
  height = max(8, 0.3* length(unique(hm_df$lr_y_id))),
  bg = "white"
)



cat("[OK] LR priority heatmap with sender boxes saved.\n")

# =========================================================
# Figure 2. Dumbbell plot： CON vs PRE
# =========================================================

dumb_df_wide <- prior %>%
  group_by(lr_id, group) %>%
  summarise(
    score = max(score, na.rm = TRUE),
    sender_label = first(sender_label),
    sender_lr_label = first(sender_lr_label),
    .groups = "drop"
  ) %>%
  complete(lr_id, group = groups_keep, fill = list(score = 0)) %>%
  group_by(lr_id) %>%
  mutate(
    sender_label = sender_label[!is.na(sender_label)][1],
    sender_lr_label = sender_lr_label[!is.na(sender_lr_label)][1]
  ) %>%
  ungroup() %>%
  filter(!is.na(sender_lr_label)) %>%
  pivot_wider(
    names_from = group,
    values_from = score,
    values_fill = 0
  )

# 
for (g in groups_keep) {
  if (!g %in% colnames(dumb_df_wide)) {
    dumb_df_wide[[g]] <- 0
  }
}

dumb_df_wide <- dumb_df_wide %>%
  mutate(
    delta_PRE_minus_CON = .data[["PRE"]] - .data[["CON"]],
    abs_delta = abs(delta_PRE_minus_CON),
    dominant_group = case_when(
      delta_PRE_minus_CON > 0 ~ "PRE higher",
      delta_PRE_minus_CON < 0 ~ "CON higher",
      TRUE ~ "Equal"
    )
  ) %>%
  arrange(desc(abs_delta)) %>%
  slice_head(n = top_n_dumbbell) %>%
  arrange(delta_PRE_minus_CON) %>%
  mutate(
    sender_lr_label = factor(sender_lr_label, levels = sender_lr_label)
  )

dumb_df_long <- dumb_df_wide %>%
  select(lr_id, sender_label, sender_lr_label, all_of(groups_keep), delta_PRE_minus_CON) %>%
  pivot_longer(
    cols = all_of(groups_keep),
    names_to = "group",
    values_to = "score"
  ) %>%
  mutate(group = factor(group, levels = groups_keep))

p_dumbbell <- ggplot(dumb_df_wide, aes(y = sender_lr_label)) +
  geom_segment(
    aes(
      x = .data[["CON"]],
      xend = .data[["PRE"]],
      yend = sender_lr_label
    ),
    color = "grey65",
    linewidth = 0.9
  ) +
  geom_point(
    data = dumb_df_long,
    aes(x = score, y = sender_lr_label, color = group),
    size = 3.5,
    alpha = 0.95
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey45"
  ) +
  scale_color_manual(
    values = c("CON" = "#4DBBD5", "PRE" = "#E64B35"),
    name = NULL
  ) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = "Prioritization score",
    y = NULL
  ) +
  theme_classic(
    base_size = plot_base_size,
    base_family = plot_font_family
  ) +
  theme(
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),
    axis.text.y = element_text(
      face = "bold",
      size = 13
    ),
    axis.text.x = element_text(
      face = "bold",
      size = 15
    ),
    axis.title.x = element_text(
      face = "bold",
      size = 17
    ),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    legend.position = "top",
    legend.text = element_text(
      face = "bold",
      size = 15
    )
  )

ggsave(
  file.path(plot_dir, "Fig_Dumbbell_CON_PRE_LR_priority_acinarReceiver.pdf"),
  p_dumbbell,
  device = cairo_pdf,
  width = 12,
  height = max(9, 0.45 * nrow(dumb_df_wide)),
  bg = "white"
)

ggsave(
  file.path(plot_dir, "Fig_Dumbbell_CON_PRE_LR_priority_acinarReceiver.png"),
  p_dumbbell,
  width = 12,
  height = max(9, 0.45 * nrow(dumb_df_wide)),
  dpi = 600,
  bg = "white"
)

cat("[OK] Dumbbell plot saved.\n")

# =========================================================
# Figure 3. Ligand-target heatmap
# =========================================================

stopifnot(file.exists(ligtar_csv))

lt0 <- read_csv(ligtar_csv, show_col_types = FALSE)

if (!all(c("ligand", "target") %in% colnames(lt0))) {
  stop("ligand_activities_targets_MLmasked.csv  ligand / target ， ligand-target heatmap")
}

lt <- lt0 %>%
  mutate(
    ligand = as.character(ligand),
    target = as.character(target)
  )

if (!"group" %in% colnames(lt)) {
  if ("contrast" %in% colnames(lt)) {
    lt <- lt %>%
      mutate(
        contrast = as.character(contrast),
        group = case_when(
          contrast == "CON-PRE" ~ "CON",
          contrast == "PRE-CON" ~ "PRE",
          TRUE ~ NA_character_
        )
      )
  } else {
    stop("ligand-target  group， contrast， CON/PRE")
  }
}

#  ligand-target 
evidence_col <- detect_col(
  lt,
  candidates = c(
    "regulatory_potential",
    "target_score",
    "target_weight",
    "weight",
    "pearson",
    "activity_scaled",
    "aupr_corrected",
    "aupr",
    "activity"
  ),
  required = FALSE,
  label = "ligand-target evidence column"
)

rank_col <- detect_col(
  lt,
  candidates = c(
    "rank_of_target",
    "target_rank",
    "rank"
  ),
  required = FALSE,
  label = "target rank column"
)

cat("[INFO] ligand-target evidence column used:", evidence_col, "\n")
cat("[INFO] ligand-target rank column used:", rank_col, "\n")

lt <- lt %>%
  mutate(
    group = as.character(group),
    direction_regulation = if ("direction_regulation" %in% colnames(.)) {
      as.character(direction_regulation)
    } else {
      "unknown"
    },
    evidence_raw = if (!is.na(evidence_col)) {
      as.numeric(.data[[evidence_col]])
    } else if (!is.na(rank_col)) {
      1 / as.numeric(.data[[rank_col]])
    } else {
      1
    },
    evidence_raw = ifelse(is.na(evidence_raw), 0, evidence_raw),
    evidence_signed = case_when(
      direction_regulation == "down" ~ -abs(evidence_raw),
      direction_regulation == "up" ~ abs(evidence_raw),
      TRUE ~ evidence_raw
    )
  ) %>%
  filter(group %in% groups_keep)

#  prioritization  top ligands， sender 
top_ligands <- prior %>%
  group_by(group, ligand) %>%
  summarise(
    ligand_score = max(score, na.rm = TRUE),
    sender_label = sender_label[which.max(score)][1],
    .groups = "drop"
  ) %>%
  group_by(group) %>%
  arrange(desc(ligand_score), .by_group = TRUE) %>%
  slice_head(n = top_n_ligands_per_group) %>%
  ungroup() %>%
  mutate(
    ligand_label = paste0(sender_label, " | ", clean_label(ligand))
  )

lt_plot <- lt %>%
  inner_join(
    top_ligands %>%
      select(group, ligand, ligand_label, ligand_score),
    by = c("group", "ligand")
  )

#  ligand  top targets
if (!is.na(rank_col)) {
  lt_plot <- lt_plot %>%
    mutate(target_rank_use = as.numeric(.data[[rank_col]])) %>%
    group_by(group, ligand) %>%
    arrange(target_rank_use, .by_group = TRUE) %>%
    slice_head(n = top_targets_per_ligand) %>%
    ungroup()
} else {
  lt_plot <- lt_plot %>%
    group_by(group, ligand) %>%
    arrange(desc(abs(evidence_signed)), .by_group = TRUE) %>%
    slice_head(n = top_targets_per_ligand) %>%
    ungroup()
}

lt_plot <- lt_plot %>%
  mutate(
    target_label = clean_label(target),
    ligand_label = factor(ligand_label, levels = rev(unique(ligand_label))),
    target_label = fct_reorder(target_label, abs(evidence_signed), .fun = max),
    group = factor(group, levels = groups_keep)
  )

cat("[INFO] ligand-target heatmap rows:", nrow(lt_plot), "\n")

if (nrow(lt_plot) == 0) {
  stop("ligand-target heatmap  top ligands  ligand-target  target ")
}

p_ligtar_heatmap <- ggplot(
  lt_plot,
  aes(x = target_label, y = ligand_label, fill = evidence_signed)
) +
  geom_tile(
    color = "white",
    linewidth = 0.3
  ) +
  facet_grid(
    rows = vars(group),
    scales = "free_y",
    space = "free_y"
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "grey95",
    high = "#B2182B",
    midpoint = 0,
    name = ifelse(
      !is.na(evidence_col),
      paste0("Signed\n", evidence_col),
      "Signed\nevidence"
    )
  ) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = "Predicted target genes",
    y = "Sender | ligand"
  ) +
  theme_classic(
    base_size = plot_base_size,
    base_family = plot_font_family
  ) +
  theme(
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),
    strip.background = element_rect(
      fill = "grey90",
      color = "black",
      linewidth = 0.5
    ),
    strip.text.y = element_text(
      face = "bold",
      size = 16
    ),
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      vjust = 1,
      face = "bold",
      size = 13
    ),
    axis.text.y = element_text(
      face = "bold",
      size = 13
    ),
    axis.title.x = element_text(
      face = "bold",
      size = 17
    ),
    axis.title.y = element_text(
      face = "bold",
      size = 17
    ),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 15)
  )

ggsave(
  file.path(plot_dir, "Fig_Ligand_target_heatmap_acinarReceiver.pdf"),
  p_ligtar_heatmap,
  device = cairo_pdf,
  width = max(15, 0.4 * length(unique(lt_plot$target_label))),
  height = max(9, 0.5 * length(unique(lt_plot$ligand_label))),
  bg = "white"
)

ggsave(
  file.path(plot_dir, "Fig_Ligand_target_heatmap_acinarReceiver.png"),
  p_ligtar_heatmap,
  width = max(15, 0.4 * length(unique(lt_plot$target_label))),
  height = max(9, 0.5 * length(unique(lt_plot$ligand_label))),
  dpi = 600,
  bg = "white"
)

cat("[OK] Ligand-target heatmap saved.\n")
cat("[DONE] All plots saved to: ", plot_dir, "\n")
