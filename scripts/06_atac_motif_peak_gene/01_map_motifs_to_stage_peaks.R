#!/usr/bin/env Rscript

# Script name: 01_map_motifs_to_stage_peaks.R
# Purpose: Cleaned public workflow script for the T2D single-cell multiome analysis.
# Input: Paths and analysis settings are read from configs/config.yaml.
# Output: Module-specific outputs are written under the configured results directory.
# Main steps: Load configuration, run the original analysis logic, and export reproducible outputs.
# Example command: Rscript scripts/06_atac_motif_peak_gene/01_map_motifs_to_stage_peaks.R

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
  library(qs)
  library(Seurat)
  library(Signac)
  library(Matrix)
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
})

# =========================
# 0) （acinar: CON vs PRE）
# =========================
project_dir <- "results/downstream/acinar_CON_PRE"

#  ATAC  peak （：CON_vs_PRE）
# -  significance ， lfc = log2(CON/PRE)
# ：peak, p_val_adj, lfc/avg_log2FC/avg_logFC, mean_SHAP, mean_abs_SHAP
da_csv <- file.path(project_dir, "motif", "DA_featurePeaks_acinar_CON_vs_PRE.csv")

chromvar_dir <- file.path(project_dir, "chromvar")
acinar_cv_qs  <- file.path(chromvar_dir, "acinar_CON_vs_PRE_chromvar.qs")
# ： PRE_vs_CON  chromVAR motif （ chromvar.R ）
# -> chromvar_lfc > 0  PRE_up
motif_da_csv  <- file.path(chromvar_dir, "chromvar_PRE_vs_CON_motifs.csv")

outdir_map <- file.path(project_dir, "mlpeaks_chromvar_mapping")
dir.create(outdir_map, recursive = TRUE, showWarnings = FALSE)

# （）
padj_cut <- 0.05
lfc_cut  <- 0.25

# =========================
# 1)  DA+SHAP（ CSV）
# =========================
if (!file.exists(da_csv)) {
  cand <- list.files(project_dir, pattern = "DA_featurePeaks.*\\.csv$", recursive = TRUE, full.names = TRUE)
  if (length(cand) > 0) {
    message("[WARN] da_csv not found. Auto-pick candidate: ", cand[1])
    da_csv <- cand[1]
  } else {
    stop(
      " acinar  DA_featurePeaks \n",
      " da_csv=", da_csv, "\n",
      " da_csv  acinar CON/PRE  ATACpeak（DA_featurePeaks_*_CON_vs_PRE.csv）"
    )
  }
}
da_peaks <- read_csv(da_csv, show_col_types = FALSE)

stopifnot(all(c("peak","p_val_adj") %in% colnames(da_peaks)))
stopifnot(all(c("mean_SHAP","mean_abs_SHAP") %in% colnames(da_peaks)))

#  logFC 
lfc_col <- if ("lfc" %in% colnames(da_peaks)) "lfc" else
  if ("avg_log2FC" %in% colnames(da_peaks)) "avg_log2FC" else
    if ("avg_logFC" %in% colnames(da_peaks)) "avg_logFC" else NA_character_
if (is.na(lfc_col)) stop("DA CSV  lfc / avg_log2FC / avg_logFC")

da_peaks <- da_peaks %>%
  mutate(
    peak = as.character(peak),
    p_val_adj = as.numeric(p_val_adj),
    p_val_adj = ifelse(is.na(p_val_adj) | p_val_adj <= 0, 1e-300, p_val_adj),
    lfc = as.numeric(.data[[lfc_col]]),
    mean_SHAP = as.numeric(mean_SHAP),
    mean_abs_SHAP = as.numeric(mean_abs_SHAP)
  ) %>%
  filter(!is.na(peak), peak != "") %>%
  distinct(peak, .keep_all = TRUE)

keep_peaks <- da_peaks$peak
cat("[INFO] keep_peaks from DA CSV:", length(keep_peaks), "\n")

#  significance （）；
if ("significance" %in% colnames(da_peaks)) {
  peaks_con_up <- da_peaks %>% filter(significance %in% c("Up in CON", "CON_up")) %>% pull(peak) %>% unique()
  peaks_pre_up <- da_peaks %>% filter(significance %in% c("Up in PRE", "PRE_up")) %>% pull(peak) %>% unique()
} else {
  # ：CON_vs_PRE（lfc=log2(CON/PRE)）
  # -> lfc>0  CON_up；lfc<0  PRE_up
  peaks_con_up <- da_peaks %>% filter(p_val_adj < padj_cut, lfc >=  lfc_cut) %>% pull(peak) %>% unique()
  peaks_pre_up <- da_peaks %>% filter(p_val_adj < padj_cut, lfc <= -lfc_cut) %>% pull(peak) %>% unique()
}
cat("[INFO] peaks_pre_up:", length(peaks_pre_up), "\n")
cat("[INFO] peaks_con_up:", length(peaks_con_up), "\n")

# =========================
# 2)  chromVAR  + motif DA 
# =========================
stopifnot(file.exists(acinar_cv_qs), file.exists(motif_da_csv))
acinar_cv <- qread(acinar_cv_qs)
acinar_cv <- tryCatch(JoinLayers(acinar_cv), error = function(e) acinar_cv)

# Seurat v4/v5 ：Assays()  character(assay names)  list(assay objects)
get_assay_names <- function(obj) {
  a <- tryCatch(Assays(obj), error = function(e) NULL)
  if (is.null(a)) return(character())
  if (is.character(a)) return(a)
  nm <- names(a)
  if (!is.null(nm)) return(nm)
  # ：/ Assays() ， @assays 
  tryCatch(names(obj@assays), error = function(e) character())
}

assay_names <- get_assay_names(acinar_cv)

atac_assay <- if ("ATAC" %in% assay_names) {
  "ATAC"
} else if ("peaks" %in% assay_names) {
  "peaks"
} else {
  stop(" ATAC assay assays: ", paste(assay_names, collapse = ", "))
}

chromvar_assay <- if ("chromvar" %in% assay_names) {
  "chromvar"
} else {
  stop(" chromvar assay assays: ", paste(assay_names, collapse = ", "))
}

cat("[INFO] Using assays: ATAC=", atac_assay, ", chromvar=", chromvar_assay, "\n", sep = "")

motif_da <- read_csv(motif_da_csv, show_col_types = FALSE)

lfc_col_motif <- if ("avg_log2FC" %in% colnames(motif_da)) "avg_log2FC" else
  if ("avg_logFC" %in% colnames(motif_da)) "avg_logFC" else NA_character_
if (is.na(lfc_col_motif)) stop("chromVAR motif DA CSV  avg_log2FC/avg_logFC")

motif_da2 <- motif_da %>%
  mutate(
    chromvar_lfc  = as.numeric(.data[[lfc_col_motif]]),
    chromvar_padj = as.numeric(p_val_adj),
    chromvar_padj = ifelse(is.na(chromvar_padj) | chromvar_padj <= 0, 1e-300, chromvar_padj),
    chromvar_log10padj = -log10(chromvar_padj)
  ) %>%
  dplyr::select(motif, TF, chromvar_lfc, chromvar_padj, chromvar_log10padj)

# =========================
# 3) motif -> ML peaks （/）
# =========================
#  ATAC assay  AddMotifs
motif_obj <- tryCatch(Motifs(acinar_cv, assay = atac_assay), error = function(e) NULL)
if (is.null(motif_obj)) {
  stop("acinar_cv[['", atac_assay, "']]@motifs ： motif  AddMotifs  qsave")
}

motif_mat <- GetMotifData(acinar_cv, assay = atac_assay, slot = "data")  # peaks x motifs ()

if (is.null(rownames(motif_mat)) || is.null(colnames(motif_mat))) {
  stop("GetMotifData(..., slot='data')  motif ")
}

#  peaks  motif_mat 
peaks_pre_up <- intersect(peaks_pre_up, rownames(motif_mat))
peaks_con_up <- intersect(peaks_con_up, rownames(motif_mat))
keep_peaks   <- intersect(keep_peaks,   rownames(motif_mat))

summarize_motif_on_peakset <- function(peakset, set_name) {
  if (length(peakset) == 0) {
    return(tibble(set = set_name, motif = character(0)))
  }
  msub <- motif_mat[peakset, , drop = FALSE]
  tibble(
    set = set_name,
    motif = colnames(msub),
    n_peaks_in_set = length(peakset),
    n_peaks_with_motif = as.numeric(Matrix::colSums(msub > 0)),
    frac_peaks_with_motif = n_peaks_with_motif / n_peaks_in_set
  ) %>%
    arrange(desc(frac_peaks_with_motif), desc(n_peaks_with_motif))
}

motif_on_pre <- summarize_motif_on_peakset(peaks_pre_up, "PRE_up_MLpeaks")
motif_on_con <- summarize_motif_on_peakset(peaks_con_up, "CON_up_MLpeaks")
motif_on_all <- summarize_motif_on_peakset(keep_peaks,   "ALL_MLpeaks")

motif_map_pre <- motif_on_pre %>% left_join(motif_da2, by = "motif")
motif_map_con <- motif_on_con %>% left_join(motif_da2, by = "motif")
motif_map_all <- motif_on_all %>% left_join(motif_da2, by = "motif")

write_csv(motif_map_pre, file.path(outdir_map, "motif_on_MLpeaks_PRE_up.csv"))
write_csv(motif_map_con, file.path(outdir_map, "motif_on_MLpeaks_CON_up.csv"))
write_csv(motif_map_all, file.path(outdir_map, "motif_on_MLpeaks_ALL.csv"))

cat("[INFO] Saved motif-on-MLpeaks tables to:", outdir_map, "\n")



suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(stringr)
})

# =========================
# 0)  / 
# =========================
in_csv <- file.path(outdir_map, "motif_on_MLpeaks_ALL.csv")

outdir <- outdir_map
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# 
hit_frac_thr <- 0.20     # ML peaks ：>=20%
padj_thr     <- 0.05     # /（，）
lfc_thr      <- 0.25     # （）
top_n_each   <- 5      # PRE/CON Top N 


# =========================
# 1) 
# =========================
stopifnot(file.exists(in_csv))
motif_map0 <- read_csv(in_csv, show_col_types = FALSE)

need_cols <- c("motif","n_peaks_in_set","n_peaks_with_motif","frac_peaks_with_motif",
               "TF","chromvar_lfc","chromvar_padj","chromvar_log10padj")
miss <- setdiff(need_cols, colnames(motif_map0))
if (length(miss) > 0) stop("CSV：", paste(miss, collapse = ", "))

motif_df <- motif_map0 %>%
  dplyr::mutate(
    # ： set ； set， ALL 
    set = if ("set" %in% colnames(motif_map0)) as.character(set) else "ALL_MLpeaks",
    motif = as.character(motif),
    TF = as.character(TF),
    n_peaks_in_set = as.numeric(n_peaks_in_set),
    n_peaks_with_motif = as.numeric(n_peaks_with_motif),
    frac_peaks_with_motif = as.numeric(frac_peaks_with_motif),
    chromvar_lfc  = as.numeric(chromvar_lfc),
    chromvar_padj = as.numeric(chromvar_padj),
    chromvar_padj = ifelse(is.na(chromvar_padj) | chromvar_padj <= 0, 1e-300, chromvar_padj),
    chromvar_log10padj = -log10(chromvar_padj),
    label = ifelse(!is.na(TF) & TF != "", TF, motif)
  )

#  CSV  set， ALL
if ("set" %in% colnames(motif_map0)) {
  motif_df <- motif_df %>% dplyr::filter(str_detect(set, regex("ALL", ignore_case = TRUE)))
}

# （ motif ）
motif_df <- motif_df %>% dplyr::distinct(motif, .keep_all = TRUE)

cat("[INFO] motifs in input (ALL):", nrow(motif_df), "\n")

# =========================
# 2) ML peaks ： >= 10%
# =========================
motif_hit <- motif_df %>%
  dplyr::filter(!is.na(n_peaks_in_set), !is.na(n_peaks_with_motif)) %>%
  dplyr::filter(n_peaks_in_set > 0) %>%
  dplyr::filter(n_peaks_with_motif >= hit_frac_thr * n_peaks_in_set)

cat("[INFO] motifs passing hit filter:",
    nrow(motif_hit), " (hit_frac_thr=", hit_frac_thr, ")\n")

# =========================
# 3)  PRE / CON  Top motif（）
#     chromVAR  + |lfc|>=，“”
# =========================
motif_hit_sig <- motif_hit %>%
  dplyr::filter(chromvar_padj < padj_thr, abs(chromvar_lfc) >= lfc_thr)

cat("[INFO] motifs passing hit + chromVAR sig:",
    nrow(motif_hit_sig), " (padj<", padj_thr, ", |lfc|>=", lfc_thr, ")\n")

# =========================
# 3.0) ：“”
# =========================
if (!exists("top_n_each", inherits = TRUE) || length(top_n_each) != 1 || is.na(top_n_each)) {
  top_n_each <- 5
}
top_n_each <- as.integer(top_n_each)
if (!exists("padj_thr", inherits = TRUE) || length(padj_thr) != 1 || is.na(padj_thr)) {
  padj_thr <- 0.05
}
if (!exists("lfc_thr", inherits = TRUE) || length(lfc_thr) != 1 || is.na(lfc_thr)) {
  lfc_thr <- 0.25
}

top_pre <- motif_hit_sig %>%
  dplyr::filter(chromvar_lfc > 0) %>%
  dplyr::arrange(desc(chromvar_lfc)) %>%          # （PRE：PRE_vs_CON）
  dplyr::slice_head(n = top_n_each)

top_con <- motif_hit_sig %>%
  dplyr::filter(chromvar_lfc < 0) %>%
  dplyr::arrange(chromvar_lfc) %>%               # （CON：PRE_vs_CON）
  dplyr::slice_head(n = top_n_each)

top_show <- bind_rows(top_pre, top_con) %>%
  dplyr::distinct(motif, .keep_all = TRUE)


cat("[INFO] top motifs to show:", nrow(top_show), "\n")

write_csv(motif_hit,     file.path(outdir, "motifs_hitFilter_ge10pct.csv"))
write_csv(motif_hit_sig, file.path(outdir, "motifs_hitFilter_sig_chromvar.csv"))
write_csv(top_show, file.path(outdir, paste0("topMotifs_PREpos_CONneg_top", top_n_each, "each.csv")))

# =========================
# 4)  + （ chromvar.R ）
# =========================
volcano_data <- motif_hit %>%
  dplyr::mutate(
    significance = case_when(
      chromvar_padj < padj_thr & chromvar_lfc >  lfc_thr ~ "PRE upregulated",
      chromvar_padj < padj_thr & chromvar_lfc < -lfc_thr ~ "CON upregulated",
      TRUE                                               ~ "Not significant"
    )
  )

# （）
volcano_sig <- volcano_data %>% dplyr::filter(chromvar_padj < padj_thr, abs(chromvar_lfc) >= lfc_thr)
n_PRE_up <- sum(volcano_sig$chromvar_lfc > 0, na.rm = TRUE)
n_CON_up <- sum(volcano_sig$chromvar_lfc < 0, na.rm = TRUE)

# 
xr <- range(volcano_data$chromvar_lfc, na.rm = TRUE)
yr <- range(volcano_data$chromvar_log10padj, na.rm = TRUE)
x_left  <- xr[1] - 0.4* diff(xr)
x_right <- xr[2] + 0.3 * diff(xr)
y_top   <- yr[2] - 0.2  * diff(yr)

p_volcano <- ggplot(volcano_data, aes(x = chromvar_lfc, y = chromvar_log10padj)) +
  geom_point(aes(color = significance), alpha = 0.70, size = 2.4) +

  # “ Top 20”
  geom_text_repel(
    data = top_show,
    aes(x = chromvar_lfc, y = chromvar_log10padj, label = label),
    size = 5.2,
    max.overlaps = 100,
    box.padding = 0.55,
    point.padding = 0.25,
    segment.color = "gray40",
    fontface = "bold"
  ) +

  geom_hline(yintercept = -log10(padj_thr), linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed", color = "gray40") +

  annotate(
    "text", x = x_left, y = y_top,
    label = paste0(n_CON_up, " motifs"),
    color = "#82b487", size = 6, fontface = "bold", hjust = 0
  ) +
  annotate(
    "text", x = x_right, y = y_top,
    label = paste0(n_PRE_up, " motifs"),
    color = "#F0A780", size = 6, fontface = "bold", hjust = 1
  ) +

  scale_color_manual(values = c(
    "PRE upregulated" = "#F0A780",
    "CON upregulated" = "#82b487",
    "Not significant" = "gray75"
  )) +

  labs(
    title = "ML-peaks-supported motif activity (PRE vs CON)",
    x = expression(bold("log"[2]*"(Fold Change) ")),
    y = expression(bold("-log"[10]*"(p.adjust)")),
    color = NULL
  ) +
  theme_classic(base_size = 16) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = "bottom",
    legend.text = element_text(size = 15, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.text.x  = element_text(size = 12, face = "bold"),
    axis.text.y  = element_text(size = 12, face = "bold")
  )

out_png <- file.path(outdir, paste0("volcano_MLpeaksSupported_chromVAR_top", top_n_each, "each.png"))
ggsave(out_png, p_volcano, width = 7.6, height = 7.6, dpi = 500, bg = "white")
ggsave(file.path(outdir, paste0("volcano_MLpeaksSupported_chromVAR_top", top_n_each, "each.pdf")),
       p_volcano, width = 7.6, height = 7.6, dpi = 500, bg = "white")
# ，
write_csv(volcano_data, file.path(outdir, "volcano_data_MLpeaksSupported.csv"))

cat("[INFO] Done.\n")
cat("[INFO] Outdir:", outdir, "\n")
cat("[INFO] Volcano:", out_png, "\n")


# =========================
# 5) MA plot（ chromVAR deviation）
# =========================
# ： A= chromVAR deviation（PRE/CON ）， M=chromvar_lfc（PRE_vs_CON）
# ： motif_hit（ ML-peaks hit-filter）， TF 

top_n_each_side_ma <- 5   # PRE_up / CON_up  TF（ 5~15）

if (!exists("top_n_each_side_ma", inherits = TRUE) || length(top_n_each_side_ma) != 1 || is.na(top_n_each_side_ma)) {
  top_n_each_side_ma <- 5
}
top_n_each_side_ma <- as.integer(top_n_each_side_ma)

get_assay_data_compat <- function(obj, assay, layer = "data") {
  out <- tryCatch(
    SeuratObject::GetAssayData(obj, assay = assay, layer = layer),
    error = function(e) NULL
  )
  if (!is.null(out)) return(out)
  SeuratObject::GetAssayData(obj, assay = assay, slot = layer)
}

clean_tf_label <- function(x) {
  if (length(x) == 0) return(character())
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- trimws(x)
  x <- gsub("\\s+", " ", x)
  x <- gsub("\\(.*$", "", x)
  x <- gsub("::.*$", "", x)
  x <- gsub("_.*$", "", x)
  x <- trimws(x)
  x
}

# ----  chromvar assay / ident  ----
Seurat::DefaultAssay(acinar_cv) <- chromvar_assay

meta_cols <- colnames(acinar_cv@meta.data)
group_col <- if ("group" %in% meta_cols) {
  "group"
} else if ("condition" %in% meta_cols) {
  "condition"
} else if ("orig.ident" %in% meta_cols) {
  "orig.ident"
} else {
  NA_character_
}

if (is.na(group_col)) {
  stop(
    "acinar_cv （group/condition/orig.ident）\n",
    " meta columns: ", paste(meta_cols, collapse = ", ")
  )
}

grp_raw <- as.character(acinar_cv@meta.data[[group_col]])
if (all(is.na(grp_raw))) {
  stop("acinar_cv@meta.data[[", group_col, "]]  NA， PRE/CON  MA ")
}

lvl <- unique(na.omit(grp_raw))
lvl_upper <- toupper(lvl)

pick_level <- function(target) {
  target_u <- toupper(target)
  if (target_u %in% lvl_upper) return(lvl[match(target_u, lvl_upper)])
  hit <- which(grepl(target_u, lvl_upper, fixed = TRUE))
  if (length(hit) >= 1) return(lvl[hit[1]])
  NA_character_
}

level_pre <- pick_level("PRE")
level_con <- pick_level("CON")
if (is.na(level_pre) || is.na(level_con)) {
  stop(
    " '", group_col, "'  PRE/CON\n",
    " levels: ", paste(lvl, collapse = ", "), "\n",
    " meta.data  PRE  CON（）"
  )
}

Seurat::Idents(acinar_cv) <- factor(acinar_cv@meta.data[[group_col]], levels = c(level_pre, level_con))

# ----  chromVAR deviation  A ----
mat_chr <- get_assay_data_compat(acinar_cv, assay = chromvar_assay, layer = "data")
pre_cells <- Seurat::WhichCells(acinar_cv, idents = level_pre)
con_cells <- Seurat::WhichCells(acinar_cv, idents = level_con)

if (length(pre_cells) == 0 || length(con_cells) == 0) {
  stop(
    "：PRE=", length(pre_cells), ", CON=", length(con_cells), "\n",
    " Idents(acinar_cv)  meta.data[[", group_col, "]]"
  )
}

mean_pre <- Matrix::rowMeans(mat_chr[, pre_cells, drop = FALSE])
mean_con <- Matrix::rowMeans(mat_chr[, con_cells, drop = FALSE])

ma_df0 <- tibble(
  motif = names(mean_pre),
  A = (mean_pre + mean_con) / 2
)

# ---- （ motif_hit：ML-peaks hit-filter  motifs）----
stats_df <- motif_hit %>%
  dplyr::transmute(
    motif,
    lfc = chromvar_lfc,
    p_val_adj = ifelse(is.na(chromvar_padj), 1, chromvar_padj),
    TF = TF
  ) %>%
  dplyr::distinct(motif, .keep_all = TRUE)

ma_df <- ma_df0 %>%
  dplyr::left_join(stats_df, by = "motif") %>%
  dplyr::mutate(
    TF_pretty = clean_tf_label(TF),
    TF_label = TF_pretty,
    group = dplyr::case_when(
      p_val_adj < padj_thr & lfc >  lfc_thr  ~ "up",   # PRE_up (PRE_vs_CON)
      p_val_adj < padj_thr & lfc < -lfc_thr  ~ "down", # CON_up (PRE_vs_CON)
      TRUE                                   ~ "none"
    ),
    group = factor(group, levels = c("up", "none", "down"))
  )

cat(
  "[INFO] MA TF mapped:",
  sum(!is.na(ma_df$TF_label) & ma_df$TF_label != ""), "/", nrow(ma_df), "\n"
)

write_csv(ma_df, file.path(outdir, "MA_data_chromVAR_MLpeaksSupported.csv"))

# ----  top TF（， TF ）----
ma_sig <- ma_df %>% dplyr::filter(p_val_adj < padj_thr)
n_up   <- sum(ma_sig$group == "up", na.rm = TRUE)
n_down <- sum(ma_sig$group == "down", na.rm = TRUE)

top_up <- ma_sig %>%
  dplyr::filter(group == "up", !is.na(TF_label), TF_label != "") %>%
  dplyr::arrange(desc(lfc)) %>%
  dplyr::distinct(TF_label, .keep_all = TRUE) %>%
  dplyr::slice_head(n = top_n_each_side_ma)

top_dn <- ma_sig %>%
  dplyr::filter(group == "down", !is.na(TF_label), TF_label != "") %>%
  dplyr::arrange(lfc) %>%
  dplyr::distinct(TF_label, .keep_all = TRUE) %>%
  dplyr::slice_head(n = top_n_each_side_ma)

sig <- dplyr::bind_rows(top_up, top_dn) %>% dplyr::distinct(motif, .keep_all = TRUE)

write_csv(top_up, file.path(outdir, paste0("MA_topTF_PRE_up_top", top_n_each_side_ma, ".csv")))
write_csv(top_dn, file.path(outdir, paste0("MA_topTF_CON_up_top", top_n_each_side_ma, ".csv")))



mycol_barbie <- c(
  up   = "#F0A780",
  none = "#b0b0b0",
  down = "#82b487"
)

p_ma <- ggplot() +
  geom_point(
    data = ma_df %>% dplyr::filter(group == "none"),
    aes(x = A, y = lfc, color = group),
    size = 1.6, alpha = 0.70
  ) +
  geom_point(
    data = ma_df %>% dplyr::filter(group == "up"),
    aes(x = A, y = lfc, color = group),
    size = 1.8, alpha = 0.75
  ) +
  geom_point(
    data = ma_df %>% dplyr::filter(group == "down"),
    aes(x = A, y = lfc, color = group),
    size = 1.8, alpha = 0.75
  ) +
  scale_color_manual(values = mycol_barbie, name = NULL) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  labs(
    title = NULL,
    x = "A: Mean chromVAR deviation",
    y = "M: avg_log2FC"
  ) +
  theme_classic(base_size = 16) +
  theme(
    text = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    legend.position = "bottom",
    legend.text = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold")
  )

# ----  ----
xr2 <- range(ma_df$A, na.rm = TRUE)
yr2 <- range(ma_df$lfc, na.rm = TRUE)
x_right2 <- xr2[2] - 0.7 * diff(xr2)
y_top2    <- yr2[2] - 0.3 * diff(yr2)
y_bottom2 <- yr2[1] + 50* diff(yr2)
if (is.finite(y_bottom2) && y_bottom2 > 0) y_bottom2 <- -0.1

x_inset <- 0.15 * diff(xr2)
x_target_right <- xr2[2] - x_inset
x_target_left  <- xr2[1] + x_inset

p_ma_labeled <- p_ma +
  annotate(
    "text",
    x = x_right2,
    y = y_top2,
    label = paste0(n_up, " motifs"),
    color = mycol_barbie[["up"]],
    size = 5.5, fontface = "bold",
    hjust = 1
  ) +
  annotate(
    "text",
    x = x_right2,
    y = -0.5,
    label = paste0(n_down, " motifs"),
    color = mycol_barbie[["down"]],
    size = 5.5, fontface = "bold",
    hjust = 1
  ) +
  geom_point(
    data = sig,
    aes(x = A, y = lfc, color = group),
    size = 3.2, alpha = 0.20
  )

if (nrow(top_up) > 0) {
  p_ma_labeled <- p_ma_labeled +
    ggrepel::geom_text_repel(
      data = top_up,
      aes(x = A, y = lfc, label = TF_label),
      seed = 233,
      size = 3.8,
      color = "black",
      fontface = "bold",
      min.segment.length = 0,
      force = 2,
      force_pull = 2,
      box.padding = 0.15,
      point.padding = 0.15,
      max.overlaps = Inf,
      segment.linetype = 3,
      segment.color = "black",
      segment.alpha = 0.5,
      nudge_x = x_target_right - top_up$A,
      direction = "y",
      hjust = 0
    )
}

if (nrow(top_dn) > 0) {
  p_ma_labeled <- p_ma_labeled +
    ggrepel::geom_text_repel(
      data = top_dn,
      aes(x = A, y = lfc, label = TF_label),
      seed = 233,
      size = 3.8,
      color = "black",
      fontface = "bold",
      min.segment.length = 0,
      force = 2,
      force_pull = 2,
      box.padding = 0.15,
      point.padding = 0.15,
      max.overlaps = Inf,
      segment.linetype = 3,
      segment.color = "black",
      segment.alpha = 0.5,
      nudge_x = x_target_left - top_dn$A,
      direction = "y",
      hjust = 1
    )
}

out_ma_png <- file.path(outdir, paste0("chromvar_MA_PRE_vs_CON_barbie_labeled_TF_top", top_n_each_side_ma, "each.png"))
out_ma_pdf <- file.path(outdir, paste0("chromvar_MA_PRE_vs_CON_barbie_labeled_TF_top", top_n_each_side_ma, "each.pdf"))
ggsave(out_ma_png, p_ma_labeled, width = 7, height = 6.8, dpi = 500, bg = "white")
ggsave(out_ma_pdf, p_ma_labeled, width = 6.2, height = 6.2, dpi = 500, bg = "white")
cat("[INFO] Saved MA plot: ", out_ma_png, "\n", sep = "")
