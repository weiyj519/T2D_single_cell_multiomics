#!/usr/bin/env Rscript

# Script name: 00_run_chromvar_motif_activity.R
# Purpose: Cleaned public workflow script for the T2D single-cell multiome analysis.
# Input: Paths and analysis settings are read from configs/config.yaml.
# Output: Module-specific outputs are written under the configured results directory.
# Main steps: Load configuration, run the original analysis logic, and export reproducible outputs.
# Example command: Rscript scripts/06_atac_motif_peak_gene/00_run_chromvar_motif_activity.R

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
# chromVAR on acinar cells: CON vs PRE (hg19)
# Input: data/processed/combined_final_celltyped.qs
# Output: motif activity (deviation) + differential motifs
# =========================================================

Sys.setenv(
  OMP_NUM_THREADS = 4,
  OPENBLAS_NUM_THREADS = 4,
  MKL_NUM_THREADS = 4,
  VECLIB_MAXIMUM_THREADS = 4,
  NUMEXPR_NUM_THREADS = 4
)

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(tibble)
  library(ggplot2)

  library(chromVAR)
  library(motifmatchr)
  library(JASPAR2020)
  library(TFBSTools)
  library(BSgenome.Hsapiens.UCSC.hg19)

  library(GenomeInfoDb)
})

qs_path <- seurat_object_path
outdir  <- "results/downstream/acinar_CON_PRE/chromvar"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("[INFO] Loading:", qs_path, "\n")
obj <- qread(qs_path)

# ---  ---
assays_obj <- Assays(obj)
assay_names <- if (is.character(assays_obj)) assays_obj else names(assays_obj)
atac_assay <- if ("ATAC" %in% assay_names) {
  "ATAC"
} else if ("peaks" %in% assay_names) {
  "peaks"
} else {
  stop(sprintf("Cannot find ATAC assay. Available assays: %s", paste(assay_names, collapse = ", ")))
}
stopifnot(all(c("celltype", "group") %in% colnames(obj@meta.data)))

# Seurat v5:  layers， GetAssayData/RunChromVAR  counts
obj <- JoinLayers(obj)

# ---  acinar cell， CON vs PRE ---
obj$celltype <- as.character(obj$celltype)
obj$group    <- as.character(obj$group)

# ---  acinar cell， CON vs PRE ---
acinar_celltype_candidates <- c(
  "Acinar cell", "Acinar_cell", "acinar cell", "acinar_cell",
  "ACINAR cell", "ACINAR_cell"
)

acinar <- subset(
  obj,
  subset = celltype %in% acinar_celltype_candidates & group %in% c("CON", "PRE")
)
DefaultAssay(acinar) <- atac_assay

cat("[INFO] acinar cells:", ncol(acinar), "\n")
print(table(acinar$group))

if (ncol(acinar) == 0) {
  stop(
    "No cells found for acinar CON/PRE. ",
    "Tried celltype candidates: ", paste(acinar_celltype_candidates, collapse = ", "), ". ",
    "Check obj@meta.data$celltype and obj@meta.data$group values."
  )
}

# --- ： peak （ GL000192.1 not found） ---
gr <- granges(acinar[[atac_assay]])

#  GRanges  names（ names(gr_std) ，subset(features=...) ）
if (is.null(names(gr)) || length(names(gr)) == 0 || all(is.na(names(gr)))) {
  names(gr) <- rownames(acinar[[atac_assay]])
}

try({
  seqlevelsStyle(gr) <- "UCSC"  # chr1/chr2...
}, silent = TRUE)

# 1) 
gr_std <- keepStandardChromosomes(gr, pruning.mode = "coarse")

# 2)  hg19 BSgenome  seqlevels
keep_seq <- intersect(
  GenomeInfoDb::seqlevels(gr_std),
  GenomeInfoDb::seqlevels(BSgenome.Hsapiens.UCSC.hg19)
)
gr_std <- GenomeInfoDb::keepSeqlevels(gr_std, keep_seq, pruning.mode = "coarse")

keep_peaks <- names(gr_std)

# 3) ： GRanges names ， peak （）
if (is.null(keep_peaks) || length(keep_peaks) == 0) {
  peak_names <- rownames(acinar[[atac_assay]])
  keep_peaks <- peak_names[grepl("^chr([0-9]+|X|Y|M)[:_-]", peak_names)]
}

# 4) ： keep_peaks  ATAC assay  feature 
features_all <- rownames(acinar[[atac_assay]])
keep_peaks_raw <- keep_peaks
keep_peaks <- intersect(keep_peaks_raw, features_all)

if (length(keep_peaks) == 0 && length(keep_peaks_raw) > 0) {
  #  dash/colon/underscore 
  norm <- function(x) gsub(":", "-", gsub("_", "-", as.character(x)))
  features_norm <- norm(features_all)
  keep_norm <- norm(keep_peaks_raw)
  idx <- match(keep_norm, features_norm)
  idx <- idx[!is.na(idx)]
  keep_peaks <- unique(features_all[idx])
}

if (length(keep_peaks) == 0) {
  stop(
    sprintf(
      "Peakfeature ATAC assay\nfeatures: %s\nkeep_peaks_raw: %s",
      paste(head(features_all, 5), collapse = ", "),
      paste(head(keep_peaks_raw, 5), collapse = ", ")
    )
  )
}

cat(sprintf(
  "[INFO] Peaks before: %d | after filtering to hg19 standard: %d\n",
  length(gr), length(keep_peaks)
))

# ：Seurat::subset(features=...)  assay  feature
#  multi-assay （ RNA/ATAC），peak  RNA assay 
# “None of the features provided found in this assay” ATAC assay
acinar[[atac_assay]] <- subset(acinar[[atac_assay]], features = keep_peaks)

# ---  motif（JASPAR2020 human CORE；hg19）---
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = "vertebrates", species = 9606, all_versions = FALSE)
)

cat("[INFO] Adding motifs (JASPAR2020 CORE human) ...\n")
acinar <- AddMotifs(
  object = acinar,
  genome = BSgenome.Hsapiens.UCSC.hg19,
  pfm    = pfm
)

# ---  motif(ID) -> TF(name) （）---
# ：chromVAR/Signac  motif feature  JASPAR ID（ MA0139.1）
#  MAxxxx， ID  TF 
pfm_map <- tibble::tibble(
  motif = as.character(TFBSTools::ID(pfm)),
  TF = as.character(TFBSTools::name(pfm))
) %>%
  dplyr::filter(!is.na(motif), motif != "") %>%
  dplyr::distinct(motif, .keep_all = TRUE)

#  TF ：
# - FOSL2::JUNB -> FOSL2
# - JUN(var.2)  -> JUN
clean_tf_label <- function(x) {
  x <- as.character(x)
  x <- gsub("\\(var\\.[^\\)]*\\)", "", x)
  x <- sub("^(.*?)::.*$", "\\1", x)
  trimws(x)
}
 
# ---  chromVAR deviation ---
cat("[INFO] Running chromVAR ...\n")

#  Signac  RunChromVAR() （ slot/new.assay.name）
# ， unused argument 
assays_before <- Assays(acinar)
rc_args <- names(formals(Signac::RunChromVAR))
rc_call <- list(
  object = acinar,
  genome = BSgenome.Hsapiens.UCSC.hg19,
  assay  = atac_assay
)
if ("slot" %in% rc_args) {
  rc_call$slot <- "counts"
}
if ("new.assay.name" %in% rc_args) {
  rc_call$`new.assay.name` <- "chromvar"
}

acinar <- do.call(Signac::RunChromVAR, rc_call)

#  new.assay.name ， assay  chromvar
assays_after <- Assays(acinar)
if (!"chromvar" %in% assays_after) {
  new_assay <- setdiff(assays_after, assays_before)
  if (length(new_assay) == 1) {
    rename_map <- setNames(list("chromvar"), new_assay)
    acinar <- do.call(Seurat::RenameAssays, c(list(acinar), rename_map))
  }
}

print(Assays(acinar))


#  chromvar assay （）
qsave(acinar, file.path(outdir, "acinar_CON_vs_PRE_chromvar.qs"))
cat("[INFO] Saved:", file.path(outdir, "acinar_CON_vs_PRE_chromvar.qs"), "\n")

# ---  motif ：PRE vs CON（ logFC>0  PRE ） ---
DefaultAssay(acinar) <- "chromvar"
Idents(acinar) <- factor(acinar$group, levels = c("CON", "PRE"))

da <- FindMarkers(
  object = acinar,
  ident.1 = "PRE",
  ident.2 = "CON",
  test.use = "wilcox",
  min.pct = 0.05,
  logfc.threshold = 0
) %>% rownames_to_column("motif")

#  Seurat  logFC 
lfc_col <- dplyr::case_when(
  "avg_log2FC" %in% colnames(da) ~ "avg_log2FC",
  "avg_logFC"  %in% colnames(da) ~ "avg_logFC",
  TRUE ~ NA_character_
)
if (is.na(lfc_col)) stop("Cannot find avg_log2FC/avg_logFC in FindMarkers result.")

#  TF ： pfm_map（JASPAR ID -> TF name）， motif.names
da <- da %>%
  dplyr::left_join(pfm_map, by = "motif")

motif_names <- tryCatch(
  GetMotifData(acinar, assay = atac_assay, slot = "motif.names"),
  error = function(e) NULL
)
if (!is.null(motif_names)) {
  tf2 <- as.character(motif_names[da$motif])
  da$TF <- dplyr::coalesce(da$TF, tf2)
}

# ： TF ， pfm_map  match （ join  NA）
da$TF <- dplyr::coalesce(da$TF, pfm_map$TF[match(da$motif, pfm_map$motif)])

cat(sprintf(
  "[INFO] Motif->TF mapped: %d/%d (%.1f%%)\n",
  sum(!is.na(da$TF) & da$TF != ""),
  nrow(da),
  100 * mean(!is.na(da$TF) & da$TF != "")
))

da <- da %>% arrange(p_val_adj, p_val)
readr::write_csv(da, file.path(outdir, "chromvar_PRE_vs_CON_motifs.csv"))
cat("[INFO] Saved:", file.path(outdir, "chromvar_PRE_vs_CON_motifs.csv"), "\n")



# =========================================================
# Volcano plot for differential chromVAR motif activity
# (PRE vs CON, acinar cells)
# Input: da (FindMarkers result) already computed above
# Output: chromvar_volcano_PRE_vs_CON.png
# =========================================================

suppressPackageStartupMessages({
  library(ggrepel)
})

# ---- 0) ： ----
padj_thr <- 0.05
lfc_thr  <- 0.5          # avg_log2FC 
top_n_each_side <- 5    # /

# ---- 1)  ----
volcano_data <- da %>%
  mutate(
    lfc = .data[[lfc_col]],
    p_val_adj = ifelse(is.na(p_val_adj), 1, p_val_adj),
    log10_padj = -log10(pmax(p_val_adj, 1e-300)),
    TF_pretty = clean_tf_label(TF),
    label = ifelse(is.na(TF_pretty) | TF_pretty == "", motif, TF_pretty),
    significance = case_when(
      p_val_adj < padj_thr & lfc >  lfc_thr  ~ "PRE upregulated",
      p_val_adj < padj_thr & lfc < -lfc_thr  ~ "CON upregulated",
      TRUE                                   ~ "Not significant"
    )
  )

# ---- 2)  &  top motifs ----
volcano_sig <- volcano_data %>% filter(p_val_adj < padj_thr)

n_PRE_up <- sum(volcano_sig$significance == "PRE upregulated", na.rm = TRUE)
n_CON_up <- sum(volcano_sig$significance == "CON upregulated", na.rm = TRUE)

top_up <- volcano_sig %>%
  filter(lfc > 0) %>%
  arrange(desc(lfc)) %>%
  slice_head(n = top_n_each_side)

top_dn <- volcano_sig %>%
  filter(lfc < 0) %>%
  arrange(lfc) %>%
  slice_head(n = top_n_each_side)

top_motifs <- bind_rows(top_up, top_dn) %>%
  distinct(motif, .keep_all = TRUE)

# ---- 3) （）----
xr <- range(volcano_data$lfc, na.rm = TRUE)
yr <- range(volcano_data$log10_padj, na.rm = TRUE)
x_left  <- xr[1] + 0.02 * diff(xr)
x_right <- xr[2] - 0.02 * diff(xr)
y_top   <- yr[2] - 0.1 * diff(yr)

# ---- 4)  ----
p_volcano <- ggplot(volcano_data, aes(x = lfc, y = log10_padj)) +
  geom_point(aes(color = significance), alpha = 0.65, size = 2.5) +

  #  top motifs（，）
  geom_text_repel(
    data = top_motifs,
    aes(label = label),
    size = 5,
    max.overlaps = 30,
    box.padding = 0.5,
    point.padding = 0.2,
    segment.color = "gray40",
    fontface = "bold"
  ) +

  geom_hline(yintercept = -log10(padj_thr), linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed", color = "gray40") +

  annotate(
    "text",
    x = x_left,
    y = y_top,
    label = paste0(n_CON_up, " motifs"),
    color = "#4DBBD5",
    size = 6,
    fontface = "bold",
    hjust = 0
  ) +
  annotate(
    "text",
    x = x_right,
    y = y_top,
    label = paste0(n_PRE_up, " motifs"),
    color = "#E64B35",
    size = 6,
    fontface = "bold",
    hjust = 1
  ) +

  scale_color_manual(
    values = c(
      "CON upregulated" = "#4DBBD5",
      "PRE upregulated" = "#E64B35",
      "Not significant" = "gray75"
    )
  ) +

  labs(
    title = "Differential Motif Activity between PRE and CON (Acinar cells, chromVAR)",
    x = expression("Average log"[2]*"(Fold Change)"),
    y = expression("-log"[10]*"(p.adjust)"),
    color = NULL
  ) +
  theme_classic(base_size = 16) +
  theme(
    text = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = "top",
    legend.text = element_text(size = 16, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 13, face = "bold"),
    axis.text.x  = element_text(size = 13, face = "bold"),
    axis.text.y  = element_text(size = 13, face = "bold")
  )

print(p_volcano)

out_png <- file.path(outdir, "chromvar_volcano_PRE_vs_CON.png")
ggsave(out_png, p_volcano, width = 10, height = 8, dpi = 500, bg = "white")
cat("[INFO] Saved:", out_png, "\n")

# ：，
readr::write_csv(volcano_data, file.path(outdir, "chromvar_volcano_data.csv"))




# =========================================================
# MA plot (Barbie palette) for chromVAR differential motifs
# (PRE vs CON, acinar cells)
# =========================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
})

# ---- （）----
top_n_each_side_ma <- 5   # up/down  TF（ 5~15）
# padj_thr / lfc_thr / lfc_col / outdir / clean_tf_label / da / acinar 

# ----  assay  ident  ----
DefaultAssay(acinar) <- "chromvar"
Idents(acinar) <- factor(acinar$group, levels = c("CON", "PRE"))

# ---- 1)  chromVAR deviation  A （CON/PRE ）----
mat_chr <- SeuratObject::GetAssayData(acinar, assay = "chromvar", layer = "data")  # features x cells

con_cells <- Seurat::WhichCells(acinar, idents = "CON")
pre_cells <- Seurat::WhichCells(acinar, idents = "PRE")

mean_con <- Matrix::rowMeans(mat_chr[, con_cells, drop = FALSE])
mean_pre <- Matrix::rowMeans(mat_chr[, pre_cells, drop = FALSE])

ma_df <- tibble(
  motif = names(mean_pre),
  A = (mean_con + mean_pre) / 2
)

# ---- 2)  FindMarkers （lfc/padj/TF）----
stats_df <- da %>%
  transmute(
    motif,
    lfc = .data[[lfc_col]],
    p_val_adj = ifelse(is.na(p_val_adj), 1, p_val_adj),
    TF = TF
  ) %>%
  distinct(motif, .keep_all = TRUE)

ma_df <- ma_df %>%
  left_join(stats_df, by = "motif") %>%
  mutate(
    TF_pretty = clean_tf_label(TF),
    TF_label  = TF_pretty,   # ： TF（ motif/MAxxxx）

    group = case_when(
      p_val_adj < padj_thr & lfc >  lfc_thr  ~ "up",
      p_val_adj < padj_thr & lfc < -lfc_thr  ~ "down",
      TRUE                                   ~ "none"
    ),
    group = factor(group, levels = c("up", "none", "down"))
  )

cat("[INFO] MA TF mapped:",
    sum(!is.na(ma_df$TF_label) & ma_df$TF_label != ""), "/", nrow(ma_df), "\n")

# ---- 3)  &  top TF（， TF ）----
ma_sig <- ma_df %>% filter(p_val_adj < padj_thr)

n_up   <- sum(ma_sig$group == "up",   na.rm = TRUE)
n_down <- sum(ma_sig$group == "down", na.rm = TRUE)

top_up <- ma_sig %>%
  filter(group == "up", !is.na(TF_label), TF_label != "") %>%
  arrange(desc(lfc)) %>%                            # up：lfc 
  distinct(TF_label, .keep_all = TRUE) %>%          # TF 
  slice_head(n = top_n_each_side_ma)

top_dn <- ma_sig %>%
  filter(group == "down", !is.na(TF_label), TF_label != "") %>%
  arrange(lfc) %>%                                  # down：lfc 
  distinct(TF_label, .keep_all = TRUE) %>%
  slice_head(n = top_n_each_side_ma)

sig <- bind_rows(top_up, top_dn)

# ---- 4) Barbie ：up / none / down ----
mycol_barbie <- c(
  up   = "#F147B1",
  none = "#fbd1e5",
  down = "#399AEB"
)

# ---- 5) MA （，）----
p_ma <- ggplot() +
  geom_point(
    data = ma_df %>% filter(group == "none"),
    aes(x = A, y = lfc, color = group),
    size = 1.6, alpha = 0.70
  ) +
  geom_point(
    data = ma_df %>% filter(group == "up"),
    aes(x = A, y = lfc, color = group),
    size = 1.8, alpha = 0.75
  ) +
  geom_point(
    data = ma_df %>% filter(group == "down"),
    aes(x = A, y = lfc, color = group),
    size = 1.8, alpha = 0.75
  ) +
  scale_color_manual(values = mycol_barbie, name = NULL) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  labs(
    title = "MA Plot of Differential Motif Activity",
    x = "A: Mean chromVAR deviation",
    y = "M: avg_log2FC"
  ) +
  theme_classic(base_size = 16) +
  theme(
    text = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    legend.position = "top",
    legend.text = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold")
  )

# ---- 6) （，）----
# ---- 6)  ----
xr <- range(ma_df$A, na.rm = TRUE)
yr <- range(ma_df$lfc, na.rm = TRUE)

x_left  <- xr[1] + 0.02 * diff(xr)
x_right <- xr[2] - 0.02 * diff(xr)

y_top    <- yr[2] - 0.04 * diff(yr)          # （ up ）
y_bottom <- yr[1] + 0.08 * diff(yr)          # （ down ）
if (y_bottom > 0) y_bottom <- -0.1           # ：

# ---- 7)  TF （up down ）----
# ---- 3) ： x， y ----
x_rng  <- xr
# ： x “”，
x_inset <- 0.15 * diff(x_rng)   # （：0.02~0.08）
x_target_right <- x_rng[2] - x_inset
x_target_left  <- x_rng[1] + x_inset

p_ma_labeled <- p_ma +
  annotate(
    "text",
    x = x_right,
    y = y_top,
    label = paste0(n_up, " motifs"),
    color = mycol_barbie[["up"]],
    size = 6, fontface = "bold",
    hjust = 1
  ) +
  annotate(
    "text",
    x = x_right,            # ✅  up  x
    y = y_bottom,           # ✅ （）
    label = paste0(n_down, " motifs"),
    color = mycol_barbie[["down"]],
    size = 6, fontface = "bold",
    hjust = 1
  ) +
  # 
  geom_point(
    data = sig,
    aes(x = A, y = lfc, color = group),
    size = 3.2, alpha = 0.20
  ) +
  # up ：
  geom_text_repel(
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
  ) +
  # down ：
  geom_text_repel(
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

print(p_ma_labeled)


# ---- 8)  ----
out_ma_png <- file.path(outdir, "chromvar_MA_PRE_vs_CON_barbie_labeled_TF.png")
ggsave(out_ma_png, p_ma_labeled, width = 6, height = 6, dpi = 500, bg = "white")
cat("[INFO] Saved:", out_ma_png, "\n")

# ： MA 
readr::write_csv(ma_df, file.path(outdir, "chromvar_MA_data.csv"))