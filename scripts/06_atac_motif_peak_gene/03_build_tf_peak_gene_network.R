#!/usr/bin/env Rscript

# Script name: 03_build_tf_peak_gene_network.R
# Purpose: Cleaned public workflow script for the T2D single-cell multiome analysis.
# Input: Paths and analysis settings are read from configs/config.yaml.
# Output: Module-specific outputs are written under the configured results directory.
# Main steps: Load configuration, run the original analysis logic, and export reproducible outputs.
# Example command: Rscript scripts/06_atac_motif_peak_gene/03_build_tf_peak_gene_network.R

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

############################################################
## netlink_acinar_CON_PRE_MLpeaks_TFmotif_peak_gene.R
## Build TF(motif) -> ML peaks -> linked genes network
## (acinar cells, CON vs PRE)
############################################################

# （）：
# 1)  ATAC  peak “ML peaks”（ peaks）， peak （CON_up/PRE_up）
# 2)  chromVAR  motif  + motif （peak×motif） TF(motif)
# 3) “TF(motif)  peak” peak（TF -> peak）
# 4)  LinkPeaks  peak-gene  peak  gene（peak -> gene）
# 5)  TF -> peak -> gene ， strength
# 6)  RNA DE ，“” TF -> peak -> gene （strength_rna）
# 7)  CSV，（CON_up  PRE_up）

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(Matrix)
  library(ggplot2)
  library(ggalluvial)
})

# =========================
# 0) Paths ()
# =========================
# ：/， Rscript 
# 
project_dir <- "results/downstream/acinar_CON_PRE"
celltype <- "acinar"

# ：/：
# - peaks:    featurePeaks（ CON_vs_PRE）
# - motifs:   chromVAR motif  PRE_vs_CON（ peaks/RNA ）
# - RNA:      FindMarkers(ident.1=CON, ident.2=PRE)， CON_vs_PRE
contrast_peak  <- "CON_vs_PRE"   # DA peaks + chromVAR 
contrast_motif <- "PRE_vs_CON"   # chromVAR motif DA 
contrast_rna   <- "CON_vs_PRE"   # RNA DE 

# （ TF/peak/gene ）
con_dir <- "CON_up"
pre_dir <- "PRE_up"

da_csv <- file.path(project_dir, "motif", sprintf("DA_featurePeaks_%s_%s.csv", celltype, contrast_peak), fsep = "/")

# chromVAR  motif DA 
cv_qs        <- file.path(project_dir, "chromvar", sprintf("%s_%s_chromvar.qs", celltype, contrast_peak), fsep = "/")
motif_da_csv <- file.path(project_dir, "chromvar", sprintf("chromvar_%s_motifs.csv", contrast_motif), fsep = "/")

#  MLpeaks-supported + chromVAR  motif （ map_chromvar_to_MLpeaks.R ）
# ：，“ motif”（ topN ）
mlpeaks_supported_motif_csv <- file.path(project_dir, "mlpeaks_chromvar_mapping", "motifs_hitFilter_sig_chromvar.csv", fsep = "/")

# LinkPeaks links
links_best_csv <- file.path(project_dir, "linkpeak", "Acinar_cell__CON_vs_PRE__LinkPeaks_hg19__featurePeaks_allGenes_SIG_fdr0.1.csv", fsep = "/")

# RNA DE（acinar ）
rna_deseq_csv <- file.path(project_dir, "GO_RNA", sprintf("DE_%s_%s_RNA.csv", celltype, contrast_rna), fsep = "/")

outdir <- file.path(project_dir, "netlink_mlpeak", fsep = "/")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# =========================
# 1) Parameters（）
# =========================
# coverage_frac/padj_cut_motif/top_n_each_dir “ motif/TF ”
# direction_match  TF(motif)  peak （，）
# keep_pos_links  peak->gene （ enhancer ；）
# max_edges_plot  TF-peak-gene triplet（）
coverage_frac     <- 0.10  # motif  ML peaks （10%）
padj_cut_motif    <- 0.05
top_n_each_dir    <- 50    # PRE/CONTop N（ chromvar_lfc ）
direction_match   <- TRUE  #  motif  peak （ TRUE）
dedup_by_TF_gene  <- TRUE  # /TF motif（ TRUE）
keep_pos_links    <- TRUE  #  peak-gene link（ enhancer ）
max_edges_plot    <- 10000    #  triplet（）
link_fdr_cut      <- 0.30  #  links  pass_sig ， FDR 

# RNA （）： RNA 
# ： RNA DE ， p_val_adj  1（）
#  RNA  gene， require_rna_padj  TRUE
require_rna_padj  <- FALSE
rna_padj_cut      <- 0.05



# =========================
# 2) Helpers
# =========================
#  peak ： peak  ':' '_' ， '-'， join/intersect 
norm_peak <- function(x) {
  x %>%
    gsub(":", "-", .) %>%
    gsub("_", "-", .) %>%
    gsub("\\s+", "", .)
}

#  "FOSL1::JUNB" / "NEUROG2(var.2)" “TF”
# ：（dedup_by_TF_gene） TF_gene ， TF  motif 
simplify_tf_gene <- function(tf) {
  tf2 <- tf %>% gsub("\\(var\\.[^\\)]*\\)", "", .)
  tf2 <- str_split_fixed(tf2, "::", 2)[,1]
  tf2 <- trimws(tf2)
  tf2
}

# Seurat/SeuratObject ：Assays(obj)  list(assay objects)
get_assay_names <- function(obj) {
  a <- tryCatch(Assays(obj), error = function(e) NULL)
  if (is.null(a)) return(names(obj@assays))
  if (is.character(a)) return(a)
  nm <- names(a)
  if (!is.null(nm)) return(nm)
  names(obj@assays)
}

detect_atac_assay <- function(obj) {
  assay_names <- get_assay_names(obj)
  prefer <- intersect(c("ATAC", "peaks", "Peaks"), assay_names)
  candidates <- unique(c(prefer, assay_names))
  for (nm in candidates) {
    ok <- tryCatch(!is.null(obj[[nm]]@motifs), error = function(e) FALSE)
    if (isTRUE(ok)) return(nm)
  }
  NA_character_
}

# =========================
# 3) Read ML peaks (DA + SHAP)
# =========================
# ： ATAC  peak ， peaks（keep_peaks）， peak 
stopifnot(file.exists(da_csv))
da_peaks <- read_csv(da_csv, show_col_types = FALSE)

if (!"peak" %in% colnames(da_peaks)) {
  stop("DA CSV  peak  peak")
}

da_peaks <- da_peaks %>%
  mutate(peak = norm_peak(peak))

#  logFC 
#  logFC  lfc  avg_log2FC，
lfc_col_peak <- if ("lfc" %in% colnames(da_peaks)) "lfc" else
  if ("avg_log2FC" %in% colnames(da_peaks)) "avg_log2FC" else NA_character_
if (is.na(lfc_col_peak)) stop("DA peaks CSV  lfc/avg_log2FC")

da_peaks <- da_peaks %>%
  mutate(peak_lfc = as.numeric(.data[[lfc_col_peak]]))

# peak（ significance；）
# peak_direction  direction_match： TF(motif)  peak 
da_peaks_dir <- da_peaks %>%
  mutate(
    peak_direction = case_when(
      "significance" %in% colnames(da_peaks) & grepl("up in\\s*con", tolower(significance)) ~ con_dir,
      "significance" %in% colnames(da_peaks) & grepl("up in\\s*pre", tolower(significance)) ~ pre_dir,
      #  significance ： featurePeaks （avg_log2FC  PRE - CON）
      peak_lfc > 0 ~ pre_dir,
      peak_lfc < 0 ~ con_dir,
      TRUE ~ "not_sig"
    )
  ) %>%
  filter(peak_direction != "not_sig")

keep_peaks <- unique(da_peaks_dir$peak)
cat("[INFO] ML peaks used:", length(keep_peaks), "\n")

# =========================
# 4) Read chromVAR object + motif DA table
# =========================
# ：
# -  chromVAR/ATAC  Seurat （cv_obj）， motif  motif_mat（peaks × motifs）
# -  motif  motif_da， chromvar_lfc/padj， motif 
stopifnot(file.exists(cv_qs), file.exists(motif_da_csv))
cv_obj <- qread(cv_qs)
cv_obj <- tryCatch(JoinLayers(cv_obj), error = function(e) cv_obj)

assay_names <- get_assay_names(cv_obj)
atac_assay <- detect_atac_assay(cv_obj)
if (is.na(atac_assay)) {
  stop(" motifs  ATAC/peaks assay assays: ", paste(assay_names, collapse = ", "))
}
if (!"chromvar" %in% assay_names) {
  stop(" 'chromvar' assay assays: ", paste(assay_names, collapse = ", "))
}

motif_da <- read_csv(motif_da_csv, show_col_types = FALSE)
# motif  logFC ， avg_log2FC/avg_logFC
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
  dplyr::select(motif, TF, chromvar_lfc, chromvar_padj, chromvar_log10padj) %>%
  distinct(motif, .keep_all = TRUE)

# motif matrix（peaks x motifs）
# motif_mat：=peak，=motif， motif /（ >0 “”）
if (is.null(cv_obj[[atac_assay]]@motifs)) {
  stop("cv_obj[[atac_assay]]@motifs ： motif （ AddMotifs ）")
}
motif_mat <- GetMotifData(cv_obj, assay = atac_assay, slot = "data")
if (is.null(rownames(motif_mat)) || is.null(colnames(motif_mat))) {
  stop("motif_mat  rownames/colnames")
}

# =========================
# 5) Intersect peaks + compute motif coverage on ML peaks
# =========================
# ： ML peaks  motif （ ML peaks  motif）
# ：coverage_frac “ peaks” motif，
keep_peaks <- intersect(keep_peaks, rownames(motif_mat))
cat("[INFO] ML peaks found in motif matrix:", length(keep_peaks), "\n")
if (length(keep_peaks) == 0) stop("keep_peaks  motif_mat ：peak")

msub_all <- motif_mat[keep_peaks, , drop = FALSE]

motif_cov <- tibble(
  set = "ALL_MLpeaks",
  motif = colnames(msub_all),
  n_peaks_in_set = length(keep_peaks),
  n_peaks_with_motif = as.numeric(Matrix::colSums(msub_all > 0)),
  frac_peaks_with_motif = n_peaks_with_motif / n_peaks_in_set
) %>%
  left_join(motif_da2, by = "motif") %>%
  filter(!is.na(TF)) %>%
  mutate(TF_gene = simplify_tf_gene(TF))

write_csv(motif_cov, file.path(outdir, "motif_coverage_on_ALL_MLpeaks.csv"))

# =========================
# 6) Select motifs for interpretation (coverage + chromVAR sig + TopN each direction)
# =========================
# ： motif/TF 
# -  MLpeaks-supported  motif （“ motif”）
# - ： stop（ coverage/topN ）
if (!file.exists(mlpeaks_supported_motif_csv)) {
  stop(
    " MLpeaks-supported sig motif CSV：", mlpeaks_supported_motif_csv, "\n",
    "（ mlpeaks_chromvar_mapping /map_chromvar_to_MLpeaks.R），/"
  )
}

cat("[INFO] Using ALL motifs from MLpeaks-supported sig CSV: ", mlpeaks_supported_motif_csv, "\n", sep = "")

motifs_raw <- read_csv(mlpeaks_supported_motif_csv, show_col_types = FALSE)
need_cols <- c("motif", "TF", "chromvar_lfc", "chromvar_padj")
miss <- setdiff(need_cols, colnames(motifs_raw))
if (length(miss) > 0) stop("MLpeaks-supported motif CSV : ", paste(miss, collapse = ", "))

motifs_sel <- motifs_raw %>%
  transmute(
    motif = as.character(motif),
    TF = as.character(TF),
    chromvar_lfc = as.numeric(chromvar_lfc),
    chromvar_padj = as.numeric(chromvar_padj),
    chromvar_log10padj = if ("chromvar_log10padj" %in% colnames(motifs_raw)) as.numeric(chromvar_log10padj) else NA_real_,
    n_peaks_with_motif = if ("n_peaks_with_motif" %in% colnames(motifs_raw)) as.numeric(n_peaks_with_motif) else NA_real_,
    frac_peaks_with_motif = if ("frac_peaks_with_motif" %in% colnames(motifs_raw)) as.numeric(frac_peaks_with_motif) else NA_real_
  ) %>%
  filter(!is.na(motif), motif != "", !is.na(TF), TF != "") %>%
  distinct(motif, .keep_all = TRUE)

n_before <- nrow(motifs_sel)
motifs_sel <- motifs_sel %>%
  filter(motif %in% colnames(motif_mat))
n_dropped <- n_before - nrow(motifs_sel)
if (n_dropped > 0) {
  cat("[WARN] Dropped ", n_dropped, " motifs not found in motif_mat columns.\n", sep = "")
}

motifs_sel <- motifs_sel %>%
  mutate(
    chromvar_padj = ifelse(is.na(chromvar_padj) | chromvar_padj <= 0, 1e-300, chromvar_padj),
    chromvar_log10padj = ifelse(is.finite(chromvar_log10padj), chromvar_log10padj, -log10(chromvar_padj)),
    TF_gene = simplify_tf_gene(TF),
    #  chromVAR motif  PRE_vs_CON：avg_log2FC > 0  PRE_up
    TF_direction = ifelse(chromvar_lfc > 0, pre_dir, con_dir)
  )

#  MLpeaks-supported CSV  coverage ， motif_cov （）
if (all(is.na(motifs_sel$frac_peaks_with_motif))) {
  motifs_sel <- motifs_sel %>%
    left_join(
      motif_cov %>% dplyr::select(motif, n_peaks_with_motif, frac_peaks_with_motif),
      by = "motif",
      suffix = c("", ".cov")
    ) %>%
    mutate(
      n_peaks_with_motif = ifelse(is.na(n_peaks_with_motif), n_peaks_with_motif.cov, n_peaks_with_motif),
      frac_peaks_with_motif = ifelse(is.na(frac_peaks_with_motif), frac_peaks_with_motif.cov, frac_peaks_with_motif)
    ) %>%
    dplyr::select(-ends_with(".cov"))
}

write_csv(motifs_sel, file.path(outdir, "selected_motifs_from_MLpeaksSupported_ALL.csv"))
cat("[INFO] Selected motifs (ALL): ", nrow(motifs_sel), " (unique TF_gene: ", n_distinct(motifs_sel$TF_gene), ")\n", sep = "")

# =========================
# 7) Map selected motifs -> peaks（）
#     netlink  motif_peak_mapping ， motifs  summary（）
# =========================
# ：“ motif/TF”“ peaks  motif”
# ：motif_mat ， summary(msub_sel) ， peak/motif 
motif_ids_sel <- motifs_sel$motif
msub_sel <- motif_mat[keep_peaks, motif_ids_sel, drop = FALSE]

mp_df <- as.data.frame(summary(msub_sel))
colnames(mp_df) <- c("peak_idx", "motif_idx", "present")

peak_names_all  <- rownames(msub_sel)
motif_names_all <- colnames(msub_sel)

motif_peak_mapping_sel <- mp_df %>%
  filter(present > 0) %>%
  mutate(
    peak  = peak_names_all[peak_idx],
    motif = motif_names_all[motif_idx]
  ) %>%
  dplyr::select(motif, peak) %>%
  distinct()

#  TF  + chromVAR 
#  TF_to_MLpeaks： TF(motif)->peak ，：
# - TF/motif  chromvar_lfc/padj（motif ）
# - peak /lfc/padj/SHAP（peak /）
TF_to_MLpeaks <- motifs_sel %>%
  dplyr::select(motif, TF, TF_gene, chromvar_lfc, chromvar_padj, chromvar_log10padj, TF_direction,
         n_peaks_with_motif, frac_peaks_with_motif) %>%
  left_join(motif_peak_mapping_sel, by = "motif") %>%
  filter(!is.na(peak)) %>%
  inner_join(
    da_peaks_dir %>%
      dplyr::select(peak, peak_direction, peak_lfc, peak_padj = p_val_adj, mean_abs_SHAP),
    by = "peak"
  )

if (direction_match) {
  # ：“TF ”“peak ” TF->peak 
  TF_to_MLpeaks <- TF_to_MLpeaks %>%
    filter(TF_direction == peak_direction)
}

TF_to_MLpeaks <- TF_to_MLpeaks %>%
  arrange(TF_gene, desc(abs(peak_lfc)))

write_csv(TF_to_MLpeaks, file.path(outdir, "TF_to_MLpeaks_edges.csv"))
cat("[INFO] TF->MLpeaks edges:", nrow(TF_to_MLpeaks),
    " | TF:", n_distinct(TF_to_MLpeaks$TF_gene),
    " | peaks:", n_distinct(TF_to_MLpeaks$peak), "\n")

# =========================
# 8) ML peaks -> genes (LinkPeaks)
# =========================
# ： peak  gene
#  links_best_csv  Signac::LinkPeaks ： peakgene(score)(FDR/pass_sig)
stopifnot(file.exists(links_best_csv))
links_best <- read_csv(links_best_csv, show_col_types = FALSE)

# 
#  Peak/peak_id ，
if (!"peak" %in% colnames(links_best)) {
  peak_col <- grep("peak", tolower(colnames(links_best)), value = TRUE)[1]
  if (is.na(peak_col)) stop("links_best  peak ")
  links_best <- links_best %>% rename(peak = all_of(peak_col))
}
links_best <- links_best %>% mutate(peak = norm_peak(peak))

gene_col <- intersect(c("gene","gene_name","symbol","SYMBOL"), colnames(links_best))[1]
if (is.na(gene_col)) stop("links_best  gene/gene_name/symbol ，")

if (!"score" %in% colnames(links_best)) {
  stop("links_best  score （LinkPeaks ）")
}
if (!"fdr" %in% colnames(links_best)) {
  warning("links_best  fdr ； FDR ")
  links_best$fdr <- NA_real_
}
if (!"pass_sig" %in% colnames(links_best)) {
  #  pass_sig，（）
  links_best <- links_best %>%
    mutate(pass_sig = ifelse(!is.na(fdr) & fdr < link_fdr_cut, TRUE, FALSE))
}

links_best2 <- links_best %>%
  transmute(
    peak,
    gene = .data[[gene_col]],
    link_score = as.numeric(score),
    link_fdr   = as.numeric(fdr),
    pass_sig
  ) %>%
  filter(!is.na(gene))

if (keep_pos_links) {
  # ：peak  ↑ -> gene  ↑（ enhancer ）
  links_best2 <- links_best2 %>% filter(link_score > 0)
}

# TF -> peak -> gene
#  TF->peak  peak->gene  peak  join，
# strength：“”（），/
TF_peak_gene <- TF_to_MLpeaks %>%
  inner_join(links_best2 %>% filter(pass_sig), by = "peak") %>%
  mutate(
    strength = abs(chromvar_lfc * peak_lfc * link_score)
  ) %>%
  arrange(desc(strength))

write_csv(TF_peak_gene, file.path(outdir, "TF_peak_gene_triplets.csv"))
cat("[INFO] TF->peak->gene triplets:", nrow(TF_peak_gene),
    " | genes:", n_distinct(TF_peak_gene$gene), "\n")

# Cytoscape edge lists（）
#  Cytoscape：
edges_TF_peak <- TF_to_MLpeaks %>%
  transmute(source = TF_gene, target = peak, edge_type = "TF_motif_to_peak",
            weight = abs(chromvar_lfc))

edges_peak_gene <- TF_peak_gene %>%
  transmute(source = peak, target = gene, edge_type = "peak_to_gene",
            weight = abs(link_score))

write_csv(edges_TF_peak,  file.path(outdir, "edges_TF_to_peak.csv"))
write_csv(edges_peak_gene, file.path(outdir, "edges_peak_to_gene.csv"))


# =========================
# 8.5) RNA DESeq2 direction-consistent filtering
# =========================
# ：“”——TF(motif)peakgene(RNA)
# ： triplets  gene  RNA DE ， RNA 
stopifnot(file.exists(rna_deseq_csv))
rna_de <- readr::read_csv(rna_deseq_csv, show_col_types = FALSE)

need_cols <- c("gene", "avg_log2FC", "p_val_adj")
miss <- setdiff(need_cols, colnames(rna_de))
if (length(miss) > 0) stop("DESeq2 CSV : ", paste(miss, collapse = ", "))

norm_gene <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  #  Ensembl ID ：ENSG... .1 -> ENSG...
  x <- gsub("\\.\\d+$", "", x)  #  Ensembl version：ENSG... .1
  toupper(x)
}

rna_de2 <- rna_de %>%
  transmute(
    gene_raw = as.character(gene),
    gene_norm = norm_gene(gene_raw),
    rna_lfc  = as.numeric(avg_log2FC),
    rna_padj = as.numeric(p_val_adj),
    gene_direction = dplyr::case_when(
      is.na(rna_lfc) ~ NA_character_,
      # RNA DE  PRE vs CON（ident.1=PRE, ident.2=CON）：avg_log2FC > 0  PRE_up
      rna_lfc > 0    ~ pre_dir,
      rna_lfc < 0    ~ con_dir,
      TRUE           ~ NA_character_
    )
  ) %>%
  filter(!is.na(gene_direction))

#  RNA  triplets
# ： gene_norm join，/
TF_peak_gene_rna <- TF_peak_gene %>%
  mutate(gene_norm = norm_gene(gene)) %>%
  left_join(rna_de2, by = "gene_norm")

# ： gene_direction == TF_direction
# ： rna_padj “”（ RNA DE ）
#  RNA ， require_rna_padj  TRUE
TF_peak_gene_consistent <- TF_peak_gene_rna %>%
  filter(!is.na(rna_lfc)) %>%
  filter(gene_direction == TF_direction)

if (require_rna_padj) {
  TF_peak_gene_consistent <- TF_peak_gene_consistent %>%
    filter(!is.na(rna_padj), rna_padj < rna_padj_cut)
}

cat("[INFO] Triplets raw:", nrow(TF_peak_gene), "\n")
cat("[INFO] Triplets after RNA-direction filter:", nrow(TF_peak_gene_consistent), "\n")

# （）“”：peak_direction  gene_direction 
# ： peak  gene 
# /， enforce_peak_gene_match  FALSE
enforce_peak_gene_match <- TRUE
if (enforce_peak_gene_match && "peak_direction" %in% colnames(TF_peak_gene_consistent)) {
  TF_peak_gene_consistent <- TF_peak_gene_consistent %>%
    filter(peak_direction == gene_direction)
  cat("[INFO] Triplets after peak-gene direction filter:", nrow(TF_peak_gene_consistent), "\n")
}

# （/）
# strength_rna  RNA ： motif peak peak-gene linkRNA 
TF_peak_gene_consistent <- TF_peak_gene_consistent %>%
  mutate(strength_rna = abs(chromvar_lfc * peak_lfc * link_score * rna_lfc)) %>%
  arrange(desc(strength_rna))

write_csv(TF_peak_gene_consistent, file.path(outdir, "TF_peak_gene_triplets_RNA_consistent.csv"))

# （）， Cytoscape
edges_peak_gene_consistent <- TF_peak_gene_consistent %>%
  transmute(source = peak, target = gene, edge_type = "peak_to_gene", weight = abs(link_score))
write_csv(edges_peak_gene_consistent, file.path(outdir, "edges_peak_to_gene_RNA_consistent.csv"))









