#!/usr/bin/env Rscript

# Script name: 02_run_peak_to_gene_links.R
# Purpose: Cleaned public workflow script for the T2D single-cell multiome analysis.
# Input: Paths and analysis settings are read from configs/config.yaml.
# Output: Module-specific outputs are written under the configured results directory.
# Main steps: Load configuration, run the original analysis logic, and export reproducible outputs.
# Example command: Rscript scripts/06_atac_motif_peak_gene/02_run_peak_to_gene_links.R

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
  library(dplyr)
  library(readr)
  library(stringr)

  library(GenomeInfoDb)
  library(GenomicRanges)

  #  hg19； hg38，“hg38 ”
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(BSgenome.Hsapiens.UCSC.hg19)
})

# =====================  =====================
obj_path  <- seurat_object_path
celltype_use <- c("Acinar cell", "acinar cell", "Acinar_cell", "acinar_cell")

peaks_csv <- "results/model_outputs/atac/Acinar_cell/PRE_vs_T2D_shap_marker_genes.csv"
genes_csv <- "results/stage_features/stage_genes/Acinar_cell__PRE_vs_T2D__test_shap_geq0p6.csv"

outdir   <- "results/downstream/atac_motif_peak_gene/linkpeaks/acinar_cell/PRE_vs_T2D/FEATURE_FROM_SHAP"
max_dist <- 500000      # +/- 500 kb
cor_thr  <- 0.10        # correlation > 0.x
use_fdr  <- TRUE        # TRUE:  BH-FDR；FALSE:  p
p_thr    <- 0.05
fdr_thr  <- 0.05
keep_direction <- "both" # "positive"|"negative"|"both"

min_cells <- 10
n_sample  <- 200
# =======================================================

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

norm_peak <- function(x) {
  x <- str_trim(x)
  x <- gsub("\\s+", "", x)
  x <- gsub(":", "-", x)
  x <- gsub("_", "-", x)
  x
}

# ============ 1)  SHAP  -> feature peaks / genes ============
peaks_tbl <- read_csv(peaks_csv, show_col_types = FALSE) %>%
  rename_with(~tolower(.x)) %>%
  dplyr::rename(peak = dplyr::any_of(c("peak"))) %>%
  mutate(
    peak = as.character(peak),
    peak_norm = norm_peak(peak),
    mean_shap = as.numeric(dplyr::coalesce(!!sym("mean_shap"), !!sym("mean_shap"))),
    mean_abs_shap = as.numeric(dplyr::coalesce(!!sym("mean_abs_shap"), !!sym("mean_abs_shap")))
  ) %>%
  dplyr::select(peak, peak_norm, everything()) %>%
  distinct(peak_norm, .keep_all = TRUE)

genes_tbl <- read_csv(genes_csv, show_col_types = FALSE) %>%
  rename_with(~tolower(.x)) %>%
  dplyr::rename(gene = dplyr::any_of(c("gene"))) %>%
  mutate(
    gene = as.character(gene),
    mean_shap = as.numeric(dplyr::coalesce(!!sym("mean_shap"), !!sym("mean_shap"))),
    mean_abs_shap = as.numeric(dplyr::coalesce(!!sym("mean_abs_shap"), !!sym("mean_abs_shap")))
  ) %>%
  dplyr::select(gene, everything()) %>%
  distinct(gene, .keep_all = TRUE)

feature_peaks_norm <- peaks_tbl$peak_norm
feature_genes      <- genes_tbl$gene

cat("[INFO] Input feature peaks:", length(feature_peaks_norm), "\n")
cat("[INFO] Input feature genes:", length(feature_genes), "\n")

# ============ 2)  ============
obj_all <- qread(obj_path)
stopifnot(all(c("RNA","ATAC") %in% names(obj_all@assays)))
stopifnot("celltype" %in% colnames(obj_all@meta.data))

obj <- subset(obj_all, subset = celltype %in% celltype_use)
if (ncol(obj) == 0) {
  ct <- unique(as.character(obj_all$celltype))
  hint <- sort(unique(ct[grepl("cinar", ct, ignore.case = TRUE)]))
  stop(
    " celltype（celltype_use=",
    paste(celltype_use, collapse = ", "),
    ")； acinar : ",
    if (length(hint) == 0) "<none>" else paste(hint, collapse = ", ")
  )
}
stopifnot(ncol(obj) > 50)

# ============ 3) peaks / genes  ============
DefaultAssay(obj) <- "ATAC"
peak_raw_all <- rownames(obj[["ATAC"]])
peak_norm_all <- norm_peak(peak_raw_all)
peak_map <- tibble(peak_raw = peak_raw_all, peak_norm = peak_norm_all) %>%
  distinct(peak_norm, .keep_all = TRUE)

peaks_keep_norm <- intersect(feature_peaks_norm, peak_map$peak_norm)
peaks_keep_raw  <- peak_map$peak_raw[match(peaks_keep_norm, peak_map$peak_norm)]

DefaultAssay(obj) <- "RNA"
genes_keep <- intersect(feature_genes, rownames(obj[["RNA"]]))

cat("[INFO] Peaks matched in object:", length(peaks_keep_raw), "/", length(feature_peaks_norm), "\n")
cat("[INFO] Genes matched in object:", length(genes_keep), "/", length(feature_genes), "\n")
if (length(peaks_keep_raw) < 5) stop(" peaks ： peak （chr-start-end / chr:start-end）")
if (length(genes_keep) < 5) stop(" genes ： gene symbol  RNA rownames ")

# ============ 4)  Feature Peaks（） ============
DefaultAssay(obj) <- "ATAC"
# ：Seurat::subset(object, features=...)  assay  feature 
#  features  ATAC peaks， RNA assay，
# "None of the features provided found in this assay"
#  assay(ATAC)  features 
obj[["ATAC"]] <- subset(obj[["ATAC"]], features = peaks_keep_raw)
stopifnot(all(peaks_keep_raw %in% rownames(obj[["ATAC"]])))

# ============ 5)  RNA data  ============
DefaultAssay(obj) <- "RNA"
if (!"data" %in% names(obj@assays$RNA@layers)) {
  obj <- NormalizeData(obj, verbose = FALSE)
}

# ============ 6)  ATAC  & RegionStats（LinkPeaks ） ============
DefaultAssay(obj) <- "ATAC"

# RegionStats/MatchRegionStats  BSgenome  seqlevels
#  peaks  genome  contigs， NA 
bs_chroms <- GenomeInfoDb::seqlevels(BSgenome.Hsapiens.UCSC.hg19)
peak_ranges <- obj[["ATAC"]]@ranges
keep_chr <- as.character(GenomicRanges::seqnames(peak_ranges)) %in% bs_chroms
if (any(!keep_chr)) {
  drop_n <- sum(!keep_chr)
  keep_features <- rownames(obj[["ATAC"]])[keep_chr]
  obj[["ATAC"]] <- subset(obj[["ATAC"]], features = keep_features)
  cat("[WARN] Dropped ", drop_n, " peaks on seqlevels not in BSgenome hg19\n", sep = "")
}

if (is.null(Annotation(obj)) || length(Annotation(obj)) == 0) {
  tx <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene, columns = c("tx_id","tx_name"))
  tx_map <- AnnotationDbi::select(
    TxDb.Hsapiens.UCSC.hg19.knownGene,
    keys = as.character(mcols(tx)$tx_id),
    keytype = "TXID",
    columns = c("TXID","GENEID")
  )
  m <- match(as.character(mcols(tx)$tx_id), as.character(tx_map$TXID))
  mcols(tx)$gene_id <- as.character(tx_map$GENEID[m])

  entrez <- unique(na.omit(mcols(tx)$gene_id))
  sym_map <- mapIds(org.Hs.eg.db, keys = entrez, keytype = "ENTREZID", column = "SYMBOL")
  mcols(tx)$gene_name <- unname(sym_map[mcols(tx)$gene_id])

  mcols(tx)$gene_biotype <- "protein_coding"
  mcols(tx)$type <- "transcript"
  seqlevelsStyle(tx) <- "UCSC"

  present_chroms <- intersect(seqlevels(tx), unique(as.character(seqnames(obj[["ATAC"]]@ranges))))
  tx <- keepSeqlevels(tx, present_chroms, pruning.mode = "coarse")
  Annotation(obj) <- tx
}

obj <- RegionStats(obj, genome = BSgenome.Hsapiens.UCSC.hg19, verbose = TRUE)

# / RegionStats  sequence.length； MatchRegionStats 
meta_feat <- obj[["ATAC"]]@meta.features
if (!"sequence.length" %in% colnames(meta_feat)) {
  obj[["ATAC"]]@meta.features$sequence.length <- as.numeric(GenomicRanges::width(obj[["ATAC"]]@ranges))
  meta_feat <- obj[["ATAC"]]@meta.features
}

# LinkPeaks  MatchRegionStats  gene  background peaks
#  n_sample ， Signac ：
#   "Error in bg.coef[...] : subscript out of bounds"
#  RegionStats  meta.features  n_sample 
bg_cols <- c("GC.percent", "sequence.length")
if (all(bg_cols %in% colnames(meta_feat))) {
  n_bg <- sum(stats::complete.cases(meta_feat[, bg_cols, drop = FALSE]))
  if (is.finite(n_bg) && n_bg > 0) {
    # ， gene  peaks 
    n_sample_eff <- min(n_sample, max(20, n_bg - 1))
    if (n_sample_eff < n_sample) {
      cat("[WARN] n_sample=", n_sample, " > available background peaks=", n_bg,
          "; set n_sample=", n_sample_eff, "\n", sep = "")
    }
    n_sample <- n_sample_eff
  }
} else {
  cat("[WARN] meta.features missing required columns for MatchRegionStats: ",
      paste(setdiff(bg_cols, colnames(meta_feat)), collapse = ", "),
      " ; consider checking RegionStats/genome build\n", sep = "")
}

# ============ 7)  Feature Genes  LinkPeaks（±500kb） ============
run_linkpeaks <- function(n_sample_use) {
  LinkPeaks(
    object = obj,
    peak.assay = "ATAC",
    expression.assay = "RNA",
    peak.slot = "data",
    expression.slot = "data",
    method = "pearson",
    distance = max_dist,
    min.cells = min_cells,
    n_sample = n_sample_use,
    pvalue_cutoff = 1,     # ，
    #  keep_direction /， LinkPeaks 
    score_cutoff = if (keep_direction %in% c("negative", "both")) -1 else 0,
    gene.id = FALSE,
    verbose = TRUE
  )
}

obj <- tryCatch(
  run_linkpeaks(n_sample),
  error = function(e) {
    msg <- conditionMessage(e)
    if (grepl("subscript out of bounds", msg, fixed = TRUE)) {
      n2 <- max(50, floor(n_sample * 0.75))
      n3 <- max(30, floor(n_sample * 0.50))
      cat("[WARN] LinkPeaks hit known background-sampling bounds error; retry with smaller n_sample: ",
          n2, " then ", n3, "\n", sep = "")
      obj2 <- tryCatch(run_linkpeaks(n2), error = function(e2) e2)
      if (inherits(obj2, "error")) {
        obj3 <- tryCatch(run_linkpeaks(n3), error = function(e3) e3)
        if (inherits(obj3, "error")) stop(obj3)
        return(obj3)
      }
      return(obj2)
    }
    stop(e)
  }
)

# ============ 8)  links  ============
links_df <- as.data.frame(Signac::Links(obj))

# Links(obj)  GRanges（ peak），/ peak/gene/score/pvalue 
# ： peak  seqnames/start/end ； gene/score/pvalue  links
if (!"peak" %in% colnames(links_df) && all(c("seqnames", "start", "end") %in% colnames(links_df))) {
  links_df$peak <- paste0(links_df$seqnames, "-", links_df$start, "-", links_df$end)
}

need_cols <- c("peak", "gene", "score", "pvalue")
if (!all(need_cols %in% colnames(links_df))) {
  miss <- setdiff(need_cols, colnames(links_df))
  cat("[WARN] Links(obj) missing columns: ", paste(miss, collapse = ", "),
      "; likely no links were found under current settings. Writing empty outputs.\n", sep = "")
  links_df <- tibble::tibble(
    peak = character(),
    gene = character(),
    score = numeric(),
    pvalue = numeric()
  )
}

links_df$peak_norm <- norm_peak(links_df$peak)

# ：Feature Peaks（ peaks） +  ML genes（ A）
links_df <- links_df %>%
  filter(peak_norm %in% peaks_keep_norm, gene %in% genes_keep)

# FDR
links_df$fdr <- p.adjust(links_df$pvalue, method = "BH")

# /
pass_cor <- switch(
  keep_direction,
  positive = links_df$score >= cor_thr,
  negative = links_df$score <= -cor_thr,
  both     = abs(links_df$score) >= cor_thr
)
pass_p <- if (use_fdr) (links_df$fdr <= fdr_thr) else (links_df$pvalue <= p_thr)

links_df$pass_sig <- pass_cor & pass_p

#  SHAP  join ，“”
links_df <- links_df %>%
  left_join(
    peaks_tbl %>% dplyr::select(peak_norm, peak_mean_shap = mean_shap, peak_mean_abs_shap = mean_abs_shap),
    by = "peak_norm"
  ) %>%
  left_join(
    genes_tbl %>% dplyr::select(gene, gene_mean_shap = mean_shap, gene_mean_abs_shap = mean_abs_shap),
    by = "gene"
  )

cat("[INFO] links total (feature×feature):", nrow(links_df), "\n")
cat("[INFO] links pass_sig:", sum(links_df$pass_sig, na.rm = TRUE), "\n")

write.csv(links_df,
          file.path(outdir, "feature_links_ALL.csv"),
          row.names = FALSE)

write.csv(links_df %>% filter(pass_sig),
          file.path(outdir, "feature_links_SIG.csv"),
          row.names = FALSE)


qsave(obj, file.path(outdir, "obj_after_feature_linkpeaks.qs"), preset = "high")