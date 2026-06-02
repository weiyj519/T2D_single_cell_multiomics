#!/usr/bin/env Rscript

# Script name: 01_filter_atac_peak_matrix.R
# Purpose: Cleaned public workflow script for the T2D single-cell multiome analysis.
# Input: Paths and analysis settings are read from configs/config.yaml.
# Output: Module-specific outputs are written under the configured results directory.
# Main steps: Load configuration, run the original analysis logic, and export reproducible outputs.
# Example command: Rscript scripts/01_feature_matrix_preparation/01_filter_atac_peak_matrix.R

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
  library(Matrix)
  library(arrow)
  library(data.table)
})

qs_path         <- seurat_object_path
atac_csv        <- file.path(processed_data_dir, "atac_peaks_autosomes.csv")

# 1.  Seurat 
combined <- qread(qs_path)

#  ATAC assay 
stopifnot("ATAC" %in% names(combined@assays))

DefaultAssay(combined) <- "ATAC"
combined <- JoinLayers(combined)  #  counts/data

# 2.  ATAC  "data"  "counts" 
# ， log-normalized ，
atac_mat <- GetAssayData(combined, assay = "ATAC", layer = "data")
#  TF-IDF  log transform， slot = "data"

cat(" ATAC :", nrow(atac_mat), " :", ncol(atac_mat), "\n")

# 3. 
peak_names <- rownames(atac_mat)  #  "chr1:100-200"

peak_df <- data.frame(
  peak = peak_names,
  stringsAsFactors = FALSE
)

# 

peak_df <- peak_df %>%
  tidyr::separate(
    col   = "peak",
    into  = c("chr", "start", "end"),
    sep   = "-",
    remove = FALSE,
    fill  = "right",
    extra = "drop"
  )

# 4. （chr1–chr22）
autosome_chr <- paste0("chr", 1:22)

# 5.  peaks
#   - ：chr1–chr22
#   - ：chrX / chrY / chrM （ GL.. ）
keep_peaks <- peak_df$peak[peak_df$chr %in% autosome_chr]

cat(": ", length(keep_peaks), "\n")

# 6.  peak
atac_mat_keep <- atac_mat[keep_peaks, , drop = FALSE]

cat(" ATAC : ", nrow(atac_mat_keep), " peaks x ",
    ncol(atac_mat_keep), " cells\n")

## 7. “” peaks（）
## -----------------------------------------------------
library(Matrix)

# （TF-IDF > 0 ）
peak_nonzero_frac <- Matrix::rowSums(atac_mat_keep != 0) / ncol(atac_mat_keep)

# ： ≥1% 
min_frac <- 0.01
keep_by_frac <- peak_nonzero_frac >= min_frac

cat(" peaks :", nrow(atac_mat_keep), "\n")
cat(" ≥", min_frac, "  peaks :", sum(keep_by_frac), "\n")

atac_mat_freq <- atac_mat_keep[keep_by_frac, , drop = FALSE]


## 8. “” peaks（）
## -----------------------------------------------------
# ：Var(X) = E[X^2] - (E[X])^2
peak_mean    <- Matrix::rowMeans(atac_mat_freq)
peak_mean_sq <- Matrix::rowMeans(atac_mat_freq^2)
peak_var     <- peak_mean_sq - peak_mean^2

#  peaks， 50%  30%
var_quantile <- 0.3   # 0.5 =  50%， 0.7 
var_cutoff   <- as.numeric(quantile(peak_var, probs = var_quantile, na.rm = TRUE))

keep_by_var <- peak_var >= var_cutoff

cat(" peaks :", nrow(atac_mat_freq), "\n")
cat("", (1 - var_quantile) * 100, "%  peaks :",
    sum(keep_by_var), "\n")

atac_mat_final <- atac_mat_freq[keep_by_var, , drop = FALSE]

cat(" ATAC : ", 
    nrow(atac_mat_final), " peaks x ", 
    ncol(atac_mat_final), " cells\n")

## 7.  cell × peak （ atac_mat_final）
## -----------------------------------------------------------
atac_mat_final_t <- Matrix::t(atac_mat_final)

# ，（）
rm(atac_mat, atac_mat_keep, atac_mat_freq, atac_mat_final); gc()

## 8.  meta 
## -----------------------------------------------------------
meta_cols <- c("orig.ident", "seurat_clusters", "celltype", "group")
stopifnot(all(meta_cols %in% colnames(combined@meta.data)))

meta_df <- combined@meta.data[, meta_cols, drop = FALSE] %>%
  tibble::rownames_to_column(var = "cell_id")

#  meta 
meta_df$cell_id         <- as.character(meta_df$cell_id)
meta_df$orig.ident      <- as.character(meta_df$orig.ident)
meta_df$seurat_clusters <- as.character(meta_df$seurat_clusters)
meta_df$celltype        <- as.character(meta_df$celltype)
meta_df$group           <- as.character(meta_df$group)

## 9.  ATAC  data.frame， cell_id
## -----------------------------------------------------------
atac_df <- atac_mat_final_t %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "cell_id")

atac_df$cell_id <- as.character(atac_df$cell_id)

## 10.  meta + ATAC（cell × peak TF-IDF ）
## -----------------------------------------------------------
df_atac_export <- dplyr::left_join(meta_df, atac_df, by = "cell_id")

cat(": ", nrow(df_atac_export), " cells x ",
    ncol(df_atac_export), " columns\n")

## 11.  CSV（ / ）
## -----------------------------------------------------------
fwrite(df_atac_export, atac_csv)
cat(" ATAC peak : ", atac_csv, "\n")
