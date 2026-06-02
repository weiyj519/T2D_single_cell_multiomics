#!/usr/bin/env Rscript

# Script name: 02_run_atac_great_enrichment.R
# Purpose: Cleaned public workflow script for the T2D single-cell multiome analysis.
# Input: Paths and analysis settings are read from configs/config.yaml.
# Output: Module-specific outputs are written under the configured results directory.
# Main steps: Load configuration, run the original analysis logic, and export reproducible outputs.
# Example command: Rscript scripts/05_functional_enrichment/02_run_atac_great_enrichment.R

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
  library(readr)
  library(dplyr)
  library(stringr)
  library(GenomicRanges)
  library(IRanges)
  library(rGREAT)
  library(Seurat)
})

# ----------  ----------
qs_path <- seurat_object_path  # （ATAC assay）
fg_csv  <- "results/model_outputs/atac/Acinar_cell/boruta_gene_stability_CON_vs_PRE.csv"
outdir  <- "results/downstream/acinar_CON_PRE/great"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---------- ：peak -> GRanges ----------
# ：chr1-100-200 / chr1:100-200 / chr1_100_200

library(stringr)
library(GenomicRanges)
library(IRanges)

peaks_to_gr <- function(peaks, keep_standard = TRUE, convert_ucsc = TRUE) {
  peaks <- as.character(peaks)
  peaks <- gsub("\\s+", "", peaks)
  peaks <- gsub("_", "-", peaks)
  peaks <- gsub(":", "-", peaks)

  #  seqname：^(anything)-(start)-(end)$
  m <- str_match(peaks, "^(.+?)-(\\d+)-(\\d+)$")
  if (any(is.na(m[, 1]))) {
    bad <- peaks[is.na(m[, 1])][1]
    stop(paste0("peak，：", bad))
  }

  seq <- m[, 2]
  start <- as.integer(m[, 3])
  end   <- as.integer(m[, 4])

  if (convert_ucsc) {
    # 1 / X / Y / M / MT -> chr1 / chrX / chrY / chrM
    seq <- ifelse(seq %in% c("MT", "Mt", "mt"), "M", seq)
    seq <- ifelse(str_detect(seq, "^(\\d+|X|Y|M)$"), paste0("chr", seq), seq)

    # hg19  contig：GL000192.1 -> chrUn_gl000192
    seq <- ifelse(
      str_detect(seq, "^GL\\d+\\.\\d+$"),
      paste0("chrUn_gl", str_replace(seq, "^GL(\\d+)\\..*$", "\\1")),
      seq
    )
  }

  gr <- GRanges(seqnames = seq, ranges = IRanges(start = start, end = end))

  # GREAT 
  if (keep_standard) {
    std <- paste0("chr", c(1:22, "X", "Y", "M"))
    std_present <- intersect(std, GenomeInfoDb::seqlevels(gr))
    gr <- GenomeInfoDb::keepSeqlevels(gr, std_present, pruning.mode = "coarse")

  }

  gr
}

# ---------- A) ：ATAC peaks ----------
combined <- qread(qs_path)
DefaultAssay(combined) <- "ATAC"
bg_peaks <- rownames(combined[["ATAC"]])
bg_gr <- peaks_to_gr(bg_peaks)

cat("[INFO] Background peaks:", length(bg_gr), "\n")


# ---------- B) ： peaks（fold_freq >= 0.6） ----------
fg <- read_csv(fg_csv, show_col_types = FALSE)

fg <- fg %>% dplyr::rename_with(tolower)

# boruta_gene_stability_*.csv  gene  peak； peak 
peak_col <- intersect(c("peak", "gene"), colnames(fg))[1]
if (is.na(peak_col)) {
  stop("Foreground CSV missing `peak`/`gene` column. Columns: ", paste(colnames(fg), collapse = ", "))
}

#  fold_freq  >=0.6 （）；
if ("fold_freq" %in% colnames(fg)) {
  fg <- fg %>% filter(as.numeric(.data$fold_freq) >= 0.6)
}

# peak：peak  gene(peak)
library(readr)
library(dplyr)

fg_peaks <- fg %>%
  pull(dplyr::all_of(peak_col)) %>%
  unique()

# “ATAC peaks”（/）
fg_peaks_norm <- gsub(":", "-", gsub("_", "-", fg_peaks))
bg_peaks_norm <- gsub(":", "-", gsub("_", "-", bg_peaks))
fg_keep_norm  <- intersect(fg_peaks_norm, bg_peaks_norm)

cat("[INFO] Foreground peaks (>=0.6) matched in ATAC:", length(fg_keep_norm), "\n")
if (length(fg_keep_norm) < 10) stop("，peak/hg19")

fg_gr <- peaks_to_gr(fg_keep_norm)

# ---------- C)  local GREAT ----------
# tss_source： "GREAT:hg19"（GREATTSS）， "hg19" :contentReference[oaicite:2]{index=2}
tss_source <- "GREAT:hg19"

# GO
res_bp <- great(fg_gr, gene_sets = "GO:BP", tss_source = tss_source, background = bg_gr, cores = 4)
res_cc <- great(fg_gr, gene_sets = "GO:CC", tss_source = tss_source, background = bg_gr, cores = 4)
res_mf <- great(fg_gr, gene_sets = "GO:MF", tss_source = tss_source, background = bg_gr, cores = 4)

# KEGG（MSigDB：C2:CP:KEGG，human） :contentReference[oaicite:3]{index=3}
res_kegg <- great(fg_gr, gene_sets = "msigdb:C2:CP:KEGG", tss_source = tss_source, background = bg_gr, cores = 4)

# ---------- D)  ----------
tb_bp   <- getEnrichmentTable(res_bp)
tb_cc   <- getEnrichmentTable(res_cc)
tb_mf   <- getEnrichmentTable(res_mf)
tb_kegg <- getEnrichmentTable(res_kegg)

write_csv(tb_bp,   file.path(outdir, "GO_BP_localGREAT.csv"))

write_csv(tb_cc,   file.path(outdir, "GO_CC_localGREAT.csv"))
write_csv(tb_mf,   file.path(outdir, "GO_MF_localGREAT.csv"))
write_csv(tb_kegg, file.path(outdir, "KEGG_localGREAT.csv"))

# ， region-gene associations
saveRDS(
  list(
    fg_csv = fg_csv,
    fg_n   = length(fg_gr),
    bg_n   = length(bg_gr),
    tss_source = tss_source,
    res = list(BP=res_bp, CC=res_cc, MF=res_mf, KEGG=res_kegg)
  ),
  file.path(outdir, "localGREAT_results.rds")
)

cat("[DONE] Results saved in:", outdir, "\n")




library(dplyr)
library(tibble)
library(ggprism)
library(readr)
library(stringr)
library(ggplot2)
outdir <- "results/downstream/acinar_CON_PRE/great"

files <- c(
  BP   = file.path(outdir, "GO_BP_localGREAT.csv"),
  CC   = file.path(outdir, "GO_CC_localGREAT.csv"),
  MF   = file.path(outdir, "GO_MF_localGREAT.csv"),
  KEGG = file.path(outdir, "KEGG_localGREAT.csv")
)

# ：
#   - "p_adjust"        GREAT binomial ()
#   - "p_adjust_hyper"  hypergeometric ()
P_COL <- "p_value"

# 
pal <- c(
  BP   = "#e6a0c4",
  MF   = "#c6cdf7",
  CC   = "#d8a499",
  KEGG = "#7294d4"
)

# ========= 1) ： =========
ok <- sapply(files, file.exists)
if (!all(ok)) {
  stop("：\n", paste(names(files)[!ok], files[!ok], sep=" -> ", collapse="\n"))
}

# ========= 2)  rGREAT ： GO/KEGG  =========
read_res_great <- function(path, onto, p_col = "p_adjust") {
  #  KEGG/GO ： 0  tibble， bind_rows/
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }

  if (isTRUE(file.size(path) == 0)) {
    return(tibble(
      id = character(),
      Description = character(),
      ONTOLOGY = character(),
      padj = numeric()
    ))
  }

  df <- read_csv(path, show_col_types = FALSE)
  if (nrow(df) == 0) {
    return(tibble(
      id = character(),
      Description = character(),
      ONTOLOGY = character(),
      padj = numeric()
    ))
  }

  # （）
  if ("id" %in% names(df)) {
    df <- df %>% filter(!is.na(id), id != "id")
  }

  #  id 
  if (!("id" %in% names(df))) stop("Missing column id in: ", path)
  if (!(p_col %in% names(df))) stop("Missing column ", p_col, " in: ", path)

  #  Description：
  # GO  description；KEGG  id
  if ("description" %in% names(df)) {
    df <- df %>% mutate(Description = as.character(description))
  } else {
    df <- df %>% mutate(Description = as.character(id))
  }

  #  KEGG 
  if (onto == "KEGG") {
    df <- df %>%
      mutate(
        Description = str_replace(Description, "^KEGG_", ""),
        Description = str_replace_all(Description, "_", " "),
        Description = str_squish(Description)
      )
  }

  # ：padj（ p.adjust ）
  df <- df %>%
    mutate(
      ONTOLOGY = onto,
      padj = suppressWarnings(as.numeric(.data[[p_col]]))
    )

  df
}

GO_BP <- read_res_great(files["BP"],   "BP",   P_COL)
GO_CC <- read_res_great(files["CC"],   "CC",   P_COL)
GO_MF <- read_res_great(files["MF"],   "MF",   P_COL)
KEGG  <- read_res_great(files["KEGG"], "KEGG", P_COL)

all_pathways <- bind_rows(GO_BP, GO_CC, GO_MF, KEGG)

# ========= 3)  Top（GO Top5 + KEGG Top5）“” =========
n_go   <- 5
n_kegg <- 5
SIG_CUTOFF <- 0.05

plot_df <- all_pathways %>%
  mutate(
    padj = suppressWarnings(as.numeric(padj)),
    padj = ifelse(is.na(padj) | padj <= 0, 1e-300, padj),
    Description = ifelse(is.na(Description) | Description == "", as.character(id), Description),
    DB = ifelse(ONTOLOGY == "KEGG", "KEGG", "GO")
  )

# ；
plot_df_sig <- plot_df %>%
  filter(!is.na(padj), padj < SIG_CUTOFF)

# rGREAT  clusterProfiler  Count ； region/gene  tie-breaker
COUNT_COL <- intersect(
  c("observed_region_hits", "observed_gene_hits", "gene_set_size"),
  names(plot_df)
)[1]

arrange_with_tie <- function(df) {
  if (!is.na(COUNT_COL) && COUNT_COL %in% names(df)) {
    df %>% arrange(padj, desc(.data[[COUNT_COL]]))
  } else {
    df %>% arrange(padj)
  }
}

go_top <- plot_df %>%
  filter(DB == "GO") %>%
  group_by(ONTOLOGY) %>%
  do(arrange_with_tie(.)) %>%
  slice_head(n = n_go) %>%
  ungroup()

kegg_top <- plot_df %>%
  filter(DB == "KEGG") %>%
  arrange_with_tie() %>%
  slice_head(n = n_kegg)

go_top <- plot_df_sig %>%
  filter(DB == "GO") %>%
  group_by(ONTOLOGY) %>%
  do(arrange_with_tie(.)) %>%
  slice_head(n = n_go) %>%
  ungroup()

kegg_top <- plot_df_sig %>%
  filter(DB == "KEGG") %>%
  arrange_with_tie() %>%
  slice_head(n = n_kegg)

plot_df2 <- bind_rows(go_top, kegg_top) %>%
  mutate(
    Description_single = str_squish(Description),
    score = -log10(padj)
  ) %>%
  mutate(ONTOLOGY = as.character(ONTOLOGY))

if (nrow(plot_df2) == 0) {
  stop(" padj < ", SIG_CUTOFF, " （ ", P_COL, "）")
}

onto_order <- c("BP", "CC", "MF", "KEGG")
onto_present <- intersect(onto_order, unique(plot_df2$ONTOLOGY))
plot_df2 <- plot_df2 %>%
  mutate(ONTOLOGY = factor(ONTOLOGY, levels = onto_present)) %>%
  group_by(ONTOLOGY) %>%
  arrange(desc(score), .by_group = TRUE) %>%
  ungroup() %>%
  arrange(ONTOLOGY, desc(score))

# （）： ONTOLOGY ， score 
plot_df2$Description_single <- factor(
  plot_df2$Description_single,
  levels = rev(plot_df2$Description_single)
)

score_breaks <- pretty(c(0, max(plot_df2$score, na.rm = TRUE)), n = 5)

p_all <- ggplot(plot_df2, aes(x = score, y = Description_single)) +
  geom_col(aes(fill = ONTOLOGY), width = 0.7) +
  geom_text(
    aes(label = Description_single, x = 0.02),
    hjust = 0, size = 4.5, color = "black", fontface = "bold"
  ) +
  scale_fill_manual(values = pal) +
  scale_x_continuous(
    breaks = score_breaks,
    labels = score_breaks,
    expand = c(0, 0),
    limits = c(0, max(plot_df2$score, na.rm = TRUE) * 1.15)
  ) +
  labs(
    title = NULL,
    x = expression(-log[10](p_value)),
    y = NULL
  ) +
  theme_classic() +
  theme(
    text = element_text(face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.x  = element_text(size = 12, face = "bold"),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_line(color = "black"),
    legend.title = element_blank(),
    legend.position = "right",
    panel.grid = element_blank()
  )

ggsave(file.path(outdir, "pathway_enrichment_insideLabel_localGREAT.pdf"),
       p_all, width = 10, height = 6)
ggsave(file.path(outdir, "pathway_enrichment_insideLabel_localGREAT.png"),
       p_all, width = 10, height = 5, dpi = 500)
