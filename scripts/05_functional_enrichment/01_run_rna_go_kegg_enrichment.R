#!/usr/bin/env Rscript

# Script name: 01_run_rna_go_kegg_enrichment.R
# Purpose: Cleaned public workflow script for the T2D single-cell multiome analysis.
# Input: Paths and analysis settings are read from configs/config.yaml.
# Output: Module-specific outputs are written under the configured results directory.
# Main steps: Load configuration, run the original analysis logic, and export reproducible outputs.
# Example command: Rscript scripts/05_functional_enrichment/01_run_rna_go_kegg_enrichment.R

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

# =========  =========
Sys.setenv(OMP_NUM_THREADS = 4, OPENBLAS_NUM_THREADS = 4, MKL_NUM_THREADS = 4,
           VECLIB_MAXIMUM_THREADS = 4, NUMEXPR_NUM_THREADS = 4)

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(dplyr)
  library(stringr)
  library(Matrix)
  library(clusterProfiler)
  library(msigdbr)
  library(org.Hs.eg.db)   # ； org.Mm.eg.db  species  "Mus musculus"
  library(enrichplot)
  library(ggplot2)
  library(readr)
})

# =========  =========
input_qs <- seurat_object_path
species <- "Homo sapiens"                 #  "Mus musculus"
target_celltype <- "Acinar cell"          # （universe）：acinar cell 
outdir <- "results/downstream/acinar_CON_PRE/GO_RNA"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ：（CSV， gene）
stable_genes_csv <- "results/stage_features/stage_genes/Acinar_cell__CON_vs_PRE__test_shap_geq0p6.csv"

# （universe）：（acinar cell）（counts > 0 in any cell； RNA assay）

# ========= （）=========
cat("[INFO] Loading:", input_qs, "\n")
combined <- qread(input_qs)

if (!"celltype" %in% colnames(combined@meta.data)) {
  stop("meta.data  celltype ： acinar cell （universe）")
}

# ========= （acinar cell ）=========
bg_assay <- if ("RNA" %in% names(combined@assays)) "RNA" else DefaultAssay(combined)
acinar_celltype_candidates <- c("Acinar cell", "Acinar_cell", "acinar cell", "acinar_cell", "ACINAR cell", "ACINAR_cell")
cells_bg <- rownames(combined@meta.data)[combined@meta.data$celltype %in% acinar_celltype_candidates]
cells_bg <- cells_bg[nzchar(cells_bg)]
stopifnot(length(cells_bg) > 0)

counts_bg <- tryCatch(
  GetAssayData(combined, assay = bg_assay, layer = "counts"),
  error = function(e) GetAssayData(combined, assay = bg_assay, slot = "counts")
)

counts_bg <- counts_bg[, cells_bg, drop = FALSE]

expressed_idx <- Matrix::rowSums(counts_bg > 0) > 0
background_genes <- rownames(counts_bg)[expressed_idx] |> as.character() |> unique()
background_genes <- background_genes[nzchar(background_genes)]

cat(sprintf(
  "[INFO] （%s ；assay=%s；cells=%d）: n=%d\n",
  target_celltype, bg_assay, length(cells_bg), length(background_genes)
))

# ========= （）=========
stopifnot(file.exists(stable_genes_csv))
stable_tbl <- readr::read_csv(stable_genes_csv, show_col_types = FALSE)
stopifnot("gene" %in% colnames(stable_tbl))
stable_symbols <- unique(stable_tbl$gene)
stable_symbols <- stable_symbols[nzchar(stable_symbols)]
cat(sprintf("[INFO]  n=%d : %s\n", length(stable_symbols), stable_genes_csv))

# stable_symbols <- intersect(stable_symbols, background_genes)
# cat(sprintf("[INFO] （）: %d\n", length(stable_symbols)))
# stopifnot(length(stable_symbols) >= 5, length(background_genes) >= 100)

# ========= SYMBOL -> ENTREZ =========
orgdb <- if (species == "Mus musculus") {
  suppressPackageStartupMessages(library(org.Mm.eg.db))
  org.Mm.eg.db
} else {
  org.Hs.eg.db
}

map_to_entrez <- function(symbols, orgdb) {
  suppressMessages({
    mapped <- AnnotationDbi::select(
      orgdb, keys = symbols, keytype = "SYMBOL", columns = c("ENTREZID")
    )
  })
  mapped <- mapped[!is.na(mapped$ENTREZID), ]
  unique(mapped$ENTREZID)
}
gene_entrez <- map_to_entrez(stable_symbols, orgdb)
bg_entrez   <- map_to_entrez(background_genes, orgdb)


cat(sprintf("[INFO]  ENTREZ：=%d, =%d\n", length(gene_entrez), length(bg_entrez)))
stopifnot(length(gene_entrez) >= 5, length(bg_entrez) >= 100)

# =========  =========
save_enrich <- function(enrich_res, name, outdir) {
  df <- as.data.frame(enrich_res)
  if (is.null(df) || nrow(df) == 0) {
    message(paste("[WARN]", name, ""))
    return(invisible(NULL))
  }

  readr::write_csv(df, file.path(outdir, paste0(name, ".csv")))

  p1 <- dotplot(enrich_res, showCategory = 10) +
    ggtitle(name) +
    theme(
      axis.title.x = element_text(size = 20, face = "bold"),
      axis.title.y = element_text(size = 20, face = "bold"),
      axis.text.x  = element_text(size = 18, face = "bold"),
      axis.text.y  = element_text(size = 18, face = "bold"),
      plot.title   = element_text(size = 18, face = "bold"),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text  = element_text(size = 12, face = "bold")
    )
  ggsave(file.path(outdir, paste0(name, "_dotplot.png")), p1, width = 10, height =6, dpi = 600)

  p2 <- barplot(enrich_res, showCategory = 10) +
    ggtitle(name) +
    theme(
      axis.title.x = element_text(size = 20, face = "bold"),
      axis.title.y = element_text(size = 20, face = "bold"),
      axis.text.x  = element_text(size = 18, face = "bold"),
      axis.text.y  = element_text(size = 18, face = "bold"),
      plot.title   = element_text(size = 18, face = "bold"),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text  = element_text(size = 12, face = "bold")
    )
  ggsave(file.path(outdir, paste0(name, "_barplot.png")), p2, width = 10, height = 6, dpi = 600)
}

# ========= GO （BP/MF/CC）=========
ego_bp <- enrichGO(
  gene          = gene_entrez,
  universe      = bg_entrez,
  OrgDb         = orgdb,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)
save_enrich(ego_bp, "GO_BP", outdir)

ego_mf <- enrichGO(
  gene          = gene_entrez,
  universe      = bg_entrez,
  OrgDb         = orgdb,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)
save_enrich(ego_mf, "GO_MF", outdir)

ego_cc <- enrichGO(
  gene          = gene_entrez,
  universe      = bg_entrez,
  OrgDb         = orgdb,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.99,
  qvalueCutoff  = 0.99,
  readable      = TRUE
)
save_enrich(ego_cc, "GO_CC", outdir)

# ========= KEGG =========
organism <- ifelse(species == "Homo sapiens", "hsa", "mmu")
ekegg <- enrichKEGG(
  gene          = gene_entrez,
  universe      = bg_entrez,
  organism      = organism,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.99,
  qvalueCutoff  = 0.99
)
if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
  ekegg <- setReadable(ekegg, OrgDb = orgdb, keyType = "ENTREZID")
}
save_enrich(ekegg, "KEGG", outdir)

cat("[DONE] ORA ：", normalizePath(outdir), "\n")
cat("- GO_BP.csv / GO_MF.csv / GO_CC.csv\n- KEGG.csv\n")
cat("：*_dotplot.png, *_barplot.png（Top 20）\n")





# =========  ORA data.frame =========
library(dplyr)
library(tibble)
library(ggprism)

read_res <- function(path, onto) {
  if (!file.exists(path)) return(NULL)
  df <- readr::read_csv(path, show_col_types = FALSE)
  if (nrow(df) == 0) return(NULL)
  df$ONTOLOGY <- onto
  df
}

GO_BP  <- read_res(file.path(outdir, "GO_BP.csv"),  "BP")
GO_MF  <- read_res(file.path(outdir, "GO_MF.csv"),  "MF")

GO_CC  <- read_res(file.path(outdir, "GO_CC.csv"),  "CC")

KEGG   <- read_res(file.path(outdir, "KEGG.csv"),   "KEGG")

all_pathways <- bind_rows(GO_BP, GO_MF, GO_CC, KEGG)

# =========  Top（GO Top5 + KEGG Top5） =========
library(dplyr)
library(stringr)
library(ggplot2)

n_go   <- 5
n_kegg <- 5

plot_df <- all_pathways %>%
  mutate(
    pvalue = suppressWarnings(as.numeric(pvalue)),
    pvalue = ifelse(is.na(pvalue) | pvalue <= 0, 1e-300, pvalue),
    DB = ifelse(ONTOLOGY == "KEGG", "KEGG", "GO")
  ) %>%
  filter(pvalue < 0.05)

if (nrow(plot_df) == 0) {
  stop("No enriched pathways with pvalue < 0.05; skip GO/KEGG Top plot.")
}
go_top <- plot_df %>%
  filter(DB == "GO") %>%
  group_by(ONTOLOGY) %>%
  arrange(pvalue, desc(Count), .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup()


kegg_top <- plot_df %>%
  filter(DB == "KEGG") %>%
  arrange(pvalue, desc(Count)) %>%
  slice_head(n = n_kegg)

plot_df2 <- bind_rows(go_top, kegg_top) %>%
  mutate(
    score = -log10(pvalue),
    ONTOLOGY = factor(ONTOLOGY, levels = c("KEGG","BP","CC","MF"))
  ) %>%
  arrange(ONTOLOGY, desc(score))   # 

# 1) 
plot_df2 <- plot_df2 %>%
  mutate(
    Description_single = str_squish(Description),
    score = -log10(pvalue)
  )

# 2) （）
plot_df2$Description_single <- factor(
  plot_df2$Description_single,
  levels = rev(unique(plot_df2$Description_single))  # unique  levels 
)

plot_df2$Description_single <- factor(plot_df2$Description_single,
                                     levels = rev(plot_df2$Description_single))
# （BP/MF/CC/KEGG）
my_colors <- c(
  "BP"   = "#e6a0c4",
  "MF"   = "#c6cdf7",
  "CC"   = "#d8a499",
  "KEGG" = "#7294d4"
)

score_breaks <- pretty(c(0, max(plot_df2$score, na.rm = TRUE)), n = 5)

# Linux  Arial（ PDF ）， PDF 
# grid.Call(C_textBounds, ...)； sans  cairo_pdf 
plot_base_family <- "sans"

p_all <- ggplot(plot_df2, aes(x = score, y = Description_single)) +
  geom_col(aes(fill = ONTOLOGY), width = 0.7) +
  # （score，0）
  geom_text(aes(label = as.character(Description_single), x = 0.02),
            hjust = 0, size = 4.5, color = "black", family = plot_base_family, fontface = "bold") +
  scale_fill_manual(values = my_colors) +
  scale_x_continuous(
    breaks = score_breaks,
    labels = score_breaks,
    expand = c(0, 0),
    limits = c(0, max(plot_df2$score, na.rm = TRUE) * 1.15)
  ) +
  labs(
    title = NULL,
    x = expression(-log[10](pvalue)),
    y = NULL
  ) +
  theme_classic(base_family = plot_base_family) +
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




ggsave(file.path(outdir, "GO_KEGG_Top5_pvalue_insideLabel.png"),
       p_all, width = 7.2, height = 6, dpi = 500)

out_pdf <- file.path(outdir, "GO_KEGG_Top5_pvalue_insideLabel.pdf")
pdf_device <- if (exists("cairo_pdf", where = asNamespace("grDevices"), inherits = FALSE)) {
  grDevices::cairo_pdf
} else {
  "pdf"
}

if (identical(pdf_device, "pdf")) {
  ggsave(out_pdf, p_all, width = 7.2, height = 6, device = "pdf", useDingbats = FALSE)
} else {
  ggsave(out_pdf, p_all, width = 7.2, height = 6, device = pdf_device)
}



