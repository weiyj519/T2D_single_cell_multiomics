#!/usr/bin/env Rscript

# Script name: 00_run_rna_de_for_enrichment.R
# Purpose: Cleaned public workflow script for the T2D single-cell multiome analysis.
# Input: Paths and analysis settings are read from configs/config.yaml.
# Output: Module-specific outputs are written under the configured results directory.
# Main steps: Load configuration, run the original analysis logic, and export reproducible outputs.
# Example command: Rscript scripts/05_functional_enrichment/00_run_rna_de_for_enrichment.R

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
  library(dplyr)
  library(ggplot2)
})

#----------------------------
# 0) 
#----------------------------
qs_path <- seurat_object_path
combined <- qread(qs_path)

stopifnot(all(c("celltype", "group") %in% colnames(combined@meta.data)))

#----------------------------
# 1)  acinar cell， CON/PRE
#----------------------------
acinar_celltype_candidates <- c("Acinar cell", "Acinar_cell", "acinar cell", "acinar_cell")
acinar <- subset(
  combined,
  subset = celltype %in% acinar_celltype_candidates & group %in% c("CON", "PRE")
)

# ，
acinar$group <- droplevels(factor(acinar$group, levels = c("CON", "PRE")))

cat("acinar cell cells:", ncol(acinar), "\n")
print(table(acinar$group))

if (ncol(acinar) == 0) {
  stop(
    "No cells matched. Please check combined@meta.data$celltype values. ",
    "Tried: ", paste(acinar_celltype_candidates, collapse = ", ")
  )
}

#----------------------------
# 2) RNA （ RNA， integrated）
#    ： v5 ， JoinLayers
#----------------------------
DefaultAssay(acinar) <- "RNA"
acinar <- JoinLayers(acinar)

#  NormalizeData ； data layer
acinar <- NormalizeData(acinar, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)


#----------------------------
# 3)  group：CON vs PRE
#----------------------------
Idents(acinar) <- acinar$group

#----------------------------
# 4) （CON vs PRE）
#     test.use = "wilcox"；latent.vars （）
#----------------------------
# （，）
candidate_latent <- c("nCount_RNA", "percent.mt", "orig.ident")
latent_use <- candidate_latent[candidate_latent %in% colnames(acinar@meta.data)]
cat("latent vars used:", paste(latent_use, collapse = ", "), "\n")

find_markers_data <- function(object, ...) {
  # Seurat v5:  layer=； slot=
  tryCatch(
    FindMarkers(object = object, layer = "data", ...),
    error = function(e) FindMarkers(object = object, slot = "data", ...)
  )
}

de_con_vs_pre <- find_markers_data(
  object          = acinar,
  ident.1         = "PRE",
  ident.2         = "CON",
  assay           = "RNA",
  test.use        = "wilcox",
  logfc.threshold = 0.25,  #  log2FC 
  min.pct         = 0.1,
  latent.vars     = latent_use
)

# 
de_con_vs_pre <- de_con_vs_pre %>%
  tibble::rownames_to_column("gene") %>%
  arrange(p_val_adj, desc(avg_log2FC))

# 
out_dir <- "results/downstream/acinar_CON_PRE/GO_RNA"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_csv <- file.path(out_dir, "DE_acinar_CON_vs_PRE_RNA.csv")
write.csv(de_con_vs_pre, out_csv, row.names = FALSE)
cat("Saved:", out_csv, "\n")

#----------------------------
# 5)  Top genes
#----------------------------
sig <- de_con_vs_pre %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.25)

cat("Significant genes (padj<0.05 & log2FC>0.25):", nrow(sig), "\n")
print(head(sig, 20))

#----------------------------
# 6) ： + 
#----------------------------
volcano <- de_con_vs_pre %>%
  mutate(
    neglog10_padj = -log10(p_val_adj + 1e-300),
    sig = (p_val_adj < 0.05 & avg_log2FC > 0.25)
  )

p_vol <- ggplot(volcano, aes(x = avg_log2FC, y = neglog10_padj)) +
  geom_point(aes(color = sig), size = 0.8, alpha = 0.8) +
  scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "#D62728")) +
  theme_classic(base_size = 12) +
  labs(
    title = "acinar cell RNA DE: CON vs PRE",
    x = "avg_log2FC (PRE / CON)",
    y = "-log10(adj p-value)"
  )

ggsave(file.path(out_dir, "volcano_acinar_CON_vs_PRE.png"), p_vol, width = 6, height = 5, dpi = 300)



