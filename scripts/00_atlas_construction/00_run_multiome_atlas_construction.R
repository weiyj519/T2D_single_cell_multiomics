#!/usr/bin/env Rscript

# Script name: 00_run_multiome_atlas_construction.R
# Purpose: Cleaned public workflow script for the T2D single-cell multiome analysis.
# Input: Paths and analysis settings are read from configs/config.yaml.
# Output: Module-specific outputs are written under the configured results directory.
# Main steps: Load configuration, run the original analysis logic, and export reproducible outputs.
# Example command: Rscript scripts/00_atlas_construction/00_run_multiome_atlas_construction.R

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
## GSE200044 multiome （ QC +  doublet）
############################################################

library(Seurat)
library(dplyr)
library(stringr)
library(purrr)
library(patchwork)
library(Signac)
library(ggplot2)
library(harmony)
library(DoubletFinder)
library(Matrix)
library(qs)
library(ensembldb)
library(EnsDb.Hsapiens.v75)
library(RColorBrewer)
library(future)
options(future.globals.maxSize = 300 * 1024^3)  # 300 GB
plan("multicore", workers = 4)  
#-----------------------------
# 0. 
#-----------------------------
data_path <- raw_data_dir
sample_dirs <- list.dirs(data_path, full.names = TRUE, recursive = FALSE)

#-----------------------------
# 1.  10x multiome 
#-----------------------------
read_10x_multiome <- function(sample_dir) {
  counts <- Read10X(data.dir = sample_dir)
  
  # RNA
  seurat_obj <- CreateSeuratObject(
    counts = counts$`Gene Expression`,
    assay  = "RNA"
  )
  
  # ATAC（ peaks ， fragments）
  seurat_obj[["ATAC"]] <- CreateChromatinAssay(
    counts = counts$Peaks,
    sep    = c(":", "-"),
    genome = "hg19"   
  )
  
  # meta 
  seurat_obj$orig.ident <- basename(sample_dir)
  group <- ifelse(grepl("CON", sample_dir), "CON",
                  ifelse(grepl("PRE", sample_dir), "PRE", "T2D"))
  seurat_obj$group <- group
  
  return(seurat_obj)
}

#-----------------------------
# 2. 
#-----------------------------
objs    <- lapply(sample_dirs, read_10x_multiome)
samples <- basename(sample_dirs)
names(objs) <- samples

combined <- merge(objs[[1]], y = objs[-1], add.cell.ids = samples)
qsave(combined, "combined_multiome_raw_merged.qs")

combined <- qread("combined_multiome_raw_merged.qs")
# ----------------------------------------------------------
# 3.  QC（RNA + ATAC）+ DoubletFinder  doublet
# ----------------------------------------------------------

# 3.1 RNA  QC 
DefaultAssay(combined) <- "RNA"
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")

# 3.2 ATAC  QC （ peak ）
DefaultAssay(combined) <- "ATAC"
atac_counts <- GetAssayData(combined, assay = "ATAC", slot = "counts")
combined$nCount_ATAC   <- Matrix::colSums(atac_counts)
combined$nFeature_ATAC <- Matrix::colSums(atac_counts > 0)

# 3.3 RNA + ATAC （）
combined <- subset(
  combined,
  subset =
    nCount_RNA    >= 1000  &          # RNA UMI  1000
    nFeature_RNA  >= 500   &          # 
    nFeature_RNA  <= 7500  &          # （ doublet）
    percent.mt    <  10    &          #  <10%
    nCount_ATAC   >= 1000  &          # ATAC 
    nFeature_ATAC >= 1000             # ATAC peaks 
)

# 3.4 DoubletFinder  doublet（ RNA）
DefaultAssay(combined) <- "RNA"

combined <- JoinLayers(combined)  # ， counts/data
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, nfeatures = 2000)
combined <- ScaleData(combined)
combined <- RunPCA(combined, npcs = 30)

# paramSweep  pK
sweep.res.list <- paramSweep(combined, PCs = 1:30, sct = FALSE)
sweep.stats    <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn          <- find.pK(sweep.stats)

best.pK <- bcmvn %>%
  dplyr::filter(BCmetric == max(BCmetric)) %>%
  dplyr::slice(1) %>%
  dplyr::pull(pK) %>%
  as.character() %>%
  as.numeric()

#  doublet （7.5%，）
doublet.rate <- 0.075
nExp_poi     <- round(doublet.rate * ncol(combined))

combined <- doubletFinder(
  seu        = combined,
  PCs        = 1:30,
  pN         = 0.25,
  pK         = best.pK,
  nExp       = nExp_poi,
  sct        = FALSE
)

df.col <- grep("DF.classification", colnames(combined@meta.data), value = TRUE)
combined <- combined[, combined@meta.data[[df.col]] == "Singlet"]

# QC +  doublet 
qsave(combined, "combined_multiome_qc_doublet.qs")
dim(combined)
table(combined$group)

combined <- qread("combined_multiome_qc_doublet.qs")
#----------------------------------------------------------
# 4. RNA （）
#----------------------------------------------------------
split_list <- SplitObject(combined, split.by = "orig.ident")
split_list <- map(split_list, NormalizeData)
split_list <- map(split_list, FindVariableFeatures,
                  selection.method = "vst", nfeatures = 2000)
anchors   <- FindIntegrationAnchors(object.list = split_list, dims = 1:30)
combined  <- IntegrateData(anchorset = anchors, dims = 1:30)
qsave(combined, "combined_integrated.qs")

#----------------------------------------------------------
# 5. RNA  + UMAP
#----------------------------------------------------------
DefaultAssay(combined) <- "integrated"
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined,
                    reduction      = "pca",
                    dims           = 1:30,
                    reduction.name = "rna.umap")
qsave(combined, "combined_rna_dr.qs")
#----------------------------------------------------------
# 6. ATAC （LSI + Harmony）
#----------------------------------------------------------
DefaultAssay(combined) <- "ATAC"
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = "q0")
combined <- RunSVD(combined)

lsi.embeddings <- Embeddings(combined, reduction = "lsi")
harmony.embeddings <- HarmonyMatrix(
  data_mat  = lsi.embeddings,
  meta_data = combined@meta.data,
  vars_use  = "orig.ident",
  do_pca    = FALSE
)

combined[["harmony"]] <- CreateDimReducObject(
  embeddings = harmony.embeddings,
  key       = "harmony_",
  assay     = "ATAC"
)

combined <- RunUMAP(
  combined,
  reduction      = "harmony",
  dims           = 1:30,
  reduction.name = "atac.umap"
)
qsave(combined, "combined_atac_harmony.qs")
#----------------------------------------------------------
# 7. WNN  +  + WNN-UMAP
#----------------------------------------------------------
combined <- FindMultiModalNeighbors(
  combined,
  reduction.list = list("pca", "harmony"),
  dims.list      = list(1:30, 1:30)
)

combined <- FindClusters(combined, graph.name = "wsnn", resolution = 0.1)

combined <- RunUMAP(
  combined,
  nn.name        = "weighted.nn",
  reduction.name = "wnn.umap"
)

qsave(combined, "combined_multiome_seurat_final.qs")
combined <- qread("combined_multiome_seurat_final.qs")

Idents(combined)[1:10]          #  cluster
table(Idents(combined))         #  cluster 
head(combined@meta.data)        # （ grouporig.ident ）


############################################################
## Step ：（ < 50 ）
############################################################

#  cluster 
cluster_sizes <- table(Idents(combined))
print(cluster_sizes)

# （， < 50 ）
small_clusters <- names(cluster_sizes[cluster_sizes < 50])
cat("Small clusters detected:", small_clusters, "\n")

# 
combined <- subset(
  combined,
  idents = setdiff(levels(combined), small_clusters)
)

#  levels
combined$seurat_clusters <- droplevels(combined$seurat_clusters)

# 
qsave(combined, "combined_noTinyClusters.qs")
cat("Saved cleaned object: combined_noTinyClusters.qs\n")

############################################################
## Step ： clean 
############################################################
cat("Cluster sizes after cleaning:\n")
print(table(Idents(combined)))



############################################################
## Step ： GENCODE v19 / Ensembl75  protein_coding, marker（ cluster ）
############################################################

DefaultAssay(combined) <- "RNA"
combined <- JoinLayers(combined)
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined, features = rownames(combined))

edb <- EnsDb.Hsapiens.v75

annot <- genes(
  edb,
  columns = c("gene_name", "gene_biotype")
) %>% as.data.frame()

rna_features <- rownames(combined@assays$RNA)

gene_info <- data.frame(gene = rna_features) %>%
  left_join(annot, by = c("gene" = "gene_name"))

protein_coding_genes <- gene_info$gene[gene_info$gene_biotype == "protein_coding"]
protein_coding_genes <- unique(na.omit(protein_coding_genes))
protein_coding_genes <- protein_coding_genes[protein_coding_genes %in% rna_features]

cat("Number of protein_coding genes used for marker detection:",
    length(protein_coding_genes), "\n")

##  protein_coding  marker，
markers_clean <- FindAllMarkers(
  combined,
  group.by        = "seurat_clusters",
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25,
  features        = protein_coding_genes   # ★ ：
)

write.csv(markers_clean, "markers_protein_coding_only.csv", row.names = FALSE)
write.csv(gene_info, "RNA_features_with_biotype_gencode19.csv", row.names = FALSE)



# 1. 
cluster_to_celltype <- c(
  "0" = "α cell",
  "1" = "β cell",
  "2" = "Acinar cell",
  "3" = "δ cell",
  "4" = "α cell",
  "5" = "α cell",
  "6" = "PP cell",
  "7" = "Stellate cell",
  "8" = "Ductal cell",
  "9" = "Immune cell"
)

# 2.  meta data
Idents(combined) <- combined$seurat_clusters
celltype_vec <- cluster_to_celltype[as.character(combined$seurat_clusters)]
names(celltype_vec) <- colnames(combined)          #  Cells(combined)
combined$celltype <- celltype_vec


# 3. （）
celltype_levels <- c(
  "α cell", "β cell", "Acinar cell", "δ cell",
  "PP cell", "Stellate cell", "Ductal cell", "Immune cell"
)
combined$celltype <- factor(combined$celltype, levels = celltype_levels)

#  celltype
Idents(combined) <- combined$celltype



suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(grid)
  library(RColorBrewer)
})
combined <- qread(seurat_object_path)
#  atac.umap， harmony  UMAP
if (!"atac.umap" %in% Reductions(combined)) {
  message(" 'atac.umap'， ATAC Harmony embeddings  UMAP...")
  DefaultAssay(combined) <- "ATAC"
  combined <- RunUMAP(
    combined,
    reduction      = "harmony",
    dims           = 1:30,
    reduction.name = "atac.umap"
  )
}

stopifnot("celltype" %in% colnames(combined@meta.data))
stopifnot("rna.umap" %in% Reductions(combined))
stopifnot("atac.umap" %in% Reductions(combined))
stopifnot("wnn.umap" %in% Reductions(combined))


# 1. （ RColorBrewer::Pastel1）
pal <- setNames(
  brewer.pal(9, "Paired")[seq_along(celltype_levels)],
  celltype_levels
)


# 2. （RNA/ATAC/WNN ）
plot_family <- "Arial"

make_celltype_umap <- function(seu_obj, reduction_name) {
  umap_mat <- Embeddings(seu_obj, reduction = reduction_name)[, 1:2]
  colnames(umap_mat) <- c("UMAP1", "UMAP2")

  plot_df <- as.data.frame(umap_mat) %>%
    mutate(celltype = seu_obj$celltype)

  label_df <- plot_df %>%
    group_by(celltype) %>%
    summarise(
      UMAP1 = median(UMAP1, na.rm = TRUE),
      UMAP2 = median(UMAP2, na.rm = TRUE),
      .groups = "drop"
    )

  point_layer <- if (requireNamespace("ggrastr", quietly = TRUE)) {
    ggrastr::geom_point_rast(size = 0.7, alpha = 0.95, raster.dpi = 600)
  } else {
    geom_point(size = 0.7, alpha = 0.95)
  }

  ggplot(plot_df, aes(UMAP1, UMAP2, colour = celltype)) +
    point_layer +
    geom_text(
      data = label_df,
      aes(x = UMAP1, y = UMAP2, label = celltype),
      inherit.aes = FALSE,
      size = 5.2,
      fontface = "bold",
      family = plot_family,
      colour = "black"
    ) +
    coord_fixed(clip = "off") +
    scale_colour_manual(values = pal, name = "Cell type", drop = FALSE) +
    theme_void(base_size = 18) +
    theme(
      text = element_text(family = plot_family),
      plot.margin = margin(20, 20, 20, 25),
      legend.position = "none",
      legend.title = element_text(face = "bold", size = 11, family = plot_family),
      legend.text = element_text(size = 10, family = plot_family),
      plot.title = element_blank()
    )
}

gg_wnn_umap <- make_celltype_umap(combined, "wnn.umap")
gg_rna_umap <- make_celltype_umap(combined, "rna.umap")
gg_atac_umap <- make_celltype_umap(combined, "atac.umap")

print(gg_wnn_umap)
ggsave("results/figures/atlas/umap_celltype.pdf", gg_wnn_umap,
  width = 5.5, height = 5.5, bg = "white",
  device = cairo_pdf, family = "Arial")

ggsave("results/figures/atlas/rna_umap_celltype_style.pdf", gg_rna_umap,
       width = 5.5, height = 5.5, bg = "white",
       device = cairo_pdf, family = "Arial")

ggsave("results/figures/atlas/atac_umap_celltype_style.pdf", gg_atac_umap,
       width = 7.5, height = 7.5, bg = "white",
       device = cairo_pdf, family = "Arial")



############################################################
## Step ： CON / PRE / T2D 
## ： +  + 
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(scales)
  library(RColorBrewer)
})

#-----------------------------
# 1. 
#-----------------------------
meta_df <- combined@meta.data %>%
  as.data.frame() %>%
  dplyr::select(group, celltype) %>%
  dplyr::mutate(
    group = factor(as.character(group), levels = c("CON", "PRE", "T2D")),
    celltype = factor(as.character(celltype), levels = names(pal))
  ) %>%
  dplyr::filter(!is.na(group), !is.na(celltype))

#  group × celltype ，
all_comb <- expand.grid(
  group = levels(meta_df$group),
  celltype = levels(meta_df$celltype),
  stringsAsFactors = FALSE
) %>%
  mutate(
    group = factor(group, levels = levels(meta_df$group)),
    celltype = factor(celltype, levels = levels(meta_df$celltype))
  )

#-----------------------------
# 2.  group  celltype 
#-----------------------------
group_cell_df <- meta_df %>%
  count(group, celltype, name = "n") %>%
  right_join(all_comb, by = c("group", "celltype")) %>%
  mutate(n = ifelse(is.na(n), 0, n)) %>%
  group_by(group) %>%
  mutate(
    group_total = sum(n),
    percentage = ifelse(group_total > 0, n / group_total, 0)
  ) %>%
  ungroup()

group_total_df <- group_cell_df %>%
  group_by(group) %>%
  summarise(total_n = sum(n), .groups = "drop")

#-----------------------------
# 3.  celltype （）
#-----------------------------
celltype_sum <- meta_df %>%
  count(celltype, name = "total_n") %>%
  right_join(
    data.frame(celltype = factor(levels(meta_df$celltype), levels = levels(meta_df$celltype))),
    by = "celltype"
  ) %>%
  mutate(total_n = ifelse(is.na(total_n), 0, total_n))




#-----------------------------
# 4. 
#-----------------------------
dotplot_main <- ggplot(
  group_cell_df,
  aes(x = group, y = celltype, size = percentage, color = celltype)
) +
  geom_point(alpha = 0.85) +
  scale_color_manual(values = pal, limits = names(pal), drop = FALSE) +
  scale_size_continuous(
    range = c(2, 12),
    labels = percent_format(accuracy = 1)
  ) +
  labs(
    x = NULL,
    y = NULL,
    size = "Cell proportion",
    color = "Cell type"
  ) +
  theme_classic(base_size = 17) +
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.position = "none",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13),
    plot.margin = margin(5, 5, 5, 5)
  )

dotplot_main <- rasterize_points(dotplot_main, dpi = 500)

#-----------------------------
# 5. （ group ）
#-----------------------------
barplot_top <- ggplot(
  group_cell_df,
  aes(x = group, y = n, fill = celltype)
) +
  geom_col(width = 0.75) +
  geom_text(
    data = group_total_df,
    aes(x = group, y = total_n, label = as.character(total_n)),
    inherit.aes = FALSE,
    vjust = -0.3,
    size = 4.8,
    color = "black"
  ) +
  scale_fill_manual(values = pal, limits = names(pal), drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Number of cells") +
  theme_classic(base_size = 17) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    plot.margin = margin(5, 5, 2, 5)
  )

#-----------------------------
# 6. （ celltype ）
#-----------------------------
barplot_right <- ggplot(
  celltype_sum,
  aes(x = total_n, y = celltype, fill = celltype)
) +
  geom_col(width = 0.75) +
  geom_text(
    aes(label = ifelse(total_n > 0, as.character(total_n), "")),
    hjust = -0.1,
    size = 4.8,
    color = "black"
  ) +
  scale_fill_manual(values = pal, limits = names(pal), drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Number of cells", y = NULL) +
  theme_classic(base_size = 17) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    plot.margin = margin(5, 5, 5, 2)
  )

#-----------------------------
# 7. 
#-----------------------------
top_panel <- barplot_top + plot_spacer() +
  plot_layout(widths = c(4.8, 1.4))

bottom_panel <- dotplot_main + barplot_right +
  plot_layout(widths = c(4.8, 1.4), guides = "collect")

final_prop_plot <- top_panel / bottom_panel +
  plot_layout(heights = c(1.2, 4.2)) &
  theme(legend.position = "right")

# 
print(final_prop_plot)

#-----------------------------
# 8. 
#-----------------------------
ggsave(
  filename = "results/figures/atlas/celltype_group_dotplot_main.pdf",
  plot = dotplot_main,
  width = 7, height = 4.5,device = cairo_pdf, family = "Arial"
)

ggsave(
  filename = "results/figures/atlas/celltype_group_barplot_top.pdf",
  plot = barplot_top,
  width = 7, height = 2.5,device = cairo_pdf, family = "Arial"
)

ggsave(
  filename = "results/figures/atlas/celltype_group_barplot_right.pdf",
  plot = barplot_right,
  width = 6, height = 4.5,device = cairo_pdf, family = "Arial"
)

ggsave(
  filename = "results/figures/atlas/celltype_group_composition_combined.pdf",
  plot = final_prop_plot,
  width = 11, height = 9,device = cairo_pdf, family = "Arial"
)

ggsave(
  filename = "results/figures/atlas/celltype_group_composition_combined.png",
  plot = final_prop_plot,
  width = 11, height = 9, dpi = 500, bg = "white"
)





############################################################
## Step：DotPlot  celltype  marker （RNA）
## ：dotplot_marker_genes.png
############################################################

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(RColorBrewer)
})

# ====== （ wnn.R  celltype ）======

stopifnot(all(c("celltype", "orig.ident") %in% colnames(combined@meta.data)))
combined$celltype <- droplevels(combined$celltype)

# ======  RNA  data layer（Seurat v5 ）======
DefaultAssay(combined) <- "RNA"
combined <- JoinLayers(combined)

if (!"data" %in% Layers(combined[["RNA"]])) {
  combined <- NormalizeData(combined, verbose = FALSE)
}

# ====== celltype （ levels）======
celltype_levels <- c(
  "α cell", "β cell", "δ cell", "PP cell",
  "Acinar cell", "Ductal cell", "Stellate cell", "Immune cell"
)
#  level，； level
celltype_levels_use <- intersect(celltype_levels, levels(combined$celltype))
combined$celltype <- factor(combined$celltype, levels = celltype_levels_use)



# ====== marker （ 8 ；）======
marker_list <- list(
  
  "α cell"        = c("GCG"),
  "β cell"        = c("INS"),
  "δ cell"        = c("SST"),
  "PP cell"       = c("PPY"),
  "Acinar cell"   = c("REG1A"),
  "Ductal cell"   = c("CFTR"),
  "Stellate cell" = c("COL6A2"),
  "Immune cell"   = c("PTPRC")
)

marker_genes <- unique(unlist(marker_list))
marker_genes <- intersect(marker_genes, rownames(combined))

if (length(marker_genes) == 0) {
  stop(" marker ，（/ID ）")
}

# ====== DotPlot （）======
dot_colours <- colorRampPalette(brewer.pal(9, "Blues"))(7)

dot_plot <- DotPlot(
  combined,
  features  = marker_genes,
  group.by  = "celltype",
  dot.scale = 7,
  scale     = TRUE,
  cols      = brewer.pal(9, "Blues")[c(2, 9)]
) +
  scale_colour_gradientn(
    colours = dot_colours,
    name    = "Avg. expression (scaled)"
  ) +
  guides(
    size = guide_legend(
      title = "Pct. expressing",
      override.aes = list(colour = "grey50")
    )
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 17) +
  theme(
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.7),
    axis.text.x       = element_text(angle = 40, hjust = 1, vjust = 1, size =14),
    axis.text.y       = element_text(face = "bold", size = 16),
    plot.title        = element_blank(),
    legend.title      = element_text(size = 16, face = "bold"),
    legend.text       = element_text(size = 14),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    strip.placement   = "outside",
    axis.line         = element_blank(),
    panel.spacing     = unit(1, "mm")
  )

out_pdf <- file.path("results/figures/atlas/dotplot_marker_genes.pdf")
ggsave(
  filename = out_pdf,
  plot     = dot_plot,
  width    = 7,
  height   = 6,
  device   = cairo_pdf,
  family   = "Arial",
  bg       = "white"
)
cat("Saved:", out_pdf, "\n")

