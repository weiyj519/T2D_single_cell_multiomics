#!/usr/bin/env Rscript

# Script name: 00_run_de_method_benchmark.R
# Purpose: Compare differential expression strategies before downstream interpretation.
# Input: Seurat object and protein-coding RNA feature matrix from configs/config.yaml.
# Output: DE benchmark summaries and method-specific Wilcoxon/DESeq2 result tables.
# Main steps: Load the atlas, run single-cell-level Wilcoxon and donor-level pseudobulk DESeq2, filter protein-coding DE genes, and summarize counts.
# Example command: Rscript scripts/04_de_method_benchmark/00_run_de_method_benchmark.R

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

required_packages <- c(
  "qs", "Seurat", "dplyr", "tidyr", "stringr", "Matrix",
  "DESeq2", "ggplot2", "readr", "tibble", "arrow"
)
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ", paste(missing_packages, collapse = ", "),
    "\nPlease install them before running this script."
  )
}

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(Matrix)
  library(DESeq2)
  library(ggplot2)
  library(readr)
  library(arrow)
})

input_qs <- seurat_object_path
outdir <- "results/de_method_benchmark"
protein_coding_parquet <- rna_expression_parquet
fc_threshold <- 0.25
padj_threshold <- 0.05
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("Input object: ", input_qs)
message("Output directory: ", outdir)
message("Protein-coding source: ", protein_coding_parquet)
message("Output filter: abs(log2FC) > ", fc_threshold,
        ", adjusted P value < ", padj_threshold,
        ", protein-coding genes only")

comparisons <- list(
  CON_vs_PRE = list(case = "PRE", ref = "CON"),
  PRE_vs_T2D = list(case = "T2D", ref = "PRE"),
  CON_vs_T2D = list(case = "T2D", ref = "CON")
)

celltypes_use <- c("β cell", "α cell", "Acinar cell")

if (!file.exists(protein_coding_parquet)) {
  stop("Protein-coding parquet does not exist: ", protein_coding_parquet)
}

metadata_cols <- c("cell_id", "orig.ident", "seurat_clusters", "celltype", "group")
protein_coding_table <- arrow::read_parquet(protein_coding_parquet, as_data_frame = FALSE)
protein_coding_genes <- setdiff(names(protein_coding_table), metadata_cols)
rm(protein_coding_table)
protein_coding_genes <- unique(protein_coding_genes[!is.na(protein_coding_genes) & nzchar(protein_coding_genes)])

if (length(protein_coding_genes) == 0) {
  stop("No protein-coding genes were found in: ", protein_coding_parquet)
}
message("Protein-coding gene whitelist size: ", length(protein_coding_genes))

old_csv <- list.files(outdir, pattern = "\\.csv$", full.names = TRUE)
if (length(old_csv) > 0) {
  message("Removing previous CSV outputs from: ", outdir)
  unlink(old_csv)
}

safe_name <- function(x) {
  x %>%
    str_replace_all("β", "Beta") %>%
    str_replace_all("α", "Alpha") %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("^_+|_+$", "")
}

collapse_messages <- function(x) {
  x <- unique(x[!is.na(x) & nzchar(x)])
  if (length(x) == 0) {
    return("")
  }
  paste(x, collapse = " | ")
}

add_warning <- function(warnings, txt) {
  c(warnings, conditionMessage(txt))
}

find_first_col <- function(cols, candidates) {
  hit <- candidates[candidates %in% cols]
  if (length(hit) == 0) {
    return(NA_character_)
  }
  hit[[1]]
}

get_assay_layer <- function(object, assay = "RNA", layer = "counts", slot = NULL) {
  if (is.null(slot)) {
    slot <- layer
  }
  out <- tryCatch(
    GetAssayData(object, assay = assay, layer = layer),
    error = function(e) {
      GetAssayData(object, assay = assay, slot = slot)
    }
  )
  out
}

get_assay_names <- function(object) {
  assay_names <- tryCatch(names(object@assays), error = function(e) character())
  if (length(assay_names) == 0) {
    assay_names <- tryCatch(as.character(SeuratObject::Assays(object)), error = function(e) character())
  }
  assay_names
}

assay_has_data_layer <- function(object, assay = "RNA") {
  assay_obj <- object[[assay]]
  layer_names <- tryCatch(Layers(assay_obj), error = function(e) NULL)
  if (!is.null(layer_names)) {
    return("data" %in% layer_names)
  }
  has_slot <- tryCatch({
    dat <- GetAssayData(object, assay = assay, slot = "data")
    nrow(dat) > 0 && ncol(dat) > 0
  }, error = function(e) FALSE)
  has_slot
}

filter_coding_degs <- function(res, fc_candidates, padj_candidates, method_name) {
  if (!"gene" %in% colnames(res)) {
    stop(method_name, " output did not contain a gene column.")
  }

  fc_col <- find_first_col(colnames(res), fc_candidates)
  padj_col <- find_first_col(colnames(res), padj_candidates)

  if (is.na(fc_col)) {
    stop(method_name, " output did not contain a log2FC column.")
  }
  if (is.na(padj_col)) {
    stop(method_name, " output did not contain an adjusted P value column.")
  }

  filtered <- res %>%
    mutate(
      .fc = suppressWarnings(as.numeric(.data[[fc_col]])),
      .padj = suppressWarnings(as.numeric(.data[[padj_col]]))
    ) %>%
    filter(
      gene %in% protein_coding_genes,
      !is.na(.fc),
      !is.na(.padj),
      abs(.fc) > fc_threshold,
      .padj < padj_threshold
    ) %>%
    select(-.fc, -.padj)

  list(
    result = filtered,
    total = nrow(filtered),
    up = sum(filtered[[fc_col]] > 0, na.rm = TRUE),
    down = sum(filtered[[fc_col]] < 0, na.rm = TRUE),
    fc_col = fc_col,
    padj_col = padj_col,
    n_input_genes = nrow(res)
  )
}

run_one_combination <- function(seu, celltype_i, comparison_name, case_group, ref_group) {
  messages <- character()
  status <- "success"
  ct_safe <- safe_name(celltype_i)

  wilcox_file <- file.path(outdir, paste0("wilcoxon_", ct_safe, "_", comparison_name, ".csv"))
  pb_file <- file.path(outdir, paste0("pseudobulk_DESeq2_", ct_safe, "_", comparison_name, ".csv"))

  empty_row <- tibble::tibble(
    celltype = celltype_i,
    comparison = comparison_name,
    ref_group = ref_group,
    case_group = case_group,
    n_cells_ref = NA_integer_,
    n_cells_case = NA_integer_,
    n_donors_ref = NA_integer_,
    n_donors_case = NA_integer_,
    wilcoxon_n_deg_total = NA_integer_,
    wilcoxon_n_deg_up = NA_integer_,
    wilcoxon_n_deg_down = NA_integer_,
    pseudobulk_n_deg_total = NA_integer_,
    pseudobulk_n_deg_up = NA_integer_,
    pseudobulk_n_deg_down = NA_integer_,
    wilcoxon_result_file = NA_character_,
    pseudobulk_result_file = NA_character_,
    status = NA_character_,
    message = NA_character_
  )

  tryCatch({
    md <- seu@meta.data
    cells_use <- rownames(md)[md$celltype == celltype_i & md$group %in% c(ref_group, case_group)]

    if (length(cells_use) == 0) {
      empty_row$status <- "skipped"
      empty_row$message <- "No cells found for this celltype and comparison."
      return(empty_row)
    }

    sub <- subset(seu, cells = cells_use)
    sub_md <- sub@meta.data %>%
      tibble::rownames_to_column("cell") %>%
      mutate(
        group = as.character(group),
        orig.ident = as.character(orig.ident)
      )

    n_cells_ref <- sum(sub_md$group == ref_group, na.rm = TRUE)
    n_cells_case <- sum(sub_md$group == case_group, na.rm = TRUE)
    n_donors_ref <- n_distinct(sub_md$orig.ident[sub_md$group == ref_group])
    n_donors_case <- n_distinct(sub_md$orig.ident[sub_md$group == case_group])

    out <- empty_row %>%
      mutate(
        n_cells_ref = n_cells_ref,
        n_cells_case = n_cells_case,
        n_donors_ref = n_donors_ref,
        n_donors_case = n_donors_case
      )

    if (n_cells_ref < 20 || n_cells_case < 20) {
      status <- "skipped"
      messages <- c(
        messages,
        paste0("Wilcoxon skipped because n_cells_ref=", n_cells_ref,
               " and n_cells_case=", n_cells_case, "; both must be >= 20.")
      )
    } else {
      wilcox_warnings <- character()
      wilcox_res <- withCallingHandlers(
        FindMarkers(
          object = sub,
          assay = "RNA",
          group.by = "group",
          ident.1 = case_group,
          ident.2 = ref_group,
          test.use = "wilcox",
          min.pct = 0.05,
          logfc.threshold = 0
        ),
        warning = function(w) {
          wilcox_warnings <<- add_warning(wilcox_warnings, w)
          invokeRestart("muffleWarning")
        }
      )

      wilcox_res <- wilcox_res %>%
        as.data.frame() %>%
        tibble::rownames_to_column("gene") %>%
        as_tibble()

      filtered_wilcox <- filter_coding_degs(
        wilcox_res,
        fc_candidates = c("avg_log2FC", "avg_logFC"),
        padj_candidates = c("p_val_adj"),
        method_name = "Wilcoxon"
      )

      write_csv(filtered_wilcox$result, wilcox_file)
      out <- out %>%
        mutate(
          wilcoxon_n_deg_total = filtered_wilcox$total,
          wilcoxon_n_deg_up = filtered_wilcox$up,
          wilcoxon_n_deg_down = filtered_wilcox$down,
          wilcoxon_result_file = wilcox_file
        )
      messages <- c(
        messages,
        paste0("Wilcoxon input genes: ", filtered_wilcox$n_input_genes,
               "; filtered coding DEGs: ", filtered_wilcox$total)
      )
      messages <- c(messages, paste0("Wilcoxon warnings: ", collapse_messages(wilcox_warnings)))
    }

    if (n_donors_ref < 2 || n_donors_case < 2) {
      status <- "skipped"
      messages <- c(
        messages,
        paste0("Pseudobulk DESeq2 skipped because n_donors_ref=", n_donors_ref,
               " and n_donors_case=", n_donors_case, "; both must be >= 2.")
      )
    } else {
      pb_warnings <- character()
      pb_res <- withCallingHandlers({
        counts <- get_assay_layer(sub, assay = "RNA", layer = "counts", slot = "counts")
        counts <- counts[, sub_md$cell, drop = FALSE]

        pb_meta <- sub_md %>%
          select(cell, orig.ident, group) %>%
          mutate(
            donor = orig.ident,
            pseudobulk_id = paste(orig.ident, group, sep = "_")
          )

        pb_levels <- unique(pb_meta$pseudobulk_id)
        pb_factor <- factor(pb_meta$pseudobulk_id, levels = pb_levels)
        aggregation_matrix <- sparse.model.matrix(~ 0 + pb_factor)
        colnames(aggregation_matrix) <- levels(pb_factor)

        pb_counts <- counts %*% aggregation_matrix
        pb_counts <- as.matrix(pb_counts)
        storage.mode(pb_counts) <- "integer"

        col_data <- pb_meta %>%
          distinct(pseudobulk_id, donor, group) %>%
          arrange(match(pseudobulk_id, colnames(pb_counts))) %>%
          as.data.frame()
        rownames(col_data) <- col_data$pseudobulk_id
        col_data$group <- factor(col_data$group, levels = c(ref_group, case_group))

        pb_counts <- pb_counts[, rownames(col_data), drop = FALSE]
        keep <- rowSums(pb_counts) >= 10
        pb_counts <- pb_counts[keep, , drop = FALSE]

        if (nrow(pb_counts) == 0) {
          stop("No genes remained after pseudobulk rowSums(counts) >= 10 filtering.")
        }

        dds <- DESeqDataSetFromMatrix(
          countData = pb_counts,
          colData = col_data,
          design = ~ group
        )
        dds$group <- relevel(dds$group, ref = ref_group)
        dds <- DESeq(dds)

        res <- results(dds, contrast = c("group", case_group, ref_group))
        res <- as.data.frame(res) %>%
          tibble::rownames_to_column("gene") %>%
          as_tibble() %>%
          arrange(padj)
        res
      }, warning = function(w) {
        pb_warnings <<- add_warning(pb_warnings, w)
        invokeRestart("muffleWarning")
      })

      filtered_pb <- filter_coding_degs(
        pb_res,
        fc_candidates = c("log2FoldChange"),
        padj_candidates = c("padj"),
        method_name = "Pseudobulk DESeq2"
      )

      write_csv(filtered_pb$result, pb_file)
      out <- out %>%
        mutate(
          pseudobulk_n_deg_total = filtered_pb$total,
          pseudobulk_n_deg_up = filtered_pb$up,
          pseudobulk_n_deg_down = filtered_pb$down,
          pseudobulk_result_file = pb_file
        )
      messages <- c(
        messages,
        paste0("Pseudobulk DESeq2 input genes: ", filtered_pb$n_input_genes,
               "; filtered coding DEGs: ", filtered_pb$total)
      )
      messages <- c(messages, paste0("Pseudobulk DESeq2 warnings: ", collapse_messages(pb_warnings)))
    }

    out$status <- status
    out$message <- collapse_messages(messages)
    out
  }, error = function(e) {
    empty_row$status <- "error"
    empty_row$message <- conditionMessage(e)
    empty_row
  })
}

message("Reading Seurat object...")
seu <- qread(input_qs)

required_meta <- c("group", "celltype", "orig.ident")
missing_meta <- setdiff(required_meta, colnames(seu@meta.data))
if (length(missing_meta) > 0) {
  stop("Missing required meta.data columns: ", paste(missing_meta, collapse = ", "))
}

assay_names <- get_assay_names(seu)
message("\nAssays found: ", paste(assay_names, collapse = ", "))
if (!"RNA" %in% assay_names) {
  stop("RNA assay was not found in the Seurat object.")
}

message("\nBasic group distribution:")
print(table(seu@meta.data$group, useNA = "ifany"))

message("\nBasic celltype distribution:")
print(table(seu@meta.data$celltype, useNA = "ifany"))

message("\nBasic orig.ident distribution:")
print(table(seu@meta.data$orig.ident, useNA = "ifany"))

DefaultAssay(seu) <- "RNA"

message("\nChecking RNA assay layers...")
rna_layers_before <- tryCatch(Layers(seu[["RNA"]]), error = function(e) NULL)
if (!is.null(rna_layers_before)) {
  message("RNA layers before JoinLayers: ", paste(rna_layers_before, collapse = ", "))
  if (length(rna_layers_before) > 1) {
    message("Joining RNA assay layers for Seurat v5 object...")
    seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
  }
}

if (!assay_has_data_layer(seu, assay = "RNA")) {
  message("RNA data layer/slot not found; running NormalizeData...")
  seu <- NormalizeData(seu, assay = "RNA", verbose = FALSE)
} else {
  message("RNA data layer/slot found; skipping NormalizeData.")
}

message("\nRunning DE benchmark...")
summary_list <- list()
i <- 1L
for (celltype_i in celltypes_use) {
  for (comparison_name in names(comparisons)) {
    case_group <- comparisons[[comparison_name]]$case
    ref_group <- comparisons[[comparison_name]]$ref
    message("Processing: ", celltype_i, " | ", comparison_name,
            " (case=", case_group, ", ref=", ref_group, ")")
    summary_list[[i]] <- run_one_combination(
      seu = seu,
      celltype_i = celltype_i,
      comparison_name = comparison_name,
      case_group = case_group,
      ref_group = ref_group
    )
    i <- i + 1L
  }
}

summary_tbl <- bind_rows(summary_list) %>%
  select(
    celltype,
    comparison,
    ref_group,
    case_group,
    n_cells_ref,
    n_cells_case,
    n_donors_ref,
    n_donors_case,
    wilcoxon_n_deg_total,
    wilcoxon_n_deg_up,
    wilcoxon_n_deg_down,
    pseudobulk_n_deg_total,
    pseudobulk_n_deg_up,
    pseudobulk_n_deg_down,
    wilcoxon_result_file,
    pseudobulk_result_file,
    status,
    message
  )

summary_file <- file.path(outdir, "de_method_gene_count_summary.csv")
write_csv(summary_tbl, summary_file)

plot_tbl <- summary_tbl %>%
  select(celltype, comparison, wilcoxon_n_deg_total, pseudobulk_n_deg_total) %>%
  pivot_longer(
    cols = c(wilcoxon_n_deg_total, pseudobulk_n_deg_total),
    names_to = "method",
    values_to = "deg_number"
  ) %>%
  mutate(
    method = recode(
      method,
      wilcoxon_n_deg_total = "Wilcoxon",
      pseudobulk_n_deg_total = "pseudobulk DESeq2"
    ),
    deg_number = replace_na(as.numeric(deg_number), 0),
    label_x = paste(celltype, comparison, sep = " | "),
    label_x = factor(label_x, levels = unique(label_x)),
    comparison = factor(comparison, levels = names(comparisons)),
    celltype = factor(celltype, levels = celltypes_use)
  )

barplot_file_pdf <- file.path(outdir, "de_method_gene_count_barplot.pdf")
barplot_file_png <- file.path(outdir, "de_method_gene_count_barplot.png")
dotplot_file_pdf <- file.path(outdir, "de_method_gene_count_dotplot.pdf")
dotplot_file_png <- file.path(outdir, "de_method_gene_count_dotplot.png")

barplot <- ggplot(plot_tbl, aes(x = label_x, y = deg_number + 1, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_y_log10(name = "Filtered coding DEG number + 1 (log10 scale)") +
  labs(x = "Celltype | comparison", fill = "Method") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave(barplot_file_pdf, barplot, width = 11, height = 6)
ggsave(barplot_file_png, barplot, width = 11, height = 6, dpi = 300)

dotplot <- ggplot(plot_tbl, aes(x = comparison, y = celltype, size = deg_number, color = method)) +
  geom_point(alpha = 0.85) +
  scale_size_continuous(range = c(1.5, 10), name = "Filtered coding DEG number") +
  labs(x = "Comparison", y = "Celltype", color = "Method") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(dotplot_file_pdf, dotplot, width = 8, height = 4.5)
ggsave(dotplot_file_png, dotplot, width = 8, height = 4.5, dpi = 300)

message("\nSummary table:")
print(summary_tbl, n = Inf, width = Inf)

message("\nSaved outputs:")
message("Summary CSV: ", summary_file)
message("Barplot PDF: ", barplot_file_pdf)
message("Barplot PNG: ", barplot_file_png)
message("Dotplot PDF: ", dotplot_file_pdf)
message("Dotplot PNG: ", dotplot_file_png)
message("Per-comparison result CSV files are saved in: ", outdir)
