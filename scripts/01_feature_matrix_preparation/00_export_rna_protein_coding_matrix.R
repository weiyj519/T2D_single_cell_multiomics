#!/usr/bin/env Rscript

# Script name: 00_export_rna_protein_coding_matrix.R
# Purpose: Cleaned public workflow script for the T2D single-cell multiome analysis.
# Input: Paths and analysis settings are read from configs/config.yaml.
# Output: Module-specific outputs are written under the configured results directory.
# Main steps: Load configuration, run the original analysis logic, and export reproducible outputs.
# Example command: Rscript scripts/01_feature_matrix_preparation/00_export_rna_protein_coding_matrix.R

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
## Export filtered RNA expression matrix for ML
##
## Source:
##   data/processed/combined_final_celltyped.qs
##
## Output:
##   cell metadata + RNA expression columns in parquet format
##
## Gene filter:
##   keep protein-coding autosomal genes only, then remove
##   MT-, RPS*, and RPL* genes.
############################################################

parse_args <- function(args) {
  opts <- list(
    input_qs = seurat_object_path,
    output_parquet = rna_expression_parquet,
    assay = "RNA",
    layer = "data",
    annotation_tsv = "",
    grch = "37"
  )

  if (length(args) == 0) {
    return(opts)
  }

  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]

    if (arg %in% c("-h", "--help")) {
      cat(
        "Usage:\n",
        "  Rscript export_expr_protein_coding_filtered.R [options]\n\n",
        "Options:\n",
        "  --input_qs PATH          Seurat qs file\n",
        "  --output_parquet PATH    Output parquet file\n",
        "  --assay NAME             Assay name, default RNA\n",
        "  --layer NAME             RNA layer/slot, default data\n",
        "  --annotation_tsv PATH    Optional gene annotation TSV/CSV\n",
        "  --grch 37|38             Ensembl genome build for biomaRt, default 37\n\n",
        "Annotation file columns can include hgnc_symbol/gene_name/gene, ",
        "ensembl_gene_id, gene_biotype/biotype, and chromosome_name/chromosome.\n",
        sep = ""
      )
      quit(save = "no", status = 0)
    }

    if (grepl("^--[^=]+=", arg)) {
      key <- sub("^--([^=]+)=.*$", "\\1", arg)
      val <- sub("^--[^=]+=", "", arg)
    } else if (grepl("^--", arg)) {
      key <- sub("^--", "", arg)
      if (i == length(args)) {
        stop("Missing value for argument: ", arg)
      }
      i <- i + 1
      val <- args[[i]]
    } else {
      stop("Unknown argument: ", arg)
    }

    key <- gsub("-", "_", key)
    if (!key %in% names(opts)) {
      stop("Unknown option: --", key)
    }
    opts[[key]] <- val
    i <- i + 1
  }

  opts
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))

required_packages <- c("qs", "Seurat", "Matrix", "tibble", "arrow")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ", paste(missing_packages, collapse = ", "),
    "\nPlease run this script in the same R environment used for wnn.R."
  )
}

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(Matrix)
  library(tibble)
  library(arrow)
})

normalize_chr <- function(x) {
  x <- as.character(x)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  x <- toupper(x)
  x
}

standardize_annotation <- function(ann) {
  names(ann) <- tolower(names(ann))

  pick_col <- function(candidates) {
    hit <- candidates[candidates %in% names(ann)]
    if (length(hit) == 0) {
      return(rep(NA_character_, nrow(ann)))
    }
    as.character(ann[[hit[[1]]]])
  }

  data.frame(
    hgnc_symbol = pick_col(c("hgnc_symbol", "gene_name", "external_gene_name", "symbol", "gene")),
    ensembl_gene_id = pick_col(c("ensembl_gene_id", "gene_id", "ensembl_id")),
    gene_biotype = pick_col(c("gene_biotype", "gene_type", "biotype")),
    chromosome_name = pick_col(c("chromosome_name", "chromosome", "chr", "seqnames")),
    stringsAsFactors = FALSE
  )
}

read_annotation_file <- function(path) {
  message("Reading annotation file: ", path)
  if (!file.exists(path)) {
    stop("Annotation file does not exist: ", path)
  }

  ann <- tryCatch(
    utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) {
      utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    }
  )

  standardize_annotation(ann)
}

query_biomart_annotation <- function(genes, grch = "37") {
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop(
      "Package biomaRt is required when --annotation_tsv is not provided.\n",
      "Either install biomaRt or pass a local annotation file with ",
      "--annotation_tsv."
    )
  }

  message("Querying Ensembl biomaRt annotation, GRCh", grch, "...")
  mart <- tryCatch(
    {
      if (grch == "37") {
        biomaRt::useEnsembl(
          biomart = "genes",
          dataset = "hsapiens_gene_ensembl",
          GRCh = 37
        )
      } else {
        biomaRt::useEnsembl(
          biomart = "genes",
          dataset = "hsapiens_gene_ensembl"
        )
      }
    },
    error = function(e) {
      if (grch == "37") {
        biomaRt::useMart(
          biomart = "ENSEMBL_MART_ENSEMBL",
          dataset = "hsapiens_gene_ensembl",
          host = "https://grch37.ensembl.org"
        )
      } else {
        stop(e)
      }
    }
  )

  attributes_use <- c(
    "hgnc_symbol",
    "ensembl_gene_id",
    "gene_biotype",
    "chromosome_name"
  )

  ann_symbol <- biomaRt::getBM(
    attributes = attributes_use,
    filters = "hgnc_symbol",
    values = genes,
    mart = mart
  )

  ensembl_like <- unique(genes[grepl("^ENSG", genes)])
  if (length(ensembl_like) > 0) {
    ann_ensembl <- biomaRt::getBM(
      attributes = attributes_use,
      filters = "ensembl_gene_id",
      values = ensembl_like,
      mart = mart
    )
    ann_symbol <- rbind(ann_symbol, ann_ensembl)
  }

  standardize_annotation(unique(ann_symbol))
}

build_gene_filter_table <- function(genes, ann) {
  ann <- ann[
    !is.na(ann$gene_biotype) &
      !is.na(ann$chromosome_name) &
      (nzchar(ann$hgnc_symbol) | nzchar(ann$ensembl_gene_id)),
    ,
    drop = FALSE
  ]

  ann$chromosome_norm <- normalize_chr(ann$chromosome_name)
  ann$match_key_symbol <- ann$hgnc_symbol
  ann$match_key_ensembl <- ann$ensembl_gene_id

  ann_by_gene <- lapply(genes, function(gene) {
    hit <- ann[
      ann$match_key_symbol == gene | ann$match_key_ensembl == gene,
      ,
      drop = FALSE
    ]

    symbol <- gene
    if (nrow(hit) > 0 && any(nzchar(hit$hgnc_symbol))) {
      symbol <- hit$hgnc_symbol[nzchar(hit$hgnc_symbol)][[1]]
    }

    is_protein_coding_autosome <- any(
      hit$gene_biotype == "protein_coding" &
        hit$chromosome_norm %in% as.character(1:22),
      na.rm = TRUE
    )

    data.frame(
      gene = gene,
      hgnc_symbol = symbol,
      matched_annotation = nrow(hit) > 0,
      protein_coding_autosome = is_protein_coding_autosome,
      is_mt_gene = grepl("^MT-", symbol, ignore.case = FALSE),
      is_rps_rpl_gene = grepl("^(RPS|RPL)", symbol, ignore.case = FALSE),
      chromosomes = if (nrow(hit) > 0) {
        paste(sort(unique(hit$chromosome_norm)), collapse = ";")
      } else {
        ""
      },
      biotypes = if (nrow(hit) > 0) {
        paste(sort(unique(hit$gene_biotype)), collapse = ";")
      } else {
        ""
      },
      stringsAsFactors = FALSE
    )
  })

  gene_filter <- do.call(rbind, ann_by_gene)
  gene_filter$keep <- with(
    gene_filter,
    protein_coding_autosome & !is_mt_gene & !is_rps_rpl_gene
  )
  gene_filter
}

get_assay_matrix <- function(seu, assay, layer) {
  DefaultAssay(seu) <- assay
  seu <- tryCatch(JoinLayers(seu), error = function(e) seu)

  layer_names <- tryCatch(Layers(seu[[assay]]), error = function(e) character())
  if (layer == "data" && length(layer_names) > 0 && !"data" %in% layer_names) {
    message("RNA data layer not found; running NormalizeData() first.")
    seu <- NormalizeData(seu, assay = assay, verbose = FALSE)
  }

  mat <- tryCatch(
    GetAssayData(seu, assay = assay, layer = layer),
    error = function(e) {
      GetAssayData(seu, assay = assay, slot = layer)
    }
  )

  list(object = seu, matrix = mat)
}

message("Input qs: ", opts$input_qs)
message("Output parquet: ", opts$output_parquet)
message("Assay/layer: ", opts$assay, "/", opts$layer)

if (!file.exists(opts$input_qs)) {
  stop("Input qs does not exist: ", opts$input_qs)
}
dir.create(dirname(opts$output_parquet), showWarnings = FALSE, recursive = TRUE)

combined <- qread(opts$input_qs)
stopifnot(opts$assay %in% names(combined@assays))

required_meta_cols <- c("orig.ident", "seurat_clusters", "celltype", "group")
missing_meta_cols <- setdiff(required_meta_cols, colnames(combined@meta.data))
if (length(missing_meta_cols) > 0) {
  stop("Missing metadata columns: ", paste(missing_meta_cols, collapse = ", "))
}

mat_result <- get_assay_matrix(combined, opts$assay, opts$layer)
combined <- mat_result$object
expr_mat <- mat_result$matrix

message(
  "Raw expression matrix: ",
  nrow(expr_mat), " genes x ", ncol(expr_mat), " cells"
)

genes <- rownames(expr_mat)
ann <- if (nzchar(opts$annotation_tsv)) {
  read_annotation_file(opts$annotation_tsv)
} else {
  query_biomart_annotation(genes, opts$grch)
}

gene_filter <- build_gene_filter_table(genes, ann)
keep_genes <- gene_filter$gene[gene_filter$keep]

if (length(keep_genes) == 0) {
  stop("No genes remained after protein-coding/autosome/MT/RPS/RPL filtering.")
}

message("Annotated genes: ", sum(gene_filter$matched_annotation), " / ", nrow(gene_filter))
message("Protein-coding autosomal genes: ", sum(gene_filter$protein_coding_autosome))
message("Dropped MT- genes: ", sum(gene_filter$is_mt_gene))
message("Dropped RPS/RPL genes: ", sum(gene_filter$is_rps_rpl_gene))
message("Final kept genes: ", length(keep_genes))

gene_filter_path <- sub("\\.parquet$", "_gene_filter_summary.csv", opts$output_parquet)
if (identical(gene_filter_path, opts$output_parquet)) {
  gene_filter_path <- paste0(opts$output_parquet, "_gene_filter_summary.csv")
}
utils::write.csv(gene_filter, gene_filter_path, row.names = FALSE)
message("Wrote gene filter summary: ", gene_filter_path)

expr_mat_keep <- expr_mat[keep_genes, , drop = FALSE]

meta_df <- combined@meta.data[, required_meta_cols, drop = FALSE] |>
  tibble::rownames_to_column(var = "cell_id")

meta_df$cell_id <- as.character(meta_df$cell_id)
meta_df$orig.ident <- as.character(meta_df$orig.ident)
meta_df$seurat_clusters <- as.character(meta_df$seurat_clusters)
meta_df$celltype <- as.character(meta_df$celltype)
meta_df$group <- as.character(meta_df$group)

message("Converting filtered matrix to cell x gene table...")
expr_df <- as.data.frame(as.matrix(Matrix::t(expr_mat_keep)), check.names = FALSE)
expr_df <- tibble::rownames_to_column(expr_df, var = "cell_id")
expr_df$cell_id <- as.character(expr_df$cell_id)

export_df <- merge(meta_df, expr_df, by = "cell_id", all.x = TRUE, sort = FALSE)
export_df <- export_df[, c("cell_id", required_meta_cols, keep_genes), drop = FALSE]

message(
  "Export table: ",
  nrow(export_df), " cells x ", ncol(export_df), " columns"
)

arrow::write_parquet(export_df, opts$output_parquet)
message("Done: ", opts$output_parquet)
