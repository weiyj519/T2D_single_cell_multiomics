#!/usr/bin/env python3

"""
Script name: 04_make_table9_candidate_motif_peak_gene_triplets.py
Purpose: Cleaned public workflow script for the T2D single-cell multiome analysis.
Input: Paths and analysis settings are read from configs/config.yaml.
Output: Module-specific outputs are written under the configured results directory.
Main steps: Load configuration, run the original analysis logic, and export reproducible outputs.
Example command: python scripts/09_figures_and_tables/04_make_table9_candidate_motif_peak_gene_triplets.py
"""

from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src" / "utils"))
from config import ensure_dir, load_config, resolve_path

CONFIG = load_config()
RAW_DATA_DIR = resolve_path(CONFIG, "paths.raw_data_dir")
PROCESSED_DATA_DIR = resolve_path(CONFIG, "paths.processed_data_dir")
SEURAT_OBJECT_PATH = resolve_path(CONFIG, "paths.seurat_object_path")
RNA_EXPRESSION_PARQUET = resolve_path(CONFIG, "paths.RNA_EXPRESSION_PARQUET")
ATAC_PEAK_PARQUET = resolve_path(CONFIG, "paths.ATAC_PEAK_PARQUET")
RESULT_DIR = ensure_dir(resolve_path(CONFIG, "paths.result_dir"))
FIGURE_DIR = ensure_dir(resolve_path(CONFIG, "paths.figure_dir"))
SUPPLEMENTARY_TABLE_DIR = ensure_dir(resolve_path(CONFIG, "paths.supplementary_table_dir"))
MODEL_RESULT_DIR = ensure_dir(resolve_path(CONFIG, "paths.model_result_dir"))
RNA_MODEL_RESULT_DIR = ensure_dir(resolve_path(CONFIG, "paths.rna_model_result_dir"))
ATAC_MODEL_RESULT_DIR = ensure_dir(resolve_path(CONFIG, "paths.atac_model_result_dir"))
STAGE_GENE_DIR = ensure_dir(resolve_path(CONFIG, "paths.stage_gene_dir"))
STAGE_PEAK_DIR = ensure_dir(resolve_path(CONFIG, "paths.stage_peak_dir"))
DOWNSTREAM_RESULT_DIR = ensure_dir(resolve_path(CONFIG, "paths.downstream_result_dir"))
ATAC_REGULATORY_RESULT_DIR = ensure_dir(resolve_path(CONFIG, "paths.atac_regulatory_result_dir"))
COMMUNICATION_RESULT_DIR = ensure_dir(resolve_path(CONFIG, "paths.communication_result_dir"))
PSEUDOTIME_RESULT_DIR = ensure_dir(resolve_path(CONFIG, "paths.pseudotime_result_dir"))

"""
Make Supplementary Table 9:
Candidate motif-peak-gene regulatory triplets.

Run with Python 3, for example:
    python \
        scripts/09_figures_and_tables/make_table9_candidate_motif_peak_gene_triplets.py

Outputs:
    results/supplementary_tables/Table9_candidate_motif_peak_gene_regulatory_triplets.csv
    results/supplementary_tables/Table9_candidate_motif_peak_gene_regulatory_triplets.xlsx
    results/supplementary_tables/Table9_candidate_motif_peak_gene_regulatory_triplets_summary.csv
"""

import glob
import os
import re

import numpy as np
import pandas as pd


BASE_DIR = "results/downstream"
OUT_DIR = "results/supplementary_tables"

OUT_CSV = os.path.join(
    OUT_DIR, "Table9_candidate_motif_peak_gene_regulatory_triplets.csv"
)
OUT_XLSX = os.path.join(
    OUT_DIR, "Table9_candidate_motif_peak_gene_regulatory_triplets.xlsx"
)
OUT_SUMMARY = os.path.join(
    OUT_DIR, "Table9_candidate_motif_peak_gene_regulatory_triplets_summary.csv"
)

TRIPLET_FILENAME = "TF_peak_gene_triplets_RNA_consistent_abs_rnaLFC0.25_padj0.05.csv"


# Directory aliases are intentionally explicit so that glob search is robust
# to alpha/beta Greek-letter naming while still limiting the table to the six
# requested cell-type/comparison combinations.
INPUT_SPECS = [
    {
        "cell_type": "Acinar cell",
        "comparison": "CON vs PRE",
        "dir_aliases": ["acinar_CON_PRE"],
    },
    {
        "cell_type": "Acinar cell",
        "comparison": "PRE vs T2D",
        "dir_aliases": ["acinar_PRE_T2D"],
    },
    {
        "cell_type": "α cell",
        "comparison": "CON vs PRE",
        "dir_aliases": ["alpha_CON_PRE", "α_CON_PRE"],
    },
    {
        "cell_type": "α cell",
        "comparison": "PRE vs T2D",
        "dir_aliases": ["alpha_PRE_T2D", "α_PRE_T2D"],
    },
    {
        "cell_type": "β cell",
        "comparison": "CON vs PRE",
        "dir_aliases": ["beta_CON_PRE", "β_CON_PRE"],
    },
    {
        "cell_type": "β cell",
        "comparison": "PRE vs T2D",
        "dir_aliases": ["beta_PRE_T2D", "β_PRE_T2D"],
    },
]


OUTPUT_COLUMNS = [
    "cell_type",
    "comparison",
    "motif_id",
    "motif_name",
    "TF_gene",
    "chromVAR_log2FC",
    "chromVAR_adjusted_p_value",
    "chromVAR_log10_adjusted_p_value",
    "TF_direction",
    "n_peaks_with_motif",
    "fraction_peaks_with_motif",
    "peak_id",
    "peak_direction",
    "peak_log2FC",
    "peak_adjusted_p_value",
    "mean_abs_SHAP",
    "linked_gene",
    "link_score",
    "link_FDR",
    "RNA_log2FC",
    "RNA_adjusted_p_value",
    "gene_direction",
    "pass_significance_filter",
    "consistency_status",
    "combined_score",
    "source_file",
]


SUMMARY_COLUMNS = [
    "cell_type",
    "comparison",
    "source_file",
    "file_exists",
    "n_rows_raw",
    "n_rows_output",
    "n_unique_motifs",
    "n_unique_peaks",
    "n_unique_genes",
    "n_unique_triplets",
    "min_chromVAR_adjusted_p_value",
    "min_link_FDR",
    "min_RNA_adjusted_p_value",
    "max_abs_chromVAR_log2FC",
    "max_abs_peak_log2FC",
    "max_abs_RNA_log2FC",
    "max_combined_score",
]


NUMERIC_COLUMNS = [
    "chromVAR_log2FC",
    "chromVAR_adjusted_p_value",
    "chromVAR_log10_adjusted_p_value",
    "n_peaks_with_motif",
    "fraction_peaks_with_motif",
    "peak_log2FC",
    "peak_adjusted_p_value",
    "mean_abs_SHAP",
    "link_score",
    "link_FDR",
    "RNA_log2FC",
    "RNA_adjusted_p_value",
    "combined_score",
]


CELL_TYPE_ORDER = ["β cell", "α cell", "Acinar cell"]
COMPARISON_ORDER = ["CON vs PRE", "PRE vs T2D"]


def find_source_file(spec):
    """
    Find the triplet CSV for one cell-type/comparison combination.

    Returns the first sorted match. If both Latin and Greek alias directories
    exist for one combination, this deterministic choice avoids double-counting.
    """
    matches = []
    for alias in spec["dir_aliases"]:
        pattern = os.path.join(BASE_DIR, alias, "netlink", TRIPLET_FILENAME)
        matches.extend(glob.glob(pattern))

    matches = sorted(set(matches))
    if len(matches) > 1:
        print("[WARNING] Multiple matching files found; using the first one:")
        for match in matches:
            print(f"          {match}")
    return matches[0] if matches else ""


def get_column(df, candidates):
    """Return the first matching column name using case-insensitive matching."""
    lookup = {str(col).strip().lower(): col for col in df.columns}
    for candidate in candidates:
        key = candidate.strip().lower()
        if key in lookup:
            return lookup[key]
    return None


def is_missing(value):
    """Treat NaN, None, and blank-like strings as missing."""
    if pd.isna(value):
        return True
    text = str(value).strip()
    return text == "" or text.lower() in {"nan", "none", "na", "null"}


def clean_text_series(series):
    """Strip whitespace and convert blank-like strings to pandas NA."""
    if not isinstance(series, pd.Series):
        series = pd.Series(series)
    cleaned = series.astype("string").str.strip()
    return cleaned.replace(
        {
            "": pd.NA,
            "nan": pd.NA,
            "NaN": pd.NA,
            "None": pd.NA,
            "NA": pd.NA,
            "null": pd.NA,
            "NULL": pd.NA,
        }
    )


def normalize_peak_id(value):
    """
    Normalize peak strings to chr-start-end.

    Supported examples:
        chr1:1000-2000 -> chr1-1000-2000
        chr1_1000_2000 -> chr1-1000-2000
        chr1-1000-2000 -> chr1-1000-2000
    """
    if is_missing(value):
        return pd.NA

    peak = str(value).strip()

    colon_match = re.match(r"^(.+):(\d+)-(\d+)$", peak)
    if colon_match:
        chrom, start, end = colon_match.groups()
        return f"{chrom}-{start}-{end}"

    underscore_match = re.match(r"^(.+)_(\d+)_(\d+)$", peak)
    if underscore_match:
        chrom, start, end = underscore_match.groups()
        return f"{chrom}-{start}-{end}"

    dash_match = re.match(r"^(.+)-(\d+)-(\d+)$", peak)
    if dash_match:
        chrom, start, end = dash_match.groups()
        return f"{chrom}-{start}-{end}"

    return peak


def column_or_na(df, column_name):
    """Return a source column or an all-NA Series with the same index."""
    if column_name is None:
        return pd.Series(pd.NA, index=df.index)
    return df[column_name]


def first_nonmissing_series(*series_list):
    """Coalesce multiple Series from left to right."""
    if not series_list:
        return pd.Series(dtype="object")

    result = series_list[0].copy()
    for series in series_list[1:]:
        missing = result.isna() | result.astype("string").str.strip().isin(
            ["", "nan", "NaN", "None", "NA", "null", "NULL"]
        )
        result = result.where(~missing, series)
    return result


def to_numeric_if_present(df, columns):
    """Convert selected columns to numeric, coercing invalid values to NaN."""
    for col in columns:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def empty_summary_row(spec, source_file="", file_exists=False, n_rows_raw=0):
    """Create a summary row for missing or unusable inputs."""
    return {
        "cell_type": spec["cell_type"],
        "comparison": spec["comparison"],
        "source_file": source_file,
        "file_exists": file_exists,
        "n_rows_raw": n_rows_raw,
        "n_rows_output": 0,
        "n_unique_motifs": 0,
        "n_unique_peaks": 0,
        "n_unique_genes": 0,
        "n_unique_triplets": 0,
        "min_chromVAR_adjusted_p_value": np.nan,
        "min_link_FDR": np.nan,
        "min_RNA_adjusted_p_value": np.nan,
        "max_abs_chromVAR_log2FC": np.nan,
        "max_abs_peak_log2FC": np.nan,
        "max_abs_RNA_log2FC": np.nan,
        "max_combined_score": np.nan,
    }


def sort_and_deduplicate(table):
    """
    Sort and remove exact biological duplicate triplets.

    For duplicates, rows with non-missing combined_score rank first by the
    maximum combined_score. If combined_score is missing, link_FDR ascending
    determines the kept row.
    """
    if table.empty:
        return pd.DataFrame(columns=OUTPUT_COLUMNS)

    table = table.copy()
    table["_cell_order"] = pd.Categorical(
        table["cell_type"], categories=CELL_TYPE_ORDER, ordered=True
    )
    table["_comparison_order"] = pd.Categorical(
        table["comparison"], categories=COMPARISON_ORDER, ordered=True
    )
    table["_has_combined_score"] = table["combined_score"].notna().astype(int)

    table = table.sort_values(
        by=[
            "_cell_order",
            "_comparison_order",
            "motif_id",
            "peak_id",
            "linked_gene",
            "_has_combined_score",
            "combined_score",
            "link_FDR",
            "chromVAR_adjusted_p_value",
        ],
        ascending=[True, True, True, True, True, False, False, True, True],
        na_position="last",
    )

    table = table.drop_duplicates(
        subset=["cell_type", "comparison", "motif_id", "peak_id", "linked_gene"],
        keep="first",
    )

    table = table.sort_values(
        by=[
            "_cell_order",
            "_comparison_order",
            "combined_score",
            "link_FDR",
            "chromVAR_adjusted_p_value",
        ],
        ascending=[True, True, False, True, True],
        na_position="last",
    )

    return table.drop(
        columns=["_cell_order", "_comparison_order", "_has_combined_score"]
    )[OUTPUT_COLUMNS]


def process_one_file(spec):
    """Read, standardize, and summarize one triplet CSV."""
    label = f"{spec['cell_type']} | {spec['comparison']}"
    source_file = find_source_file(spec)

    print(f"[START] {label}")
    if source_file:
        print(f"        {source_file}")
    else:
        searched = ", ".join(spec["dir_aliases"])
        print(f"[WARNING] Missing input file for {label}; searched aliases: {searched}")
        return pd.DataFrame(columns=OUTPUT_COLUMNS), empty_summary_row(spec)

    if not os.path.exists(source_file):
        print(f"[WARNING] Missing input file: {source_file}")
        return (
            pd.DataFrame(columns=OUTPUT_COLUMNS),
            empty_summary_row(spec, source_file=source_file),
        )

    df = pd.read_csv(source_file)
    n_rows_raw = len(df)
    print(f"        raw rows: {n_rows_raw:,}")

    motif_col = get_column(df, ["motif"])
    tf_col = get_column(df, ["TF"])
    tf_gene_col = get_column(df, ["TF_gene"])
    peak_col = get_column(df, ["peak"])
    gene_col = get_column(df, ["gene"])
    strength_rna_col = get_column(df, ["strength_rna"])
    strength_col = get_column(df, ["strength"])

    motif_id = clean_text_series(column_or_na(df, motif_col))
    tf_name = clean_text_series(column_or_na(df, tf_col))
    tf_gene = clean_text_series(column_or_na(df, tf_gene_col))
    motif_name = clean_text_series(first_nonmissing_series(tf_name, tf_gene, motif_id))

    peak_id = column_or_na(df, peak_col).map(normalize_peak_id)
    linked_gene = clean_text_series(column_or_na(df, gene_col))

    combined_score = first_nonmissing_series(
        column_or_na(df, strength_rna_col), column_or_na(df, strength_col)
    )

    out = pd.DataFrame(
        {
            "cell_type": spec["cell_type"],
            "comparison": spec["comparison"],
            "motif_id": motif_id,
            "motif_name": motif_name,
            "TF_gene": tf_gene,
            "chromVAR_log2FC": column_or_na(df, get_column(df, ["chromvar_lfc"])),
            "chromVAR_adjusted_p_value": column_or_na(
                df, get_column(df, ["chromvar_padj"])
            ),
            "chromVAR_log10_adjusted_p_value": column_or_na(
                df, get_column(df, ["chromvar_log10padj"])
            ),
            "TF_direction": clean_text_series(
                column_or_na(df, get_column(df, ["TF_direction"]))
            ),
            "n_peaks_with_motif": column_or_na(
                df, get_column(df, ["n_peaks_with_motif"])
            ),
            "fraction_peaks_with_motif": column_or_na(
                df, get_column(df, ["frac_peaks_with_motif"])
            ),
            "peak_id": clean_text_series(peak_id),
            "peak_direction": clean_text_series(
                column_or_na(df, get_column(df, ["peak_direction"]))
            ),
            "peak_log2FC": column_or_na(df, get_column(df, ["peak_lfc"])),
            "peak_adjusted_p_value": column_or_na(df, get_column(df, ["peak_padj"])),
            "mean_abs_SHAP": column_or_na(df, get_column(df, ["mean_abs_SHAP"])),
            "linked_gene": linked_gene,
            "link_score": column_or_na(df, get_column(df, ["link_score"])),
            "link_FDR": column_or_na(df, get_column(df, ["link_fdr"])),
            "RNA_log2FC": column_or_na(df, get_column(df, ["rna_lfc"])),
            "RNA_adjusted_p_value": column_or_na(df, get_column(df, ["rna_padj"])),
            "gene_direction": clean_text_series(
                column_or_na(df, get_column(df, ["gene_direction"]))
            ),
            "pass_significance_filter": column_or_na(df, get_column(df, ["pass_sig"])),
            "consistency_status": "direction_consistent",
            "combined_score": combined_score,
            "source_file": source_file,
        }
    )

    out = to_numeric_if_present(out, NUMERIC_COLUMNS)

    before_drop = len(out)
    out = out.dropna(subset=["peak_id", "motif_id", "linked_gene"]).copy()
    print(f"        dropped missing peak/motif/gene rows: {before_drop - len(out):,}")

    out = sort_and_deduplicate(out)

    triplet_keys = ["motif_id", "peak_id", "linked_gene"]
    summary = {
        "cell_type": spec["cell_type"],
        "comparison": spec["comparison"],
        "source_file": source_file,
        "file_exists": True,
        "n_rows_raw": n_rows_raw,
        "n_rows_output": len(out),
        "n_unique_motifs": out["motif_id"].nunique(dropna=True),
        "n_unique_peaks": out["peak_id"].nunique(dropna=True),
        "n_unique_genes": out["linked_gene"].nunique(dropna=True),
        "n_unique_triplets": out[triplet_keys].drop_duplicates().shape[0],
        "min_chromVAR_adjusted_p_value": out["chromVAR_adjusted_p_value"].min(
            skipna=True
        ),
        "min_link_FDR": out["link_FDR"].min(skipna=True),
        "min_RNA_adjusted_p_value": out["RNA_adjusted_p_value"].min(skipna=True),
        "max_abs_chromVAR_log2FC": out["chromVAR_log2FC"].abs().max(skipna=True),
        "max_abs_peak_log2FC": out["peak_log2FC"].abs().max(skipna=True),
        "max_abs_RNA_log2FC": out["RNA_log2FC"].abs().max(skipna=True),
        "max_combined_score": out["combined_score"].max(skipna=True),
    }

    print(f"[DONE]  output rows: {len(out):,}")
    return out[OUTPUT_COLUMNS], summary


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    tables = []
    summaries = []

    for spec in INPUT_SPECS:
        table, summary = process_one_file(spec)
        tables.append(table)
        summaries.append(summary)

    if tables:
        combined = pd.concat(tables, ignore_index=True)
    else:
        combined = pd.DataFrame(columns=OUTPUT_COLUMNS)

    combined = sort_and_deduplicate(combined)
    summary_df = pd.DataFrame(summaries, columns=SUMMARY_COLUMNS)

    combined.to_csv(OUT_CSV, index=False)
    combined.to_excel(OUT_XLSX, index=False, engine="openpyxl")
    summary_df.to_csv(OUT_SUMMARY, index=False)

    print("\n[FINISHED] Table 9 candidate motif-peak-gene regulatory triplets")
    print(f"  CSV:     {OUT_CSV}")
    print(f"  XLSX:    {OUT_XLSX}")
    print(f"  Summary: {OUT_SUMMARY}")
    print(f"  Total output rows: {len(combined):,}")


if __name__ == "__main__":
    main()
