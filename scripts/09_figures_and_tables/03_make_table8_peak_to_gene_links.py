#!/usr/bin/env python3

"""
Script name: 03_make_table8_peak_to_gene_links.py
Purpose: Cleaned public workflow script for the T2D single-cell multiome analysis.
Input: Paths and analysis settings are read from configs/config.yaml.
Output: Module-specific outputs are written under the configured results directory.
Main steps: Load configuration, run the original analysis logic, and export reproducible outputs.
Example command: python scripts/09_figures_and_tables/03_make_table8_peak_to_gene_links.py
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
Make Supplementary Table 8: Peak-to-gene links.

Run with the xgb310 conda environment, for example:
    python \
        scripts/09_figures_and_tables/make_table8_peak_to_gene_links.py

Outputs:
    results/supplementary_tables/Table8_peak_to_gene_links.csv
    results/supplementary_tables/Table8_peak_to_gene_links.xlsx
    results/supplementary_tables/Table8_peak_to_gene_links_summary.csv
"""

import os
import re

import numpy as np
import pandas as pd


OUT_DIR = "results/supplementary_tables"
OUT_CSV = os.path.join(OUT_DIR, "Table8_peak_to_gene_links.csv")
OUT_XLSX = os.path.join(OUT_DIR, "Table8_peak_to_gene_links.xlsx")
OUT_SUMMARY = os.path.join(OUT_DIR, "Table8_peak_to_gene_links_summary.csv")


INPUT_SPECS = [
    {
        "cell_type": "Acinar cell",
        "comparison": "PRE vs T2D",
        "source_file": "results/downstream/acinar_PRE_T2D/linkpeak/acinar_LinkPeaks_hg19_featurePeaks_allGenes_SIG_fdr0.1.csv",
    },
    {
        "cell_type": "Acinar cell",
        "comparison": "CON vs PRE",
        "source_file": "results/downstream/acinar_CON_PRE/linkpeak/Acinar_cell__CON_vs_PRE__LinkPeaks_hg19__featurePeaks_allGenes_SIG_fdr0.1.csv",
    },
    {
        "cell_type": "α cell",
        "comparison": "CON vs PRE",
        "source_file": "results/downstream/alpha_CON_PRE/linkpeak/CON_vs_PRE_boruta_peaks_ALLLINKS_FDR0.10_ABS_SCORE0.00.csv",
    },
    {
        "cell_type": "α cell",
        "comparison": "PRE vs T2D",
        "source_file": "results/downstream/alpha_PRE_T2D/linkpeak/PRE_vs_T2D_boruta_peaks_ALLLINKS_FDR0.10_ABS_SCORE0.00.csv",
    },
    {
        "cell_type": "β cell",
        "comparison": "CON vs PRE",
        "source_file": "results/downstream/beta_CON_PRE/linkpeak/β_cell__CON_vs_PRE__LinkPeaks_hg19__featurePeaks_allGenes_SIG_fdr0.1.csv",
    },
    {
        "cell_type": "β cell",
        "comparison": "PRE vs T2D",
        "source_file": "results/downstream/beta_PRE_T2D/GO_ATAC/PRE_vs_T2D_boruta_peaks_ALLLINKS_FDR0.1_SCORE0_POS.csv",
    },
]


OUTPUT_COLUMNS = [
    "cell_type",
    "comparison",
    "peak_id",
    "linked_gene",
    "link_score",
    "link_p_value",
    "link_FDR",
]


FDR_THRESHOLD = 0.1


SUMMARY_COLUMNS = [
    "cell_type",
    "comparison",
    "source_file",
    "file_exists",
    "n_rows_raw",
    "n_rows_output",
    "n_unique_peaks",
    "n_unique_genes",
    "min_link_FDR",
    "max_abs_link_score",
]


def get_column(df, candidates):
    """Return the first matching column name, using case-insensitive matching."""
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


def normalize_peak_id(value):
    """
    Normalize peak strings to chr-start-end.

    Supported examples:
        chr1:1000-2000 -> chr1-1000-2000
        chr1_1000_2000 -> chr1-1000-2000
        chr1-1000-2000 -> chr1-1000-2000
    """
    if is_missing(value):
        return np.nan

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


def coord_peak_id(row, seq_col, start_col, end_col):
    """Build chr-start-end from coordinate columns when peak is unavailable."""
    if seq_col is None or start_col is None or end_col is None:
        return np.nan

    chrom = row.get(seq_col)
    start = row.get(start_col)
    end = row.get(end_col)
    if is_missing(chrom) or is_missing(start) or is_missing(end):
        return np.nan

    return f"{str(chrom).strip()}-{format_coordinate(start)}-{format_coordinate(end)}"


def format_coordinate(value):
    """Format numeric-looking coordinates without trailing .0."""
    if is_missing(value):
        return ""
    try:
        number = float(value)
        if np.isfinite(number) and number.is_integer():
            return str(int(number))
    except (TypeError, ValueError):
        pass
    return str(value).strip()


def clean_text_series(series):
    """Strip whitespace and convert blank-like strings to NaN."""
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


def to_numeric_if_present(df, columns):
    """Convert columns to numeric when present; invalid values become NaN."""
    for col in columns:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def empty_summary_row(spec, file_exists, n_rows_raw=0):
    """Create a summary row for missing or unusable inputs."""
    return {
        "cell_type": spec["cell_type"],
        "comparison": spec["comparison"],
        "source_file": spec["source_file"],
        "file_exists": file_exists,
        "n_rows_raw": n_rows_raw,
        "n_rows_output": 0,
        "n_unique_peaks": 0,
        "n_unique_genes": 0,
        "min_link_FDR": np.nan,
        "max_abs_link_score": np.nan,
    }


def process_one_file(spec):
    """Read and standardize one LinkPeaks result file."""
    source_file = spec["source_file"]
    label = f"{spec['cell_type']} | {spec['comparison']}"
    print(f"[START] {label}")
    print(f"        {source_file}")

    if not os.path.exists(source_file):
        print(f"[WARNING] Missing input file: {source_file}")
        return pd.DataFrame(columns=OUTPUT_COLUMNS), empty_summary_row(spec, False)

    df = pd.read_csv(source_file)
    n_rows_raw = len(df)
    print(f"        raw rows: {n_rows_raw:,}")

    # Locate columns with tolerant matching. These cover the observed LinkPeaks exports.
    peak_col = get_column(df, ["peak", "peak_id", "peak_norm"])
    seq_col = get_column(df, ["seqnames", "seqname", "chrom", "chr", "chromosome"])
    start_col = get_column(df, ["start", "peak_start"])
    end_col = get_column(df, ["end", "peak_end"])
    gene_col = get_column(df, ["gene", "linked_gene", "gene_name", "symbol"])
    score_col = get_column(df, ["score", "link_score"])
    pvalue_col = get_column(df, ["pvalue", "p_value", "p.val", "pval", "link_p_value"])
    fdr_col = get_column(df, ["fdr", "padj", "qvalue", "q_value", "link_fdr"])

    # Best-effort numeric conversion for source columns, including SHAP columns
    # that are not exported in the final table but may be useful during debugging.
    numeric_source_cols = [
        score_col,
        pvalue_col,
        fdr_col,
        get_column(df, ["mean_SHAP", "featurePeak_mean_shap", "peak_mean_shap"]),
        get_column(df, ["mean_abs_SHAP", "featurePeak_mean_abs_shap", "peak_mean_abs_shap"]),
    ]
    numeric_source_cols = [col for col in numeric_source_cols if col is not None]
    df = to_numeric_if_present(df, numeric_source_cols)

    if peak_col is not None:
        peak_id = df[peak_col].map(normalize_peak_id)
    else:
        peak_id = pd.Series(np.nan, index=df.index)

    # Fill missing peak IDs from genomic coordinate columns.
    missing_peak = peak_id.isna() | peak_id.astype("string").str.strip().isin(["", "nan", "None"])
    if missing_peak.any():
        coord_ids = df.apply(lambda row: coord_peak_id(row, seq_col, start_col, end_col), axis=1)
        peak_id = peak_id.where(~missing_peak, coord_ids)

    out = pd.DataFrame(
        {
            "cell_type": spec["cell_type"],
            "comparison": spec["comparison"],
            "peak_id": clean_text_series(peak_id),
            "linked_gene": clean_text_series(df[gene_col]) if gene_col is not None else pd.NA,
            "link_score": df[score_col] if score_col is not None else np.nan,
            "link_p_value": df[pvalue_col] if pvalue_col is not None else np.nan,
            "link_FDR": df[fdr_col] if fdr_col is not None else np.nan,
        }
    )

    out = to_numeric_if_present(out, ["link_score", "link_p_value", "link_FDR"])

    before_drop = len(out)
    out = out.dropna(subset=["peak_id", "linked_gene"]).copy()
    print(f"        dropped missing peak/gene rows: {before_drop - len(out):,}")

    # Enforce the global FDR threshold for inclusion in Table 8.
    before_fdr = len(out)
    out = out[out["link_FDR"].notna() & (out["link_FDR"] < FDR_THRESHOLD)].copy()
    print(f"        dropped by link_FDR < {FDR_THRESHOLD}: {before_fdr - len(out):,}")

    # Deduplicate within this source. Global deduplication is also applied after
    # combining all files, but doing it here makes per-file summaries match output.
    out["_abs_link_score"] = out["link_score"].abs()
    out = out.sort_values(
        by=["cell_type", "comparison", "peak_id", "linked_gene", "link_FDR", "_abs_link_score"],
        ascending=[True, True, True, True, True, False],
        na_position="last",
    )
    out = out.drop_duplicates(
        subset=["cell_type", "comparison", "peak_id", "linked_gene"],
        keep="first",
    )
    out = out.drop(columns=["_abs_link_score"])

    summary = {
        "cell_type": spec["cell_type"],
        "comparison": spec["comparison"],
        "source_file": source_file,
        "file_exists": True,
        "n_rows_raw": n_rows_raw,
        "n_rows_output": len(out),
        "n_unique_peaks": out["peak_id"].nunique(dropna=True),
        "n_unique_genes": out["linked_gene"].nunique(dropna=True),
        "min_link_FDR": out["link_FDR"].min(skipna=True),
        "max_abs_link_score": out["link_score"].abs().max(skipna=True),
    }

    print(f"[DONE]  output rows: {len(out):,}")
    return out[OUTPUT_COLUMNS], summary


def sort_and_deduplicate(all_links):
    """Apply final duplicate handling and requested table ordering."""
    if all_links.empty:
        return pd.DataFrame(columns=OUTPUT_COLUMNS)

    all_links = all_links.copy()
    all_links["_abs_link_score"] = all_links["link_score"].abs()
    all_links = all_links.sort_values(
        by=["cell_type", "comparison", "link_FDR", "_abs_link_score"],
        ascending=[True, True, True, False],
        na_position="last",
    )
    all_links = all_links.drop_duplicates(
        subset=["cell_type", "comparison", "peak_id", "linked_gene"],
        keep="first",
    )
    all_links = all_links.sort_values(
        by=["cell_type", "comparison", "link_FDR", "_abs_link_score"],
        ascending=[True, True, True, False],
        na_position="last",
    )
    all_links = all_links.drop(columns=["_abs_link_score"])
    return all_links[OUTPUT_COLUMNS]


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

    print("\n[FINISHED] Table 8 peak-to-gene links")
    print(f"  CSV:     {OUT_CSV}")
    print(f"  XLSX:    {OUT_XLSX}")
    print(f"  Summary: {OUT_SUMMARY}")
    print(f"  Total output rows: {len(combined):,}")


if __name__ == "__main__":
    main()
