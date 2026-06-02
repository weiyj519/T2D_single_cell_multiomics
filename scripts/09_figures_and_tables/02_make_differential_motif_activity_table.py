#!/usr/bin/env python3

"""
Script name: 02_make_differential_motif_activity_table.py
Purpose: Cleaned public workflow script for the T2D single-cell multiome analysis.
Input: Paths and analysis settings are read from configs/config.yaml.
Output: Module-specific outputs are written under the configured results directory.
Main steps: Load configuration, run the original analysis logic, and export reproducible outputs.
Example command: python scripts/09_figures_and_tables/02_make_differential_motif_activity_table.py
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

# -*- coding: utf-8 -*-

"""
Generate supplementary table:
Differential motif activity based on STAGE peaks.

Input chromVAR motif activity files are searched under:
    results/downstream

Expected input columns:
    set,motif,n_peaks_in_set,n_peaks_with_motif,frac_peaks_with_motif,
    TF,chromvar_lfc,chromvar_padj,chromvar_log10padj,label

Outputs:
    results/model_outputs/atac/supplementary_tables/
        Differential_motif_activity_based_on_STAGE_peaks.csv
        Differential_motif_activity_based_on_STAGE_peaks.xlsx
        Differential_motif_activity_based_on_STAGE_peaks_summary.csv
"""

import glob
import os
import re

import pandas as pd


# =========================================================
# User settings
# =========================================================

INPUT_ROOT = "results/downstream"
OUT_DIR = "results/model_outputs/atac/supplementary_tables"
OUT_PREFIX = "Differential_motif_activity_based_on_STAGE_peaks"

MOTIF_FILENAME = "motifs_hitFilter_sig_chromvar.csv"
MAPPING_DIR = "mlpeaks_chromvar_mapping"

CELL_TYPES = {
    "beta": r"(beta|β)",
    "alpha": r"(alpha|α)",
    "acinar": r"acinar",
}

COMPARISONS = {
    "CON vs PRE": ("CON", "PRE"),
    "PRE vs T2D": ("PRE", "T2D"),
}

OUTPUT_COLUMNS = [
    "cell_type",
    "comparison",
    "motif_id",
    "motif_name",
    "n_STAGE_peaks",
    "n_peaks_with_motif",
    "fraction_peaks_with_motif",
    "chromVAR_log2FC",
    "adjusted_p_value",
    "direction",
]

REQUIRED_COLUMNS = [
    "motif",
    "n_peaks_in_set",
    "n_peaks_with_motif",
    "frac_peaks_with_motif",
    "chromvar_lfc",
    "chromvar_padj",
]

NUMERIC_COLUMNS = [
    "n_STAGE_peaks",
    "n_peaks_with_motif",
    "fraction_peaks_with_motif",
    "chromVAR_log2FC",
    "adjusted_p_value",
]


# =========================================================
# Helper functions
# =========================================================

def print_step(message: str) -> None:
    """Print progress with immediate flushing for terminal runs."""
    print(message, flush=True)


def comparison_regex(first_group: str, second_group: str) -> str:
    """Match comparison names such as CON_PRE, CON-vs-PRE, or con.pre."""
    return rf"{re.escape(first_group)}.*{re.escape(second_group)}"


def get_candidate_files() -> list[str]:
    """
    Find all chromVAR motif CSVs below INPUT_ROOT with glob.

    The later regex filtering is case-insensitive, which makes the search more
    robust than relying on a case-sensitive Unix glob for beta/alpha/acinar.
    """
    pattern = os.path.join(INPUT_ROOT, "*", MAPPING_DIR, MOTIF_FILENAME)
    files = sorted(glob.glob(pattern))
    return files


def find_files_for_combo(
    all_files: list[str],
    cell_type: str,
    cell_pattern: str,
    comparison: str,
    first_group: str,
    second_group: str,
) -> list[str]:
    """Return files whose run directory matches the requested cell/comparison."""
    combo_pattern = comparison_regex(first_group, second_group)
    matches = []

    for path in all_files:
        run_dir = os.path.basename(os.path.dirname(os.path.dirname(path)))
        has_cell = re.search(cell_pattern, run_dir, flags=re.IGNORECASE) is not None
        has_comparison = re.search(combo_pattern, run_dir, flags=re.IGNORECASE) is not None
        if has_cell and has_comparison:
            matches.append(path)

    if not matches:
        example_pattern = os.path.join(
            INPUT_ROOT,
            f"*{cell_type}*{first_group}*{second_group}*",
            MAPPING_DIR,
            MOTIF_FILENAME,
        )
        print_step(
            f"[WARN] No file found for {cell_type} | {comparison}. "
            f"Example glob: {example_pattern}"
        )
    elif len(matches) > 1:
        print_step(
            f"[WARN] Found {len(matches)} files for {cell_type} | {comparison}; "
            "all matched files will be included."
        )
        for path in matches:
            print_step(f"       - {path}")

    return matches


def require_columns(df: pd.DataFrame, path: str) -> None:
    """Validate required input columns before transforming."""
    missing = [col for col in REQUIRED_COLUMNS if col not in df.columns]
    if "TF" not in df.columns and "label" not in df.columns:
        missing.append("TF or label")

    if missing:
        raise ValueError(f"{path} is missing required column(s): {', '.join(missing)}")


def choose_motif_name(df: pd.DataFrame) -> pd.Series:
    """Use TF first; if TF is missing/blank, fall back to label."""
    if "TF" in df.columns:
        motif_name = df["TF"].astype("string")
    else:
        motif_name = pd.Series(pd.NA, index=df.index, dtype="string")

    if "label" in df.columns:
        label = df["label"].astype("string")
        motif_name = motif_name.mask(motif_name.isna() | (motif_name.str.strip() == ""), label)

    return motif_name


def direction_from_lfc(lfc: float, comparison: str) -> str:
    """
    Assign activity direction.

    chromvar_lfc is assumed to be the latter group relative to the former group.
    """
    if pd.isna(lfc):
        return "unknown"

    if lfc > 0:
        if comparison == "CON vs PRE":
            return "PRE_up"
        if comparison == "PRE vs T2D":
            return "T2D_up"
    elif lfc < 0:
        if comparison == "CON vs PRE":
            return "CON_up"
        if comparison == "PRE vs T2D":
            return "PRE_up"

    return "no_change"


def load_and_transform_file(path: str, cell_type: str, comparison: str) -> pd.DataFrame:
    """Read one motif CSV and return standardized output columns."""
    print_step(f"[INFO] Reading: {path}")
    df = pd.read_csv(path)
    require_columns(df, path)

    out = pd.DataFrame(
        {
            "cell_type": cell_type,
            "comparison": comparison,
            "motif_id": df["motif"],
            "motif_name": choose_motif_name(df),
            "n_STAGE_peaks": df["n_peaks_in_set"],
            "n_peaks_with_motif": df["n_peaks_with_motif"],
            "fraction_peaks_with_motif": df["frac_peaks_with_motif"],
            "chromVAR_log2FC": df["chromvar_lfc"],
            "adjusted_p_value": df["chromvar_padj"],
        }
    )

    for col in NUMERIC_COLUMNS:
        out[col] = pd.to_numeric(out[col], errors="coerce")

    out["direction"] = out["chromVAR_log2FC"].apply(
        lambda value: direction_from_lfc(value, comparison)
    )

    return out[OUTPUT_COLUMNS]


def build_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Create per-cell/per-comparison summary statistics."""
    summary = (
        df.groupby(["cell_type", "comparison"], as_index=False)
        .agg(
            n_rows=("motif_id", "size"),
            n_unique_motifs=("motif_id", pd.Series.nunique),
            n_STAGE_peaks=("n_STAGE_peaks", "max"),
            min_adjusted_p_value=("adjusted_p_value", "min"),
            max_abs_chromVAR_log2FC=("chromVAR_log2FC", lambda x: x.abs().max()),
        )
        .sort_values(["cell_type", "comparison"], ascending=[True, True])
        .reset_index(drop=True)
    )
    return summary


def write_outputs(df: pd.DataFrame, summary: pd.DataFrame) -> None:
    """Write CSV, XLSX, and summary CSV."""
    os.makedirs(OUT_DIR, exist_ok=True)

    out_csv = os.path.join(OUT_DIR, f"{OUT_PREFIX}.csv")
    out_xlsx = os.path.join(OUT_DIR, f"{OUT_PREFIX}.xlsx")
    summary_csv = os.path.join(OUT_DIR, f"{OUT_PREFIX}_summary.csv")

    print_step(f"[INFO] Writing CSV: {out_csv}")
    df.to_csv(out_csv, index=False)

    print_step(f"[INFO] Writing XLSX: {out_xlsx}")
    df.to_excel(out_xlsx, index=False)

    print_step(f"[INFO] Writing summary CSV: {summary_csv}")
    summary.to_csv(summary_csv, index=False)


# =========================================================
# Main
# =========================================================

def main() -> None:
    print_step("[INFO] Starting supplementary motif activity table generation.")
    print_step(f"[INFO] Input root: {INPUT_ROOT}")
    print_step(f"[INFO] Output directory: {OUT_DIR}")

    all_files = get_candidate_files()
    print_step(f"[INFO] Found {len(all_files)} chromVAR motif CSV candidate file(s).")

    tables = []
    for cell_type, cell_pattern in CELL_TYPES.items():
        for comparison, (first_group, second_group) in COMPARISONS.items():
            print_step(f"[INFO] Searching for {cell_type} | {comparison}")
            matched_files = find_files_for_combo(
                all_files=all_files,
                cell_type=cell_type,
                cell_pattern=cell_pattern,
                comparison=comparison,
                first_group=first_group,
                second_group=second_group,
            )

            for path in matched_files:
                try:
                    table = load_and_transform_file(path, cell_type, comparison)
                except Exception as exc:
                    print_step(f"[WARN] Skipping file due to error: {path}")
                    print_step(f"       {type(exc).__name__}: {exc}")
                    continue

                print_step(f"[INFO] Loaded {len(table)} row(s) for {cell_type} | {comparison}")
                tables.append(table)

    if not tables:
        raise RuntimeError("No valid input files were loaded. No outputs were written.")

    final_df = pd.concat(tables, ignore_index=True)
    final_df = (
        final_df.sort_values(
            ["cell_type", "comparison", "adjusted_p_value"],
            ascending=[True, True, True],
            na_position="last",
        )
        .reset_index(drop=True)
    )

    summary = build_summary(final_df)

    print_step(f"[INFO] Final table rows: {len(final_df)}")
    print_step(f"[INFO] Summary rows: {len(summary)}")

    write_outputs(final_df, summary)
    print_step("[INFO] Done.")


if __name__ == "__main__":
    main()
