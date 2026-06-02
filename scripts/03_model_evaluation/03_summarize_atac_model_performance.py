#!/usr/bin/env python3

"""
Script name: 03_summarize_atac_model_performance.py
Purpose: Cleaned public workflow script for the T2D single-cell multiome analysis.
Input: Paths and analysis settings are read from configs/config.yaml.
Output: Module-specific outputs are written under the configured results directory.
Main steps: Load configuration, run the original analysis logic, and export reproducible outputs.
Example command: python scripts/03_model_evaluation/03_summarize_atac_model_performance.py
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
Generate Supplementary Table 4 for ATAC machine-learning model performance.

Input:
    results/model_outputs/atac/results__*.pkl
    data/processed/atac_peaks_autosomes.parquet

Output:
    results/supplementary_tables/
        Supplementary_Table_4_ATAC_machine_learning_model_performance.csv
        Supplementary_Table_4_ATAC_machine_learning_model_performance.xlsx

Column definitions:
    n_features_raw:
        Total number of peak features before training-set detection filtering.

    n_features_input:
        Number of peak features retained after the training-set detection filter.

    n_features_selected:
        Number of features actually used in the final model.
        This comes from the selected_features field in the pkl.
        If no feature reached fold_freq >= 0.6 during training, the model may have used
        fallback top-ranked features, and these are still counted here.

    n_features_stable:
        Number of stable BorutaShap features with fold_freq >= 0.6.
        These are the strict stable STAGE peaks.

Important:
    Positive class is the latter condition:
        CON vs PRE  -> PRE
        PRE vs T2D  -> T2D
        CON vs T2D  -> T2D
"""

import os
import re
import gc
import glob
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from sklearn.model_selection import train_test_split


# =========================================================
# 1. User settings
# =========================================================

GLOBAL_SEED = 42

MODALITY = "ATAC"

PARQUET_PATH = str(ATAC_PEAK_PARQUET)

PKL_GLOB = "results/model_outputs/atac/results__*.pkl"

OUT_DIR = "results/supplementary_tables"

OUT_PREFIX = "Supplementary_Table_4_ATAC_machine_learning_model_performance"

FREQ_THRESHOLD = 0.6

META_COLS = {
    "cell_id",
    "orig.ident",
    "seurat_clusters",
    "celltype",
    "group",
    "sample_id",
    "condition",
}

PAIR_ORDER = [
    ("CON", "PRE"),
    ("PRE", "T2D"),
    ("CON", "T2D"),
]


# =========================================================
# 2. Helper functions
# =========================================================

def normalize_celltype_name(x: str) -> str:
    """
    Convert cell type name to pkl-compatible stem.

    Examples:
        Acinar cell -> Acinar_cell
        β cell      -> β_cell
        δ cell      -> δ_cell
    """
    return str(x).replace(" ", "_").replace("/", "_")


def get_celltype_stem_from_pkl(pkl_path: str) -> str:
    """
    Convert pkl filename to cell type stem.

    Examples:
        results__Acinar_cell.pkl -> Acinar_cell
        results__β_cell.pkl      -> β_cell
        results__δ_cell.pkl      -> δ_cell
    """
    name = Path(pkl_path).name
    name = name.replace("results__", "")
    name = name.replace(".pkl", "")
    name = re.sub(r"__test.*$", "", name)
    return name


def ensure_condition_and_sample_id(df: pd.DataFrame) -> pd.DataFrame:
    """
    Parse sample_id and condition from orig.ident if needed.

    Expected orig.ident pattern:
        A1CON, B2PRE, C3T2D, etc.
    """
    df = df.copy()

    if "condition" not in df.columns or df["condition"].isna().all():
        if "orig.ident" not in df.columns:
            raise ValueError("No condition column and no orig.ident column found.")

        extracted = df["orig.ident"].astype(str).str.extract(r"([A-Z]\d+)(CON|PRE|T2D)")
        df["sample_id"] = extracted[0]
        df["condition"] = extracted[1]

    if "sample_id" not in df.columns or df["sample_id"].isna().all():
        if "orig.ident" in df.columns:
            extracted = df["orig.ident"].astype(str).str.extract(r"([A-Z]\d+)(CON|PRE|T2D)")
            df["sample_id"] = extracted[0]

    return df


def build_stem_to_celltype_map(parquet_path: str) -> dict:
    """
    Build mapping from normalized pkl stem to original cell type name.

    Example:
        Acinar_cell -> Acinar cell
        β_cell      -> β cell
    """
    df_ct = pd.read_parquet(parquet_path, columns=["celltype"])
    celltypes = sorted(df_ct["celltype"].dropna().unique().tolist())
    return {normalize_celltype_name(ct): ct for ct in celltypes}


def load_celltype_data(parquet_path: str, target_celltype: str) -> pd.DataFrame:
    """
    Read one cell type from parquet.
    pyarrow filtering is used first to reduce memory usage.
    """
    filters = [("celltype", "=", target_celltype)]

    try:
        table = pq.read_table(parquet_path, filters=filters)
        df = table.to_pandas()
    except Exception as e:
        print(f"[WARN] Filtered reading failed for {target_celltype}: {e}")
        print("[WARN] Falling back to full parquet reading.")
        df = pd.read_parquet(parquet_path)
        df = df[df["celltype"] == target_celltype].copy()

    if df.empty:
        return df

    df = ensure_condition_and_sample_id(df)
    return df


def recompute_input_feature_number(
    data_cell: pd.DataFrame,
    negative_class: str,
    positive_class: str,
) -> dict:
    """
    Recalculate:
        n_cells_group1
        n_cells_group2
        n_features_raw
        n_features_input

    This follows the original ATAC model logic:
        1. Pairwise subset
        2. 70/30 stratified train/test split
        3. Training-set detection filter:
           features detected in >= floor(0.05 * n_train) + 1 cells

    Here, n_features_input means the number of peak features entering BorutaShap
    after training-set detection filtering.
    """
    data_pair = data_cell[data_cell["condition"].isin([negative_class, positive_class])].copy()

    out = {
        "n_cells_group1": int((data_pair["condition"] == negative_class).sum()),
        "n_cells_group2": int((data_pair["condition"] == positive_class).sum()),
        "n_features_raw": np.nan,
        "n_features_input": np.nan,
    }

    if data_pair.empty or data_pair["condition"].nunique() < 2:
        return out

    y = np.where(data_pair["condition"] == positive_class, 1, 0)

    feature_cols = [c for c in data_pair.columns if c not in META_COLS]

    out["n_features_raw"] = int(len(feature_cols))

    if len(feature_cols) == 0:
        return out

    X = data_pair[feature_cols].astype(np.float32)

    try:
        X_train, _, _, _ = train_test_split(
            X,
            y,
            test_size=0.3,
            stratify=y,
            random_state=GLOBAL_SEED,
        )
    except ValueError as e:
        print(f"[WARN] train_test_split failed for {negative_class} vs {positive_class}: {e}")
        return out

    n_train = X_train.shape[0]
    min_cells = max(5, int(np.floor(n_train * 0.05)) + 1)

    keep_mask = (X_train > 0).sum(axis=0) >= min_cells
    out["n_features_input"] = int(keep_mask.sum())

    return out


def get_feature_selection_info(res: dict, freq_threshold: float = 0.6) -> dict:
    """
    Return both:
        n_features_selected:
            Actual number of features used in the final model.
            This is extracted from selected_features in the pkl.

        n_features_stable:
            Number of features satisfying fold_freq >= freq_threshold.
            These are strict stable BorutaShap features.

    This solves the situation where:
        n_features_stable = 0
    but:
        n_features_selected > 0

    In that case, the model used fallback top-ranked features, but no strict stable
    feature should be counted as a stable STAGE peak.
    """
    stability_df = res.get("boruta_stability", None)

    stable_features = []

    if isinstance(stability_df, pd.DataFrame):
        if "gene" in stability_df.columns and "fold_freq" in stability_df.columns:
            stable_features = stability_df.loc[
                stability_df["fold_freq"] >= freq_threshold,
                "gene"
            ].tolist()

    model_features = res.get("selected_features", None)

    # Compatibility with possible posthoc pkl naming
    if model_features is None:
        model_features = res.get("selected_features_geq06", [])

    if model_features is None:
        model_features = []

    model_features = list(model_features)

    return {
        "stable_features": stable_features,
        "model_features": model_features,
        "n_features_stable": int(len(stable_features)),
        "n_features_selected": int(len(model_features)),
    }


def get_auc(res: dict):
    """
    Extract test AUC.
    """
    if "roc_auc" in res:
        return res["roc_auc"]

    if "test_auc_geq06" in res:
        return res["test_auc_geq06"]

    return np.nan


def get_classification_report(res: dict):
    """
    Extract classification report.
    """
    if "classification_report" in res:
        return res["classification_report"]

    if "test_classification_report_geq06" in res:
        return res["test_classification_report_geq06"]

    return {}


def get_confusion_matrix(res: dict):
    """
    Extract confusion matrix.
    """
    if "confusion_matrix" in res:
        return res["confusion_matrix"]

    if "test_confusion_matrix_geq06" in res:
        return res["test_confusion_matrix_geq06"]

    return None


def extract_report_metrics(report: dict, negative_class: str, positive_class: str) -> dict:
    """
    Extract positive-class metrics from sklearn classification_report output_dict.

    In the original pkl, report keys are usually condition names:
        CON, PRE, T2D

    If the report uses numeric labels:
        0, 1
    then 1 is treated as the positive class.
    """
    if not isinstance(report, dict):
        report = {}

    positive_key = None
    negative_key = None

    if positive_class in report:
        positive_key = positive_class
    elif "1" in report:
        positive_key = "1"
    elif 1 in report:
        positive_key = 1

    if negative_class in report:
        negative_key = negative_class
    elif "0" in report:
        negative_key = "0"
    elif 0 in report:
        negative_key = 0

    def get_metric(class_key, metric_name):
        if class_key is None:
            return np.nan
        try:
            return float(report[class_key].get(metric_name, np.nan))
        except Exception:
            return np.nan

    try:
        accuracy = float(report.get("accuracy", np.nan))
    except Exception:
        accuracy = np.nan

    out = {
        "accuracy": accuracy,
        "precision": get_metric(positive_key, "precision"),
        "recall": get_metric(positive_key, "recall"),
        "F1_score": get_metric(positive_key, "f1-score"),
        "test_n_group1": get_metric(negative_key, "support"),
        "test_n_group2": get_metric(positive_key, "support"),
        "macro_F1": np.nan,
        "weighted_F1": np.nan,
    }

    if "macro avg" in report:
        try:
            out["macro_F1"] = float(report["macro avg"].get("f1-score", np.nan))
        except Exception:
            pass

    if "weighted avg" in report:
        try:
            out["weighted_F1"] = float(report["weighted avg"].get("f1-score", np.nan))
        except Exception:
            pass

    return out


def extract_sensitivity_specificity(cm) -> dict:
    """
    Extract sensitivity and specificity from binary confusion matrix.

    sklearn confusion_matrix with labels [0, 1]:

        TN  FP
        FN  TP

    sensitivity = TP / (TP + FN)
    specificity = TN / (TN + FP)
    """
    out = {
        "sensitivity": np.nan,
        "specificity": np.nan,
        "TN": np.nan,
        "FP": np.nan,
        "FN": np.nan,
        "TP": np.nan,
    }

    if cm is None:
        return out

    cm = np.asarray(cm)

    if cm.shape != (2, 2):
        return out

    tn, fp, fn, tp = cm.ravel()

    out["TN"] = int(tn)
    out["FP"] = int(fp)
    out["FN"] = int(fn)
    out["TP"] = int(tp)

    out["sensitivity"] = float(tp / (tp + fn)) if (tp + fn) > 0 else np.nan
    out["specificity"] = float(tn / (tn + fp)) if (tn + fp) > 0 else np.nan

    return out


def comparison_label(negative_class: str, positive_class: str) -> str:
    return f"{negative_class} vs {positive_class}"


# =========================================================
# 3. Main
# =========================================================

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    if not os.path.exists(PARQUET_PATH):
        raise FileNotFoundError(f"Parquet file not found: {PARQUET_PATH}")

    pkl_files = sorted(glob.glob(PKL_GLOB))

    if len(pkl_files) == 0:
        raise FileNotFoundError(f"No pkl files found: {PKL_GLOB}")

    print(f"[INFO] Found {len(pkl_files)} pkl files.")
    print(f"[INFO] Parquet: {PARQUET_PATH}")
    print(f"[INFO] PKL glob: {PKL_GLOB}")
    print(f"[INFO] Output dir: {OUT_DIR}")

    stem_to_celltype = build_stem_to_celltype_map(PARQUET_PATH)

    all_rows = []

    for pkl_path in pkl_files:
        stem = get_celltype_stem_from_pkl(pkl_path)

        if stem not in stem_to_celltype:
            print(f"[WARN] Cannot map pkl stem to celltype: {stem}")
            print(f"       Available examples: {list(stem_to_celltype.keys())[:10]}")
            continue

        cell_type = stem_to_celltype[stem]

        print(f"\n========== {cell_type} ==========")
        print(f"[INFO] PKL: {pkl_path}")

        try:
            with open(pkl_path, "rb") as f:
                results = pickle.load(f)
        except Exception as e:
            print(f"[WARN] Failed to load pkl: {pkl_path}; {e}")
            continue

        data_cell = load_celltype_data(PARQUET_PATH, cell_type)

        if data_cell.empty:
            print(f"[WARN] Empty data for cell type: {cell_type}")
            continue

        for negative_class, positive_class in PAIR_ORDER:
            key = (negative_class, positive_class)

            if key not in results:
                print(f"  [SKIP] Missing pair: {negative_class} vs {positive_class}")
                continue

            res = results[key]

            count_info = recompute_input_feature_number(
                data_cell=data_cell,
                negative_class=negative_class,
                positive_class=positive_class,
            )

            feature_info = get_feature_selection_info(
                res=res,
                freq_threshold=FREQ_THRESHOLD,
            )

            auc = get_auc(res)
            report = get_classification_report(res)
            cm = get_confusion_matrix(res)

            report_metrics = extract_report_metrics(
                report=report,
                negative_class=negative_class,
                positive_class=positive_class,
            )

            cm_metrics = extract_sensitivity_specificity(cm)

            row = {
                "modality": MODALITY,
                "cell_type": cell_type,
                "comparison": comparison_label(negative_class, positive_class),
                "group1": negative_class,
                "group2": positive_class,
                "positive_class": positive_class,

                # All cells in this comparison
                "n_cells_group1": count_info["n_cells_group1"],
                "n_cells_group2": count_info["n_cells_group2"],

                # Test-set cells from classification_report
                "test_n_group1": report_metrics["test_n_group1"],
                "test_n_group2": report_metrics["test_n_group2"],

                # Feature counts
                "n_features_raw": count_info["n_features_raw"],
                "n_features_input": count_info["n_features_input"],
                "n_features_selected": feature_info["n_features_selected"],
                "n_features_stable": feature_info["n_features_stable"],

                # Performance metrics
                "AUC": float(auc) if pd.notna(auc) else np.nan,
                "accuracy": report_metrics["accuracy"],
                "precision": report_metrics["precision"],
                "recall": report_metrics["recall"],
                "F1_score": report_metrics["F1_score"],
                "sensitivity": cm_metrics["sensitivity"],
                "specificity": cm_metrics["specificity"],

                # Optional diagnostic columns
                "macro_F1": report_metrics["macro_F1"],
                "weighted_F1": report_metrics["weighted_F1"],
                "TN": cm_metrics["TN"],
                "FP": cm_metrics["FP"],
                "FN": cm_metrics["FN"],
                "TP": cm_metrics["TP"],

                # Reproducibility
                "boruta_frequency_rule": f"fold_freq >= {FREQ_THRESHOLD}",
                "pkl_path": pkl_path,
            }

            all_rows.append(row)

            print(
                f"  {negative_class} vs {positive_class}: "
                f"AUC={row['AUC']:.4f}, "
                f"n_input={row['n_features_input']}, "
                f"n_selected={row['n_features_selected']}, "
                f"n_stable={row['n_features_stable']}"
            )

        del data_cell
        gc.collect()

    if len(all_rows) == 0:
        raise RuntimeError("No rows generated. Please check pkl files and celltype names.")

    df = pd.DataFrame(all_rows)

    pair_order = {
        "CON vs PRE": 0,
        "PRE vs T2D": 1,
        "CON vs T2D": 2,
    }

    celltype_order = {
        "β cell": 0,
        "α cell": 1,
        "Acinar cell": 2,
        "δ cell": 3,
        "PP cell": 4,
        "Ductal cell": 5,
        "Stellate cell": 6,
        "Immune cell": 7,
    }

    df["_celltype_order"] = df["cell_type"].map(celltype_order).fillna(99)
    df["_pair_order"] = df["comparison"].map(pair_order).fillna(99)

    df = df.sort_values(
        ["_celltype_order", "cell_type", "_pair_order"]
    ).drop(columns=["_celltype_order", "_pair_order"])

    metric_cols = [
        "AUC",
        "accuracy",
        "precision",
        "recall",
        "F1_score",
        "sensitivity",
        "specificity",
        "macro_F1",
        "weighted_F1",
    ]

    for col in metric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").round(4)

    ordered_cols = [
        "modality",
        "cell_type",
        "comparison",
        "group1",
        "group2",
        "positive_class",

        "n_cells_group1",
        "n_cells_group2",
        "test_n_group1",
        "test_n_group2",

        "n_features_raw",
        "n_features_input",
        "n_features_selected",
        "n_features_stable",

        "AUC",
        "accuracy",
        "precision",
        "recall",
        "F1_score",
        "sensitivity",
        "specificity",
        "macro_F1",
        "weighted_F1",

        "TN",
        "FP",
        "FN",
        "TP",

        "boruta_frequency_rule",
        "pkl_path",
    ]

    existing_cols = [c for c in ordered_cols if c in df.columns]
    remaining_cols = [c for c in df.columns if c not in existing_cols]
    df = df[existing_cols + remaining_cols]

    csv_path = os.path.join(OUT_DIR, f"{OUT_PREFIX}.csv")
    xlsx_path = os.path.join(OUT_DIR, f"{OUT_PREFIX}.xlsx")

    df.to_csv(csv_path, index=False)

    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
        df.to_excel(writer, sheet_name="ATAC_model_performance", index=False)

    print("\n========== Done ==========")
    print(f"Saved CSV : {csv_path}")
    print(f"Saved XLSX: {xlsx_path}")
    print(f"Total rows: {df.shape[0]}")

    expected_rows = len(pkl_files) * 3
    print(f"Expected rows if every pkl contains all three comparisons: {expected_rows}")


if __name__ == "__main__":
    main()