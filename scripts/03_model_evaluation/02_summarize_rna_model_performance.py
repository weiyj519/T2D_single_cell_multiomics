#!/usr/bin/env python3

"""
Script name: 02_summarize_rna_model_performance.py
Purpose: Cleaned public workflow script for the T2D single-cell multiome analysis.
Input: Paths and analysis settings are read from configs/config.yaml.
Output: Module-specific outputs are written under the configured results directory.
Main steps: Load configuration, run the original analysis logic, and export reproducible outputs.
Example command: python scripts/03_model_evaluation/02_summarize_rna_model_performance.py
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
Generate Supplementary Table 4. Machine-learning model performance.

This script summarizes pairwise XGBoost/BorutaShap model performance from saved pkl files.

Compatible with:
1. original pkl:
   results__Acinar_cell.pkl

2. posthoc >=0.6 pkl:
   results__Acinar_cell__test_geq0p6.pkl

Output:
- Supplementary_Table_4_machine_learning_model_performance.csv
- Supplementary_Table_4_machine_learning_model_performance.xlsx
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
# 1. 
# =========================================================

GLOBAL_SEED = 42

DATASETS = [
    {
        "modality": "RNA",
        "parquet_path": str(RNA_EXPRESSION_PARQUET),
        "pkl_glob": "results/stage_features/stage_genes/results__*__test_geq0p6.pkl",
    },
]

OUT_DIR = "results/supplementary_tables"
OUT_PREFIX = "Supplementary_Table_4_machine_learning_model_performance"


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
# 2. 
# =========================================================

def normalize_celltype_name(x: str) -> str:
    """Make cell type name comparable with pkl filename stem."""
    return str(x).replace(" ", "_").replace("/", "_")


def get_celltype_stem_from_pkl(pkl_path: str) -> str:
    """
    Convert:
    results__Acinar_cell.pkl -> Acinar_cell
    results__Acinar_cell__test_geq0p6.pkl -> Acinar_cell
    """
    name = Path(pkl_path).name
    name = name.replace("results__", "")
    name = name.replace(".pkl", "")
    name = re.sub(r"__test.*$", "", name)
    return name


def ensure_condition_and_sample_id(df: pd.DataFrame) -> pd.DataFrame:
    """Parse sample_id and condition from orig.ident if not already present."""
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


def load_celltype_data(parquet_path: str, target_celltype: str) -> pd.DataFrame:
    """Read one cell type from parquet."""
    filters = [("celltype", "=", target_celltype)]

    try:
        table = pq.read_table(parquet_path, filters=filters)
        df = table.to_pandas()
    except Exception as e:
        print(f"[WARN] pyarrow filter failed for {target_celltype}: {e}")
        print("[WARN] Falling back to full parquet read.")
        df = pd.read_parquet(parquet_path)
        df = df[df["celltype"] == target_celltype].copy()

    if df.empty:
        return df

    df = ensure_condition_and_sample_id(df)
    return df


def build_stem_to_celltype_map(parquet_path: str) -> dict:
    """Build mapping from normalized cell type stem to original cell type name."""
    df_ct = pd.read_parquet(parquet_path, columns=["celltype"])
    celltypes = sorted(df_ct["celltype"].dropna().unique().tolist())
    return {normalize_celltype_name(ct): ct for ct in celltypes}


def recompute_input_feature_number(
    data_cell: pd.DataFrame,
    negative_class: str,
    positive_class: str,
) -> dict:
    """
    Recompute:
    - total cells in each group
    - raw feature number
    - feature number after >5% detection filtering in the training set

    This follows the logic in your original model script:
    train_test_split(test_size=0.3, stratify=y, random_state=42)
    features detected in >= floor(0.05*n_train)+1 cells are retained.
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

    feature_cols = [
        c for c in data_pair.columns
        if c not in META_COLS
    ]

    out["n_features_raw"] = int(len(feature_cols))

    if len(feature_cols) == 0:
        return out

    X = data_pair[feature_cols].astype(np.float32)

    try:
        X_train, _, y_train, _ = train_test_split(
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
    min_expr_pct = 0.05
    min_cells = max(5, int(np.floor(n_train * min_expr_pct)) + 1)

    keep_mask = (X_train > 0).sum(axis=0) >= min_cells
    out["n_features_input"] = int(keep_mask.sum())

    return out


def get_result_field(res: dict, posthoc_key: str, original_key: str, default=np.nan):
    """Prefer posthoc field if available, otherwise original field."""
    if isinstance(res, dict):
        if posthoc_key in res:
            return res[posthoc_key]
        if original_key in res:
            return res[original_key]
    return default


def get_classification_report(res: dict):
    return get_result_field(
        res,
        posthoc_key="test_classification_report_geq06",
        original_key="classification_report",
        default={},
    )


def get_confusion_matrix(res: dict):
    return get_result_field(
        res,
        posthoc_key="test_confusion_matrix_geq06",
        original_key="confusion_matrix",
        default=None,
    )


def get_auc(res: dict):
    return get_result_field(
        res,
        posthoc_key="test_auc_geq06",
        original_key="roc_auc",
        default=np.nan,
    )


def get_selected_features(res: dict):
    feats = get_result_field(
        res,
        posthoc_key="selected_features_geq06",
        original_key="selected_features",
        default=[],
    )
    if feats is None:
        return []
    return list(feats)


def extract_report_metrics(report: dict, negative_class: str, positive_class: str) -> dict:
    """
    Extract positive-class metrics.

    In original pkl:
        report keys are usually class names, e.g. "CON", "PRE".

    In posthoc pkl:
        report keys may be "0", "1".
        Here, "1" corresponds to the positive class, i.e. the latter condition.
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
        out["macro_F1"] = float(report["macro avg"].get("f1-score", np.nan))
    if "weighted avg" in report:
        out["weighted_F1"] = float(report["weighted avg"].get("f1-score", np.nan))

    return out


def extract_sensitivity_specificity(cm) -> dict:
    """
    confusion_matrix default order is [0, 1].
    For binary classification:
        TN FP
        FN TP

    sensitivity = TP / (TP + FN), positive-class recall
    specificity = TN / (TN + FP), negative-class recall
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
# 3. 
# =========================================================

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    all_rows = []

    for ds in DATASETS:
        modality = ds["modality"]
        parquet_path = ds["parquet_path"]
        pkl_glob = ds["pkl_glob"]

        print(f"\n========== Processing modality: {modality} ==========")
        print(f"Parquet: {parquet_path}")
        print(f"PKL glob: {pkl_glob}")

        if not os.path.exists(parquet_path):
            print(f"[SKIP] Parquet file not found: {parquet_path}")
            continue

        pkl_files = sorted(glob.glob(pkl_glob))
        if len(pkl_files) == 0:
            print(f"[SKIP] No pkl files found: {pkl_glob}")
            continue

        stem_to_celltype = build_stem_to_celltype_map(parquet_path)

        for pkl_path in pkl_files:
            stem = get_celltype_stem_from_pkl(pkl_path)

            if stem not in stem_to_celltype:
                print(f"[WARN] Cannot map pkl stem to celltype: {stem}")
                print(f"       Available examples: {list(stem_to_celltype.keys())[:5]}")
                continue

            cell_type = stem_to_celltype[stem]
            print(f"\n[Cell type] {cell_type}")
            print(f"PKL: {pkl_path}")

            try:
                with open(pkl_path, "rb") as f:
                    results = pickle.load(f)
            except Exception as e:
                print(f"[WARN] Failed to load pkl: {pkl_path}; {e}")
                continue

            data_cell = load_celltype_data(parquet_path, cell_type)

            if data_cell.empty:
                print(f"[WARN] Empty data for cell type: {cell_type}")
                continue

            for negative_class, positive_class in PAIR_ORDER:
                key = (negative_class, positive_class)

                if key not in results:
                    print(f"  [SKIP] Missing pair in pkl: {key}")
                    continue

                res = results[key]

                count_info = recompute_input_feature_number(
                    data_cell=data_cell,
                    negative_class=negative_class,
                    positive_class=positive_class,
                )

                auc = get_auc(res)
                selected_features = get_selected_features(res)
                report = get_classification_report(res)
                cm = get_confusion_matrix(res)

                report_metrics = extract_report_metrics(
                    report=report,
                    negative_class=negative_class,
                    positive_class=positive_class,
                )

                cm_metrics = extract_sensitivity_specificity(cm)

                row = {
                    "modality": modality,
                    "cell_type": cell_type,
                    "comparison": comparison_label(negative_class, positive_class),
                    "group1": negative_class,
                    "group2": positive_class,
                    "positive_class": positive_class,

                    # full dataset cells in this pair
                    "n_cells_group1": count_info["n_cells_group1"],
                    "n_cells_group2": count_info["n_cells_group2"],

                    # optional but useful
                    "test_n_group1": report_metrics["test_n_group1"],
                    "test_n_group2": report_metrics["test_n_group2"],

                    # feature numbers
                    "n_features_raw": count_info["n_features_raw"],
                    "n_features_input": count_info["n_features_input"],
                    "n_features_selected": int(len(selected_features)),

                    # performance
                    "AUC": float(auc) if pd.notna(auc) else np.nan,
                    "accuracy": report_metrics["accuracy"],
                    "precision": report_metrics["precision"],
                    "recall": report_metrics["recall"],
                    "F1_score": report_metrics["F1_score"],
                    "sensitivity": cm_metrics["sensitivity"],
                    "specificity": cm_metrics["specificity"],

                    # optional diagnostic columns
                    "macro_F1": report_metrics["macro_F1"],
                    "weighted_F1": report_metrics["weighted_F1"],
                    "TN": cm_metrics["TN"],
                    "FP": cm_metrics["FP"],
                    "FN": cm_metrics["FN"],
                    "TP": cm_metrics["TP"],
                    "pkl_path": pkl_path,
                }

                all_rows.append(row)

                print(
                    f"  {negative_class} vs {positive_class}: "
                    f"AUC={row['AUC']:.4f}, "
                    f"n_input={row['n_features_input']}, "
                    f"n_selected={row['n_features_selected']}"
                )

            del data_cell
            gc.collect()

    if len(all_rows) == 0:
        raise RuntimeError("No rows were generated. Please check paths and pkl structure.")

    df = pd.DataFrame(all_rows)

    # Sort table
    modality_order = {"RNA": 0, "ATAC": 1}
    pair_order = {
        "CON vs PRE": 0,
        "PRE vs T2D": 1,
        "CON vs T2D": 2,
    }

    df["_modality_order"] = df["modality"].map(modality_order).fillna(99)
    df["_pair_order"] = df["comparison"].map(pair_order).fillna(99)

    df = df.sort_values(
        by=["_modality_order", "cell_type", "_pair_order"]
    ).drop(columns=["_modality_order", "_pair_order"])

    # Round numeric metrics
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

    # Reorder main columns
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
        "pkl_path",
    ]

    existing_cols = [c for c in ordered_cols if c in df.columns]
    remaining_cols = [c for c in df.columns if c not in existing_cols]
    df = df[existing_cols + remaining_cols]

    csv_path = os.path.join(OUT_DIR, f"{OUT_PREFIX}.csv")
    xlsx_path = os.path.join(OUT_DIR, f"{OUT_PREFIX}.xlsx")

    df.to_csv(csv_path, index=False)

    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
        df.to_excel(writer, sheet_name="model_performance", index=False)

    print("\n========== Done ==========")
    print(f"Saved CSV : {csv_path}")
    print(f"Saved XLSX: {xlsx_path}")
    print(f"Total rows: {df.shape[0]}")

    expected_rows = len(DATASETS) * 8 * 3
    print(f"Expected rows if all 8 cell types × 3 comparisons × modalities are present: {expected_rows}")


if __name__ == "__main__":
    main()