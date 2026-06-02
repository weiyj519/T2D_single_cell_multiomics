#!/usr/bin/env python3

"""
Script name: 01_train_atac_xgboost_borutashap.py
Purpose: Cleaned public workflow script for the T2D single-cell multiome analysis.
Input: Paths and analysis settings are read from configs/config.yaml.
Output: Module-specific outputs are written under the configured results directory.
Main steps: Load configuration, run the original analysis logic, and export reproducible outputs.
Example command: python scripts/02_ml_feature_selection/01_train_atac_xgboost_borutashap.py
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

import os
import gc
import time
import pickle
import warnings
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import shap
import pyarrow.parquet as pq

import xgboost as xgb



from xgboost import XGBClassifier
from sklearn.model_selection import (
    StratifiedGroupKFold,
    StratifiedKFold,
    GroupShuffleSplit,
    train_test_split,
)
from sklearn.metrics import (
    roc_auc_score,
    roc_curve,
    classification_report,
    confusion_matrix,
    f1_score,
)
from sklearn.feature_selection import SelectKBest, f_classif

from BorutaShap import BorutaShap  # pip install BorutaShap

warnings.filterwarnings("ignore")
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

# 
GLOBAL_SEED = 42
rng_global = np.random.RandomState(GLOBAL_SEED)

import numpy as np
from xgboost import XGBClassifier

class XGBClassifierForBoruta(XGBClassifier):
    """
     BorutaShap  importance_measure="gini"  feature_importances_
    ：xgboost  feature_importances_  property， set
     property（ setter）
    """
    def __init__(self, importance_type="gain", **kwargs):
        super().__init__(importance_type=importance_type, **kwargs)
        self._boruta_feature_importances = None

    @property
    def feature_importances_(self):
        # BorutaShap  check_model ；
        if self._boruta_feature_importances is None:
            return np.array([], dtype=float)
        return self._boruta_feature_importances

    @feature_importances_.setter
    def feature_importances_(self, value):
        self._boruta_feature_importances = np.asarray(value, dtype=float)

    def fit(self, X, y, **kwargs):
        super().fit(X, y, **kwargs)

        booster = self.get_booster()

        #  DataFrame 
        if hasattr(X, "columns"):
            feature_names = list(X.columns)
        else:
            feature_names = booster.feature_names
            if feature_names is None:
                feature_names = [f"f{i}" for i in range(self.n_features_in_)]

        score = booster.get_score(importance_type="gain")
        imp = np.array([score.get(fn, 0.0) for fn in feature_names], dtype=float)

        # ：， 0 
        s = imp.sum()
        if s > 0:
            imp = imp / s

        self.feature_importances_ = imp
        return self

# =========================================================
# 
# =========================================================
def get_xgb_params_by_count(
    n_cells: int, scale_pos_weight: float = 1.0, random_state: int = GLOBAL_SEED
) -> dict:
    """
     XGBoost （ early stopping）
    - n_cells >= 2000: Rich Mode
    - n_cells < 2000 : Sparse Mode
    """
    base = dict(
        objective="binary:logistic",
        random_state=random_state,
        n_jobs=20,
        n_estimators=500,      # ，
        scale_pos_weight=scale_pos_weight,
        tree_method="hist",
        eval_metric="auc",
    )

    if n_cells >= 2000:
        # ：
        print(f"  []  'Rich Mode' (n_cells={n_cells} >= 2000)")
        base.update(
            dict(
                learning_rate=0.04,
                max_depth=3,
                min_child_weight=3,
                subsample=0.8,
                colsample_bytree=0.7,
                reg_alpha=3.0,
                reg_lambda=6.0,
            )
        )
    else:
        # ：，
        print(f"  []  'Sparse Mode' (n_cells={n_cells} < 2000)")
        base.update(
            dict(
                learning_rate=0.03,
                max_depth=3,
                min_child_weight=3,
                subsample=0.7,
                colsample_bytree=0.6,
                reg_alpha=2.0,
                reg_lambda=4.0,
            )
        )

    return base


def build_xgb(params_extra=None, scale_pos_weight=1.0, random_state=GLOBAL_SEED):
    """
    （ Rich Mode），
    """
    base = get_xgb_params_by_count(10000, scale_pos_weight, random_state)
    if params_extra:
        base.update(params_extra)
    return XGBClassifier(**base)


def load_celltype_data(parquet_path: str, target_celltype: str) -> pd.DataFrame:
    """
     Parquet ， celltype 
     orig.ident  sample_id  condition
    """
    filters = [("celltype", "=", target_celltype)]
    try:
        table = pq.read_table(parquet_path, filters=filters)
        df = table.to_pandas()
    except Exception as e:
        print(f"，: {e}")
        df = pd.read_parquet(parquet_path)
        df = df[df["celltype"] == target_celltype].copy()

    if df.empty:
        return df

    pattern = r"([A-Z]\d+)(CON|PRE|T2D)"
    extracted = df["orig.ident"].str.extract(pattern)
    df["sample_id"] = extracted[0]
    df["condition"] = extracted[1]

    return df


def tune_xgb_hyperparams_cv(
    X_train,
    y_train,
    groups_train,
    base_params,
    n_splits: int = 5,
):
    """ StratifiedGroupKFold ， AUC """

    print("\n=====  (StratifiedGroupKFold, 5-fold) =====")
    y_train = np.asarray(y_train)
    groups_train = np.asarray(groups_train)

    if len(np.unique(y_train)) < 2:
        print("⚠️ ， base_params")
        return base_params

    X_arr = (
        X_train
        if isinstance(X_train, np.ndarray)
        else X_train.astype(np.float32).values
    )

    # ： donor ， K 
    cv = StratifiedKFold(
        n_splits=n_splits, shuffle=True, random_state=GLOBAL_SEED
    )

    param_grid = {
        "max_depth": [3, 4],
        "learning_rate": [0.03, 0.06],
        "subsample": [0.7, 0.9],
        "colsample_bytree": [0.6, 0.8],
        "min_child_weight": [1, 3],
    }
    from itertools import product

    keys = list(param_grid.keys())
    value_lists = [param_grid[k] for k in keys]
    param_list = []
    for values in product(*value_lists):
        p = {k: v for k, v in zip(keys, values)}
        param_list.append(p)

    print(f": {len(param_list)}")

    best_auc = -np.inf
    best_params = base_params.copy()

    for i, p_extra in enumerate(param_list, 1):
        params = base_params.copy()
        params.update(p_extra)
        params["n_estimators"] = 600

        fold_aucs = []
        print(f"  [ {i}/{len(param_list)}] : {p_extra}")

        for fold, (tr_idx, va_idx) in enumerate(
            cv.split(X_arr, y_train), 1
        ):
            model = XGBClassifier(**params)
            model.fit(X_arr[tr_idx], y_train[tr_idx])

            proba = model.predict_proba(X_arr[va_idx])[
                :, np.where(model.classes_ == 1)[0][0]
            ]
            try:
                auc = roc_auc_score(y_train[va_idx], proba)
            except ValueError:
                print(f"    ->  {fold} ，")
                continue
            fold_aucs.append(auc)
            print(f"    ->  {fold}  AUC = {auc:.4f}")

        if len(fold_aucs) == 0:
            print("    -> ，")
            continue

        mean_auc = float(np.mean(fold_aucs))
        print(f"    ->  AUC = {mean_auc:.4f}")

        if mean_auc > best_auc:
            best_auc = mean_auc
            best_params = params.copy()

    print("\n=====  =====")
    if best_auc == -np.inf:
        print("⚠️  AUC， base_params")
        return base_params

    print(f" AUC = {best_auc:.4f}")
    print(":")
    for k, v in best_params.items():
        print(f"  {k}: {v}")

    return best_params


# =========================================================
# ：（ donor ）
# =========================================================
class PairwiseXGBSHAPCellAnalyzer:
    def __init__(
        self,
        data_path,
        cell_types=None,
        results_dir="pairwise_xgb_boruta_results",
    ):
        self.data_path = data_path

        if cell_types is None:
            df_ct = pd.read_parquet(self.data_path, columns=["celltype"])
            names = sorted(df_ct["celltype"].dropna().unique().tolist())
        else:
            names = sorted(set(list(cell_types)))

        self.target_cell_type_names = names
        self.output_dir = results_dir
        os.makedirs(self.output_dir, exist_ok=True)

        #  pairwise 
        self.pair_classifiers = [
            ("CON", "PRE"),
            ("PRE", "T2D"),
            ("CON", "T2D"),
        ]

        print("\n=====  70%/30% /， donor =====\n")

    def analyze_all_cell_types(self):
        for cell_type_name in self.target_cell_type_names:
            print(f" {cell_type_name} ")
            data_cell = load_celltype_data(self.data_path, cell_type_name)
            if data_cell is None or data_cell.shape[0] == 0:
                print(f" {cell_type_name} ，")
                continue

            print(f"{cell_type_name} : {len(data_cell)}")

            if data_cell.shape[0] < 10:
                print(f" {cell_type_name} ，！")
                del data_cell
                gc.collect()
                continue

            outdir = os.path.join(
                self.output_dir,
                cell_type_name.replace(" ", "_").replace("/", "_"),
            )
            os.makedirs(outdir, exist_ok=True)
            os.makedirs(os.path.join(outdir, "shap"), exist_ok=True)

            analyzer = PairwiseXGBSHAPAnalyzerCore(
                data_cell,
                outdir,
                self.pair_classifiers,
                cell_type_name,
            )
            analyzer.run_all()
            with open(
                f"{self.output_dir}/results__{cell_type_name.replace(' ', '_')}.pkl",
                "wb",
            ) as f:
                pickle.dump(analyzer.results, f)
            print(
                f" {cell_type_name}  "
                f"{self.output_dir}/results__{cell_type_name.replace(' ', '_')}.pkl"
            )

            del analyzer
            del data_cell
            gc.collect()

        print("\n！")


# =========================================================
# ： pairwise XGB + Boruta
# =========================================================
class PairwiseXGBSHAPAnalyzerCore:
    def __init__(
        self,
        data,
        output_dir,
        pair_classifiers,
        cell_type_name,
    ):
        self.data = data
        self.output_dir = output_dir
        self.cell_type_name = cell_type_name
        self.pair_classifiers = pair_classifiers

        self.models = {}
        self.marker_genes = {}
        self.shap_values = {}
        self.results = {}

    # -----------------------------
    # ： pair (CON vs PRE )， 70%/30% 
    # -----------------------------
    def prepare_data_for_pair(self, negative_class, positive_class):
        """
        -  pair 
        -  y (positive_class  1)
        -  70%  / 30% （stratify=y）
        - groups_train  sample_id， GroupCV / Boruta 
        """
        data_pair = self.data[
            self.data["condition"].isin([negative_class, positive_class])
        ].copy()
        if data_pair["condition"].nunique() < 2:
            return None

        y = np.where(data_pair["condition"] == positive_class, 1, 0)

        feature_cols = [
            col
            for col in data_pair.columns
            if col
            not in [
                "cell_id",
                "orig.ident",
                "seurat_clusters",
                "celltype",
                "group",
                "sample_id",
                "condition",
            ]
        ]
        X = data_pair[feature_cols].astype(np.float32)

        #  70/30 ， stratify 
        try:
            X_train, X_test, y_train, y_test = train_test_split(
                X,
                y,
                test_size=0.3,
                stratify=y,
                random_state=GLOBAL_SEED,
            )
        except ValueError as e:
            print(
                f"  ✗ {negative_class} vs {positive_class}: "
                f" 70/30 stratify （）：{e}"
            )
            return None

        #  donor ， GroupCV / Boruta
        groups_train = data_pair.loc[X_train.index, "sample_id"].values

        sample_info = {
            "train": sorted(
                data_pair.loc[X_train.index, "sample_id"].unique().tolist()
            ),
            "test": sorted(
                data_pair.loc[X_test.index, "sample_id"].unique().tolist()
            ),
        }
        print(
            f"  :  {X_train.shape[0]},  {X_test.shape[0]}, "
            f" donor  {len(sample_info['train'])},  donor  {len(sample_info['test'])}"
        )
        feature_names = X.columns.tolist()
        class_names = np.array([negative_class, positive_class])

        return (
            X_train,
            X_test,
            y_train,
            y_test,
            class_names,
            sample_info,
            feature_names,
            groups_train,
        )

    # -----------------------------
    # Boruta ：，
    # -----------------------------
    def run_boruta_stability(
        self,
        X_train_filtered: pd.DataFrame,
        y_train: np.ndarray,
        groups_train: np.ndarray,
        anova_top_k: int = 3000,    # ，
        outer_folds: int = 5,
        boruta_n_trials: int = 20,
        freq_threshold: float = 0.6,
    ):
        """ 5  BorutaShap，，>=freq_threshold """

        y_train = np.asarray(y_train)
        X_df = X_train_filtered
        genes_all = X_df.columns.to_list()
        p = X_df.shape[1]
        print(f"  [Boruta] 5-fold  BorutaShap， {p} ")

        # Boruta ： K ， donor 
        cv = StratifiedKFold(
            n_splits=outer_folds, shuffle=True, random_state=GLOBAL_SEED
        )

        fold_selected = []

        for fold, (tr_idx, va_idx) in enumerate(
            cv.split(X_df.values, y_train), 1
        ):
            print(f"  [Boruta]  {fold}/{outer_folds} ")

            X_tr = X_df.iloc[tr_idx]
            y_tr = y_train[tr_idx]

            pos_fold = (y_tr == 1).sum()
            neg_fold = (y_tr == 0).sum()
            scale_fold = (
                neg_fold / max(pos_fold, 1) if pos_fold > 0 else 1.0
            )

            boruta_params = get_xgb_params_by_count(
                n_cells=len(y_tr),
                scale_pos_weight=scale_fold,
                random_state=GLOBAL_SEED + fold,
            )
            boruta_params.pop("early_stopping_rounds", None)
            boruta_params["n_estimators"] = 200

         

            selector = BorutaShap(
                model = XGBClassifierForBoruta(**boruta_params),
                importance_measure="gini",
                classification=True,
            )

            selector.fit(
                X=X_tr,
                y=y_tr,
                n_trials=boruta_n_trials,
                sample=False,
                verbose=False,
            )

            features_to_remove = set(selector.features_to_remove)
            selected_genes_fold = [
                g for g in genes_all if g not in features_to_remove
            ]
            print(
                f"    ->  Boruta  = {len(selected_genes_fold)}"
            )
            fold_selected.append(set(selected_genes_fold))

        # 
        gene_counts = defaultdict(int)
        for s in fold_selected:
            for g in s:
                gene_counts[g] += 1

        stab_list = []
        for g in genes_all:
            c = gene_counts.get(g, 0)
            freq = c / float(outer_folds)
            stab_list.append((g, c, freq))

        stability_df = pd.DataFrame(
            stab_list, columns=["gene", "fold_count", "fold_freq"]
        ).sort_values("fold_freq", ascending=False)

        selected_genes = stability_df.loc[
            stability_df["fold_freq"] >= freq_threshold, "gene"
        ].tolist()

        if len(selected_genes) == 0:
            print(
                f"  [Boruta] ： {freq_threshold}， 100 "
            )
            top_k = min(100, stability_df.shape[0])
            selected_genes = stability_df.head(top_k)["gene"].tolist()

        print(
            f"  [Boruta] ，(>={freq_threshold}) = {len(selected_genes)}"
        )

        return selected_genes, stability_df, outer_folds

    # -----------------------------
    # ： pairwise  + 
    # -----------------------------
    def run_all(self):
        print(
            f"\n=== {self.output_dir} "
            f" + Boruta  +  + SHAP  "
            f"( 70/30 /) ==="
        )

        os.makedirs(self.output_dir, exist_ok=True)
        shap_dir = os.path.join(self.output_dir, "shap")
        os.makedirs(shap_dir, exist_ok=True)

        roc_curves = []

        for negative_class, positive_class in self.pair_classifiers:
            print(f"\n==== {negative_class}(0) vs {positive_class}(1) ====")
            prep = self.prepare_data_for_pair(negative_class, positive_class)
            if prep is None:
                print("  ，")
                continue

            (
                X_train,
                X_test,
                y_train,
                y_test,
                classes,
                sample_info,
                feature_names,
                groups_train,
            ) = prep

            print(
                f": {X_train.shape[0]}, "
                f": {X_test.shape[0]}"
            )

            # ========= 1. ：>5%  =========
            n_cells = X_train.shape[0]
            min_expr_pct = 0.05
            #  >5%：floor(0.05*n)+1，5
            min_cells = max(5, int(np.floor(n_cells * min_expr_pct)) + 1)
            print(
                f":  > {min_expr_pct*100:.1f}%  "
                f"( >= {min_cells} )"
            )

            genes_keep_mask = (X_train > 0).sum(axis=0) >= min_cells
            genes_keep = X_train.columns[genes_keep_mask]

            X_train_filtered = X_train[genes_keep]
            X_test_filtered = X_test[genes_keep]
            feature_names_filtered = genes_keep.tolist()

            print(
                f": {X_train_filtered.shape[1]} / "
                f"{X_train.shape[1]}"
            )

            print(
                f": {classes[0]}=0(), {classes[1]}=1()"
            )
            print(":", np.unique(y_train, return_counts=True))
            print(":", np.unique(y_test, return_counts=True))

            pos_count = np.sum(y_train == 1)
            neg_count = np.sum(y_train == 0)
            scale_weight = (
                neg_count / max(pos_count, 1) if pos_count > 0 else 1.0
            )

            # ========= 2. ：5  BorutaShap  =========
            stable_genes, stability_df, n_folds_boruta = self.run_boruta_stability(
                X_train_filtered=X_train_filtered,
                y_train=y_train,
                groups_train=groups_train,
                anova_top_k=3000,
                outer_folds=5,
                boruta_n_trials=20,
                freq_threshold=0.6,
            )

            # 
            stab_path = os.path.join(
                self.output_dir,
                f"boruta_gene_stability_{negative_class}_vs_{positive_class}.csv",
            )
            stability_df.to_csv(stab_path, index=False)
            print(f"Boruta : {stab_path}")

            if len(stable_genes) == 0:
                print("  ⚠️ ， pair")
                continue

            selected_gene_names = stable_genes
            print(f": {len(selected_gene_names)}")

            X_train_final_full = X_train_filtered[selected_gene_names].values
            X_test_final = X_test_filtered[selected_gene_names].values

            # ========= 3.  5  StratifiedGroupKFold  XGBoost （） =========
            base_params_final = get_xgb_params_by_count(
                n_cells=X_train.shape[0],
                scale_pos_weight=scale_weight,
                random_state=GLOBAL_SEED,
            )
            tuned_params = tune_xgb_hyperparams_cv(
                X_train=X_train_final_full,
                y_train=y_train,
                groups_train=groups_train,
                base_params=base_params_final,
                n_splits=5,
            )

            # ========= 4.  F1 ， 0.5 =========
            best_thresh = 0.5

            # ========= 5. ： tuned_params  =========
            final_model = XGBClassifier(**tuned_params)
            final_model.fit(X_train_final_full, y_train)

            pos_class_idx = np.where(final_model.classes_ == 1)[0]
            if len(pos_class_idx) == 0:
                raise ValueError(" classes_ (1)")
            pos_col = pos_class_idx[0]

            # ========= 6. （） =========
            proba_train_all = final_model.predict_proba(
                X_train_final_full
            )[:, pos_col]
            auc_train_all = roc_auc_score(y_train, proba_train_all)
            print(f"AUC: {auc_train_all:.4f}")

            # ========= 7. （30% ） =========
            y_pred_proba_full = final_model.predict_proba(X_test_final)
            y_pred_proba = y_pred_proba_full[:, pos_col]
            roc_auc = roc_auc_score(y_test, y_pred_proba)

            print(f" classes_: {final_model.classes_}")
            print(f"AUC: {roc_auc:.4f}")

            # ========= 8. SHAP （） =========
            try:
                X_train_final_df = pd.DataFrame(
                    X_train_final_full, columns=selected_gene_names
                )
                explainer = shap.TreeExplainer(
                    final_model, feature_perturbation="tree_path_dependent"
                )
                shap_values = explainer(X_train_final_df)

                mean_shap = shap_values.values.mean(axis=0)
                mean_abs_shap = np.abs(shap_values.values).mean(axis=0)

                shap_df = pd.DataFrame(
                    {
                        "gene": selected_gene_names,
                        "mean_SHAP": mean_shap,
                        "mean_abs_SHAP": mean_abs_shap,
                    }
                ).sort_values("mean_abs_SHAP", ascending=False)

                shap_table_path = os.path.join(
                    self.output_dir,
                    f"{negative_class}_vs_{positive_class}_shap_marker_genes.csv",
                )
                shap_df.to_csv(shap_table_path, index=False)
                print(f"SHAP: {shap_table_path}")

                topn = min(20, len(selected_gene_names))
                plt.figure(figsize=(8, 6))
                shap.plots.bar(shap_values, max_display=topn, show=False)
                plt.title(
                    f"{negative_class}(0) vs {positive_class}(1) - SHAP Top{topn}"
                )
                plt.tight_layout()
                out_path = os.path.join(
                    self.output_dir,
                    f"shap/{negative_class}_vs_{positive_class}_shap_bar.png",
                )
                plt.savefig(out_path, dpi=300)
                plt.close()

                plt.figure(figsize=(8, 6))
                shap.summary_plot(
                    shap_values.values,
                    X_train_final_df,
                    feature_names=selected_gene_names,
                    max_display=topn,
                    show=False,
                )
                plt.title(
                    f"{negative_class}(0) vs {positive_class}(1) - SHAP Summary"
                )
                plt.tight_layout()
                out_path2 = os.path.join(
                    self.output_dir,
                    f"shap/{negative_class}_vs_{positive_class}_shap_summary.png",
                )
                plt.savefig(out_path2, dpi=300)
                plt.close()

            except Exception as e:
                import traceback

                print(f"SHAP: {e}")
                tb = traceback.format_exc()
                print(tb)

                # ：xgboost  base_score  "[5E-1]" ，
                # shap.TreeExplainer  base_score  float 
                # ： "base_score" （），
                #  traceback  shap  xgboost loader 
                msg = str(e)
                is_tree_loader_tb = (
                    "shap/explainers/_tree.py" in tb
                    and ("XGBTreeModelLoader" in tb or "learner_model_param" in tb)
                )
                looks_like_bracketed_float = (
                    "could not convert string to float" in msg and "[" in msg and "]" in msg
                )
                if (
                    ("could not convert string to float" in msg and ("base_score" in msg or is_tree_loader_tb))
                    or looks_like_bracketed_float
                ):
                    try:
                        print(" XGBoost pred_contribs  SHAP (fallback)...")

                        booster = final_model.get_booster()
                        dtrain = xgb.DMatrix(
                            X_train_final_full,
                            feature_names=list(selected_gene_names),
                        )
                        contribs = booster.predict(dtrain, pred_contribs=True)

                        # contribs: (n_samples, n_features + 1),  bias/base_value
                        values = contribs[:, :-1]
                        base_values = contribs[:, -1]

                        shap_values_fb = shap.Explanation(
                            values=values,
                            base_values=base_values,
                            data=X_train_final_df,
                            feature_names=list(selected_gene_names),
                        )

                        mean_shap = values.mean(axis=0)
                        mean_abs_shap = np.abs(values).mean(axis=0)
                        shap_df = pd.DataFrame(
                            {
                                "gene": selected_gene_names,
                                "mean_SHAP": mean_shap,
                                "mean_abs_SHAP": mean_abs_shap,
                            }
                        ).sort_values("mean_abs_SHAP", ascending=False)

                        shap_table_path = os.path.join(
                            self.output_dir,
                            f"{negative_class}_vs_{positive_class}_shap_marker_genes.csv",
                        )
                        shap_df.to_csv(shap_table_path, index=False)
                        print(f"SHAP( fallback ): {shap_table_path}")

                        topn = min(20, len(selected_gene_names))

                        plt.figure(figsize=(8, 6))
                        shap.plots.bar(shap_values_fb, max_display=topn, show=False)
                        plt.title(
                            f"{negative_class}(0) vs {positive_class}(1) - SHAP Top{topn}"
                        )
                        plt.tight_layout()
                        out_path = os.path.join(
                            self.output_dir,
                            f"shap/{negative_class}_vs_{positive_class}_shap_bar.png",
                        )
                        plt.savefig(out_path, dpi=300)
                        plt.close()

                        plt.figure(figsize=(8, 6))
                        shap.summary_plot(
                            values,
                            X_train_final_df,
                            feature_names=selected_gene_names,
                            max_display=topn,
                            show=False,
                        )
                        plt.title(
                            f"{negative_class}(0) vs {positive_class}(1) - SHAP Summary"
                        )
                        plt.tight_layout()
                        out_path2 = os.path.join(
                            self.output_dir,
                            f"shap/{negative_class}_vs_{positive_class}_shap_summary.png",
                        )
                        plt.savefig(out_path2, dpi=300)
                        plt.close()

                        print("SHAP fallback ")
                    except Exception as e2:
                        print(f"SHAP fallback : {e2}")
                        traceback.print_exc()

            # ========= 9. （ 30% ） =========
            y_pred_default = (y_pred_proba >= best_thresh).astype(int)
            best_f1 = f1_score(y_test, y_pred_default, zero_division=0)
            print(f" {best_thresh:.4f} :")
            print(
                classification_report(
                    y_test, y_pred_default, target_names=classes
                )
            )
            cm = confusion_matrix(y_test, y_pred_default)
            plt.figure(figsize=(5, 4))
            sns.heatmap(
                cm,
                annot=True,
                fmt="d",
                cmap="Blues",
                xticklabels=classes,
                yticklabels=classes,
            )
            plt.xlabel("Predicted")
            plt.ylabel("Actual")
            plt.title(
                f"Confusion Matrix: {negative_class}(0) vs {positive_class}(1)"
            )
            plt.tight_layout()
            save_path = (
                f"{self.output_dir}/confmat_{negative_class}_vs_{positive_class}.png"
            )
            plt.savefig(save_path, dpi=300)
            plt.close()

            # ========= 10. ROC  =========
            fpr, tpr, _ = roc_curve(y_test, y_pred_proba)
            plt.figure(figsize=(6, 5))
            plt.plot(fpr, tpr, "b", label=f"AUC = {roc_auc:.4f}")
            plt.plot([0, 1], [0, 1], "k--", label="Random")
            plt.xlabel("False Positive Rate")
            plt.ylabel("True Positive Rate")
            plt.title(
                f"ROC: {self.cell_type_name} "
                f"({negative_class}(0) vs {positive_class}(1))"
            )
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            save_path = (
                f"{self.output_dir}/roc_{negative_class}_vs_{positive_class}.png"
            )
            plt.savefig(save_path, dpi=300)
            plt.close()
            roc_curves.append(
                (
                    fpr,
                    tpr,
                    roc_auc,
                    f"{negative_class}(0) vs {positive_class}(1)",
                )
            )

            # ========= 11.  =========
            self.results[(negative_class, positive_class)] = {
                "y_test": y_test,
                "y_pred": y_pred_default,
                "y_pred_proba": y_pred_proba,
                "roc_auc": roc_auc,
                "confusion_matrix": cm,
                "classification_report": classification_report(
                    y_test, y_pred_default, target_names=classes, output_dict=True
                ),
                "samples": sample_info,
                "selected_features": list(selected_gene_names),
                "boruta_stability": stability_df,
                "boruta_n_folds": n_folds_boruta,
                "boruta_freq_threshold": 0.6,
                "final_params": final_model.get_params(),
                "best_threshold": float(best_thresh),
                "best_f1": float(best_f1),
                "class_labels": {
                    "negative": negative_class,
                    "positive": positive_class,
                },
                "model_classes": final_model.classes_.tolist(),
            }

            sample_info_df = pd.DataFrame(
                {
                    "sample_type": ["train"] * len(sample_info["train"])
                    + ["test"] * len(sample_info["test"]),
                    "sample_id": sample_info["train"] + sample_info["test"],
                }
            )
            sample_info_df.to_csv(
                f"{self.output_dir}/sample_split_{negative_class}_vs_{positive_class}.csv",
                index=False,
            )

        # ========= 12.  ROC  =========
        if len(roc_curves) > 0:
            plt.figure(figsize=(7, 6))
            colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]
            for i, (fpr_i, tpr_i, auc_i, label_i) in enumerate(roc_curves):
                c = colors[i % len(colors)]
                plt.plot(
                    fpr_i,
                    tpr_i,
                    linestyle="-",
                    color=c,
                    linewidth=2,
                    label=f"{label_i} (AUC={auc_i:.3f})",
                )
            plt.plot([0, 1], [0, 1], "k--", linewidth=1, label="Random")
            plt.xlim([0.0, 1.0])
            plt.ylim([0.0, 1.05])
            plt.xlabel("False Positive Rate")
            plt.ylabel("True Positive Rate")
            plt.title(f"ROC comparison: {self.cell_type_name}")
            plt.legend(loc="lower right", frameon=True)
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            save_path_all = os.path.join(self.output_dir, "roc_all_pairs.png")
            plt.savefig(save_path_all, dpi=300)
            plt.close()
            print(f"\n pairwise ！")


# =========================================================
# 
# =========================================================
if __name__ == "__main__":
    cell_types = ["Acinar cell", "Ductal cell", "Immune cell", "Stellate cell","PP cell","δ cell"]
    analyzer = PairwiseXGBSHAPCellAnalyzer(
        data_path=str(ATAC_PEAK_PARQUET),
        cell_types=cell_types,
    )
    analyzer.analyze_all_cell_types()
    print(
        "\n！ pairwise  "
        "pairwise_xgb_boruta_results "
    )
