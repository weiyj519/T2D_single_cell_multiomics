#!/usr/bin/env python3

"""
Script name: 00_compute_scriskcell_scores.py
Purpose: Cleaned public workflow script for the T2D single-cell multiome analysis.
Input: Paths and analysis settings are read from configs/config.yaml.
Output: Module-specific outputs are written under the configured results directory.
Main steps: Load configuration, run the original analysis logic, and export reproducible outputs.
Example command: python scripts/08_pseudotime/00_compute_scriskcell_scores.py
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

import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
import scRiskCell
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu


def remove_extreme_outliers_iqr(df, cols, iqr_k=1.5):
    """ IQR ，"""
    out = df.copy()
    mask = pd.Series(True, index=out.index)
    for col in cols:
        s = pd.to_numeric(out[col], errors="coerce")
        q1 = s.quantile(0.25)
        q3 = s.quantile(0.75)
        iqr = q3 - q1
        lower = q1 - iqr_k * iqr
        upper = q3 + iqr_k * iqr
        mask &= s.between(lower, upper)
    return out.loc[mask].copy()


CATEGORY_ORDER = ["CON", "PRE", "T2D"]
CATEGORY_PAIRS = [("CON", "PRE"), ("PRE", "T2D"), ("CON", "T2D")]


def p_adjust_bh(pvalues):
    """Benjamini-Hochberg FDR correction without requiring statsmodels."""
    pvalues = np.asarray(pvalues, dtype=float)
    adjusted = np.full(pvalues.shape, np.nan, dtype=float)
    valid = ~np.isnan(pvalues)
    if not valid.any():
        return adjusted

    valid_p = pvalues[valid]
    order = np.argsort(valid_p)
    ranked = valid_p[order]
    n = len(ranked)
    ranked_adjusted = ranked * n / np.arange(1, n + 1)
    ranked_adjusted = np.minimum.accumulate(ranked_adjusted[::-1])[::-1]
    ranked_adjusted = np.clip(ranked_adjusted, 0, 1)

    valid_adjusted = np.empty_like(ranked_adjusted)
    valid_adjusted[order] = ranked_adjusted
    adjusted[valid] = valid_adjusted
    return adjusted


def significance_stars(pvalue):
    if pd.isna(pvalue):
        return "NA"
    if pvalue < 1e-4:
        return "****"
    if pvalue < 1e-3:
        return "***"
    if pvalue < 1e-2:
        return "**"
    if pvalue < 5e-2:
        return "*"
    return "ns"


def pairwise_mannwhitney_tests(df, value_col, pairs=CATEGORY_PAIRS):
    """Pairwise two-sided Mann-Whitney U tests with BH-adjusted q values."""
    rows = []
    for group1, group2 in pairs:
        values1 = pd.to_numeric(
            df.loc[df["Category"] == group1, value_col], errors="coerce"
        ).dropna()
        values2 = pd.to_numeric(
            df.loc[df["Category"] == group2, value_col], errors="coerce"
        ).dropna()
        if len(values1) == 0 or len(values2) == 0:
            stat = np.nan
            pvalue = np.nan
        else:
            try:
                result = mannwhitneyu(
                    values1,
                    values2,
                    alternative="two-sided",
                    method="asymptotic",
                )
            except TypeError:
                result = mannwhitneyu(values1, values2, alternative="two-sided")
            stat = result.statistic
            pvalue = result.pvalue

        rows.append({
            "group1": group1,
            "group2": group2,
            "n1": len(values1),
            "n2": len(values2),
            "median1": values1.median() if len(values1) else np.nan,
            "median2": values2.median() if len(values2) else np.nan,
            "mannwhitney_u": stat,
            "p_value": pvalue,
        })

    stats_df = pd.DataFrame(rows)
    stats_df["q_value_bh"] = p_adjust_bh(stats_df["p_value"].to_numpy())
    stats_df["significance"] = stats_df["q_value_bh"].map(significance_stars)
    return stats_df


def add_significance_annotations(ax, stats_df, order, value_col, data):
    """Draw significance brackets above category plots."""
    present = set(order)
    stats_df = stats_df[
        stats_df["group1"].isin(present) & stats_df["group2"].isin(present)
    ].copy()
    if stats_df.empty:
        return

    values = pd.to_numeric(data[value_col], errors="coerce").dropna()
    if values.empty:
        return

    y_min = values.min()
    y_max = values.max()
    y_range = y_max - y_min
    if y_range == 0:
        y_range = max(abs(y_max), 1.0)

    bracket_height = y_range * 0.025
    step = y_range * 0.12
    start = y_max + y_range * 0.08
    positions = {category: idx for idx, category in enumerate(order)}

    for i, (_, row) in enumerate(stats_df.iterrows()):
        x1 = positions[row["group1"]]
        x2 = positions[row["group2"]]
        y = start + i * step
        label = row["significance"]
        ax.plot(
            [x1, x1, x2, x2],
            [y, y + bracket_height, y + bracket_height, y],
            color="black",
            linewidth=1.2,
        )
        ax.text(
            (x1 + x2) / 2,
            y + bracket_height,
            label,
            ha="center",
            va="bottom",
            color="black",
            fontsize=11,
        )

    ax.set_ylim(y_min - y_range * 0.05, start + len(stats_df) * step + y_range * 0.08)


def plot_category_with_significance(
    df,
    value_col,
    ylabel,
    colors,
    save_path,
    kind="box",
    stats_df=None,
    y_limits=None,
):
    """Plot category violin/boxplot with pairwise significance annotations."""
    order = [c for c in CATEGORY_ORDER if c in set(df["Category"])]
    order += sorted(set(df["Category"]) - set(order))
    palette = dict(zip(CATEGORY_ORDER, colors))

    width = max(5, 3 + len(order) * 1)
    plt.figure(figsize=(width, 5))
    ax = plt.gca()

    if kind == "violin":
        sns.violinplot(
            x="Category",
            y=value_col,
            data=df,
            order=order,
            hue="Category",
            hue_order=order,
            palette=palette,
            width=0.65,
            cut=0,
            inner="box",
            legend=False,
            ax=ax,
        )
    else:
        sns.boxplot(
            x="Category",
            y=value_col,
            data=df,
            order=order,
            hue="Category",
            hue_order=order,
            palette=palette,
            width=0.65,
            showfliers=False,
            legend=False,
            ax=ax,
        )

    ax.set_xlabel("Category")
    ax.set_ylabel(ylabel)
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)

    if y_limits is not None:
        ax.set_ylim(*y_limits)

    if stats_df is None:
        stats_df = pairwise_mannwhitney_tests(df, value_col)
    add_significance_annotations(ax, stats_df, order, value_col, df)

    plt.tight_layout()
    plt.savefig(save_path, bbox_inches="tight")
    plt.close()
    return stats_df

# =========================
# 0. 
# =========================
infile = str(RNA_EXPRESSION_PARQUET)
df = pd.read_parquet(infile)

print(":", df.shape)
print("celltype:")
print(df["celltype"].value_counts())
print("group:")
print(df["group"].value_counts())

# =========================
# 1.  Acinar 
# =========================
acinar_df = df[df["celltype"] == "Acinar cell"].copy()

print("\nAcinar:", acinar_df.shape)
print("Acinar:")
print(acinar_df["group"].value_counts())

# =========================
# 2.  scRiskCell 
# =========================
acinar_df = acinar_df.rename(columns={
    "cell_id": "Cell_id",
    "orig.ident": "Donor",
    "group": "Category"
})

# 
label_map = {"CON": 0, "PRE": 1, "T2D": 2}
acinar_df["Label"] = acinar_df["Category"].map(label_map)

# 
if acinar_df["Label"].isna().any():
    bad = acinar_df.loc[acinar_df["Label"].isna(), "Category"].unique()
    raise ValueError(f"CategoryLabel: {bad}")

# =========================
# 3. 
# =========================

meta_cols = ["Cell_id", "Donor", "Category", "Label", "celltype", "seurat_clusters"]
gene_cols = [c for c in acinar_df.columns if c not in meta_cols]

print("\n:", len(gene_cols))
print("10:", gene_cols[:10])

# ，0
acinar_df[gene_cols] = acinar_df[gene_cols].apply(pd.to_numeric, errors="coerce").fillna(0)

# =========================
# 4.  Acinar
# =========================
scaler = StandardScaler(with_mean=True, with_std=True)
acinar_df[gene_cols] = scaler.fit_transform(acinar_df[gene_cols])

#  scRiskCell 
sc_df = acinar_df[gene_cols + ["Donor", "Category", "Label", "Cell_id"]].copy()

print("\nscRiskCell:", sc_df.shape)
print(sc_df.iloc[:5, -4:])

# 
print("\n:")
print("global min:", sc_df[gene_cols].min().min())
print("global max:", sc_df[gene_cols].max().max())
print("fraction negative:", (sc_df[gene_cols] < 0).sum().sum() / sc_df[gene_cols].size)

# 
sc_df.to_parquet("results/downstream/riskcell/acinar/acinar_scRiskCell_input.parquet", index=False)


threshold_params ={"ratio": 0.70}
# =========================
# 5.  scRiskCell：CON/T2D ，PRE 
# =========================
pca_endpoint_df, pca_intermediate_df, model, index_df, threshold, risk_df, donor_risk_ratio = scRiskCell.scRiskCell(
    df=sc_df,
    category_group1=["CON", "T2D"],
    category_group2=["PRE"],
    n_components=20,
    get_threshold_params=threshold_params
)

print("\nDisease index:")
print(index_df.head())

print("\nthreshold =", threshold)

print("\nRisk:")
print(risk_df["Risk"].value_counts(dropna=False))

print("\ndonor:")
print(donor_risk_ratio)

# 
index_df.to_csv("results/downstream/riskcell/acinar/acinar_index_df.tsv", sep="\t", index=False)
risk_df.to_csv("results/downstream/riskcell/acinar/acinar_risk_df.tsv", sep="\t", index=False)
donor_risk_ratio.to_csv("results/downstream/riskcell/acinar/acinar_donor_risk_ratio.tsv", sep="\t", index=False)


print("results/downstream/riskcell/acinar/acinar_index_df.tsv")
print("results/downstream/riskcell/acinar/acinar_risk_df.tsv")
print("results/downstream/riskcell/acinar/acinar_donor_risk_ratio.tsv")





# =========================
# A. CON vs PRE
# =========================
df_cp = sc_df[sc_df["Category"].isin(["CON", "PRE"])].copy()
df_cp["Label"] = df_cp["Category"].map({"CON": 0, "PRE": 1})

_, _, model_cp, index_cp, th_cp, risk_cp, donor_cp = scRiskCell.scRiskCell(
    df=df_cp,
    category_group1=["CON", "PRE"],
    category_group2=[],
    n_components=20,
    get_threshold_params=threshold_params
)

print("CON vs PRE threshold:", th_cp)
print(risk_cp["Risk"].value_counts())

# =========================
# B. PRE vs T2D
# =========================
df_pt = sc_df[sc_df["Category"].isin(["PRE", "T2D"])].copy()
df_pt["Label"] = df_pt["Category"].map({"PRE": 0, "T2D": 1})

_, _, model_pt, index_pt, th_pt, risk_pt, donor_pt = scRiskCell.scRiskCell(
    df=df_pt,
    category_group1=["PRE", "T2D"],
    category_group2=[],
    n_components=20,
    get_threshold_params=threshold_params
)

print("PRE vs T2D threshold:", th_pt)
print(risk_pt["Risk"].value_counts())

# =========================
# C. 
# =========================
cp_keep = risk_cp[["Cell_id", "Donor", "Category", "Disease_index", "Risk"]].copy()
cp_keep = cp_keep.rename(columns={
    "Disease_index": "DI_CON_PRE",
    "Risk": "Risk_CON_PRE"
})

pt_keep = risk_pt[["Cell_id", "Donor", "Category", "Disease_index", "Risk"]].copy()
pt_keep = pt_keep.rename(columns={
    "Disease_index": "DI_PRE_T2D",
    "Risk": "Risk_PRE_T2D"
})

merged = sc_df[["Cell_id", "Donor", "Category"]].drop_duplicates().copy()
merged = merged.merge(cp_keep[["Cell_id", "DI_CON_PRE", "Risk_CON_PRE"]], on="Cell_id", how="left")
merged = merged.merge(pt_keep[["Cell_id", "DI_PRE_T2D", "Risk_PRE_T2D"]], on="Cell_id", how="left")

# =========================
# D. “”Acinar
# =========================
# CONCON-like
con_stage = merged[
    (merged["Category"] == "CON") &
    (merged["Risk_CON_PRE"] == 0)
].copy()

# PRE：CON，T2D
pre_stage = merged[
    (merged["Category"] == "PRE") &
    (merged["Risk_CON_PRE"] == 1) &
    (merged["Risk_PRE_T2D"] == 0)
].copy()

# T2DT2D-like
t2d_stage = merged[
    (merged["Category"] == "T2D") &
    (merged["Risk_PRE_T2D"] == 1)
].copy()

selected = pd.concat([con_stage, pre_stage, t2d_stage], ignore_index=True)
selected["Stage_selected"] = selected["Category"]

print("\n：")
print(selected["Stage_selected"].value_counts())

# 
selected.to_csv("results/downstream/riskcell/acinar/acinar_selected_cells_for_pseudotime.tsv", sep="\t", index=False)

print("\n：results/downstream/riskcell/acinar/acinar_selected_cells_for_pseudotime.tsv")
print(selected.head())


# 
colors3 = ["#D4D4D4", "#F4B7AD", "#90A4C4"]   # CON, PRE, T2D
colors2 = ["#D9D9D9", "#B2ADD5"]              # non-risk, risk

# （，）
index_df_plot = remove_extreme_outliers_iqr(index_df, ["Disease_index"], iqr_k=1.5)
print(f"index_df : {len(index_df) - len(index_df_plot)} / {len(index_df)}")

# Disease index 
index_category_stats = pairwise_mannwhitney_tests(index_df_plot, "Disease_index")

plot_category_with_significance(
    df=index_df_plot,
    value_col="Disease_index",
    ylabel="Disease index",
    colors=colors3,
    save_path="results/downstream/riskcell/acinar/acinar_index_violin_by_category.pdf",
    kind="violin",
    stats_df=index_category_stats,
)

plot_category_with_significance(
    df=index_df_plot,
    value_col="Disease_index",
    ylabel="Disease index",
    colors=colors3,
    save_path="results/downstream/riskcell/acinar/acinar_index_boxplot_by_category.pdf",
    kind="box",
    stats_df=index_category_stats,
)

category_significance_stats = index_category_stats.copy()

# Disease index  donor
scRiskCell.PlotIndexViolin_by_Donor(
    index_df=index_df_plot,
    colors=colors3,
    save_path="results/downstream/riskcell/acinar/acinar_index_violin_by_donor.pdf"
)

scRiskCell.PlotIndexBoxplot_by_Donor(
    index_df=index_df_plot,
    colors=colors3,
    save_path="results/downstream/riskcell/acinar/acinar_index_boxplot_by_donor.pdf"
)

# donor risk ratio
risk_ratio_category_stats = pairwise_mannwhitney_tests(donor_risk_ratio, "Risk_ratio")

plot_category_with_significance(
    df=donor_risk_ratio,
    value_col="Risk_ratio",
    ylabel="Risk cell ratio",
    colors=colors3,
    save_path="results/downstream/riskcell/acinar/acinar_riskratio_boxplot_by_category.pdf",
    kind="box",
    stats_df=risk_ratio_category_stats,
)

category_significance_stats = pd.concat(
    [
        index_category_stats.assign(metric="Disease_index"),
        risk_ratio_category_stats.assign(metric="Risk_ratio"),
    ],
    ignore_index=True,
)
category_significance_stats = category_significance_stats[
    [
        "metric",
        "group1",
        "group2",
        "n1",
        "n2",
        "median1",
        "median2",
        "mannwhitney_u",
        "p_value",
        "q_value_bh",
        "significance",
    ]
]
category_significance_stats.to_csv(
    "results/downstream/riskcell/acinar/acinar_category_significance_tests.tsv",
    sep="\t",
    index=False,
)
print("：results/downstream/riskcell/acinar/acinar_category_significance_tests.tsv")

scRiskCell.PlotRatioStackBar_by_Donor(
    donor_risk_ratio=donor_risk_ratio,
    colors=["#D9D9D9", "#B2ADD5"],
    save_path="results/downstream/riskcell/acinar/acinar_riskratio_stackbar_by_donor.pdf"
)

# ROC（）
scRiskCell.PlotROC(
    donor_risk_ratio=donor_risk_ratio,
    color="#4C78A8",
    save_path="results/downstream/riskcell/acinar/acinar_riskratio_ROC.pdf"
)







import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

index_df2 = index_df.copy()
selected_ids = set(selected["Cell_id"])
index_df2["Selected"] = index_df2["Cell_id"].isin(selected_ids)
index_df2["Selected"] = index_df2["Selected"].map({True: "Selected", False: "Unselected"})

plt.figure(figsize=(7, 5))
sns.boxplot(
    data=index_df2,
    x="Category",
    y="Disease_index",
    hue="Selected",
    showfliers=False
)
plt.ylim(-20, 20)
plt.tight_layout()
plt.savefig("results/downstream/riskcell/acinar/acinar_selected_vs_unselected_boxplot.pdf")
plt.show()




tmp = merged.copy()
tmp["Selected"] = tmp["Cell_id"].isin(selected["Cell_id"])

donor_sel = (
    tmp.groupby(["Donor", "Category"])["Selected"]
    .mean()
    .reset_index(name="Selected_ratio")
    .sort_values(["Category", "Selected_ratio"], ascending=[True, False])
)

cat_order = ["CON", "PRE", "T2D"]
cat_palette = dict(zip(cat_order, colors3))

plt.figure(figsize=(10, 5))
sns.barplot(
    data=donor_sel,
    x="Donor",
    y="Selected_ratio",
    hue="Category",
    hue_order=cat_order,
    palette=cat_palette
)
plt.xticks(rotation=45, ha="right")
plt.ylabel("Selected cell ratio")
plt.tight_layout()
plt.savefig("results/downstream/riskcell/acinar/acinar_selected_ratio_by_donor.pdf")
plt.show()