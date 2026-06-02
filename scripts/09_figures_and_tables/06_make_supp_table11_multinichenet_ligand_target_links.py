#!/usr/bin/env python3

"""
Script name: 06_make_supp_table11_multinichenet_ligand_target_links.py
Purpose: Cleaned public workflow script for the T2D single-cell multiome analysis.
Input: Paths and analysis settings are read from configs/config.yaml.
Output: Module-specific outputs are written under the configured results directory.
Main steps: Load configuration, run the original analysis logic, and export reproducible outputs.
Example command: python scripts/09_figures_and_tables/06_make_supp_table11_multinichenet_ligand_target_links.py
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
Make Table S11:
Ligand-target links supporting receiver-cell STAGE gene programs.

The script recursively searches MultiNicheNet ligand-target files under
BASE_DIR, removes rows without a concrete target gene, supplements each ligand
with the highest-priority sender/receptor pair from the matching
prioritization table when available, ranks targets and ligands, and writes CSV
and XLSX outputs.

This version intentionally uses only the Python standard library, so it can run
even in environments without pandas/openpyxl.

Run:
    python3 scripts/09_figures_and_tables/make_supp_table11_multinichenet_ligand_target_links.py

Outputs:
    results/supplementary_tables/Supplementary_Table_11_MultiNicheNet_ligand_target_links.csv
    results/supplementary_tables/Supplementary_Table_11_MultiNicheNet_ligand_target_links.xlsx
"""

import csv
import math
import posixpath
import re
import zipfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from xml.sax.saxutils import escape


# ---------------------------------------------------------------------------
# User-editable settings
# ---------------------------------------------------------------------------

BASE_DIR = Path("results/downstream")
OUT_DIR = Path("results/supplementary_tables")

LIGAND_TARGET_FILENAME = "ligand_activities_targets_MLmasked.csv"
PRIORITIZATION_FILENAME = "prioritization_group_prioritization_table_source.csv"

OUT_CSV = OUT_DIR / "Supplementary_Table_11_MultiNicheNet_ligand_target_links.csv"
OUT_XLSX = OUT_DIR / "Supplementary_Table_11_MultiNicheNet_ligand_target_links.xlsx"

# Set to None to keep all targets per ligand.
TOP_TARGETS_PER_LIGAND: int | None = 20

# Set to None to keep all ligands per receiver/comparison/target_direction_group.
TOP_LIGANDS_PER_GROUP: int | None = 50


OUTPUT_COLUMNS = [
    "receiver_cell_type",
    "comparison",
    "target_direction_group",
    "direction_regulation",
    "ligand",
    "predicted_target_gene",
    "ligand_activity",
    "ligand_activity_scaled",
    "ligand_target_weight",
    "sender_cell_type",
    "receptor",
    "interaction",
    "prioritization_score",
    "target_gene_log2FC",
    "target_gene_padj",
    "rank_within_ligand",
    "rank_within_receiver_comparison_group",
    "source_file",
]

TEXT_COLUMNS = {
    "receiver_cell_type",
    "comparison",
    "target_direction_group",
    "direction_regulation",
    "ligand",
    "predicted_target_gene",
    "sender_cell_type",
    "receptor",
    "interaction",
    "source_file",
}

P_VALUE_COLUMNS = {"target_gene_padj"}
NUMERIC_COLUMNS = set(OUTPUT_COLUMNS) - TEXT_COLUMNS

RECEIVER_ORDER = {"beta cell": 0, "alpha cell": 1, "acinar cell": 2}
COMPARISON_ORDER = {"CON vs PRE": 0, "PRE vs T2D": 1}
GROUP_ORDER = {"CON": 0, "PRE": 1, "T2D": 2, "NA": 99}


LIGAND_TARGET_ALIASES = {
    "ligand": ["ligand"],
    "ligand_activity": ["activity", "ligand_activity"],
    "contrast": ["contrast", "comparison"],
    "predicted_target_gene": ["target", "predicted_target_gene", "target_gene", "gene"],
    "ligand_target_weight": [
        "ligand_target_weight",
        "target_weight",
        "weight",
        "regulatory_potential",
    ],
    "receiver": ["receiver", "receiver_cell_type", "receiver_cell"],
    "direction_regulation": ["direction_regulation", "direction", "regulation"],
    "ligand_activity_scaled": [
        "activity_scaled",
        "ligand_activity_scaled",
        "scaled_activity",
        "max_scaled_activity",
    ],
    "target_gene_log2FC": [
        "target_gene_log2FC",
        "target_gene_logFC",
        "target_log2FC",
        "target_logFC",
        "rna_lfc",
        "rna_logfc",
        "rna_log2fc",
        "log2FC",
        "logFC",
        "avg_log2FC",
    ],
    "target_gene_padj": [
        "target_gene_padj",
        "target_gene_p_adj",
        "target_gene_adjusted_p_value",
        "target_padj",
        "target_p_adj",
        "rna_padj",
        "rna_p_adj",
        "rna_adjusted_p_value",
        "padj",
        "p_adj",
        "adj_p_val",
        "p_val_adj",
    ],
}

PRIORITIZATION_ALIASES = {
    "contrast": ["contrast", "comparison"],
    "receiver": ["receiver", "receiver_cell_type", "receiver_cell"],
    "sender_cell_type": ["sender", "sender_cell_type", "sender_cell"],
    "ligand": ["ligand"],
    "receptor": ["receptor"],
    "interaction": ["lr_interaction", "interaction", "ligand_receptor", "lr_pair"],
    "prioritization_score": [
        "prioritization_score",
        "priority_score",
        "score",
        "scaled_prioritization_score",
    ],
}


def clean_text(value: Any) -> str:
    """Return stripped text; missing-like values become an empty string."""
    if value is None:
        return ""
    text = str(value).strip()
    if text.lower() in {"", "na", "nan", "none", "null"}:
        return ""
    return text


def output_text(value: Any) -> str:
    text = clean_text(value)
    return text if text else "NA"


def normalize_for_matching(value: Any) -> str:
    text = clean_text(value).lower()
    text = text.replace("β", "beta").replace("α", "alpha")
    return re.sub(r"[^a-z0-9]+", "_", text).strip("_")


def normalize_receiver(value: Any) -> str | None:
    """Map receiver spellings to unified paper-style names."""
    text = normalize_for_matching(value)
    parts = [part for part in text.split("_") if part]

    if "beta" in parts or text in {"betacell", "beta_cell"}:
        return "beta cell"
    if "alpha" in parts or text in {"alphacell", "alpha_cell"}:
        return "alpha cell"
    if "acinar" in parts or text in {"acinarcell", "acinar_cell"}:
        return "acinar cell"
    return None


def tokens_from_text(value: Any) -> set[str]:
    text = normalize_for_matching(value)
    return {token.upper() for token in text.split("_") if token}


def normalize_comparison(*values: Any) -> str | None:
    """Standardize comparison names from path and/or contrast text."""
    tokens: set[str] = set()
    for value in values:
        tokens.update(tokens_from_text(value))

    if {"CON", "PRE"}.issubset(tokens):
        return "CON vs PRE"
    if {"PRE", "T2D"}.issubset(tokens):
        return "PRE vs T2D"
    return None


def normalize_direction(value: Any) -> str:
    text = normalize_for_matching(value)
    if text in {"up", "upregulated", "up_regulated"}:
        return "up"
    if text in {"down", "downregulated", "down_regulated"}:
        return "down"
    return clean_text(value)


def normalize_contrast(value: Any) -> str:
    """Return contrast in TOKEN-TOKEN style, e.g. CON-PRE or PRE-T2D."""
    text = normalize_for_matching(value)
    tokens = [token.upper() for token in text.split("_") if token]
    if len(tokens) >= 2:
        return f"{tokens[0]}-{tokens[1]}"
    return clean_text(value).upper()


def infer_target_direction_group(contrast: Any, direction: Any) -> str:
    """Infer which disease-stage group the target gene program supports."""
    contrast_norm = normalize_contrast(contrast)
    direction_norm = normalize_direction(direction)

    mapping = {
        ("CON-PRE", "up"): "CON",
        ("CON-PRE", "down"): "PRE",
        ("PRE-CON", "up"): "PRE",
        ("PRE-CON", "down"): "CON",
        ("PRE-T2D", "up"): "PRE",
        ("PRE-T2D", "down"): "T2D",
        ("T2D-PRE", "up"): "T2D",
        ("T2D-PRE", "down"): "PRE",
    }
    return mapping.get((contrast_norm, direction_norm), "NA")


def header_lookup(fieldnames: list[str]) -> dict[str, str]:
    return {normalize_for_matching(name): name for name in fieldnames}


def find_column(fieldnames: list[str], aliases: list[str]) -> str | None:
    lookup = header_lookup(fieldnames)
    for alias in aliases:
        key = normalize_for_matching(alias)
        if key in lookup:
            return lookup[key]
    return None


def resolve_columns(fieldnames: list[str], aliases: dict[str, list[str]]) -> dict[str, str | None]:
    """Resolve logical column names once per file for speed and robustness."""
    return {
        logical_name: find_column(fieldnames, candidate_names)
        for logical_name, candidate_names in aliases.items()
    }


def get_value(row: dict[str, Any], columns: dict[str, str | None], logical_col: str) -> str:
    source_col = columns.get(logical_col)
    if source_col is None:
        return ""
    return clean_text(row.get(source_col, ""))


def parse_float(value: Any) -> float | None:
    text = clean_text(value)
    if not text:
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    if math.isnan(number) or math.isinf(number):
        return None
    return number


def numeric_sort_value(value: Any, missing: float = -math.inf) -> float:
    number = parse_float(value)
    return number if number is not None else missing


def abs_numeric_sort_value(value: Any, missing: float = -math.inf) -> float:
    number = parse_float(value)
    return abs(number) if number is not None else missing


def path_context(path: Path) -> str:
    try:
        return str(path.relative_to(BASE_DIR))
    except ValueError:
        return str(path)


def collect_ligand_target_files() -> list[Path]:
    return sorted(BASE_DIR.rglob(LIGAND_TARGET_FILENAME))


def source_key(receiver: str, comparison: str, ligand: str) -> tuple[str, str, str]:
    return (receiver, comparison, ligand.upper())


def read_prioritization_support(multinichenet_dir: Path) -> dict[tuple[str, str, str], dict[str, str]]:
    """
    For each receiver/comparison/ligand, keep the sender-receptor interaction
    with the highest prioritization_score from the same MultiNicheNet directory.
    """
    path = multinichenet_dir / PRIORITIZATION_FILENAME
    if not path.exists():
        return {}

    path_receiver = normalize_receiver(path_context(path))
    path_comparison = normalize_comparison(path_context(path))
    support: dict[tuple[str, str, str], dict[str, str]] = {}

    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        columns = resolve_columns(fieldnames, PRIORITIZATION_ALIASES)
        contrast_col = columns.get("contrast")
        receiver_col = columns.get("receiver")

        for row in reader:
            ligand = get_value(row, columns, "ligand")
            if not ligand:
                continue

            row_receiver = normalize_receiver(row.get(receiver_col, "") if receiver_col else "")
            receiver = path_receiver or row_receiver
            if receiver is None:
                continue
            if path_receiver is not None and row_receiver is not None and row_receiver != path_receiver:
                continue

            contrast = row.get(contrast_col, "") if contrast_col else ""
            comparison = path_comparison or normalize_comparison(contrast)
            if comparison is None:
                continue

            receptor = get_value(row, columns, "receptor")
            interaction = get_value(row, columns, "interaction")
            if not interaction and ligand and receptor:
                interaction = f"{ligand}_{receptor}"

            support_row = {
                "sender_cell_type": get_value(row, columns, "sender_cell_type"),
                "receptor": receptor,
                "interaction": interaction,
                "prioritization_score": get_value(row, columns, "prioritization_score"),
            }
            key = source_key(receiver, comparison, ligand)

            if key not in support:
                support[key] = support_row
                continue

            old_score = numeric_sort_value(support[key].get("prioritization_score"))
            new_score = numeric_sort_value(support_row.get("prioritization_score"))
            if new_score > old_score:
                support[key] = support_row

    return support


def read_ligand_target_file(
    path: Path,
    support: dict[tuple[str, str, str], dict[str, str]],
) -> tuple[list[dict[str, str]], dict[str, Any]]:
    path_text = path_context(path)
    path_receiver = normalize_receiver(path_text)
    path_comparison = normalize_comparison(path_text)

    standardized_rows: list[dict[str, str]] = []
    raw_rows = 0
    kept_target_rows = 0
    skipped_target = 0
    skipped_receiver = 0
    skipped_comparison = 0
    supplemented_rows = 0

    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        columns = resolve_columns(fieldnames, LIGAND_TARGET_ALIASES)
        contrast_col = columns.get("contrast")
        receiver_col = columns.get("receiver")

        for row in reader:
            raw_rows += 1

            target = get_value(row, columns, "predicted_target_gene")
            if not target:
                skipped_target += 1
                continue
            kept_target_rows += 1

            row_receiver = normalize_receiver(row.get(receiver_col, "") if receiver_col else "")
            receiver = path_receiver or row_receiver
            if receiver is None:
                skipped_receiver += 1
                continue

            # Path naming defines the analysis receiver. The row-level receiver
            # in these target files can reflect target metadata, so it is used
            # only as a fallback when the path cannot be parsed.
            contrast = row.get(contrast_col, "") if contrast_col else ""
            comparison = path_comparison or normalize_comparison(contrast)
            if comparison is None:
                skipped_comparison += 1
                continue

            ligand = get_value(row, columns, "ligand")
            direction = normalize_direction(get_value(row, columns, "direction_regulation"))
            target_group = infer_target_direction_group(contrast, direction)

            support_row = support.get(source_key(receiver, comparison, ligand), {})
            if support_row:
                supplemented_rows += 1

            standardized_rows.append(
                {
                    "receiver_cell_type": receiver,
                    "comparison": comparison,
                    "target_direction_group": target_group,
                    "direction_regulation": direction,
                    "ligand": ligand,
                    "predicted_target_gene": target,
                    "ligand_activity": get_value(row, columns, "ligand_activity"),
                    "ligand_activity_scaled": get_value(
                        row, columns, "ligand_activity_scaled"
                    ),
                    "ligand_target_weight": get_value(
                        row, columns, "ligand_target_weight"
                    ),
                    "sender_cell_type": support_row.get("sender_cell_type", ""),
                    "receptor": support_row.get("receptor", ""),
                    "interaction": support_row.get("interaction", ""),
                    "prioritization_score": support_row.get("prioritization_score", ""),
                    "target_gene_log2FC": get_value(row, columns, "target_gene_log2FC"),
                    "target_gene_padj": get_value(row, columns, "target_gene_padj"),
                    "rank_within_ligand": "",
                    "rank_within_receiver_comparison_group": "",
                    "source_file": str(path),
                }
            )

    stats = {
        "path": str(path),
        "raw_rows": raw_rows,
        "after_target_filter": kept_target_rows,
        "kept_rows": len(standardized_rows),
        "skipped_target": skipped_target,
        "skipped_receiver": skipped_receiver,
        "skipped_comparison": skipped_comparison,
        "supplemented_rows": supplemented_rows,
    }
    return standardized_rows, stats


def dedupe_key(row: dict[str, str]) -> tuple[str, str, str, str, str]:
    return (
        row.get("receiver_cell_type", ""),
        row.get("comparison", ""),
        row.get("target_direction_group", ""),
        row.get("ligand", ""),
        row.get("predicted_target_gene", ""),
    )


def dedupe_evidence_key(row: dict[str, str]) -> tuple[float, float, float]:
    return (
        numeric_sort_value(row.get("ligand_target_weight")),
        abs_numeric_sort_value(row.get("ligand_activity_scaled")),
        numeric_sort_value(row.get("ligand_activity")),
    )


def deduplicate_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    best_by_key: dict[tuple[str, str, str, str, str], dict[str, str]] = {}
    for row in rows:
        key = dedupe_key(row)
        if key not in best_by_key or dedupe_evidence_key(row) > dedupe_evidence_key(best_by_key[key]):
            best_by_key[key] = row
    return list(best_by_key.values())


def base_sort_key(row: dict[str, str]) -> tuple[Any, ...]:
    return (
        RECEIVER_ORDER.get(row.get("receiver_cell_type", ""), 99),
        row.get("receiver_cell_type", ""),
        COMPARISON_ORDER.get(row.get("comparison", ""), 99),
        row.get("comparison", ""),
        GROUP_ORDER.get(row.get("target_direction_group", ""), 99),
        row.get("target_direction_group", ""),
        -abs_numeric_sort_value(row.get("ligand_activity_scaled")),
        -numeric_sort_value(row.get("ligand_target_weight")),
        row.get("ligand", ""),
        row.get("predicted_target_gene", ""),
    )


def apply_ligand_filter(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Keep top ligands by absolute activity per receiver/comparison/target group."""
    if TOP_LIGANDS_PER_GROUP is None:
        return rows

    rows_by_group: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        rows_by_group[
            (
                row.get("receiver_cell_type", ""),
                row.get("comparison", ""),
                row.get("target_direction_group", ""),
            )
        ].append(row)

    kept_rows: list[dict[str, str]] = []
    for group_key, group_rows in rows_by_group.items():
        best_by_ligand: dict[str, tuple[float, float]] = {}
        for row in group_rows:
            ligand = row.get("ligand", "")
            score = (
                abs_numeric_sort_value(row.get("ligand_activity_scaled")),
                numeric_sort_value(row.get("ligand_target_weight")),
            )
            if ligand not in best_by_ligand or score > best_by_ligand[ligand]:
                best_by_ligand[ligand] = score

        top_ligands = {
            ligand
            for ligand, _score in sorted(
                best_by_ligand.items(),
                key=lambda item: (-item[1][0], -item[1][1], item[0]),
            )[:TOP_LIGANDS_PER_GROUP]
        }
        kept_rows.extend(row for row in group_rows if row.get("ligand", "") in top_ligands)

    return kept_rows


def rank_and_filter_targets(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Rank targets within each ligand and optionally keep top targets."""
    rows_by_ligand: dict[tuple[str, str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        rows_by_ligand[
            (
                row.get("receiver_cell_type", ""),
                row.get("comparison", ""),
                row.get("target_direction_group", ""),
                row.get("ligand", ""),
            )
        ].append(row)

    kept_rows: list[dict[str, str]] = []
    for ligand_key in sorted(
        rows_by_ligand,
        key=lambda k: (
            RECEIVER_ORDER.get(k[0], 99),
            COMPARISON_ORDER.get(k[1], 99),
            GROUP_ORDER.get(k[2], 99),
            k,
        ),
    ):
        ligand_rows = sorted(
            rows_by_ligand[ligand_key],
            key=lambda row: (
                -numeric_sort_value(row.get("ligand_target_weight")),
                -abs_numeric_sort_value(row.get("ligand_activity_scaled")),
                -numeric_sort_value(row.get("ligand_activity")),
                row.get("predicted_target_gene", ""),
            ),
        )
        for rank, row in enumerate(ligand_rows, start=1):
            row["rank_within_ligand"] = str(rank)
        if TOP_TARGETS_PER_LIGAND is not None:
            ligand_rows = ligand_rows[:TOP_TARGETS_PER_LIGAND]
        kept_rows.extend(ligand_rows)

    return kept_rows


def add_group_ranks(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rows_by_group: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        rows_by_group[
            (
                row.get("receiver_cell_type", ""),
                row.get("comparison", ""),
                row.get("target_direction_group", ""),
            )
        ].append(row)

    ranked_rows: list[dict[str, str]] = []
    for group_key in sorted(
        rows_by_group,
        key=lambda k: (
            RECEIVER_ORDER.get(k[0], 99),
            COMPARISON_ORDER.get(k[1], 99),
            GROUP_ORDER.get(k[2], 99),
            k,
        ),
    ):
        group_rows = sorted(
            rows_by_group[group_key],
            key=lambda row: (
                -abs_numeric_sort_value(row.get("ligand_activity_scaled")),
                -numeric_sort_value(row.get("ligand_target_weight")),
                row.get("ligand", ""),
                row.get("predicted_target_gene", ""),
            ),
        )
        for rank, row in enumerate(group_rows, start=1):
            row["rank_within_receiver_comparison_group"] = str(rank)
        ranked_rows.extend(group_rows)

    return sorted(ranked_rows, key=base_sort_key)


def finalize_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rows = apply_ligand_filter(rows)
    rows = rank_and_filter_targets(rows)
    rows = add_group_ranks(rows)
    return rows


def write_csv(rows: list[dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({col: output_text(row.get(col, "")) for col in OUTPUT_COLUMNS})


def excel_column_name(index: int) -> str:
    letters = []
    while index:
        index, remainder = divmod(index - 1, 26)
        letters.append(chr(65 + remainder))
    return "".join(reversed(letters))


def xml_text(value: Any) -> str:
    return escape(str(value), {'"': "&quot;"})


def cell_xml(row_idx: int, col_idx: int, value: Any, style_id: int = 0) -> str:
    ref = f"{excel_column_name(col_idx)}{row_idx}"
    text = clean_text(value)
    style_attr = f' s="{style_id}"' if style_id else ""

    if not text:
        return f'<c r="{ref}"{style_attr}/>'

    number = parse_float(text)
    if number is not None:
        return f'<c r="{ref}"{style_attr}><v>{number:.15g}</v></c>'

    return f'<c r="{ref}" t="inlineStr"{style_attr}><is><t>{xml_text(text)}</t></is></c>'


def compute_widths(rows: list[dict[str, str]]) -> dict[str, float]:
    widths = {col: float(len(col) + 2) for col in OUTPUT_COLUMNS}
    for row in rows:
        for col in OUTPUT_COLUMNS:
            text = output_text(row.get(col, ""))
            widths[col] = max(widths[col], float(len(text) + 2))
    return widths


def worksheet_xml(rows: list[dict[str, str]], widths: dict[str, float]) -> str:
    n_rows = len(rows) + 1
    n_cols = len(OUTPUT_COLUMNS)
    last_col = excel_column_name(n_cols)

    cols_xml = ["<cols>"]
    for idx, col in enumerate(OUTPUT_COLUMNS, start=1):
        width = min(max(widths.get(col, 12.0), 8.0), 55.0)
        cols_xml.append(
            f'<col min="{idx}" max="{idx}" width="{width:.1f}" customWidth="1"/>'
        )
    cols_xml.append("</cols>")

    sheet_rows = []
    header_cells = [
        cell_xml(1, col_idx, col, style_id=1)
        for col_idx, col in enumerate(OUTPUT_COLUMNS, start=1)
    ]
    sheet_rows.append(f'<row r="1">{"".join(header_cells)}</row>')

    for row_idx, row in enumerate(rows, start=2):
        cells = []
        for col_idx, col in enumerate(OUTPUT_COLUMNS, start=1):
            value = output_text(row.get(col, ""))
            if col in P_VALUE_COLUMNS:
                style_id = 3
            elif col in NUMERIC_COLUMNS:
                style_id = 2
            else:
                style_id = 0
            cells.append(cell_xml(row_idx, col_idx, value, style_id=style_id))
        sheet_rows.append(f'<row r="{row_idx}">{"".join(cells)}</row>')

    auto_filter_ref = f"A1:{last_col}{max(n_rows, 1)}"
    return f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"
 xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
<sheetViews>
  <sheetView workbookViewId="0">
    <pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/>
    <selection pane="bottomLeft" activeCell="A2" sqref="A2"/>
  </sheetView>
</sheetViews>
{"".join(cols_xml)}
<sheetData>
{"".join(sheet_rows)}
</sheetData>
<autoFilter ref="{auto_filter_ref}"/>
<pageMargins left="0.7" right="0.7" top="0.75" bottom="0.75" header="0.3" footer="0.3"/>
</worksheet>'''


def write_xlsx(rows: list[dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    widths = compute_widths(rows)
    now = datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")

    content_types = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
  <Default Extension="xml" ContentType="application/xml"/>
  <Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>
  <Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>
  <Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>
  <Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>
  <Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>
</Types>'''

    root_rels = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>
  <Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/>
  <Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" Target="docProps/app.xml"/>
</Relationships>'''

    workbook_xml = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"
 xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
  <sheets>
    <sheet name="Table S11" sheetId="1" r:id="rId1"/>
  </sheets>
</workbook>'''

    workbook_rels = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/>
  <Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>
</Relationships>'''

    styles_xml = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
  <numFmts count="2">
    <numFmt numFmtId="164" formatCode="0.0000"/>
    <numFmt numFmtId="165" formatCode="0.00E+00"/>
  </numFmts>
  <fonts count="2">
    <font><sz val="11"/><name val="Calibri"/></font>
    <font><b/><sz val="11"/><name val="Calibri"/></font>
  </fonts>
  <fills count="2">
    <fill><patternFill patternType="none"/></fill>
    <fill><patternFill patternType="gray125"/></fill>
  </fills>
  <borders count="1">
    <border><left/><right/><top/><bottom/><diagonal/></border>
  </borders>
  <cellStyleXfs count="1">
    <xf numFmtId="0" fontId="0" fillId="0" borderId="0"/>
  </cellStyleXfs>
  <cellXfs count="4">
    <xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>
    <xf numFmtId="0" fontId="1" fillId="0" borderId="0" xfId="0" applyFont="1"/>
    <xf numFmtId="164" fontId="0" fillId="0" borderId="0" xfId="0" applyNumberFormat="1"/>
    <xf numFmtId="165" fontId="0" fillId="0" borderId="0" xfId="0" applyNumberFormat="1"/>
  </cellXfs>
  <cellStyles count="1">
    <cellStyle name="Normal" xfId="0" builtinId="0"/>
  </cellStyles>
</styleSheet>'''

    core_xml = f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties"
 xmlns:dc="http://purl.org/dc/elements/1.1/"
 xmlns:dcterms="http://purl.org/dc/terms/"
 xmlns:dcmitype="http://purl.org/dc/dcmitype/"
 xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <dc:creator>make_supp_table11_multinichenet_ligand_target_links.py</dc:creator>
  <cp:lastModifiedBy>make_supp_table11_multinichenet_ligand_target_links.py</cp:lastModifiedBy>
  <dcterms:created xsi:type="dcterms:W3CDTF">{now}</dcterms:created>
  <dcterms:modified xsi:type="dcterms:W3CDTF">{now}</dcterms:modified>
</cp:coreProperties>'''

    app_xml = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties"
 xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes">
  <Application>Python standard library</Application>
</Properties>'''

    files = {
        "[Content_Types].xml": content_types,
        "_rels/.rels": root_rels,
        "xl/workbook.xml": workbook_xml,
        "xl/_rels/workbook.xml.rels": workbook_rels,
        "xl/styles.xml": styles_xml,
        "xl/worksheets/sheet1.xml": worksheet_xml(rows, widths),
        "docProps/core.xml": core_xml,
        "docProps/app.xml": app_xml,
    }

    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as xlsx:
        for name, content in files.items():
            xlsx.writestr(posixpath.normpath(name), content)


def print_summary(
    source_files: list[Path],
    file_stats: list[dict[str, Any]],
    merged_rows: int,
    deduped_rows: int,
    final_rows: list[dict[str, str]],
) -> None:
    print("\nSummary")
    print("=======")
    print(f"Found ligand-target files: {len(source_files)}")
    print("\nPer-file rows:")
    for stats in file_stats:
        print(
            f"- {stats['path']}: raw {stats['raw_rows']}, "
            f"after target filter {stats['after_target_filter']}, "
            f"kept {stats['kept_rows']}, "
            f"skipped_target {stats['skipped_target']}, "
            f"skipped_receiver {stats['skipped_receiver']}, "
            f"skipped_comparison {stats['skipped_comparison']}, "
            f"supplemented {stats['supplemented_rows']}"
        )

    print(f"\nMerged rows after target/comparison filters: {merged_rows}")
    print(f"Rows after de-duplication: {deduped_rows}")
    print(f"Final output rows: {len(final_rows)}")

    supplemented_final = sum(
        1
        for row in final_rows
        if clean_text(row.get("sender_cell_type"))
        or clean_text(row.get("receptor"))
        or clean_text(row.get("interaction"))
    )
    print(f"Rows with supplemented sender/receptor information: {supplemented_final}")

    print("\nRows per receiver/comparison/target_direction_group:")
    counts = Counter(
        (
            row.get("receiver_cell_type", ""),
            row.get("comparison", ""),
            row.get("target_direction_group", ""),
        )
        for row in final_rows
    )
    for key, count in sorted(
        counts.items(),
        key=lambda item: (
            RECEIVER_ORDER.get(item[0][0], 99),
            COMPARISON_ORDER.get(item[0][1], 99),
            GROUP_ORDER.get(item[0][2], 99),
            item[0],
        ),
    ):
        receiver, comparison, target_group = key
        print(f"- {receiver} | {comparison} | {target_group}: {count}")

    print(f"\nWrote CSV:  {OUT_CSV}")
    print(f"Wrote XLSX: {OUT_XLSX}")


def main() -> None:
    source_files = collect_ligand_target_files()
    all_rows: list[dict[str, str]] = []
    file_stats: list[dict[str, Any]] = []

    # Prioritization support is read once per MultiNicheNet directory.
    support_cache: dict[Path, dict[tuple[str, str, str], dict[str, str]]] = {}

    for path in source_files:
        multinichenet_dir = path.parent
        if multinichenet_dir not in support_cache:
            support_cache[multinichenet_dir] = read_prioritization_support(multinichenet_dir)

        rows, stats = read_ligand_target_file(path, support_cache[multinichenet_dir])
        all_rows.extend(rows)
        file_stats.append(stats)

    deduped = deduplicate_rows(all_rows)
    final_rows = finalize_rows(deduped)

    write_csv(final_rows, OUT_CSV)
    write_xlsx(final_rows, OUT_XLSX)

    print_summary(
        source_files=source_files,
        file_stats=file_stats,
        merged_rows=len(all_rows),
        deduped_rows=len(deduped),
        final_rows=final_rows,
    )


if __name__ == "__main__":
    main()
