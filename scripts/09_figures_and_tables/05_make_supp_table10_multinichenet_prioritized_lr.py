#!/usr/bin/env python3

"""
Script name: 05_make_supp_table10_multinichenet_prioritized_lr.py
Purpose: Cleaned public workflow script for the T2D single-cell multiome analysis.
Input: Paths and analysis settings are read from configs/config.yaml.
Output: Module-specific outputs are written under the configured results directory.
Main steps: Load configuration, run the original analysis logic, and export reproducible outputs.
Example command: python scripts/09_figures_and_tables/05_make_supp_table10_multinichenet_prioritized_lr.py
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
Make Supplementary Table 10:
MultiNicheNet prioritized ligand-receptor interactions.

This script intentionally uses only the Python standard library. It can write
both CSV and a lightly formatted XLSX file without requiring pandas/openpyxl.

Run:
    python3 scripts/09_figures_and_tables/make_supp_table10_multinichenet_prioritized_lr.py

Outputs:
    results/downstream/supplementary_tables/Supplementary_Table_10_MultiNicheNet_prioritized_LR.csv
    results/downstream/supplementary_tables/Supplementary_Table_10_MultiNicheNet_prioritized_LR.xlsx
"""

import csv
import math
import os
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
SOURCE_FILENAME = "prioritization_group_prioritization_table_source.csv"

OUT_DIR = BASE_DIR / "supplementary_tables"
OUT_CSV = OUT_DIR / "Supplementary_Table_10_MultiNicheNet_prioritized_LR.csv"
OUT_XLSX = OUT_DIR / "Supplementary_Table_10_MultiNicheNet_prioritized_LR.xlsx"

# Set to None to output all rows after de-duplication.
TOP_N_PER_GROUP: int | None = 50


OUTPUT_COLUMNS = [
    "receiver_cell_type",
    "comparison",
    "top_group",
    "group",
    "sender_cell_type",
    "ligand",
    "receptor",
    "interaction",
    "ligand_logFC",
    "receptor_logFC",
    "ligand_receptor_logFC_avg",
    "ligand_p_adj",
    "receptor_p_adj",
    "ligand_activity",
    "direction_regulation",
    "ligand_activity_scaled",
    "activity_up",
    "activity_scaled_up",
    "activity_down",
    "activity_scaled_down",
    "fraction_ligand_group",
    "fraction_receptor_group",
    "pb_ligand_group",
    "pb_receptor_group",
    "prioritization_score",
    "rank_within_receiver_comparison_group",
]

TEXT_COLUMNS = {
    "receiver_cell_type",
    "comparison",
    "top_group",
    "group",
    "sender_cell_type",
    "ligand",
    "receptor",
    "interaction",
    "direction_regulation",
}

P_VALUE_COLUMNS = {"ligand_p_adj", "receptor_p_adj"}

NUMERIC_COLUMNS = set(OUTPUT_COLUMNS) - TEXT_COLUMNS

RECEIVER_ORDER = {"beta cell": 0, "alpha cell": 1, "acinar cell": 2}
COMPARISON_ORDER = {"CON vs PRE": 0, "PRE vs T2D": 1}
GROUP_ORDER = {"CON": 0, "PRE": 1, "T2D": 2}


# Candidate input columns. The first existing column in each list is used.
# The later scaled-column fallbacks are included for MultiNicheNet exports that
# already store prioritization components with scaled names.
COLUMN_ALIASES = {
    "contrast": ["contrast", "comparison"],
    "group": ["group"],
    "sender_cell_type": ["sender", "sender_cell_type", "sender_cell"],
    "receiver": ["receiver", "receiver_cell_type", "receiver_cell"],
    "ligand": ["ligand"],
    "receptor": ["receptor"],
    "interaction": ["lr_interaction", "interaction", "ligand_receptor", "lr_pair"],
    "ligand_logFC": [
        "lfc_ligand",
        "logfc_ligand",
        "ligand_logfc",
        "ligand_logFC",
        "scaled_lfc_ligand",
    ],
    "receptor_logFC": [
        "lfc_receptor",
        "logfc_receptor",
        "receptor_logfc",
        "receptor_logFC",
        "scaled_lfc_receptor",
    ],
    "ligand_receptor_logFC_avg": [
        "ligand_receptor_lfc_avg",
        "lr_lfc_avg",
        "lfc_ligand_receptor_avg",
        "scaled_ligand_receptor_lfc_avg",
    ],
    "ligand_p_adj": [
        "p_adj_ligand",
        "padj_ligand",
        "ligand_p_adj",
        "ligand_adj_p",
        "p_val_ligand_adapted",
        "scaled_p_val_ligand_adapted",
        "p_val_ligand",
    ],
    "receptor_p_adj": [
        "p_adj_receptor",
        "padj_receptor",
        "receptor_p_adj",
        "receptor_adj_p",
        "p_val_receptor_adapted",
        "scaled_p_val_receptor_adapted",
        "p_val_receptor",
    ],
    "ligand_activity": ["activity", "ligand_activity", "max_activity"],
    "direction_regulation": ["direction_regulation", "direction", "regulation"],
    "ligand_activity_scaled": [
        "activity_scaled",
        "scaled_activity",
        "ligand_activity_scaled",
        "max_scaled_activity",
    ],
    "activity_up": ["activity_up"],
    "activity_scaled_up": ["activity_scaled_up", "scaled_activity_up"],
    "activity_down": ["activity_down"],
    "activity_scaled_down": ["activity_scaled_down", "scaled_activity_down"],
    "fraction_ligand_group": [
        "fraction_ligand_group",
        "fraction_ligand",
        "fraction_expressing_ligand",
        "fraction_expressing_ligand_receptor",
    ],
    "fraction_receptor_group": [
        "fraction_receptor_group",
        "fraction_receptor",
        "fraction_expressing_receptor",
        "fraction_expressing_ligand_receptor",
    ],
    "pb_ligand_group": ["pb_ligand_group", "pb_ligand", "scaled_pb_ligand"],
    "pb_receptor_group": ["pb_receptor_group", "pb_receptor", "scaled_pb_receptor"],
    "prioritization_score": [
        "prioritization_score",
        "priority_score",
        "score",
        "scaled_prioritization_score",
    ],
    "top_group": ["top_group"],
}


def clean_text(value: Any) -> str:
    """Return stripped text; missing-like values become an empty string."""
    if value is None:
        return ""
    text = str(value).strip()
    if text.lower() in {"", "na", "nan", "none", "null"}:
        return ""
    return text


def normalize_for_matching(value: Any) -> str:
    """Normalize strings for loose matching of headers, paths, and cell names."""
    text = clean_text(value).lower()
    text = text.replace("β", "beta").replace("α", "alpha")
    return re.sub(r"[^a-z0-9]+", "_", text).strip("_")


def normalize_receiver(value: Any) -> str | None:
    """Map receiver spellings to paper-style cell type names."""
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
    """Standardize comparison names using path and/or contrast text."""
    tokens: set[str] = set()
    for value in values:
        tokens.update(tokens_from_text(value))

    if {"CON", "PRE"}.issubset(tokens):
        return "CON vs PRE"
    if {"PRE", "T2D"}.issubset(tokens):
        return "PRE vs T2D"
    return None


def normalize_group(value: Any) -> str:
    text = clean_text(value)
    upper = normalize_for_matching(text).upper()
    if upper in {"CON", "PRE", "T2D"}:
        return upper
    return text


def header_lookup(fieldnames: list[str]) -> dict[str, str]:
    """Map normalized header names to their original spelling."""
    return {normalize_for_matching(name): name for name in fieldnames}


def find_column(fieldnames: list[str], aliases: list[str]) -> str | None:
    lookup = header_lookup(fieldnames)
    for alias in aliases:
        key = normalize_for_matching(alias)
        if key in lookup:
            return lookup[key]
    return None


def resolve_columns(fieldnames: list[str]) -> dict[str, str | None]:
    """Resolve all potentially used input columns once per file."""
    return {
        logical_name: find_column(fieldnames, aliases)
        for logical_name, aliases in COLUMN_ALIASES.items()
    }


def get_value(row: dict[str, Any], columns: dict[str, str | None], output_col: str) -> str:
    col = columns.get(output_col)
    if col is None:
        return ""
    return clean_text(row.get(col, ""))


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
    """Use the relative path as extra context for receiver/comparison detection."""
    try:
        return str(path.relative_to(BASE_DIR))
    except ValueError:
        return str(path)


def collect_source_files() -> list[Path]:
    return sorted(BASE_DIR.rglob(SOURCE_FILENAME))


def read_source_file(path: Path) -> tuple[list[dict[str, str]], dict[str, Any]]:
    """Read and standardize one source CSV."""
    path_text = path_context(path)
    receiver_from_path = normalize_receiver(path_text)
    comparison_from_path = normalize_comparison(path_text)

    standardized_rows: list[dict[str, str]] = []
    skipped_receiver = 0
    skipped_comparison = 0

    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        columns = resolve_columns(fieldnames)
        contrast_col = columns.get("contrast")
        receiver_col = columns.get("receiver")

        raw_rows = 0
        for row in reader:
            raw_rows += 1

            row_receiver_raw = row.get(receiver_col, "") if receiver_col else ""
            row_receiver = normalize_receiver(row_receiver_raw)

            # The directory is expected to define the target receiver. The row
            # value is used as a fallback and also protects against mixed files.
            receiver = receiver_from_path or row_receiver
            if receiver is None:
                skipped_receiver += 1
                continue
            if row_receiver is not None and receiver_from_path is not None:
                if row_receiver != receiver_from_path:
                    skipped_receiver += 1
                    continue

            contrast_raw = row.get(contrast_col, "") if contrast_col else ""
            comparison = comparison_from_path or normalize_comparison(contrast_raw)
            if comparison is None:
                skipped_comparison += 1
                continue

            ligand = get_value(row, columns, "ligand")
            receptor = get_value(row, columns, "receptor")
            interaction = get_value(row, columns, "interaction")
            if not interaction and ligand and receptor:
                interaction = f"{ligand}_{receptor}"

            out = {
                "receiver_cell_type": receiver,
                "comparison": comparison,
                "top_group": normalize_group(get_value(row, columns, "top_group")),
                "group": normalize_group(get_value(row, columns, "group")),
                "sender_cell_type": get_value(row, columns, "sender_cell_type"),
                "ligand": ligand,
                "receptor": receptor,
                "interaction": interaction,
                "ligand_logFC": get_value(row, columns, "ligand_logFC"),
                "receptor_logFC": get_value(row, columns, "receptor_logFC"),
                "ligand_receptor_logFC_avg": get_value(
                    row, columns, "ligand_receptor_logFC_avg"
                ),
                "ligand_p_adj": get_value(row, columns, "ligand_p_adj"),
                "receptor_p_adj": get_value(row, columns, "receptor_p_adj"),
                "ligand_activity": get_value(row, columns, "ligand_activity"),
                "direction_regulation": get_value(
                    row, columns, "direction_regulation"
                ),
                "ligand_activity_scaled": get_value(
                    row, columns, "ligand_activity_scaled"
                ),
                "activity_up": get_value(row, columns, "activity_up"),
                "activity_scaled_up": get_value(row, columns, "activity_scaled_up"),
                "activity_down": get_value(row, columns, "activity_down"),
                "activity_scaled_down": get_value(
                    row, columns, "activity_scaled_down"
                ),
                "fraction_ligand_group": get_value(
                    row, columns, "fraction_ligand_group"
                ),
                "fraction_receptor_group": get_value(
                    row, columns, "fraction_receptor_group"
                ),
                "pb_ligand_group": get_value(row, columns, "pb_ligand_group"),
                "pb_receptor_group": get_value(row, columns, "pb_receptor_group"),
                "prioritization_score": get_value(
                    row, columns, "prioritization_score"
                ),
                "rank_within_receiver_comparison_group": "",
            }
            standardized_rows.append(out)

    stats = {
        "path": str(path),
        "raw_rows": raw_rows,
        "kept_rows": len(standardized_rows),
        "skipped_receiver": skipped_receiver,
        "skipped_comparison": skipped_comparison,
    }
    return standardized_rows, stats


def evidence_key(row: dict[str, str]) -> tuple[float, float, float]:
    """Sort key for keeping the strongest duplicate evidence."""
    return (
        numeric_sort_value(row.get("prioritization_score")),
        abs_numeric_sort_value(row.get("ligand_activity_scaled")),
        abs_numeric_sort_value(row.get("ligand_activity")),
    )


def deduplicate_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    dedupe_key_cols = [
        "receiver_cell_type",
        "comparison",
        "top_group",
        "group",
        "sender_cell_type",
        "ligand",
        "receptor",
        "interaction",
    ]

    best_by_key: dict[tuple[str, ...], dict[str, str]] = {}
    for row in rows:
        key = tuple(row.get(col, "") for col in dedupe_key_cols)
        if key not in best_by_key or evidence_key(row) > evidence_key(best_by_key[key]):
            best_by_key[key] = row

    return list(best_by_key.values())


def final_sort_key(row: dict[str, str]) -> tuple[Any, ...]:
    return (
        RECEIVER_ORDER.get(row.get("receiver_cell_type", ""), 99),
        row.get("receiver_cell_type", ""),
        COMPARISON_ORDER.get(row.get("comparison", ""), 99),
        row.get("comparison", ""),
        GROUP_ORDER.get(row.get("top_group", ""), 99),
        row.get("top_group", ""),
        -numeric_sort_value(row.get("prioritization_score")),
        -abs_numeric_sort_value(row.get("ligand_activity_scaled")),
        -abs_numeric_sort_value(row.get("ligand_activity")),
        row.get("sender_cell_type", ""),
        row.get("ligand", ""),
        row.get("receptor", ""),
    )


def add_ranks_and_limit(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Rank within receiver/comparison/top_group and optionally retain top N."""
    grouped: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        key = (
            row.get("receiver_cell_type", ""),
            row.get("comparison", ""),
            row.get("top_group", ""),
        )
        grouped[key].append(row)

    ranked_rows: list[dict[str, str]] = []
    for key in sorted(
        grouped,
        key=lambda k: (
            RECEIVER_ORDER.get(k[0], 99),
            COMPARISON_ORDER.get(k[1], 99),
            GROUP_ORDER.get(k[2], 99),
            k,
        ),
    ):
        group_rows = sorted(grouped[key], key=final_sort_key)
        for rank, row in enumerate(group_rows, start=1):
            row["rank_within_receiver_comparison_group"] = str(rank)
        if TOP_N_PER_GROUP is not None:
            group_rows = group_rows[:TOP_N_PER_GROUP]
        ranked_rows.extend(group_rows)

    return sorted(ranked_rows, key=final_sort_key)


def csv_value(row: dict[str, str], col: str) -> str:
    value = clean_text(row.get(col, ""))
    return value if value else "NA"


def write_csv(rows: list[dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({col: csv_value(row, col) for col in OUTPUT_COLUMNS})


def excel_column_name(index: int) -> str:
    """Convert a 1-based column index to Excel letters."""
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


def worksheet_xml(rows: list[dict[str, str]], widths: dict[str, float]) -> str:
    n_rows = len(rows) + 1
    n_cols = len(OUTPUT_COLUMNS)
    last_col = excel_column_name(n_cols)

    cols_xml = ["<cols>"]
    for idx, col in enumerate(OUTPUT_COLUMNS, start=1):
        width = min(max(widths.get(col, 12.0), 8.0), 45.0)
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
            value = csv_value(row, col)
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


def compute_widths(rows: list[dict[str, str]]) -> dict[str, float]:
    widths = {col: float(len(col) + 2) for col in OUTPUT_COLUMNS}
    for row in rows:
        for col in OUTPUT_COLUMNS:
            text = csv_value(row, col)
            widths[col] = max(widths[col], float(len(text) + 2))
    return widths


def write_xlsx(rows: list[dict[str, str]], path: Path) -> None:
    """Write a simple formatted XLSX workbook using OpenXML parts."""
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
    <sheet name="Supplementary Table 10" sheetId="1" r:id="rId1"/>
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
  <dc:creator>make_supp_table10_multinichenet_prioritized_lr.py</dc:creator>
  <cp:lastModifiedBy>make_supp_table10_multinichenet_prioritized_lr.py</cp:lastModifiedBy>
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
            # OpenXML requires POSIX separators inside the zip archive.
            xlsx.writestr(posixpath.normpath(name), content)


def print_summary(
    source_files: list[Path],
    file_stats: list[dict[str, Any]],
    merged_count: int,
    deduped_count: int,
    final_rows: list[dict[str, str]],
) -> None:
    print("\nSummary")
    print("=======")
    print(f"Found source files: {len(source_files)}")
    print("\nPer-file rows:")
    for stats in file_stats:
        print(
            f"- {stats['path']}: read {stats['raw_rows']} rows, "
            f"kept {stats['kept_rows']}, "
            f"skipped_receiver {stats['skipped_receiver']}, "
            f"skipped_comparison {stats['skipped_comparison']}"
        )

    print(f"\nMerged rows: {merged_count}")
    print(f"Rows after de-duplication: {deduped_count}")
    print(f"Final output rows: {len(final_rows)}")

    print("\nRows per receiver/comparison/top_group:")
    counts = Counter(
        (
            row.get("receiver_cell_type", ""),
            row.get("comparison", ""),
            row.get("top_group", ""),
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
        receiver, comparison, top_group = key
        print(f"- {receiver} | {comparison} | {top_group}: {count}")

    print(f"\nWrote CSV:  {OUT_CSV}")
    print(f"Wrote XLSX: {OUT_XLSX}")


def main() -> None:
    source_files = collect_source_files()
    all_rows: list[dict[str, str]] = []
    file_stats: list[dict[str, Any]] = []

    for path in source_files:
        rows, stats = read_source_file(path)
        all_rows.extend(rows)
        file_stats.append(stats)

    deduped_rows = deduplicate_rows(all_rows)
    final_rows = add_ranks_and_limit(deduped_rows)

    write_csv(final_rows, OUT_CSV)
    write_xlsx(final_rows, OUT_XLSX)

    print_summary(
        source_files=source_files,
        file_stats=file_stats,
        merged_count=len(all_rows),
        deduped_count=len(deduped_rows),
        final_rows=final_rows,
    )


if __name__ == "__main__":
    main()
