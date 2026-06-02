#!/usr/bin/env python3

"""
Script name: 07_make_supp_table12_tradeseq_pseudotime_dynamic_genes.py
Purpose: Cleaned public workflow script for the T2D single-cell multiome analysis.
Input: Paths and analysis settings are read from configs/config.yaml.
Output: Module-specific outputs are written under the configured results directory.
Main steps: Load configuration, run the original analysis logic, and export reproducible outputs.
Example command: python scripts/09_figures_and_tables/07_make_supp_table12_tradeseq_pseudotime_dynamic_genes.py
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
Make Table S12:
Pseudotime-associated dynamic genes identified by tradeSeq.

This script combines the three fixed tradeSeq associationTest significant-gene
tables for beta, alpha, and acinar receiver cell programs.

Run:
    python3 scripts/09_figures_and_tables/make_supp_table12_tradeseq_pseudotime_dynamic_genes.py

Outputs:
    results/supplementary_tables/Supplementary_Table_12_tradeSeq_pseudotime_dynamic_genes.csv
    results/supplementary_tables/Supplementary_Table_12_tradeSeq_pseudotime_dynamic_genes.xlsx
"""

import csv
import math
import posixpath
import re
import zipfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from xml.sax.saxutils import escape


# ---------------------------------------------------------------------------
# User-editable settings
# ---------------------------------------------------------------------------

OUT_DIR = Path("results/supplementary_tables")
OUT_CSV = OUT_DIR / "Supplementary_Table_12_tradeSeq_pseudotime_dynamic_genes.csv"
OUT_XLSX = OUT_DIR / "Supplementary_Table_12_tradeSeq_pseudotime_dynamic_genes.xlsx"

# The three fixed significant associationTest result files.
INPUT_SPECS = [
    {
        "receiver_cell_type": "beta cell",
        "lineage": "partition 1",
        "source_file": Path(
            "results/downstream/riskcell_monocle/beta/union_native/"
            "partition1_only/tradeSeq/associationTest_sig_q005.csv"
        ),
    },
    {
        "receiver_cell_type": "alpha cell",
        "lineage": "lineage 1",
        "source_file": Path(
            "results/downstream/riskcell_monocle/alpha_union/"
            "tradeSeq_all_lineages/associationTest_lineage1_sig_q005.csv"
        ),
    },
    {
        "receiver_cell_type": "acinar cell",
        "lineage": "pseudotime trajectory",
        "source_file": Path(
            "results/downstream/riskcell_monocle/acinar/tradeSeq/"
            "tradeSeq_associationTest_sig.csv"
        ),
    },
]


OUTPUT_COLUMNS = [
    "receiver_cell_type",
    "lineage",
    "gene",
    "wald_statistic",
    "degrees_of_freedom",
    "p_value",
    "q_value",
    "mean_logFC",
    "detection_fraction",
    "mean_count",
    "rank_within_receiver",
    "source_file",
]

TEXT_COLUMNS = {"receiver_cell_type", "lineage", "gene", "source_file"}
P_VALUE_COLUMNS = {"p_value", "q_value"}
NUMERIC_COLUMNS = set(OUTPUT_COLUMNS) - TEXT_COLUMNS

RECEIVER_ORDER = {"beta cell": 0, "alpha cell": 1, "acinar cell": 2}


COLUMN_ALIASES = {
    "gene": ["gene", "genes", "symbol"],
    "wald_statistic": ["waldStat", "wald_statistic", "wald_stat", "wald"],
    "degrees_of_freedom": ["df", "degrees_of_freedom"],
    "p_value": ["pvalue", "p_value", "p.val", "pval"],
    "q_value": ["qvalue", "qvalue_pvalue", "q_value", "padj", "p_adj", "fdr"],
    "mean_logFC": ["meanLogFC", "mean_logFC", "mean_log2FC", "logFC", "log2FC"],
    "detection_fraction": ["detect_frac", "detection_fraction", "pct_detected"],
    "mean_count": ["mean_count", "mean_counts", "avg_count", "average_count"],
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
    return re.sub(r"[^a-z0-9]+", "_", text).strip("_")


def header_lookup(fieldnames: list[str]) -> dict[str, str]:
    return {normalize_for_matching(name): name for name in fieldnames}


def find_column(fieldnames: list[str], aliases: list[str]) -> str | None:
    lookup = header_lookup(fieldnames)
    for alias in aliases:
        key = normalize_for_matching(alias)
        if key in lookup:
            return lookup[key]
    return None


def resolve_columns(fieldnames: list[str]) -> dict[str, str | None]:
    return {
        logical_name: find_column(fieldnames, aliases)
        for logical_name, aliases in COLUMN_ALIASES.items()
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


def numeric_sort_value(value: Any, missing: float = math.inf) -> float:
    number = parse_float(value)
    return number if number is not None else missing


def read_one_file(spec: dict[str, Any]) -> tuple[list[dict[str, str]], dict[str, Any]]:
    path = spec["source_file"]
    rows: list[dict[str, str]] = []
    raw_rows = 0
    skipped_missing_gene = 0

    if not path.exists():
        return [], {
            "path": str(path),
            "exists": False,
            "raw_rows": 0,
            "kept_rows": 0,
            "skipped_missing_gene": 0,
        }

    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        columns = resolve_columns(fieldnames)

        for row in reader:
            raw_rows += 1
            gene = get_value(row, columns, "gene")
            if not gene:
                skipped_missing_gene += 1
                continue

            rows.append(
                {
                    "receiver_cell_type": spec["receiver_cell_type"],
                    "lineage": spec["lineage"],
                    "gene": gene,
                    "wald_statistic": get_value(row, columns, "wald_statistic"),
                    "degrees_of_freedom": get_value(row, columns, "degrees_of_freedom"),
                    "p_value": get_value(row, columns, "p_value"),
                    "q_value": get_value(row, columns, "q_value"),
                    "mean_logFC": get_value(row, columns, "mean_logFC"),
                    "detection_fraction": get_value(row, columns, "detection_fraction"),
                    "mean_count": get_value(row, columns, "mean_count"),
                    "rank_within_receiver": "",
                    "source_file": str(path),
                }
            )

    return rows, {
        "path": str(path),
        "exists": True,
        "raw_rows": raw_rows,
        "kept_rows": len(rows),
        "skipped_missing_gene": skipped_missing_gene,
    }


def sort_key(row: dict[str, str]) -> tuple[Any, ...]:
    return (
        RECEIVER_ORDER.get(row.get("receiver_cell_type", ""), 99),
        numeric_sort_value(row.get("q_value")),
        numeric_sort_value(row.get("p_value")),
        -numeric_sort_value(row.get("wald_statistic"), missing=-math.inf),
        row.get("gene", ""),
    )


def add_ranks(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault(row["receiver_cell_type"], []).append(row)

    ranked_rows: list[dict[str, str]] = []
    for receiver in sorted(grouped, key=lambda x: (RECEIVER_ORDER.get(x, 99), x)):
        receiver_rows = sorted(grouped[receiver], key=sort_key)
        for rank, row in enumerate(receiver_rows, start=1):
            row["rank_within_receiver"] = str(rank)
        ranked_rows.extend(receiver_rows)

    return sorted(ranked_rows, key=sort_key)


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
            widths[col] = max(widths[col], float(len(output_text(row.get(col, ""))) + 2))
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
    <sheet name="Table S12" sheetId="1" r:id="rId1"/>
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
  <dc:creator>make_supp_table12_tradeseq_pseudotime_dynamic_genes.py</dc:creator>
  <cp:lastModifiedBy>make_supp_table12_tradeseq_pseudotime_dynamic_genes.py</cp:lastModifiedBy>
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


def print_summary(file_stats: list[dict[str, Any]], final_rows: list[dict[str, str]]) -> None:
    print("\nSummary")
    print("=======")
    print(f"Input files configured: {len(INPUT_SPECS)}")

    print("\nPer-file rows:")
    for stats in file_stats:
        print(
            f"- {stats['path']}: exists {stats['exists']}, "
            f"raw {stats['raw_rows']}, kept {stats['kept_rows']}, "
            f"skipped_missing_gene {stats['skipped_missing_gene']}"
        )

    print(f"\nFinal output rows: {len(final_rows)}")
    print("\nRows per receiver:")
    counts = Counter(row["receiver_cell_type"] for row in final_rows)
    for receiver, count in sorted(
        counts.items(), key=lambda item: (RECEIVER_ORDER.get(item[0], 99), item[0])
    ):
        print(f"- {receiver}: {count}")

    print(f"\nWrote CSV:  {OUT_CSV}")
    print(f"Wrote XLSX: {OUT_XLSX}")


def main() -> None:
    all_rows: list[dict[str, str]] = []
    file_stats: list[dict[str, Any]] = []

    for spec in INPUT_SPECS:
        rows, stats = read_one_file(spec)
        all_rows.extend(rows)
        file_stats.append(stats)

    final_rows = add_ranks(all_rows)
    write_csv(final_rows, OUT_CSV)
    write_xlsx(final_rows, OUT_XLSX)
    print_summary(file_stats, final_rows)


if __name__ == "__main__":
    main()
