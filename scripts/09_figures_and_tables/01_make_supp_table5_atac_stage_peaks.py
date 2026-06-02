#!/usr/bin/env python3

"""
Script name: 01_make_supp_table5_atac_stage_peaks.py
Purpose: Cleaned public workflow script for the T2D single-cell multiome analysis.
Input: Paths and analysis settings are read from configs/config.yaml.
Output: Module-specific outputs are written under the configured results directory.
Main steps: Load configuration, run the original analysis logic, and export reproducible outputs.
Example command: python scripts/09_figures_and_tables/01_make_supp_table5_atac_stage_peaks.py
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
Generate Supplementary Table 5. ATAC-derived STAGE peaks.

This table includes all ATAC-derived STAGE peaks from STAGE_peaks, joined with
their five-fold Boruta selection frequency from selection_frequency. If a STAGE
peak is absent from the matching stability table, its selection_frequency is
reported as 0.0.

Output columns:
- cell_type
- comparison
- peak
- selection_frequency
- mean_SHAP

Outputs:
- Supplementary_Table_5_ATAC_derived_STAGE_peaks.csv
- Supplementary_Table_5_ATAC_derived_STAGE_peaks.xlsx
"""

import csv
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from xml.sax.saxutils import escape


ATAC_DIR = Path("results/model_outputs/atac")
STAGE_DIR = ATAC_DIR / "STAGE_peaks"
SELECTION_FREQUENCY_DIR = ATAC_DIR / "selection_frequency"
OUT_DIR = ATAC_DIR / "supplementary_tables"
OUT_PREFIX = "Supplementary_Table_5_ATAC_derived_STAGE_peaks"

COLUMNS = [
    "cell_type",
    "comparison",
    "peak",
    "selection_frequency",
    "mean_SHAP",
]

NUMERIC_COLUMNS = {"selection_frequency", "mean_SHAP"}

COMPARISON_ORDER = {
    "CON vs PRE": 0,
    "PRE vs T2D": 1,
    "CON vs T2D": 2,
}


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def first_present_column(row: dict[str, str], candidates: tuple[str, ...]) -> str:
    for column in candidates:
        if column in row:
            return column
    raise ValueError(f"None of these columns were found: {candidates}")


def parse_stage_path(path: Path) -> tuple[str, str]:
    """
    Convert:
    Acinar_cell/CON_vs_PRE_shap_marker_genes.csv
    β_cell/CON_vs_PRE_TEST_shap_mean.csv
    to:
    ("Acinar_cell", "CON vs PRE")
    """
    cell_type = path.parent.name
    name = path.stem

    for suffix in ("_shap_marker_genes", "_TEST_shap_mean"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break

    if "_vs_" not in name:
        raise ValueError(f"Unexpected STAGE peak filename format: {path}")

    comparison = name.replace("_vs_", " vs ")
    return cell_type, comparison


def load_selection_frequency(cell_type: str, comparison: str) -> dict[str, str]:
    comparison_stem = comparison.replace(" vs ", "_vs_")
    path = (
        SELECTION_FREQUENCY_DIR
        / cell_type
        / f"boruta_gene_stability_{comparison_stem}.csv"
    )
    if not path.exists():
        raise FileNotFoundError(f"Missing selection-frequency file: {path}")

    rows = read_csv_rows(path)
    if not rows:
        return {}

    feature_col = first_present_column(rows[0], ("peak", "gene"))
    if "fold_freq" not in rows[0]:
        raise ValueError(f"{path} must contain a fold_freq column")

    return {row[feature_col]: row["fold_freq"] for row in rows}


def sort_key(row: dict[str, str]) -> tuple:
    comparison_rank = COMPARISON_ORDER.get(row["comparison"], 99)
    try:
        selection_frequency = float(row["selection_frequency"])
    except (TypeError, ValueError):
        selection_frequency = -1.0
    try:
        abs_mean_shap = abs(float(row["mean_SHAP"]))
    except (TypeError, ValueError):
        abs_mean_shap = -1.0

    return (
        row["cell_type"],
        comparison_rank,
        -selection_frequency,
        -abs_mean_shap,
        row["peak"],
    )


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def cell_ref(row_idx: int, col_idx: int) -> str:
    col_name = ""
    while col_idx:
        col_idx, remainder = divmod(col_idx - 1, 26)
        col_name = chr(65 + remainder) + col_name
    return f"{col_name}{row_idx}"


def xml_text(value: object) -> str:
    return escape("" if value is None else str(value))


def numeric_value(value: object) -> str | None:
    try:
        return format(float(value), ".15g")
    except (TypeError, ValueError):
        return None


def xlsx_value(column: str, value: object) -> object:
    if column in NUMERIC_COLUMNS:
        parsed = numeric_value(value)
        if parsed is not None:
            return float(parsed)
    return value


def write_basic_xlsx(path: Path, rows: list[dict[str, str]]) -> None:
    """Write a simple one-sheet XLSX without third-party dependencies."""
    now = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    table_rows = [COLUMNS] + [[row[col] for col in COLUMNS] for row in rows]

    sheet_rows = []
    for row_idx, row_values in enumerate(table_rows, start=1):
        cells = []
        for col_idx, value in enumerate(row_values, start=1):
            ref = cell_ref(row_idx, col_idx)
            column = COLUMNS[col_idx - 1]
            number = numeric_value(value) if row_idx > 1 and column in NUMERIC_COLUMNS else None
            if number is None:
                cells.append(
                    f'<c r="{ref}" t="inlineStr"><is><t>{xml_text(value)}</t></is></c>'
                )
            else:
                cells.append(f'<c r="{ref}"><v>{number}</v></c>')
        sheet_rows.append(f'<row r="{row_idx}">{"".join(cells)}</row>')

    worksheet_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        '<sheetData>'
        f'{"".join(sheet_rows)}'
        '</sheetData>'
        '</worksheet>'
    )

    files = {
        "[Content_Types].xml": (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
            '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
            '<Default Extension="xml" ContentType="application/xml"/>'
            '<Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>'
            '<Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>'
            '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
            '<Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
            '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>'
            '</Types>'
        ),
        "_rels/.rels": (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>'
            '<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/>'
            '<Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" Target="docProps/app.xml"/>'
            '</Relationships>'
        ),
        "docProps/app.xml": (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties" '
            'xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes">'
            '<Application>Python</Application>'
            '</Properties>'
        ),
        "docProps/core.xml": (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" '
            'xmlns:dc="http://purl.org/dc/elements/1.1/" '
            'xmlns:dcterms="http://purl.org/dc/terms/" '
            'xmlns:dcmitype="http://purl.org/dc/dcmitype/" '
            'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'
            '<dc:creator>Python</dc:creator>'
            f'<dcterms:created xsi:type="dcterms:W3CDTF">{now}</dcterms:created>'
            f'<dcterms:modified xsi:type="dcterms:W3CDTF">{now}</dcterms:modified>'
            '</cp:coreProperties>'
        ),
        "xl/workbook.xml": (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
            'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
            '<sheets>'
            '<sheet name="ATAC_STAGE_peaks" sheetId="1" r:id="rId1"/>'
            '</sheets>'
            '</workbook>'
        ),
        "xl/_rels/workbook.xml.rels": (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/>'
            '<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>'
            '</Relationships>'
        ),
        "xl/styles.xml": (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
            '<fonts count="1"><font><sz val="11"/><name val="Calibri"/></font></fonts>'
            '<fills count="1"><fill><patternFill patternType="none"/></fill></fills>'
            '<borders count="1"><border/></borders>'
            '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
            '<cellXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/></cellXfs>'
            '</styleSheet>'
        ),
        "xl/worksheets/sheet1.xml": worksheet_xml,
    }

    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for file_name, content in files.items():
            archive.writestr(file_name, content)


def write_xlsx(path: Path, rows: list[dict[str, str]]) -> str:
    try:
        from openpyxl import Workbook
    except ImportError:
        write_basic_xlsx(path, rows)
        return "stdlib"

    workbook = Workbook()
    worksheet = workbook.active
    worksheet.title = "ATAC_STAGE_peaks"
    worksheet.append(COLUMNS)
    for row in rows:
        worksheet.append([xlsx_value(col, row[col]) for col in COLUMNS])

    workbook.save(path)
    return "openpyxl"


def main() -> None:
    if not STAGE_DIR.exists():
        raise FileNotFoundError(f"Missing STAGE peaks directory: {STAGE_DIR}")
    if not SELECTION_FREQUENCY_DIR.exists():
        raise FileNotFoundError(
            f"Missing selection-frequency directory: {SELECTION_FREQUENCY_DIR}"
        )

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    all_rows = []
    missing_frequency = []

    stage_files = sorted(STAGE_DIR.glob("*/*.csv"))
    if not stage_files:
        raise RuntimeError(f"No STAGE peak CSV files found in {STAGE_DIR}")

    for stage_path in stage_files:
        cell_type, comparison = parse_stage_path(stage_path)
        frequency_by_peak = load_selection_frequency(cell_type, comparison)
        stage_rows = read_csv_rows(stage_path)

        if not stage_rows:
            continue

        feature_col = first_present_column(stage_rows[0], ("peak", "gene"))
        if "mean_SHAP" not in stage_rows[0]:
            raise ValueError(f"{stage_path} must contain a mean_SHAP column")

        for stage_row in stage_rows:
            peak = stage_row[feature_col]
            selection_frequency = frequency_by_peak.get(peak, "0.0")
            if peak not in frequency_by_peak:
                missing_frequency.append((cell_type, comparison, peak))

            all_rows.append(
                {
                    "cell_type": cell_type,
                    "comparison": comparison,
                    "peak": peak,
                    "selection_frequency": selection_frequency,
                    "mean_SHAP": stage_row["mean_SHAP"],
                }
            )

    all_rows = sorted(all_rows, key=sort_key)

    csv_path = OUT_DIR / f"{OUT_PREFIX}.csv"
    xlsx_path = OUT_DIR / f"{OUT_PREFIX}.xlsx"

    write_csv(csv_path, all_rows)
    xlsx_writer = write_xlsx(xlsx_path, all_rows)

    print("\n========== Done ==========")
    print(f"Saved CSV : {csv_path}")
    print(f"Saved XLSX: {xlsx_path} ({xlsx_writer})")
    print(f"STAGE files read: {len(stage_files)}")
    print(f"Total rows: {len(all_rows)}")
    if missing_frequency:
        print(
            "[WARN] STAGE peaks absent from stability files; "
            f"selection_frequency set to 0.0: {len(missing_frequency)}"
        )


if __name__ == "__main__":
    main()
