from __future__ import annotations

from datetime import UTC, datetime
from decimal import Decimal, InvalidOperation
from pathlib import Path
from xml.sax.saxutils import escape
from zipfile import ZIP_DEFLATED, ZipFile


def write_rows_xlsx(
    path: str | Path,
    fieldnames: list[str],
    rows: list[dict[str, str]],
    *,
    sheet_name: str = "Sheet1",
    header_labels: dict[str, str] | None = None,
    percentage_fields: set[str] | None = None,
    column_widths: dict[str, float] | None = None,
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    sheet_name = _safe_sheet_name(sheet_name)
    with ZipFile(path, "w", ZIP_DEFLATED) as archive:
        archive.writestr("[Content_Types].xml", _content_types_xml())
        archive.writestr("_rels/.rels", _root_rels_xml())
        archive.writestr("docProps/app.xml", _app_xml())
        archive.writestr("docProps/core.xml", _core_xml())
        archive.writestr("xl/workbook.xml", _workbook_xml(sheet_name))
        archive.writestr("xl/_rels/workbook.xml.rels", _workbook_rels_xml())
        archive.writestr("xl/styles.xml", _styles_xml())
        archive.writestr(
            "xl/worksheets/sheet1.xml",
            _worksheet_xml(
                fieldnames,
                rows,
                header_labels or {},
                percentage_fields or set(),
                column_widths or {},
            ),
        )


def _safe_sheet_name(value: str) -> str:
    invalid = set("[]:*?/\\")
    cleaned = "".join("_" if char in invalid else char for char in value).strip()
    return (cleaned or "Sheet1")[:31]


def _column_name(index: int) -> str:
    name = ""
    while index:
        index, remainder = divmod(index - 1, 26)
        name = chr(65 + remainder) + name
    return name


def _cell_ref(row_number: int, column_number: int) -> str:
    return f"{_column_name(column_number)}{row_number}"


def _is_integer_text(value: str) -> bool:
    if not value:
        return False
    return value.isdigit() or (value.startswith("-") and value[1:].isdigit())


def _percentage_fraction_text(value: str) -> str | None:
    if not value.endswith("%"):
        return None
    try:
        fraction = Decimal(value[:-1]) / Decimal(100)
    except InvalidOperation:
        return None
    return format(fraction, "f")


def _cell_xml(row_number: int, column_number: int, value: str, *, percentage: bool = False) -> str:
    ref = _cell_ref(row_number, column_number)
    if percentage:
        fraction = _percentage_fraction_text(value)
        if fraction is not None:
            return f'<c r="{ref}" s="1"><v>{fraction}</v></c>'
    if _is_integer_text(value):
        return f'<c r="{ref}"><v>{value}</v></c>'
    return f'<c r="{ref}" t="inlineStr"><is><t>{escape(value)}</t></is></c>'


def _worksheet_xml(
    fieldnames: list[str],
    rows: list[dict[str, str]],
    header_labels: dict[str, str],
    percentage_fields: set[str],
    column_widths: dict[str, float],
) -> str:
    max_row = len(rows) + 1
    max_col = max(1, len(fieldnames))
    dimension = f"A1:{_cell_ref(max_row, max_col)}"
    xml_rows = []

    header_cells = [
        _cell_xml(1, index, header_labels.get(name, name))
        for index, name in enumerate(fieldnames, start=1)
    ]
    xml_rows.append(f'<row r="1">{"".join(header_cells)}</row>')

    for row_number, row in enumerate(rows, start=2):
        cells = [
            _cell_xml(
                row_number,
                index,
                str(row.get(field, "")),
                percentage=field in percentage_fields,
            )
            for index, field in enumerate(fieldnames, start=1)
        ]
        xml_rows.append(f'<row r="{row_number}">{"".join(cells)}</row>')

    auto_filter = f'<autoFilter ref="{dimension}"/>' if fieldnames else ""
    columns = []
    for index, field in enumerate(fieldnames, start=1):
        width = column_widths.get(field)
        if width is not None:
            columns.append(
                f'<col min="{index}" max="{index}" width="{width:g}" customWidth="1"/>'
            )
    columns_xml = f'<cols>{"".join(columns)}</cols>' if columns else ""
    return f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
  <dimension ref="{dimension}"/>
  <sheetViews>
    <sheetView workbookViewId="0">
      <pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/>
      <selection pane="bottomLeft"/>
    </sheetView>
  </sheetViews>
  <sheetFormatPr defaultRowHeight="15"/>
  {columns_xml}
  <sheetData>
    {"".join(xml_rows)}
  </sheetData>
  {auto_filter}
</worksheet>'''


def _content_types_xml() -> str:
    return '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
  <Default Extension="xml" ContentType="application/xml"/>
  <Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>
  <Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>
  <Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>
  <Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>
  <Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>
</Types>'''


def _root_rels_xml() -> str:
    return '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>
  <Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/>
  <Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" Target="docProps/app.xml"/>
</Relationships>'''


def _workbook_xml(sheet_name: str) -> str:
    return f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
  <sheets>
    <sheet name="{escape(sheet_name)}" sheetId="1" r:id="rId1"/>
  </sheets>
</workbook>'''


def _workbook_rels_xml() -> str:
    return '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/>
  <Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>
</Relationships>'''


def _styles_xml() -> str:
    return '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
  <numFmts count="1"><numFmt numFmtId="164" formatCode="0.000000%"/></numFmts>
  <fonts count="1"><font><sz val="11"/><name val="Calibri"/></font></fonts>
  <fills count="1"><fill><patternFill patternType="none"/></fill></fills>
  <borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>
  <cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>
  <cellXfs count="2"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/><xf numFmtId="164" fontId="0" fillId="0" borderId="0" xfId="0" applyNumberFormat="1"/></cellXfs>
  <cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>
</styleSheet>'''


def _app_xml() -> str:
    return '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties" xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes">
  <Application>CPM Paired Fastq IgBLAST AIRR tsv Ver3.0</Application>
</Properties>'''


def _core_xml() -> str:
    created = datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")
    return f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:dcmitype="http://purl.org/dc/dcmitype/" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <dc:creator>CPM Paired Fastq IgBLAST AIRR tsv Ver3.0</dc:creator>
  <cp:lastModifiedBy>CPM Paired Fastq IgBLAST AIRR tsv Ver3.0</cp:lastModifiedBy>
  <dcterms:created xsi:type="dcterms:W3CDTF">{created}</dcterms:created>
  <dcterms:modified xsi:type="dcterms:W3CDTF">{created}</dcterms:modified>
</cp:coreProperties>'''
