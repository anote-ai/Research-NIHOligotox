"""
Convert the Phase 2 submission markdown documents into NIH-format-compliant PDFs.

Format rules (NIH OligoTox Challenge Announcement, "How to Enter"):
  - PDF, page size 8.5" x 11", margins >= 1"
  - Font >= 11pt Arial (Helvetica substituted here — no Arial TTF available in this
    environment; swap the FONT_* constants below if you have real Arial available, e.g.
    via Word/Google Docs export, for exact compliance)
  - Line spacing >= 1.0
  - No HHS/NIH/NCATS logos or seals (none are used here)

This is a small, purpose-built markdown subset renderer (headers, bold, bullet/numbered
lists, blockquotes, simple pipe tables, horizontal rules) — not a general markdown
converter — tuned to exactly the syntax used in the three source documents.

Usage:
    python3 phase2_submission/tools/build_pdfs.py
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import (
    SimpleDocTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
    HRFlowable,
)

FONT = "Helvetica"
FONT_BOLD = "Helvetica-Bold"
FONT_SIZE = 11
LEADING = FONT_SIZE * 1.3  # line spacing > 1.0

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

DOCS = [
    ("01_narrative_document.md", "01_narrative_document.pdf"),
    ("02_methodology_document.md", "02_methodology_document.pdf"),
    ("03_public_access_and_dissemination_plan.md", "03_public_access_and_dissemination_plan.pdf"),
]

styles = {
    "body": ParagraphStyle(
        "body", fontName=FONT, fontSize=FONT_SIZE, leading=LEADING,
        spaceAfter=8, alignment=4,  # justified
    ),
    "quote": ParagraphStyle(
        "quote", fontName=FONT, fontSize=FONT_SIZE, leading=LEADING,
        spaceAfter=8, leftIndent=18, textColor=colors.HexColor("#555555"),
    ),
    "h1": ParagraphStyle(
        "h1", fontName=FONT_BOLD, fontSize=18, leading=22, spaceBefore=4, spaceAfter=12,
    ),
    "h2": ParagraphStyle(
        "h2", fontName=FONT_BOLD, fontSize=14, leading=18, spaceBefore=14, spaceAfter=8,
    ),
    "h3": ParagraphStyle(
        "h3", fontName=FONT_BOLD, fontSize=12, leading=15, spaceBefore=10, spaceAfter=6,
    ),
    "meta": ParagraphStyle(
        "meta", fontName=FONT, fontSize=10, leading=13, textColor=colors.HexColor("#444444"),
        spaceAfter=4,
    ),
    "li": ParagraphStyle(
        "li", fontName=FONT, fontSize=FONT_SIZE, leading=LEADING,
    ),
    "ul_item": ParagraphStyle(
        "ul_item", fontName=FONT, fontSize=FONT_SIZE, leading=LEADING,
        leftIndent=18, bulletIndent=6, spaceAfter=4,
    ),
    "ol_item": ParagraphStyle(
        "ol_item", fontName=FONT, fontSize=FONT_SIZE, leading=LEADING,
        leftIndent=20, bulletIndent=6, spaceAfter=4,
    ),
}


def inline(text: str) -> str:
    """Convert a subset of inline markdown to reportlab mini-XML."""
    text = text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    text = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", text)
    text = re.sub(r"(?<!\*)\*([^*]+?)\*(?!\*)", r"<i>\1</i>", text)
    text = re.sub(r"`(.+?)`", r"<font face='Courier'>\1</font>", text)
    text = re.sub(r"\[(.+?)\]\((.+?)\)", r'<link href="\2">\1</link>', text)
    return text


def parse_table(lines: list[str]) -> Table:
    rows = []
    for line in lines:
        if re.match(r"^\|?\s*:?-{2,}:?\s*(\|\s*:?-{2,}:?\s*)+\|?$", line):
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        rows.append([Paragraph(inline(c), styles["li"]) for c in cells])
    t = Table(rows, hAlign="LEFT")
    t.setStyle(TableStyle([
        ("FONT", (0, 0), (-1, -1), FONT, 9),
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#e8e8e8")),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#999999")),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 4),
        ("RIGHTPADDING", (0, 0), (-1, -1), 4),
        ("TOPPADDING", (0, 0), (-1, -1), 3),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
    ]))
    return t


def build_flowables(md_text: str) -> list:
    flow = []
    lines = md_text.splitlines()
    i = 0
    n = len(lines)
    list_buffer: list[tuple[str, str]] = []  # (kind, text)

    def flush_list():
        nonlocal list_buffer
        if not list_buffer:
            return
        ol_counter = 0
        for kind, text in list_buffer:
            if kind == "ul":
                flow.append(Paragraph(f"•&nbsp;&nbsp;{inline(text)}", styles["ul_item"]))
            else:
                ol_counter += 1
                flow.append(Paragraph(f"{ol_counter}.&nbsp;&nbsp;{inline(text)}", styles["ol_item"]))
        flow.append(Spacer(1, 4))
        list_buffer = []

    while i < n:
        line = lines[i].rstrip("\n")
        stripped = line.strip()

        if not stripped:
            flush_list()
            i += 1
            continue

        if stripped.startswith("|"):
            flush_list()
            table_lines = []
            while i < n and lines[i].strip().startswith("|"):
                table_lines.append(lines[i])
                i += 1
            flow.append(parse_table(table_lines))
            flow.append(Spacer(1, 10))
            continue

        if re.match(r"^-{3,}$", stripped):
            flush_list()
            flow.append(HRFlowable(width="100%", thickness=0.75, color=colors.HexColor("#bbbbbb"),
                                    spaceBefore=6, spaceAfter=6))
            i += 1
            continue

        if stripped.startswith("> "):
            flush_list()
            quote_lines = []
            while i < n and lines[i].strip().startswith(">"):
                quote_lines.append(lines[i].strip().lstrip(">").strip())
                i += 1
            flow.append(Paragraph(inline(" ".join(quote_lines)), styles["quote"]))
            continue

        m = re.match(r"^(#{1,3})\s+(.*)$", stripped)
        if m:
            flush_list()
            level = len(m.group(1))
            style = {1: "h1", 2: "h2", 3: "h3"}[level]
            flow.append(Paragraph(inline(m.group(2)), styles[style]))
            i += 1
            continue

        m = re.match(r"^[-*]\s+(.*)$", stripped)
        if m:
            list_buffer.append(("ul", m.group(1)))
            i += 1
            continue

        m = re.match(r"^\d+\.\s+(.*)$", stripped)
        if m:
            list_buffer.append(("ol", m.group(1)))
            i += 1
            continue

        if stripped.startswith("**") and stripped.count("**") <= 2 and len(stripped) < 120:
            flush_list()
            flow.append(Paragraph(inline(stripped), styles["meta"]))
            i += 1
            continue

        flush_list()
        flow.append(Paragraph(inline(stripped), styles["body"]))
        i += 1

    flush_list()
    return flow


def build_pdf(md_path: Path, pdf_path: Path) -> int:
    text = md_path.read_text()
    doc = SimpleDocTemplate(
        str(pdf_path), pagesize=letter,
        leftMargin=1 * inch, rightMargin=1 * inch,
        topMargin=1 * inch, bottomMargin=1 * inch,
        title=md_path.stem,
    )
    flow = build_flowables(text)
    doc.build(flow)
    # Count pages by re-reading (reportlab doesn't expose page count post-build directly;
    # use pypdf if available, else fall back to canvas page counter written during build)
    try:
        from pypdf import PdfReader
        return len(PdfReader(str(pdf_path)).pages)
    except ImportError:
        return -1


def main() -> None:
    out_dir = ROOT / "pdf"
    out_dir.mkdir(exist_ok=True)
    for src_name, out_name in DOCS:
        src = ROOT / src_name
        out = out_dir / out_name
        pages = build_pdf(src, out)
        page_str = f"{pages} pages" if pages >= 0 else "page count unknown (pypdf not installed)"
        print(f"{src_name} -> {out.relative_to(ROOT.parent)}  ({page_str})")


if __name__ == "__main__":
    main()
