"""
Convert the Phase 2 submission markdown documents into NIH-format-compliant, professionally
formatted PDFs.

Format rules (NIH OligoTox Challenge Announcement, "How to Enter"):
  - PDF, page size 8.5" x 11", margins >= 1"
  - Font >= 11pt Arial (Helvetica substituted here — no Arial TTF available in this
    environment; swap the FONT_* constants below if you have real Arial available, e.g.
    via Word/Google Docs export, for exact compliance)
  - Line spacing >= 1.0
  - No HHS/NIH/NCATS logos or seals (none are used here; the navy accent color below is a
    generic design choice, not any agency's branding)

This is a small, purpose-built markdown subset renderer (headers, bold/italic, bullet/numbered
lists, blockquotes-as-callouts, simple pipe tables, horizontal rules) — not a general markdown
converter — tuned to exactly the syntax used in the three source documents.

Usage:
    python3 phase2_submission/tools/build_pdfs.py
"""
from __future__ import annotations

import re
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import inch
from reportlab.pdfgen import canvas as pdfcanvas
from reportlab.platypus import (
    SimpleDocTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
    HRFlowable,
    KeepTogether,
)

FONT = "Helvetica"
FONT_BOLD = "Helvetica-Bold"
FONT_ITALIC = "Helvetica-Oblique"
FONT_SIZE = 11
LEADING = FONT_SIZE * 1.35  # line spacing > 1.0

ACCENT = colors.HexColor("#1F3B57")       # dark navy — generic, not agency branding
ACCENT_LIGHT = colors.HexColor("#EEF2F6")
WARN_BG = colors.HexColor("#FBF0DA")
WARN_BORDER = colors.HexColor("#C98A1A")
RULE = colors.HexColor("#C7CFD6")
MUTED = colors.HexColor("#54606B")

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

DOCS = [
    dict(src="01_narrative_document.md", out="01_narrative_document.pdf",
         short_title="OligoTox Phase 2 — Narrative Document",
         page_limit=12, draft_footer=True),
    dict(src="02_methodology_document.md", out="02_methodology_document.pdf",
         short_title="OligoTox Phase 2 — Methodology Document",
         page_limit=5, draft_footer=True),
    dict(src="03_public_access_and_dissemination_plan.md",
         out="03_public_access_and_dissemination_plan.pdf",
         short_title="OligoTox Phase 2 — Public Access and Dissemination Plan",
         page_limit=5, draft_footer=False),
]

styles = {
    "body": ParagraphStyle(
        "body", fontName=FONT, fontSize=FONT_SIZE, leading=LEADING,
        spaceAfter=8, alignment=4,  # justified
    ),
    "h1": ParagraphStyle(
        "h1", fontName=FONT_BOLD, fontSize=19, leading=23, spaceBefore=0, spaceAfter=6,
        textColor=ACCENT,
    ),
    "h2": ParagraphStyle(
        "h2", fontName=FONT_BOLD, fontSize=13.5, leading=17, spaceBefore=16, spaceAfter=6,
        textColor=ACCENT,
    ),
    "h3": ParagraphStyle(
        "h3", fontName=FONT_BOLD, fontSize=11.5, leading=15, spaceBefore=10, spaceAfter=5,
        textColor=colors.HexColor("#2E4A66"),
    ),
    "meta_label": ParagraphStyle(
        "meta_label", fontName=FONT_BOLD, fontSize=9.5, leading=13.5, textColor=MUTED,
    ),
    "meta_value": ParagraphStyle(
        "meta_value", fontName=FONT, fontSize=9.5, leading=13.5, textColor=colors.black,
    ),
    "li": ParagraphStyle(
        "li", fontName=FONT, fontSize=FONT_SIZE, leading=LEADING,
    ),
    "ul_item": ParagraphStyle(
        "ul_item", fontName=FONT, fontSize=FONT_SIZE, leading=LEADING,
        leftIndent=18, spaceAfter=4,
    ),
    "ol_item": ParagraphStyle(
        "ol_item", fontName=FONT, fontSize=FONT_SIZE, leading=LEADING,
        leftIndent=20, spaceAfter=4,
    ),
    "callout": ParagraphStyle(
        "callout", fontName=FONT, fontSize=9.7, leading=13.2, textColor=colors.HexColor("#5C3E0A"),
    ),
    "table_cell": ParagraphStyle(
        "table_cell", fontName=FONT, fontSize=9, leading=12,
    ),
    "table_head": ParagraphStyle(
        "table_head", fontName=FONT_BOLD, fontSize=9, leading=12, textColor=colors.white,
    ),
    "footer": ParagraphStyle(
        "footer", fontName=FONT, fontSize=8, leading=10, textColor=MUTED,
    ),
}


def inline(text: str) -> str:
    """Convert a subset of inline markdown to reportlab mini-XML."""
    text = text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    text = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", text)
    text = re.sub(r"(?<!\*)\*([^*]+?)\*(?!\*)", r"<i>\1</i>", text)
    text = re.sub(r"`(.+?)`", r"<font face='Courier'>\1</font>", text)
    text = re.sub(r"\[(.+?)\]\((.+?)\)", r'<link href="\2"><u>\1</u></link>', text)
    return text


def parse_table(lines: list[str]) -> Table:
    rows = []
    for line in lines:
        if re.match(r"^\|?\s*:?-{2,}:?\s*(\|\s*:?-{2,}:?\s*)+\|?$", line):
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        rows.append(cells)
    header, *body = rows
    data = [[Paragraph(inline(c), styles["table_head"]) for c in header]]
    for r in body:
        data.append([Paragraph(inline(c), styles["table_cell"]) for c in r])
    t = Table(data, hAlign="LEFT", repeatRows=1)
    style = [
        ("BACKGROUND", (0, 0), (-1, 0), ACCENT),
        ("LINEBELOW", (0, 0), (-1, 0), 1, ACCENT),
        ("LINEBELOW", (0, -1), (-1, -1), 0.75, RULE),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 5),
        ("RIGHTPADDING", (0, 0), (-1, -1), 5),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
    ]
    for row_idx in range(1, len(data)):
        if row_idx % 2 == 0:
            style.append(("BACKGROUND", (0, row_idx), (-1, row_idx), ACCENT_LIGHT))
    t.setStyle(TableStyle(style))
    return t


def make_callout(text: str) -> Table:
    p = Paragraph(inline(text), styles["callout"])
    t = Table([[p]], colWidths=[6.5 * inch], hAlign="LEFT")
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), WARN_BG),
        ("LINEBEFORE", (0, 0), (0, -1), 3, WARN_BORDER),
        ("TOPPADDING", (0, 0), (-1, -1), 8),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 8),
        ("LEFTPADDING", (0, 0), (-1, -1), 10),
        ("RIGHTPADDING", (0, 0), (-1, -1), 10),
    ]))
    return t


def build_flowables(md_text: str) -> tuple[list, list, list[tuple[str, str]]]:
    """Returns (title_flowables, body_flowables, metadata_kv_pairs)."""
    title_flow = []
    flow = []
    lines = md_text.splitlines()
    i = 0
    n = len(lines)
    list_buffer: list[tuple[str, str]] = []
    meta_pairs: list[tuple[str, str]] = []
    title_done = False

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

        # Title (H1) rendered specially; bold "**Key:** value" lines right after it become
        # a metadata block instead of plain paragraphs.
        m1 = re.match(r"^#\s+(.*)$", stripped)
        if m1 and not title_done:
            title_flow.append(Paragraph(inline(m1.group(1)), styles["h1"]))
            title_flow.append(HRFlowable(width="100%", thickness=1.4, color=ACCENT, spaceBefore=2, spaceAfter=10))
            title_done = True
            i += 1
            # consume immediately-following **Key:** value metadata lines
            while i < n:
                meta_line = lines[i].strip()
                mm = re.match(r"^\*\*(.+?)\*\*:?\s*(.*)$", meta_line)
                if mm:
                    meta_pairs.append((mm.group(1).rstrip(":"), mm.group(2)))
                    i += 1
                elif not meta_line:
                    i += 1
                else:
                    break
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
            flow.append(HRFlowable(width="100%", thickness=0.75, color=RULE, spaceBefore=8, spaceAfter=8))
            i += 1
            continue

        if stripped.startswith("> "):
            flush_list()
            quote_lines = []
            while i < n and lines[i].strip().startswith(">"):
                quote_lines.append(lines[i].strip().lstrip(">").strip())
                i += 1
            flow.append(Spacer(1, 4))
            flow.append(make_callout(" ".join(quote_lines)))
            flow.append(Spacer(1, 8))
            continue

        m = re.match(r"^(#{1,3})\s+(.*)$", stripped)
        if m:
            flush_list()
            level = len(m.group(1))
            style = {1: "h1", 2: "h2", 3: "h3"}[level]
            heading = Paragraph(inline(m.group(2)), styles[style])
            if level == 2:
                flow.append(KeepTogether([
                    heading,
                    HRFlowable(width="100%", thickness=0.6, color=RULE, spaceBefore=2, spaceAfter=6),
                ]))
            else:
                flow.append(heading)
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

        flush_list()
        flow.append(Paragraph(inline(stripped), styles["body"]))
        i += 1

    flush_list()
    return title_flow, flow, meta_pairs


def make_metadata_box(meta_pairs: list[tuple[str, str]]) -> Table:
    rows = [[Paragraph(f"{k}:", styles["meta_label"]), Paragraph(inline(v), styles["meta_value"])]
            for k, v in meta_pairs]
    t = Table(rows, colWidths=[1.35 * inch, 5.15 * inch], hAlign="LEFT")
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), ACCENT_LIGHT),
        ("BOX", (0, 0), (-1, -1), 0.75, RULE),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 8),
        ("RIGHTPADDING", (0, 0), (-1, -1), 8),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
    ]))
    return t


def make_canvas_factory(short_title: str, draft_footer: bool):
    class NumberedCanvas(pdfcanvas.Canvas):
        def __init__(self, *args, **kwargs):
            pdfcanvas.Canvas.__init__(self, *args, **kwargs)
            self._saved_states = []

        def showPage(self):
            self._saved_states.append(dict(self.__dict__))
            self._startPage()

        def save(self):
            total = len(self._saved_states)
            for state in self._saved_states:
                self.__dict__.update(state)
                self._draw_furniture(total)
                pdfcanvas.Canvas.showPage(self)
            pdfcanvas.Canvas.save(self)

        def _draw_furniture(self, total_pages):
            page_num = self._pageNumber
            w, h = letter
            self.setStrokeColor(RULE)
            self.setLineWidth(0.5)
            self.line(1 * inch, 0.75 * inch, w - 1 * inch, 0.75 * inch)
            self.setFont(FONT, 8)
            self.setFillColor(MUTED)
            self.drawString(1 * inch, 0.58 * inch, short_title)
            self.drawRightString(w - 1 * inch, 0.58 * inch, f"Page {page_num} of {total_pages}")
            if draft_footer:
                self.setFillColor(WARN_BORDER)
                self.setFont(FONT_BOLD, 7.5)
                self.drawCentredString(w / 2, 0.58 * inch, "DRAFT — PENDING REAL EXPERIMENTAL DATA")

    return NumberedCanvas


def build_pdf(md_path: Path, pdf_path: Path, short_title: str, draft_footer: bool) -> int:
    text = md_path.read_text()
    doc = SimpleDocTemplate(
        str(pdf_path), pagesize=letter,
        leftMargin=1 * inch, rightMargin=1 * inch,
        topMargin=1 * inch, bottomMargin=1 * inch,
        title=md_path.stem,
    )
    title_flow, flow, meta_pairs = build_flowables(text)
    full_flow = list(title_flow)
    if meta_pairs:
        full_flow.append(make_metadata_box(meta_pairs))
        full_flow.append(Spacer(1, 14))
    full_flow.extend(flow)
    doc.build(full_flow, canvasmaker=make_canvas_factory(short_title, draft_footer))
    try:
        from pypdf import PdfReader
        return len(PdfReader(str(pdf_path)).pages)
    except ImportError:
        return -1


def main() -> None:
    out_dir = ROOT / "pdf"
    out_dir.mkdir(exist_ok=True)
    for spec in DOCS:
        src = ROOT / spec["src"]
        out = out_dir / spec["out"]
        pages = build_pdf(src, out, spec["short_title"], spec["draft_footer"])
        limit = spec["page_limit"]
        flag = "OK" if 0 <= pages <= limit else "CHECK LIMIT" if pages >= 0 else "?"
        print(f"{spec['src']} -> {out.relative_to(ROOT.parent)}  ({pages} pages, limit {limit}) [{flag}]")


if __name__ == "__main__":
    main()
