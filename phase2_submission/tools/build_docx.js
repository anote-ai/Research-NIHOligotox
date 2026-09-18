/**
 * Build editable .docx versions of the three Phase 2 submission documents.
 *
 * Requires the `docx` npm package (not vendored in this repo):
 *   npm install docx
 *   node phase2_submission/tools/build_docx.js
 *
 * Mirrors the layout of build_pdfs.py (title, metadata box, styled callouts,
 * section headers, tables, footer with page numbers) but produces real,
 * Word-editable .docx output using true Arial instead of the PDF's Helvetica
 * substitute. Reuses the same wrapped-line-merging fix as build_pdfs.py so
 * multi-line bullets/paragraphs in the source markdown don't get split.
 */
const fs = require("fs");
const path = require("path");
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, AlignmentType,
  Table, TableRow, TableCell, WidthType, BorderStyle, ShadingType,
  Header, Footer, PageNumber, NumberOfTotalPages, TabStopType, TabStopPosition,
  convertInchesToTwip, VerticalAlign, LevelFormat,
} = require("docx");

const ACCENT = "1F3B57";
const ACCENT_LIGHT = "EEF2F6";
const WARN_BG = "FBF0DA";
const WARN_BORDER = "C98A1A";
const RULE = "C7CFD6";
const MUTED = "54606B";
const FONT = "Arial";

// ---------- inline markdown -> TextRun[] ----------
function inlineRuns(text, baseOpts = {}) {
  const runs = [];
  let i = 0;
  const push = (t, opts) => {
    if (!t) return;
    runs.push(new TextRun({ text: t, font: FONT, ...baseOpts, ...opts }));
  };
  const re = /\*\*(.+?)\*\*|`([^`]+?)`|\*([^*]+?)\*/g;
  let last = 0, m;
  while ((m = re.exec(text)) !== null) {
    if (m.index > last) push(text.slice(last, m.index), {});
    if (m[1] !== undefined) push(m[1], { bold: true });
    else if (m[2] !== undefined) push(m[2], { font: "Courier New", size: 20 });
    else if (m[3] !== undefined) push(m[3], { italics: true });
    last = re.lastIndex;
  }
  if (last < text.length) push(text.slice(last), {});
  if (runs.length === 0) push(text, {});
  return runs;
}

// ---------- markdown block parser (mirrors tools/build_pdfs.py logic) ----------
function parseDoc(mdText) {
  const lines = mdText.split("\n");
  let i = 0;
  const n = lines.length;
  let title = "";
  const meta = [];
  const blocks = [];
  const NEW_BLOCK = /^(\|)|^-{3,}$|^>|^#{1,3}\s+|^[-*]\s+|^\d+\.\s+/;

  // title (H1)
  while (i < n && lines[i].trim() === "") i++;
  if (i < n) {
    const m1 = lines[i].trim().match(/^#\s+(.*)$/);
    if (m1) {
      title = m1[1];
      i++;
      while (i < n) {
        const mline = lines[i].trim();
        const mm = mline.match(/^\*\*(.+?)\*\*:?\s*(.*)$/);
        if (mm) { meta.push([mm[1].replace(/:$/, ""), mm[2]]); i++; }
        else if (mline === "") { i++; }
        else break;
      }
    }
  }

  while (i < n) {
    const raw = lines[i];
    const s = raw.trim();
    if (s === "") { i++; continue; }

    if (s.startsWith("|")) {
      const tlines = [];
      while (i < n && lines[i].trim().startsWith("|")) { tlines.push(lines[i].trim()); i++; }
      const rows = tlines
        .filter((l) => !/^\|?\s*:?-{2,}:?\s*(\|\s*:?-{2,}:?\s*)+\|?$/.test(l))
        .map((l) => l.replace(/^\|/, "").replace(/\|$/, "").split("|").map((c) => c.trim()));
      blocks.push({ type: "table", header: rows[0], rows: rows.slice(1) });
      continue;
    }

    if (/^-{3,}$/.test(s)) { blocks.push({ type: "hr" }); i++; continue; }

    if (s.startsWith(">")) {
      const qlines = [];
      while (i < n && lines[i].trim().startsWith(">")) {
        qlines.push(lines[i].trim().replace(/^>\s?/, ""));
        i++;
      }
      blocks.push({ type: "callout", text: qlines.join(" ") });
      continue;
    }

    const hm = s.match(/^(#{1,3})\s+(.*)$/);
    if (hm) { blocks.push({ type: "h" + hm[1].length, text: hm[2] }); i++; continue; }

    const ulm = s.match(/^[-*]\s+(.*)$/);
    if (ulm) {
      const items = [["ul", ulm[1]]];
      i++;
      while (i < n) {
        const ls = lines[i].trim();
        if (ls === "" || NEW_BLOCK.test(ls) && !ls.match(/^[-*]\s+/) && !ls.match(/^\d+\.\s+/)) break;
        const u2 = ls.match(/^[-*]\s+(.*)$/);
        const o2 = ls.match(/^\d+\.\s+(.*)$/);
        if (u2) { items.push(["ul", u2[1]]); i++; }
        else if (o2) { items.push(["ol", o2[1]]); i++; }
        else { items[items.length - 1][1] += " " + ls; i++; }
      }
      blocks.push({ type: "list", items });
      continue;
    }

    const olm = s.match(/^\d+\.\s+(.*)$/);
    if (olm) {
      const items = [["ol", olm[1]]];
      i++;
      while (i < n) {
        const ls = lines[i].trim();
        if (ls === "" || (NEW_BLOCK.test(ls) && !ls.match(/^[-*]\s+/) && !ls.match(/^\d+\.\s+/))) break;
        const u2 = ls.match(/^[-*]\s+(.*)$/);
        const o2 = ls.match(/^\d+\.\s+(.*)$/);
        if (u2) { items.push(["ul", u2[1]]); i++; }
        else if (o2) { items.push(["ol", o2[1]]); i++; }
        else { items[items.length - 1][1] += " " + ls; i++; }
      }
      blocks.push({ type: "list", items });
      continue;
    }

    // paragraph with continuation merging
    const plines = [s];
    i++;
    while (i < n) {
      const ls = lines[i].trim();
      if (ls === "" || NEW_BLOCK.test(ls)) break;
      plines.push(ls);
      i++;
    }
    blocks.push({ type: "p", text: plines.join(" ") });
  }

  return { title, meta, blocks };
}

// ---------- docx rendering ----------
const numbering = {
  config: [
    {
      reference: "bullets",
      levels: [{ level: 0, format: LevelFormat.BULLET, text: "•", alignment: AlignmentType.LEFT,
        style: { paragraph: { indent: { left: convertInchesToTwip(0.35), hanging: convertInchesToTwip(0.2) } } } }],
    },
  ],
};

function metaTable(meta) {
  return new Table({
    width: { size: 100, type: WidthType.PERCENTAGE },
    columnWidths: [2200, 8000],
    borders: allBorders(RULE, 4),
    rows: meta.map(([k, v]) =>
      new TableRow({
        children: [
          cell(k + ":", { bold: true, size: 19, color: MUTED }, 2200, ACCENT_LIGHT),
          cell(v, { size: 19 }, 8000, ACCENT_LIGHT, true),
        ],
      })
    ),
  });
}

function allBorders(color, size) {
  const b = { style: BorderStyle.SINGLE, size, color };
  return { top: b, bottom: b, left: b, right: b, insideHorizontal: b, insideVertical: b };
}

function cell(text, runOpts, width, shade, richInline) {
  return new TableCell({
    width: { size: width, type: WidthType.DXA },
    shading: shade ? { type: ShadingType.CLEAR, fill: shade } : undefined,
    verticalAlign: VerticalAlign.TOP,
    margins: { top: 80, bottom: 80, left: 120, right: 120 },
    children: [new Paragraph({ children: richInline ? inlineRuns(text, runOpts) : [new TextRun({ text, font: FONT, ...runOpts })] })],
  });
}

function calloutBlock(text) {
  return new Table({
    width: { size: 100, type: WidthType.PERCENTAGE },
    columnWidths: [10200],
    borders: {
      top: { style: BorderStyle.SINGLE, size: 4, color: WARN_BG },
      bottom: { style: BorderStyle.SINGLE, size: 4, color: WARN_BG },
      right: { style: BorderStyle.SINGLE, size: 4, color: WARN_BG },
      left: { style: BorderStyle.SINGLE, size: 24, color: WARN_BORDER },
    },
    rows: [
      new TableRow({
        children: [
          new TableCell({
            width: { size: 10200, type: WidthType.DXA },
            shading: { type: ShadingType.CLEAR, fill: WARN_BG },
            margins: { top: 160, bottom: 160, left: 200, right: 200 },
            children: [new Paragraph({ children: inlineRuns(text, { size: 19, color: "5C3E0A" }) })],
          }),
        ],
      }),
    ],
  });
}

function dataTable(header, rows) {
  const nCols = header.length;
  const total = 10200;
  const colW = Math.floor(total / nCols);
  const widths = new Array(nCols).fill(colW);
  return new Table({
    width: { size: 100, type: WidthType.PERCENTAGE },
    columnWidths: widths,
    borders: allBorders(RULE, 4),
    rows: [
      new TableRow({
        tableHeader: true,
        children: header.map((h) => cell(h, { bold: true, color: "FFFFFF", size: 18 }, colW, ACCENT, true)),
      }),
      ...rows.map((r, idx) =>
        new TableRow({
          children: r.map((c) => cell(c, { size: 18 }, colW, idx % 2 === 1 ? ACCENT_LIGHT : undefined, true)),
        })
      ),
    ],
  });
}

function headingParagraph(text, level) {
  const sizeMap = { h2: 27, h3: 23 };
  return new Paragraph({
    heading: level === "h2" ? HeadingLevel.HEADING_2 : HeadingLevel.HEADING_3,
    spacing: { before: level === "h2" ? 320 : 200, after: 120 },
    border: level === "h2" ? { bottom: { style: BorderStyle.SINGLE, size: 4, color: RULE, space: 4 } } : undefined,
    children: [new TextRun({ text, bold: true, font: FONT, size: sizeMap[level], color: level === "h2" ? ACCENT : "2E4A66" })],
  });
}

function bodyParagraph(text) {
  return new Paragraph({
    alignment: AlignmentType.JUSTIFIED,
    spacing: { after: 160, line: 300, lineRule: "atLeast" },
    children: inlineRuns(text),
  });
}

function listParagraphs(items) {
  let olCounter = 0;
  return items.map(([kind, text]) => {
    if (kind === "ul") {
      return new Paragraph({
        numbering: { reference: "bullets", level: 0 },
        spacing: { after: 100, line: 300, lineRule: "atLeast" },
        children: inlineRuns(text),
      });
    } else {
      olCounter += 1;
      return new Paragraph({
        indent: { left: convertInchesToTwip(0.35), hanging: convertInchesToTwip(0.35) },
        spacing: { after: 100, line: 300, lineRule: "atLeast" },
        tabStops: [{ type: TabStopType.LEFT, position: convertInchesToTwip(0.35) }],
        children: [new TextRun({ text: `${olCounter}.\t`, font: FONT }), ...inlineRuns(text)],
      });
    }
  });
}

function buildBody(parsed) {
  const out = [];
  out.push(new Paragraph({
    spacing: { after: 80 },
    border: { bottom: { style: BorderStyle.SINGLE, size: 8, color: ACCENT, space: 6 } },
    children: [new TextRun({ text: parsed.title, bold: true, font: FONT, size: 38, color: ACCENT })],
  }));
  if (parsed.meta.length) {
    out.push(new Paragraph({ spacing: { before: 160, after: 200 }, children: [] }));
    out.push(metaTable(parsed.meta));
  }
  out.push(new Paragraph({ spacing: { after: 120 }, children: [] }));

  for (const b of parsed.blocks) {
    if (b.type === "h2" || b.type === "h3") out.push(headingParagraph(b.text, b.type));
    else if (b.type === "p") out.push(bodyParagraph(b.text));
    else if (b.type === "list") out.push(...listParagraphs(b.items));
    else if (b.type === "callout") {
      out.push(new Paragraph({ spacing: { after: 40 }, children: [] }));
      out.push(calloutBlock(b.text));
      out.push(new Paragraph({ spacing: { after: 160 }, children: [] }));
    } else if (b.type === "table") {
      out.push(dataTable(b.header, b.rows));
      out.push(new Paragraph({ spacing: { after: 200 }, children: [] }));
    } else if (b.type === "hr") {
      out.push(new Paragraph({
        spacing: { before: 120, after: 120 },
        border: { bottom: { style: BorderStyle.SINGLE, size: 4, color: RULE } },
        children: [],
      }));
    }
  }
  return out;
}

function buildFooter(shortTitle, draftFooter) {
  const children = [];
  if (draftFooter) {
    children.push(new Paragraph({
      alignment: AlignmentType.CENTER,
      spacing: { after: 60 },
      children: [new TextRun({ text: "DRAFT — PENDING REAL EXPERIMENTAL DATA", bold: true, font: FONT, size: 15, color: WARN_BORDER })],
    }));
  }
  children.push(new Paragraph({
    border: { top: { style: BorderStyle.SINGLE, size: 4, color: RULE, space: 4 } },
    tabStops: [{ type: TabStopType.RIGHT, position: convertInchesToTwip(6.5) }],
    children: [
      new TextRun({ text: shortTitle, font: FONT, size: 16, color: MUTED }),
      new TextRun({ text: "\t", font: FONT }),
      new TextRun({ text: "Page ", font: FONT, size: 16, color: MUTED }),
      new TextRun({ children: [PageNumber.CURRENT], font: FONT, size: 16, color: MUTED }),
      new TextRun({ text: " of ", font: FONT, size: 16, color: MUTED }),
      new TextRun({ children: [PageNumber.TOTAL_PAGES], font: FONT, size: 16, color: MUTED }),
    ],
  }));
  return new Footer({ children });
}

function buildDocx(mdPath, outPath, shortTitle, draftFooter) {
  const mdText = fs.readFileSync(mdPath, "utf8");
  const parsed = parseDoc(mdText);
  const doc = new Document({
    numbering,
    styles: {
      default: { document: { run: { font: FONT, size: 22 } } }, // 11pt
    },
    sections: [
      {
        properties: {
          page: {
            size: { width: 12240, height: 15840 }, // US Letter
            margin: { top: 1440, bottom: 1440, left: 1440, right: 1440 },
          },
        },
        footers: { default: buildFooter(shortTitle, draftFooter) },
        children: buildBody(parsed),
      },
    ],
  });
  return Packer.toBuffer(doc).then((buf) => {
    fs.mkdirSync(path.dirname(outPath), { recursive: true });
    fs.writeFileSync(outPath, buf);
    console.log("wrote", outPath);
  });
}

const ROOT = "/home/user/Research-NIHOligotox/phase2_submission";
const DOCS = [
  ["01_narrative_document.md", "01_narrative_document.docx", "OligoTox Phase 2 — Narrative Document", true],
  ["02_methodology_document.md", "02_methodology_document.docx", "OligoTox Phase 2 — Methodology Document", false],
  ["03_public_access_and_dissemination_plan.md", "03_public_access_and_dissemination_plan.docx", "OligoTox Phase 2 — Public Access and Dissemination Plan", false],
];

(async () => {
  for (const [src, out, short, draft] of DOCS) {
    await buildDocx(path.join(ROOT, src), path.join(ROOT, "docx", out), short, draft);
  }
})();
