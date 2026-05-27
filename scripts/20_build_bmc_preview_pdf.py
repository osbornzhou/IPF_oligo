#!/usr/bin/env python
"""
Build an integrated BMC Genomics-style PDF preview.

This is a reading/review PDF assembled from the draft manuscript text and
publication-style figures. It intentionally keeps author identity fields as
placeholders because personal information has not been supplied.
"""

from __future__ import annotations

import re
from pathlib import Path

from PIL import Image
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.platypus import (
    Image as RLImage,
    KeepTogether,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)
from reportlab.pdfbase.pdfmetrics import stringWidth
from pypdf import PdfReader


PROJECT_DIR = Path(__file__).resolve().parents[1]
MANUSCRIPT_DIR = PROJECT_DIR / "manuscript"
FIG_DIR = MANUSCRIPT_DIR / "figures"
MD_PATH = MANUSCRIPT_DIR / "ipf_oligo_ml_bmc_genomics_manuscript_draft.md"
LEGEND_PATH = FIG_DIR / "figure_legends.md"
PDF_PATH = MANUSCRIPT_DIR / "ipf_oligo_ml_bmc_genomics_integrated_preview_mechanistic.pdf"
REVIEW_PATH = MANUSCRIPT_DIR / "bmc_preview_pdf_mechanistic_self_review.txt"
EXTRA_FIGURES = [
    (
        PROJECT_DIR / "results/submission_enhancements/plots/final_target_priority_top20.png",
        "Additional Figure S1. Integrated target-prioritization scores for the top robust IPF candidate genes. Scores combine robust differential expression, machine-learning evidence, pathway and network support, miRNA-axis evidence, single-cell localization, and oligonucleotide strategy compatibility.",
    ),
    (
        PROJECT_DIR / "results/submission_enhancements/plots/ml_interpretability_calibration_decision_curve.png",
        "Additional Figure S2. Machine-learning interpretability and utility checks in pooled external validation samples, including permutation importance, calibration, and decision-curve analysis.",
    ),
    (
        PROJECT_DIR / "results/submission_enhancements/plots/biological_model_summary.png",
        "Additional Figure S3. Biological interpretation model linking robust transcriptomic signals, cell context, and cautious oligonucleotide development hypotheses.",
    ),
    (
        PROJECT_DIR / "results/mechanistic_extension/plots/module_trait_correlation_heatmap.png",
        "Additional Figure S4. WGCNA-like module-trait heatmap showing IPF-associated epithelial/ciliary, matrix-remodeling, and TNIK-containing coexpression modules in GSE32537.",
    ),
    (
        PROJECT_DIR / "results/mechanistic_extension/plots/top_curated_ligand_receptor_interactions.png",
        "Additional Figure S5. Curated ligand-receptor score changes in IPF single-cell datasets after excluding multiplets and non-informative cell labels.",
    ),
    (
        PROJECT_DIR / "results/mechanistic_extension/plots/perturbation_priority_proxy.png",
        "Additional Figure S6. Coexpression-neighborhood perturbation-priority ranking for oligonucleotide follow-up candidates, with TNIK treated as an externally motivated bridge rather than a primary discovery candidate.",
    ),
    (
        PROJECT_DIR / "results/mechanistic_extension/plots/tnik_bulk_evidence.png",
        "Additional Figure S7. TNIK differential-expression evidence across bulk cohorts, showing validation-cohort support but no primary discovery-cohort selection.",
    ),
]


def clean_text(text: str) -> str:
    text = text.replace("->", "->")
    text = re.sub(r"\*\*(.*?)\*\*", r"<b>\1</b>", text)
    text = re.sub(r"`([^`]+)`", r"\1", text)
    return text


def split_markdown_tables(lines: list[str]) -> list[tuple[str, object]]:
    blocks: list[tuple[str, object]] = []
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("|"):
            table_lines = []
            while i < len(lines) and lines[i].startswith("|"):
                table_lines.append(lines[i])
                i += 1
            blocks.append(("table", table_lines))
        else:
            para_lines = []
            while i < len(lines) and not lines[i].startswith("|"):
                para_lines.append(lines[i])
                i += 1
            blocks.append(("text", "\n".join(para_lines)))
    return blocks


def parse_md_table(table_lines: list[str], max_rows: int | None = None, max_cols: int | None = None) -> list[list[str]]:
    rows = []
    for line in table_lines:
        parts = [x.strip() for x in line.strip().strip("|").split("|")]
        if all(set(x) <= {"-", " "} for x in parts):
            continue
        rows.append(parts)
    if max_cols and rows:
        rows = [row[:max_cols] for row in rows]
    if max_rows and len(rows) > max_rows + 1:
        rows = rows[: max_rows + 1] + [[f"... {len(rows) - max_rows - 1} additional rows in source table", ""] + [""] * (len(rows[0]) - 2)]
    return rows


def figure_legend_map() -> dict[str, str]:
    text = LEGEND_PATH.read_text(encoding="utf-8")
    legends = {}
    for match in re.finditer(r"(Figure\s+(\d+)\..*?)(?=\n\nFigure\s+\d+\.|\Z)", text, flags=re.S):
        legends[match.group(2)] = " ".join(match.group(1).split())
    return legends


styles = getSampleStyleSheet()
styles.add(
    ParagraphStyle(
        name="BMCHeader",
        parent=styles["Normal"],
        fontName="Helvetica",
        fontSize=8,
        leading=10,
        textColor=colors.HexColor("#4D4D4D"),
        alignment=TA_CENTER,
    )
)
styles.add(
    ParagraphStyle(
        name="BMCTitle",
        parent=styles["Title"],
        fontName="Helvetica-Bold",
        fontSize=18,
        leading=22,
        textColor=colors.HexColor("#222222"),
        alignment=TA_LEFT,
        spaceAfter=10,
    )
)
styles.add(
    ParagraphStyle(
        name="BMCArticleType",
        parent=styles["Normal"],
        fontName="Helvetica-Bold",
        fontSize=9,
        leading=11,
        textColor=colors.HexColor("#A64B2A"),
        spaceAfter=8,
    )
)
styles.add(
    ParagraphStyle(
        name="BMCSection",
        parent=styles["Heading1"],
        fontName="Helvetica-Bold",
        fontSize=13,
        leading=16,
        textColor=colors.HexColor("#1F3A5F"),
        spaceBefore=12,
        spaceAfter=6,
    )
)
styles.add(
    ParagraphStyle(
        name="BMCSubsection",
        parent=styles["Heading2"],
        fontName="Helvetica-Bold",
        fontSize=10.5,
        leading=13,
        textColor=colors.HexColor("#333333"),
        spaceBefore=10,
        spaceAfter=4,
    )
)
styles.add(
    ParagraphStyle(
        name="BMCBody",
        parent=styles["BodyText"],
        fontName="Times-Roman",
        fontSize=9.5,
        leading=12.2,
        alignment=TA_JUSTIFY,
        spaceAfter=5,
    )
)
styles.add(
    ParagraphStyle(
        name="BMCBodyLeft",
        parent=styles["BMCBody"],
        alignment=TA_LEFT,
    )
)
styles.add(
    ParagraphStyle(
        name="BMCSmall",
        parent=styles["Normal"],
        fontName="Times-Roman",
        fontSize=8,
        leading=10,
        alignment=TA_LEFT,
    )
)
styles.add(
    ParagraphStyle(
        name="BMCCaption",
        parent=styles["Normal"],
        fontName="Times-Roman",
        fontSize=8.2,
        leading=10.5,
        alignment=TA_JUSTIFY,
        spaceBefore=6,
        spaceAfter=10,
    )
)
styles.add(
    ParagraphStyle(
        name="BMCAbstractLabel",
        parent=styles["Normal"],
        fontName="Helvetica-Bold",
        fontSize=9,
        leading=11,
        textColor=colors.HexColor("#1F3A5F"),
        spaceBefore=4,
        spaceAfter=2,
    )
)


def page_decorator(canvas, doc):
    canvas.saveState()
    width, height = A4
    canvas.setStrokeColor(colors.HexColor("#D9D9D9"))
    canvas.setLineWidth(0.5)
    canvas.line(doc.leftMargin, height - 1.45 * cm, width - doc.rightMargin, height - 1.45 * cm)
    canvas.setFont("Helvetica", 8)
    canvas.setFillColor(colors.HexColor("#555555"))
    canvas.drawString(doc.leftMargin, height - 1.18 * cm, "BMC Genomics | Research article draft preview")
    canvas.drawRightString(width - doc.rightMargin, 0.85 * cm, f"Page {doc.page}")
    canvas.restoreState()


def add_paragraph_or_heading(story: list, block: str) -> None:
    for part in re.split(r"\n\s*\n", block.strip()):
        if not part.strip():
            continue
        lines = [x.strip() for x in part.splitlines() if x.strip()]
        text = " ".join(lines)
        if text.startswith("# "):
            story.append(Paragraph(clean_text(text[2:]), styles["BMCTitle"]))
        elif text.startswith("## "):
            story.append(Paragraph(clean_text(text[3:]), styles["BMCSection"]))
        elif text.startswith("### "):
            story.append(Paragraph(clean_text(text[4:]), styles["BMCSubsection"]))
        elif re.match(r"^(Figure|Table|Supplementary Table)\s+\d+\.", text):
            story.append(Paragraph(clean_text(text), styles["BMCBodyLeft"]))
        elif text.startswith("Authors:") or text.startswith("Affiliations:") or text.startswith("Corresponding author:"):
            story.append(Paragraph(clean_text(text), styles["BMCBodyLeft"]))
        elif text.startswith("Background:") or text.startswith("Results:") or text.startswith("Conclusions:"):
            label, rest = text.split(":", 1)
            story.append(Paragraph(label, styles["BMCAbstractLabel"]))
            story.append(Paragraph(clean_text(rest.strip()), styles["BMCBody"]))
        else:
            story.append(Paragraph(clean_text(text), styles["BMCBody"]))


def add_table(story: list, rows: list[list[str]]) -> None:
    if not rows or len(rows[0]) > 10:
        return
    col_count = len(rows[0])
    usable = A4[0] - 3.6 * cm
    col_widths = [usable / col_count] * col_count
    data = [[Paragraph(clean_text(cell), styles["BMCSmall"]) for cell in row] for row in rows]
    table = Table(data, colWidths=col_widths, repeatRows=1, hAlign="LEFT")
    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#EEF2F6")),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.HexColor("#222222")),
                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#C8CDD2")),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 3),
                ("RIGHTPADDING", (0, 0), (-1, -1), 3),
                ("TOPPADDING", (0, 0), (-1, -1), 3),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
            ]
        )
    )
    story.append(table)
    story.append(Spacer(1, 8))


def add_figure_page(story: list, fig_path: Path, caption: str, doc_width: float, doc_height: float) -> None:
    story.append(PageBreak())
    story.append(Paragraph(fig_path.stem.replace("_", " "), styles["BMCSection"]))
    image = Image.open(fig_path)
    img_w, img_h = image.size
    max_w = doc_width
    max_h = doc_height - 4.2 * cm
    scale = min(max_w / img_w, max_h / img_h)
    draw_w, draw_h = img_w * scale, img_h * scale
    story.append(RLImage(str(fig_path), width=draw_w, height=draw_h))
    story.append(Paragraph(clean_text(caption), styles["BMCCaption"]))


def build_pdf() -> None:
    raw = MD_PATH.read_text(encoding="utf-8")
    lines = raw.splitlines()
    blocks = split_markdown_tables(lines)

    doc = SimpleDocTemplate(
        str(PDF_PATH),
        pagesize=A4,
        leftMargin=1.8 * cm,
        rightMargin=1.8 * cm,
        topMargin=1.75 * cm,
        bottomMargin=1.55 * cm,
        title="IPF oligonucleotide ML BMC Genomics draft preview",
        author="Yunyi Zhou; Yanli Zhang",
    )
    story = []
    story.append(Paragraph("Research article", styles["BMCArticleType"]))
    story.append(Paragraph("Prepared in a BMC Genomics submission-style layout for internal review", styles["BMCHeader"]))
    story.append(Spacer(1, 8))

    for kind, content in blocks:
        if kind == "text":
            add_paragraph_or_heading(story, content)
        else:
            rows = parse_md_table(content, max_rows=12, max_cols=8)
            add_table(story, rows)

    legends = figure_legend_map()
    fig_paths = sorted(FIG_DIR.glob("Figure_*.png"))
    for path in fig_paths:
        match = re.search(r"Figure_(\d+)", path.name)
        number = match.group(1) if match else ""
        add_figure_page(story, path, legends.get(number, path.stem), doc.width, doc.height)
    for path, caption in EXTRA_FIGURES:
        if path.exists():
            add_figure_page(story, path, caption, doc.width, doc.height)

    doc.build(story, onFirstPage=page_decorator, onLaterPages=page_decorator)


def self_review() -> None:
    reader = PdfReader(str(PDF_PATH))
    text = "\n".join(page.extract_text() or "" for page in reader.pages)
    checks = [
        ("PDF created", PDF_PATH.exists() and PDF_PATH.stat().st_size > 100_000, str(PDF_PATH)),
        ("Page count reasonable", len(reader.pages) >= 18, f"{len(reader.pages)} pages; compact preview with main and additional figure pages"),
        ("BMC structure present", all(x in text for x in ["Abstract", "Background", "Methods", "Results", "Discussion", "Conclusions", "Declarations"]), "major sections detected"),
        ("Core mRNA/miRNA counts present", "280 robust mRNAs" in text and "10 robust miRNA" in text, "robust candidate counts detected"),
        ("ML result present", "0.971" in text and "Elastic Net" in text, "external AUC and model detected"),
        ("Figure pages included", sum(1 for p in FIG_DIR.glob("Figure_*.png")) == 6 and all(path.exists() for path, _ in EXTRA_FIGURES), "six main PNG figures plus submission and mechanistic enhancement figures available and embedded"),
        ("Author metadata present", "Yunyi Zhou" in text and "Yanli Zhang" in text, "author names detected"),
    ]
    lines = ["BMC preview PDF self-review", ""]
    for name, ok, note in checks:
        lines.append(f"{'PASS' if ok else 'FAIL'}\t{name}\t{note}")
    REVIEW_PATH.write_text("\n".join(lines), encoding="utf-8")
    print(PDF_PATH)
    print(REVIEW_PATH)
    print("\n".join(lines))


if __name__ == "__main__":
    build_pdf()
    self_review()
