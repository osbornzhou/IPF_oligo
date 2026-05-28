#!/usr/bin/env python
"""
Build editable SVG source files for all six main manuscript figures.

The submitted DOCX keeps rendered PNGs for stable journal layout. These SVG
sources are for editing: text, bars, boxes, arrows, labels, heatmap cells, and
plot annotations are written as separate SVG elements instead of flattened PNGs.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from xml.sax.saxutils import escape

import numpy as np
import pandas as pd

PROJECT_DIR = Path(__file__).resolve().parents[1]
EDIT_DIR = PROJECT_DIR / "manuscript" / "editable_figures"
SVG_DIR = EDIT_DIR / "svg"
for directory in (EDIT_DIR, SVG_DIR):
    directory.mkdir(parents=True, exist_ok=True)

COL_UP = "#B55D60"
COL_DOWN = "#4C78A8"
COL_CTRL = "#5B84A4"
COL_ACCENT = "#3D7F6F"
COL_NEUTRAL = "#B8B8B8"
COL_TEXT = "#222222"
COL_GRID = "#E6EBF2"

MANIFEST_ROWS: list[dict[str, str]] = []


def read_csv(rel_path: str) -> pd.DataFrame:
    return pd.read_csv(PROJECT_DIR / rel_path)


def clean(value: object) -> str:
    return escape("" if pd.isna(value) else str(value))


def lerp(a: float, b: float, t: float) -> float:
    return a + (b - a) * max(0.0, min(1.0, t))


def hex_to_rgb(hex_color: str) -> tuple[int, int, int]:
    h = hex_color.strip("#")
    return int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)


def rgb_to_hex(rgb: tuple[float, float, float]) -> str:
    return "#" + "".join(f"{max(0, min(255, int(round(v)))):02X}" for v in rgb)


def mix(c1: str, c2: str, t: float) -> str:
    r1, g1, b1 = hex_to_rgb(c1)
    r2, g2, b2 = hex_to_rgb(c2)
    return rgb_to_hex((lerp(r1, r2, t), lerp(g1, g2, t), lerp(b1, b2, t)))


def div_color(value: float, vmax: float) -> str:
    if not np.isfinite(value):
        return "#FFFFFF"
    vmax = max(vmax, 1e-9)
    if value >= 0:
        return mix("#F7F7F7", COL_UP, min(value / vmax, 1))
    return mix("#F7F7F7", COL_DOWN, min(abs(value) / vmax, 1))


def xmap(value: float, vmin: float, vmax: float, x0: float, x1: float) -> float:
    if vmax == vmin:
        return (x0 + x1) / 2
    return x0 + (value - vmin) / (vmax - vmin) * (x1 - x0)


def ymap(value: float, vmin: float, vmax: float, y0: float, y1: float) -> float:
    if vmax == vmin:
        return (y0 + y1) / 2
    return y1 - (value - vmin) / (vmax - vmin) * (y1 - y0)


class SVG:
    def __init__(self, width: int, height: int, title: str):
        self.width = width
        self.height = height
        self.title = title
        self.parts: list[str] = [
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
            f"<title>{clean(title)}</title>",
            '<rect x="0" y="0" width="100%" height="100%" fill="white"/>',
            '<style>text{font-family:Arial,DejaVu Sans,sans-serif;} .small{fill:#333;} .axis{stroke:#333;stroke-width:1.2;} .grid{stroke:#E6EBF2;stroke-width:1;} .panel{font-weight:700;font-size:30px;}</style>',
        ]

    def rect(self, x, y, w, h, fill="none", stroke="none", sw=1, rx=0, opacity=1, cls=""):
        self.parts.append(
            f'<rect class="{cls}" x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" rx="{rx}" ry="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}" opacity="{opacity}"/>'
        )

    def line(self, x1, y1, x2, y2, stroke="#333333", sw=1, dash="", marker_end=False):
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        marker = ' marker-end="url(#arrow)"' if marker_end else ""
        self.parts.append(f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" stroke="{stroke}" stroke-width="{sw}"{dash_attr}{marker}/>')

    def circle(self, x, y, r, fill, stroke="none", sw=0, opacity=1):
        self.parts.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{r:.2f}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}" opacity="{opacity}"/>')

    def text(self, x, y, text, size=18, weight="400", anchor="start", fill=COL_TEXT, rotate=0, leading=1.15):
        transform = f' transform="rotate({rotate:.2f} {x:.2f} {y:.2f})"' if rotate else ""
        lines = str(text).split("\n")
        self.parts.append(f'<text x="{x:.2f}" y="{y:.2f}" font-size="{size}" font-weight="{weight}" text-anchor="{anchor}" fill="{fill}"{transform}>')
        for i, line in enumerate(lines):
            dy = "0" if i == 0 else f"{size * leading:.2f}"
            self.parts.append(f'<tspan x="{x:.2f}" dy="{dy}">{clean(line)}</tspan>')
        self.parts.append("</text>")

    def path(self, d, fill="none", stroke="#333333", sw=1, opacity=1):
        self.parts.append(f'<path d="{d}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}" opacity="{opacity}"/>')

    def arrow_defs(self):
        self.parts.append(
            '<defs><marker id="arrow" markerWidth="10" markerHeight="8" refX="9" refY="4" orient="auto" markerUnits="strokeWidth"><path d="M0,0 L10,4 L0,8 Z" fill="#4A5568"/></marker></defs>'
        )

    def save(self, path: Path):
        self.parts.append("</svg>")
        path.write_text("\n".join(self.parts), encoding="utf-8")


def save_svg(svg: SVG, stem: str, note: str) -> None:
    svg.save(SVG_DIR / f"{stem}.svg")
    MANIFEST_ROWS.append(
        {
            "figure": stem.replace("_editable", ""),
            "editable_svg": f"svg/{stem}.svg",
            "editable_elements": note,
        }
    )


def make_figure1() -> None:
    svg = SVG(1800, 900, "Figure 1 editable workflow")
    svg.arrow_defs()
    svg.text(900, 90, "Boundary-tested perturbation-triage study design", 42, "700", "middle")
    boxes = [
        (90, 215, 260, 140, "Public GEO input\n8 bulk/miRNA datasets\n2 scRNA-seq datasets"),
        (500, 215, 270, 140, "Sample-label audit\nand matrix QC\n691 profiled samples"),
        (930, 215, 330, 140, "Reproducible transcriptomic\nevidence base\n280 mRNAs / 10 miRNAs"),
        (130, 570, 330, 145, "Disease-boundary\nstress tests\nfibrotic ILD state score"),
        (560, 570, 330, 145, "miRNA evidence gates\n3 exact hsa-miR-375 axes\nweak axes downgraded"),
        (990, 570, 340, 145, "Donor-aware single-cell\ncontext\npseudobulk localization"),
        (1410, 570, 330, 145, "Perturbation triage\nscreening, markers,\nrestoration and axes"),
    ]
    for x, y, w, h, label in boxes:
        svg.rect(x, y, w, h, "#F6F8FA", "#4A5568", 2.2, rx=12)
        svg.text(x + w / 2, y + 48, label, 20, "400", "middle")
    svg.line(350, 285, 500, 285, "#4A5568", 3, marker_end=True)
    svg.line(770, 285, 930, 285, "#4A5568", 3, marker_end=True)
    for x2 in [295, 725, 1160, 1575]:
        svg.line(1095, 355, x2, 570, "#4A5568", 2.2, marker_end=True)
    svg.text(900, 830, "All downstream gates use QC-passed datasets and preserve external validation cohorts for disease-state boundary testing.", 19, "400", "middle", "#444444")
    save_svg(svg, "Figure_1_workflow_editable", "editable text boxes, rounded rectangles, arrows, title, and footnote")


def draw_axis(svg: SVG, x0, y0, x1, y1, x_label="", y_label="", xticks=None, yticks=None, xlim=(0, 1), ylim=(0, 1)):
    svg.line(x0, y1, x1, y1, "#333333", 2)
    svg.line(x0, y0, x0, y1, "#333333", 2)
    if xticks:
        for tick, label in xticks:
            x = xmap(tick, xlim[0], xlim[1], x0, x1)
            svg.line(x, y1, x, y1 + 10, "#333333", 2)
            svg.text(x, y1 + 34, label, 16, "400", "middle")
    if yticks:
        for tick, label in yticks:
            y = ymap(tick, ylim[0], ylim[1], y0, y1)
            svg.line(x0 - 10, y, x0, y, "#333333", 2)
            svg.text(x0 - 16, y + 5, label, 16, "400", "end")
            svg.line(x0, y, x1, y, COL_GRID, 1)
    if x_label:
        svg.text((x0 + x1) / 2, y1 + 72, x_label, 22, "400", "middle")
    if y_label:
        svg.text(x0 - 70, (y0 + y1) / 2, y_label, 22, "400", "middle", rotate=-90)


def make_figure2() -> None:
    svg = SVG(2200, 1650, "Figure 2 editable transcriptomic evidence base")
    de_qc = read_csv("results/differential_expression/differential_expression_qc.csv")
    mrna_de = read_csv("results/differential_expression_annotated/GSE32537_IPF_vs_Control_de_results_annotated_gene_or_mirna_level.csv")
    mrna_dir = read_csv("results/robust_candidates/mrna_direction_matrix.csv")
    mirna_dir = read_csv("results/robust_candidates/mirna_direction_matrix.csv")
    robust_mrna = read_csv("results/robust_candidates/robust_mrna_candidates_strict.csv")
    robust_mirna = read_csv("results/robust_candidates/robust_mirna_candidates_strict.csv")

    svg.text(55, 60, "A", 32, "700")
    svg.text(1125, 60, "B", 32, "700")
    svg.text(55, 880, "C", 32, "700")
    svg.text(1125, 880, "D", 32, "700")

    # Panel A
    x0, y0, x1, y1 = 150, 110, 950, 720
    vals = de_qc["significant_fdr_0_05_logfc_1"].astype(float).to_numpy()
    labels = (de_qc["series_id"] + "\n" + de_qc["data_type"].str.replace("bulk ", "", regex=False)).tolist()
    draw_axis(svg, x0, y0, x1, y1, y_label="Significant features\nFDR < 0.05 and |logFC| >= 1", xticks=[], yticks=[(0, "0"), (1000, "1000"), (2000, "2000"), (3000, "3000"), (4000, "4000"), (5000, "5000")], xlim=(0, len(vals)), ylim=(0, 5500))
    svg.text((x0 + x1) / 2, 72, "Differential-expression yield", 28, "400", "middle")
    bw = (x1 - x0) / len(vals) * 0.78
    for i, (val, label, dtype) in enumerate(zip(vals, labels, de_qc["data_type"])):
        cx = x0 + (i + 0.5) * (x1 - x0) / len(vals)
        y = ymap(val, 0, 5500, y0, y1)
        svg.rect(cx - bw / 2, y, bw, y1 - y, COL_ACCENT if dtype == "miRNA" else COL_CTRL)
        svg.text(cx, y - 8, f"{int(val)}", 16, "400", "middle")
        svg.text(cx - 8, y1 + 48, label, 15, "400", "end", rotate=-50)

    # Panel B
    x0, y0, x1, y1 = 1220, 110, 2050, 720
    work = mrna_de.copy()
    work["adj.P.Val"] = pd.to_numeric(work["adj.P.Val"], errors="coerce").clip(lower=1e-300)
    work["logFC"] = pd.to_numeric(work["logFC"], errors="coerce")
    work["neglog10_fdr"] = -np.log10(work["adj.P.Val"])
    draw_axis(svg, x0, y0, x1, y1, "log2 fold change", "-log10 FDR", xticks=[(-3, "-3"), (-2, "-2"), (-1, "-1"), (0, "0"), (1, "1"), (2, "2"), (3, "3")], yticks=[(0, "0"), (10, "10"), (20, "20"), (30, "30"), (40, "40"), (50, "50")], xlim=(-3.1, 3.6), ylim=(-3, 55))
    svg.text((x0 + x1) / 2, 72, "Discovery mRNA volcano (GSE32537)", 28, "400", "middle")
    for xv in [-1, 1]:
        x = xmap(xv, -3.1, 3.6, x0, x1)
        svg.line(x, y0, x, y1, "#777777", 1.5, dash="6 5")
    ythr = ymap(-math.log10(0.05), -3, 55, y0, y1)
    svg.line(x0, ythr, x1, ythr, "#777777", 1.5, dash="6 5")
    for _, row in work.iterrows():
        lfc = row["logFC"]
        yy = row["neglog10_fdr"]
        if not np.isfinite(lfc) or not np.isfinite(yy):
            continue
        color = COL_NEUTRAL
        opacity = 0.35
        if row["adj.P.Val"] < 0.05 and lfc >= 1:
            color, opacity = COL_UP, 0.70
        elif row["adj.P.Val"] < 0.05 and lfc <= -1:
            color, opacity = COL_DOWN, 0.70
        svg.circle(xmap(lfc, -3.1, 3.6, x0, x1), ymap(yy, -3, 55, y0, y1), 3.3, color, opacity=opacity)
    for _, row in work.sort_values("adj.P.Val").head(8).iterrows():
        svg.text(xmap(row["logFC"], -3.1, 3.6, x0, x1), ymap(row["neglog10_fdr"], -3, 55, y0, y1) - 8, row.get("gene_symbol", row.get("standard_feature_id", "")), 15, "400", "middle")

    # Panel C/D heatmaps
    def heatmap(panel_x, panel_y, data, row_names, col_names, title, cell_w, cell_h, max_abs, label_size=15):
        svg.text(panel_x + cell_w * len(col_names) / 2 + 120, panel_y - 45, title, 26, "400", "middle")
        for i, row_name in enumerate(row_names):
            svg.text(panel_x - 12, panel_y + i * cell_h + cell_h * 0.62, row_name, label_size, "400", "end")
            for j, col_name in enumerate(col_names):
                val = data[i, j]
                svg.rect(panel_x + j * cell_w, panel_y + i * cell_h, cell_w, cell_h, div_color(val, max_abs), "#F1F1F1", 0.8)
                if np.isfinite(val):
                    svg.text(panel_x + j * cell_w + cell_w / 2, panel_y + i * cell_h + cell_h * 0.62, f"{val:.1f}", 13, "400", "middle")
        for j, col_name in enumerate(col_names):
            svg.text(panel_x + j * cell_w + cell_w * 0.35, panel_y + len(row_names) * cell_h + 50, col_name, 15, "400", "end", rotate=-45)
        svg.text(panel_x + cell_w * len(col_names) / 2, panel_y + len(row_names) * cell_h + 105, "Dataset", 20, "400", "middle")

    mrna_sets = ["GSE32537", "GSE110147", "GSE150910", "GSE53845", "GSE92592"]
    mirna_sets = ["GSE32538", "GSE21394", "GSE27430"]
    mrna_order = robust_mrna.head(20)["standard_feature_id"].tolist()
    mirna_order = robust_mirna.head(10)["standard_feature_id"].tolist()
    mrna_dir = mrna_dir[mrna_dir["standard_feature_id"].isin(mrna_order)].copy()
    mirna_dir = mirna_dir[mirna_dir["standard_feature_id"].isin(mirna_order)].copy()
    mrna_dir["order"] = mrna_dir["standard_feature_id"].map({g: i for i, g in enumerate(mrna_order)})
    mirna_dir["order"] = mirna_dir["standard_feature_id"].map({g: i for i, g in enumerate(mirna_order)})
    mrna_dir = mrna_dir.sort_values("order")
    mirna_dir = mirna_dir.sort_values("order")
    mat_m = np.array([[float(row.get(f"{ds}_logFC", np.nan)) for ds in mrna_sets] for _, row in mrna_dir.iterrows()])
    mat_i = np.array([[float(row.get(f"{ds}_logFC", np.nan)) for ds in mirna_sets] for _, row in mirna_dir.iterrows()])
    vmax = max(np.nanmax(np.abs(mat_m)), np.nanmax(np.abs(mat_i)), 1)
    heatmap(150, 960, mat_m, mrna_dir["standard_feature_id"].tolist(), mrna_sets, "Top robust mRNAs: cross-cohort logFC", 150, 31, vmax)
    heatmap(1240, 960, mat_i, mirna_dir["standard_feature_id"].tolist(), mirna_sets, "Robust miRNAs: cross-cohort logFC", 250, 62, vmax, 16)
    # Editable colorbar
    cbx, cby, cbw, cbh = 2045, 955, 30, 620
    for k in range(80):
        value = vmax - 2 * vmax * k / 79
        svg.rect(cbx, cby + k * cbh / 80, cbw, cbh / 80 + 0.5, div_color(value, vmax))
    svg.rect(cbx, cby, cbw, cbh, "none", "#333333", 1)
    svg.text(cbx + 60, cby + cbh / 2, "log2 fold change", 18, "400", "middle", rotate=-90)
    save_svg(svg, "Figure_2_transcriptomic_evidence_base_editable", "editable bars, volcano points, heatmap cells, values, labels, axes, and colorbar")


def make_figure3() -> None:
    svg = SVG(1800, 1050, "Figure 3 editable disease-boundary stress tests")
    svg.text(50, 70, "Figure 3. Disease-boundary stress tests for the fibrotic ILD state score", 38, "700")
    svg.text(90, 170, "A", 30, "700")
    svg.text(1020, 170, "B", 30, "700")
    svg.text(130, 170, "Fibrotic-disease-state score distribution", 26)
    svg.text(1065, 170, "Matched random-panel baseline", 26)
    scores = read_csv("results/ml_stress_tests/gse110147_excluded_ild_disease_control_scores.csv")
    panels = read_csv("results/ml_stress_tests/matched_random_discovery_feature_panel_baseline.csv")
    observed = float(panels[panels["panel_type"].eq("observed_final_panel_regularized_refit")]["mean_roc_auc"].iloc[0])
    random_vals = panels[panels["panel_type"].eq("matched_random_discovery_feature_panel")]["mean_roc_auc"].astype(float).to_numpy()
    p95 = float(np.percentile(random_vals, 95))
    groups = ["Normal control", "NSIP", "Mixed IPF-NSIP", "IPF"]
    colors = {"Normal control": "#7A8793", "NSIP": "#B07A55", "Mixed IPF-NSIP": "#B07A55", "IPF": COL_UP}
    x0, y0, x1, y1 = 110, 220, 840, 930
    draw_axis(svg, x0, y0, x1, y1, y_label="Score", yticks=[(0, "0.0"), (0.2, "0.2"), (0.4, "0.4"), (0.6, "0.6"), (0.8, "0.8"), (1.0, "1.0")], ylim=(0, 1.05))
    rng = np.random.default_rng(12)
    for i, group in enumerate(groups):
        vals = scores.loc[scores["stress_group"].eq(group), "refit_final_panel_score"].astype(float).to_numpy()
        cx = x0 + (i + 0.5) * (x1 - x0) / len(groups)
        q1, med, q3 = np.percentile(vals, [25, 50, 75])
        ymin, ymax = vals.min(), vals.max()
        svg.rect(cx - 35, ymap(q3, 0, 1.05, y0, y1), 70, ymap(q1, 0, 1.05, y0, y1) - ymap(q3, 0, 1.05, y0, y1), "white", colors[group], 3)
        svg.line(cx - 42, ymap(med, 0, 1.05, y0, y1), cx + 42, ymap(med, 0, 1.05, y0, y1), colors[group], 4)
        svg.line(cx, ymap(ymin, 0, 1.05, y0, y1), cx, ymap(ymax, 0, 1.05, y0, y1), colors[group], 2)
        for val in vals:
            svg.circle(cx + rng.normal(0, 12), ymap(val, 0, 1.05, y0, y1), 5, colors[group], opacity=0.9)
        svg.text(cx, 205, f"n={len(vals)}", 16, "400", "middle")
        svg.text(cx, 975, group.replace(" ", "\n", 1), 17, "400", "middle")
    hx0, hy0, hx1, hy1 = 1070, 240, 1700, 920
    bins = np.linspace(0.88, 0.97, 18)
    counts, edges = np.histogram(random_vals, bins=bins)
    ymax = max(counts) * 1.1
    for c, left, right in zip(counts, edges[:-1], edges[1:]):
        x = xmap(left, 0.86, 1.0, hx0, hx1)
        w = xmap(right, 0.86, 1.0, hx0, hx1) - x
        h = hy1 - ymap(c, 0, ymax, hy0, hy1)
        svg.rect(x, hy1 - h, w - 2, h, "#B9CBDB", "white", 1)
    svg.line(xmap(observed, 0.86, 1.0, hx0, hx1), hy0, xmap(observed, 0.86, 1.0, hx0, hx1), hy1, COL_UP, 4)
    svg.line(xmap(p95, 0.86, 1.0, hx0, hx1), hy0, xmap(p95, 0.86, 1.0, hx0, hx1), hy1, COL_DOWN, 4)
    for tick in [0.86, 0.90, 0.93, 0.96, 1.00]:
        x = xmap(tick, 0.86, 1.0, hx0, hx1)
        svg.text(x, hy1 + 35, f"{tick:.2f}", 16, "400", "middle")
    svg.text((hx0 + hx1) / 2, hy1 + 80, "Mean external ROC AUC", 22, "400", "middle")
    svg.text(hx0 + 15, hy0 + 40, f"Observed refit={observed:.3f}", 21, "400", fill=COL_UP)
    svg.text(hx0 + 15, hy0 + 80, f"Matched random 95th={p95:.3f}", 21, "400", fill=COL_DOWN)
    save_svg(svg, "Figure_3_disease_boundary_stress_tests_editable", "editable boxplots, points, histogram bars, reference lines, titles, labels, and annotations")


def make_figure4() -> None:
    svg = SVG(1800, 900, "Figure 4 editable miRNA evidence-gate stress tests")
    svg.text(50, 70, "Figure 4. miRNA evidence-gate stress tests", 38, "700")
    svg.text(90, 165, "A", 30, "700")
    svg.text(1070, 165, "B", 30, "700")
    svg.text(130, 165, "Robust downregulated miRNAs", 26)
    svg.text(1110, 165, "hsa-miR-375 target release-like score", 26)
    enrich = read_csv("results/mirna_program_support/robust_mirna_target_set_enrichment.csv")
    rel = read_csv("results/mirna_program_support/hsa_mir_375_target_repression_release_score.csv").iloc[0]
    work = enrich[enrich["target_source"].eq("miRTarBase_arm_recoverable")].copy()
    keep = ["hsa-miR-141", "hsa-miR-30a", "hsa-miR-30d", "hsa-miR-375", "hsa-miR-92a", "hsa-miR-423-5p", "hsa-miR-203"]
    work["order"] = work["candidate_mirna"].map({m: i for i, m in enumerate(keep)})
    work = work[work["candidate_mirna"].isin(keep)].sort_values("order")
    x0, y0, x1, y1 = 320, 210, 780, 760
    svg.line(xmap(1, 0, 1.8, x0, x1), y0 - 20, xmap(1, 0, 1.8, x0, x1), y1 + 20, "#999999", 2)
    svg.text(xmap(1, 0, 1.8, x0, x1), y1 + 55, "OR=1", 17, "400", "middle", "#666666")
    for i, (_, row) in enumerate(work.iterrows()):
        y = y0 + i * 78
        mir = row["candidate_mirna"]
        val = row["odds_ratio"]
        svg.text(95, y + 18, mir, 22)
        if pd.notna(val):
            width = xmap(float(val), 0, 1.8, x0, x1) - x0
            svg.rect(x0, y - 10, width, 42, COL_UP if mir == "hsa-miR-375" else COL_DOWN)
            svg.text(x0 + width + 12, y + 17, f"OR={float(val):.2g}; FDR={row['fdr_bh']:.0f}", 18)
        else:
            svg.text(x0 + 15, y + 17, "OR=NA; FDR=NA", 18)
    svg.text((x0 + x1) / 2, 850, "Odds ratio for target overlap", 22, "400", "middle")
    bx0, by0, bx1, by1 = 1160, 500, 1700, 640
    means = [float(rel["mean_target_discovery_logFC"]), float(rel["mean_non_target_discovery_logFC"])]
    svg.line(bx0, by1, bx1, by1, "#999999", 2)
    for i, (label, val, col) in enumerate(zip(["Targets", "Non-targets"], means, [COL_UP, "#7A8793"])):
        cx = bx0 + 130 + i * 300
        h = max(3, abs(val) * 1400)
        svg.rect(cx - 55, by1 - h, 110, h, col)
        svg.text(cx, by1 - h - 15, f"{val:.2f}", 19, "400", "middle")
        svg.text(cx, by1 + 95, label, 21, "400", "middle")
    svg.text(bx0 + 20, 250, f"Permutation p={float(rel['one_sided_permutation_p_target_greater_than_background']):.3f}", 21)
    svg.text(bx0 + 20, 300, f"n targets={int(rel['targets_in_discovery_mrna_background'])}", 21)
    svg.text((bx0 + bx1) / 2, 830, "Mean discovery mRNA logFC", 22, "400", "middle")
    save_svg(svg, "Figure_4_mirna_evidence_gate_stress_tests_editable", "editable bars, reference line, target/non-target columns, labels, titles, and annotations")


def make_figure5() -> None:
    svg = SVG(1700, 2200, "Figure 5 editable single-cell pseudobulk localization")
    svg.text(40, 55, "A", 30, "700")
    svg.text(925, 55, "B", 30, "700")
    svg.text(500, 980, "C  Donor-aware pseudobulk validation", 40, "700")
    delta = read_csv("results/single_cell_validation/single_cell_clean_ipf_control_delta_by_celltype.csv")
    top = read_csv("results/single_cell_validation/single_cell_clean_top_celltype_gene_changes.csv").head(12)
    pseudo = read_csv("results/single_cell_pseudobulk/single_cell_pseudobulk_core_candidate_summary.csv").sort_values("ipf_minus_control_log1p_norm", ascending=True)
    genes = ["SPP1", "COL1A1", "COL3A1", "POSTN", "GPX3", "COL14A1", "THY1", "TPPP3", "CD24", "ASPN"]
    celltypes = [("GSE135893", "Immune", "Macrophages"), ("GSE136831", "Myeloid", "Macrophage"), ("GSE136831", "Stromal", "Myofibroblast"), ("GSE135893", "Mesenchymal", "Myofibroblasts")]
    rows = []
    for sid, broad, fine in celltypes:
        sub = delta[(delta.series_id.eq(sid)) & (delta.broad_celltype.eq(broad)) & (delta.fine_celltype.eq(fine))]
        rows.append([float(sub.loc[sub.gene_symbol.eq(g), "ipf_minus_control_log1p_mean_norm"].iloc[0]) if (sub.gene_symbol.eq(g)).any() else np.nan for g in genes])
    vmax = 2.1
    x0, y0, cw, ch = 120, 110, 50, 105
    svg.text(x0 + cw * len(genes) / 2, 90, "Single-cell IPF-control expression differences", 20, "400", "middle")
    for i, row in enumerate(rows):
        label = f"{celltypes[i][0]}\n{celltypes[i][2]}"
        svg.text(x0 - 18, y0 + i * ch + 50, label, 14, "400", "end")
        for j, val in enumerate(row):
            svg.rect(x0 + j * cw, y0 + i * ch, cw, ch, div_color(val, vmax), "#F1F1F1", 0.8)
            if np.isfinite(val):
                svg.text(x0 + j * cw + cw / 2, y0 + i * ch + 56, f"{val:.1f}", 12, "400", "middle")
    for j, g in enumerate(genes):
        svg.text(x0 + j * cw + 10, y0 + len(rows) * ch + 40, g, 13, "400", "end", rotate=-45)
    # Panel B bars
    bx0, by0, bx1, by1 = 1030, 110, 1590, 820
    svg.line(xmap(0, -2.1, 2.3, bx0, bx1), by0, xmap(0, -2.1, 2.3, bx0, bx1), by1, "#333333", 2)
    svg.text((bx0 + bx1) / 2, 90, "Largest validated cell-type candidate shifts", 20, "400", "middle")
    plot = top.iloc[::-1].copy()
    for i, (_, row) in enumerate(plot.iterrows()):
        val = float(row["ipf_minus_control_log1p_mean_norm"])
        y = by0 + i * 55
        label = f"{str(row['series_id']).replace('GSE','')} | {row['fine_celltype']} | {row['gene_symbol']}".replace("Vascular-Endothelial_A", "Vasc-Endo A")
        svg.text(bx0 - 15, y + 16, label, 13, "400", "end")
        xz = xmap(0, -2.1, 2.3, bx0, bx1)
        xv = xmap(val, -2.1, 2.3, bx0, bx1)
        svg.rect(min(xz, xv), y - 7, abs(xv - xz), 34, COL_UP if val >= 0 else COL_DOWN)
    svg.text((bx0 + bx1) / 2, 870, "IPF - control log1p mean", 16, "400", "middle")
    # Panel C
    cx0, cy0, cx1, cy1 = 570, 1110, 1320, 2060
    svg.text(1160, 1045, "n=IPF/control donors represented", 18, "400", fill="#666666")
    svg.rect(1160, 1072, 24, 24, COL_UP)
    svg.text(1192, 1092, "IPF-increased", 18, fill="#666666")
    svg.rect(1370, 1072, 24, 24, COL_DOWN)
    svg.text(1402, 1092, "IPF-decreased", 18, fill="#666666")
    for tick in [-1, -0.5, 0, 0.5, 1, 1.5]:
        x = xmap(tick, -1.25, 1.9, cx0, cx1)
        svg.line(x, cy0, x, cy1, COL_GRID, 1.2)
        svg.text(x, cy1 + 35, f"{tick:g}", 18, "400", "middle")
    zero = xmap(0, -1.25, 1.9, cx0, cx1)
    svg.line(zero, cy0, zero, cy1, "#666666", 3)
    row_h = (cy1 - cy0) / len(pseudo)
    for i, (_, row) in enumerate(pseudo.iterrows()):
        y = cy0 + i * row_h + row_h * 0.5
        val = float(row["ipf_minus_control_log1p_norm"])
        xv = xmap(val, -1.25, 1.9, cx0, cx1)
        label = f"{row['gene_symbol']} ({row['broad_celltype']}, {row['series_id']})"
        fdr = float(row["fdr_bh"])
        fdr_label = f"{fdr:.1e}" if fdr < 1e-3 else f"{fdr:.3f}"
        svg.text(35, y + 8, label, 24)
        svg.rect(min(zero, xv), y - 15, abs(xv - zero), 30, COL_UP if val >= 0 else COL_DOWN)
        svg.text(xv + (16 if val >= 0 else -16), y + 8, f"n={int(row['ipf_donors'])}/{int(row['control_donors'])}; FDR={fdr_label}", 18, "400", "start" if val >= 0 else "end")
    svg.text((cx0 + cx1) / 2, 2150, "IPF-control donor-level pseudobulk log1p normalized expression difference", 18, "400", "middle")
    save_svg(svg, "Figure_5_single_cell_pseudobulk_editable", "editable heatmap cells, cell-level bars, pseudobulk bars, donor/FDR labels, axes, and legends")


def make_figure6() -> None:
    svg = SVG(1600, 1600, "Figure 6 editable perturbation-triage map")
    table = read_csv("results/oligonucleotide_actionability/oligonucleotide_actionability_index.csv")
    order = {"Knockdown-screening candidate": 0, "Context-dependent candidate": 1, "Restoration/pathway marker": 2, "miRNA-axis hypothesis": 3, "External TNIK bridge": 4}
    colors = {"Knockdown-screening candidate": COL_UP, "Context-dependent candidate": "#D99F3D", "Restoration/pathway marker": COL_DOWN, "miRNA-axis hypothesis": "#4FA188", "External TNIK bridge": "#7E61A8"}
    plot = table.copy()
    plot["class_order"] = plot["actionability_class"].map(order)
    plot = plot.sort_values(["class_order", "oligo_actionability_score"], ascending=[False, True])
    svg.text(800, 60, "Oligonucleotide-focused perturbation-triage map", 32, "700", "middle")
    x0, y0, x1, y1 = 275, 140, 1340, 1460
    for tick in range(0, 6):
        x = xmap(tick, 0, 5.6, x0, x1)
        svg.line(x, y0, x, y1, COL_GRID, 1.3)
        svg.text(x, y1 + 35, str(tick), 16, "400", "middle", "#555555")
    row_h = (y1 - y0) / len(plot)
    for i, (_, row) in enumerate(plot.iterrows()):
        y = y0 + i * row_h + row_h * 0.17
        score = float(row["oligo_actionability_score"])
        width = xmap(score, 0, 5.6, x0, x1) - x0
        cls = row["actionability_class"]
        svg.text(35, y + row_h * 0.40, row["gene_symbol"], 21)
        svg.rect(x0, y, width, row_h * 0.62, colors.get(cls, "#7A8793"))
        svg.text(x0 + width + 16, y + row_h * 0.40, str(cls).replace(" candidate", ""), 16)
    svg.line(x0, y1, x1, y1, "#333333", 1.5)
    svg.text((x0 + x1) / 2, 1535, "Oligonucleotide validation-planning score", 19, "400", "middle")
    save_svg(svg, "Figure_6_perturbation_triage_map_editable", "editable horizontal bars, candidate labels, class labels, axis title, gridlines, and main title")


def write_index() -> None:
    readme = """# Editable Main-Figure Sources

This folder contains editable SVG source versions for all six main figures in the BMC Genomics Research article manuscript.

- `svg/` contains the editable vector files.
- `editable_figure_manifest.csv` maps each figure to its source and lists which components are editable.
- `editable_figures_index.html` provides a quick browser preview.

The manuscript DOCX still uses rendered images for stable journal layout. Edit these SVG files in PowerPoint, Illustrator, Inkscape, Affinity Designer, or another vector editor, then export the corrected PNG/PDF for the final manuscript. Heatmaps were redrawn as individual SVG rectangles rather than raster images, so cell blocks and numeric labels can be adjusted.
"""
    (EDIT_DIR / "README.md").write_text(readme, encoding="utf-8")
    with (EDIT_DIR / "editable_figure_manifest.csv").open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["figure", "editable_svg", "editable_elements"])
        writer.writeheader()
        writer.writerows(MANIFEST_ROWS)
    cards = []
    for row in MANIFEST_ROWS:
        cards.append(f"<h2>{clean(row['figure'])}</h2><p>{clean(row['editable_elements'])}</p><object data='{clean(row['editable_svg'])}' type='image/svg+xml' style='width:100%;border:1px solid #ddd'></object>")
    html = "<!doctype html><meta charset='utf-8'><title>Editable figure preview</title><body style='font-family:Arial,sans-serif;max-width:1200px;margin:32px auto'>" + "\n".join(cards) + "</body>"
    (EDIT_DIR / "editable_figures_index.html").write_text(html, encoding="utf-8")


def main() -> None:
    make_figure1()
    make_figure2()
    make_figure3()
    make_figure4()
    make_figure5()
    make_figure6()
    write_index()
    print(EDIT_DIR)
    print(f"Editable SVG figures written: {len(MANIFEST_ROWS)}")


if __name__ == "__main__":
    main()
