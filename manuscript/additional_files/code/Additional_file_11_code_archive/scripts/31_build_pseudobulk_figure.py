#!/usr/bin/env python
"""
Build a donor-aware pseudobulk figure and merge it into Figure 6.

Panel C labels donor counts and FDR values so the pseudobulk validation is not
hidden in text or supplementary tables.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from PIL import Image, ImageDraw, ImageFont


PROJECT_DIR = Path(__file__).resolve().parents[1]
PSEUDO_DIR = PROJECT_DIR / "results" / "single_cell_pseudobulk"
FIG_DIR = PROJECT_DIR / "manuscript" / "figures"


def load_font(size: int, bold: bool = False):
    candidates = ["arialbd.ttf", "arial.ttf"] if bold else ["arial.ttf", "DejaVuSans.ttf"]
    for name in candidates:
        try:
            return ImageFont.truetype(name, size=size)
        except OSError:
            continue
    return ImageFont.load_default()


def fdr_label(value: float) -> str:
    if value < 1e-3:
        return f"FDR={value:.1e}"
    return f"FDR={value:.3f}"


def make_pseudobulk_panel() -> Path:
    summary = pd.read_csv(PSEUDO_DIR / "single_cell_pseudobulk_core_candidate_summary.csv")
    summary = summary.sort_values("ipf_minus_control_log1p_norm", ascending=True).reset_index(drop=True)

    width, height = 2300, 1650
    left, right, top, bottom = 720, 440, 240, 190
    plot_w = width - left - right
    row_h = (height - top - bottom) / len(summary)
    x_min = min(-1.25, float(summary["ipf_minus_control_log1p_norm"].min()) - 0.15)
    x_max = max(1.90, float(summary["ipf_minus_control_log1p_norm"].max()) + 0.15)

    def x_pos(value: float) -> int:
        return int(left + (value - x_min) / (x_max - x_min) * plot_w)

    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    title_font = load_font(62, bold=True)
    font = load_font(44)
    small = load_font(38)
    tiny = load_font(36)

    draw.text((left, 35), "C  Donor-aware pseudobulk validation", fill=(0, 0, 0), font=title_font)
    axis_y0, axis_y1 = top - 15, height - bottom + 10
    for tick in [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5]:
        x = x_pos(tick)
        draw.line((x, axis_y0, x, axis_y1), fill=(230, 235, 242), width=2)
        draw.text((x - 25, height - bottom + 25), f"{tick:g}", fill=(35, 35, 35), font=small)
    zero = x_pos(0)
    draw.line((zero, axis_y0, zero, axis_y1), fill=(95, 95, 95), width=4)

    for i, row in summary.iterrows():
        y = int(top + i * row_h + row_h * 0.5)
        value = float(row["ipf_minus_control_log1p_norm"])
        color = (181, 93, 96) if value >= 0 else (76, 120, 168)
        label = f"{row['gene_symbol']} ({row['broad_celltype']}, {row['series_id']})"
        donor_label = f"n={int(row['ipf_donors'])}/{int(row['control_donors'])}; {fdr_label(float(row['fdr_bh']))}"
        draw.text((30, y - 18), label, fill=(0, 0, 0), font=font)
        x0, x1 = (zero, x_pos(value)) if value >= 0 else (x_pos(value), zero)
        draw.rectangle((x0, y - 22, x1, y + 22), fill=color)
        tx = x1 + 22
        draw.text((tx, y - 18), donor_label, fill=(40, 40, 40), font=tiny)

    draw.text(
        (left, height - 55),
        "IPF-control donor-level pseudobulk log1p normalized expression difference",
        fill=(0, 0, 0),
        font=small,
    )
    legend_x = width - 650
    draw.text((legend_x, 105), "n=IPF/control donors represented", fill=(80, 80, 80), font=tiny)
    draw.rectangle((legend_x, 160, legend_x + 34, 190), fill=(181, 93, 96))
    draw.text((legend_x + 48, 153), "IPF-increased", fill=(80, 80, 80), font=tiny)
    draw.rectangle((legend_x + 330, 160, legend_x + 364, 190), fill=(76, 120, 168))
    draw.text((legend_x + 378, 153), "IPF-decreased", fill=(80, 80, 80), font=tiny)
    out = PSEUDO_DIR / "single_cell_pseudobulk_core_candidates_annotated.png"
    img.save(out)
    return out


def merge_figure6(panel_c: Path) -> None:
    original = Image.open(FIG_DIR / "Figure_6_single_cell_validation_celllevel.png").convert("RGB")
    panel = Image.open(panel_c).convert("RGB")
    target_w = original.width
    scale = target_w / panel.width
    panel = panel.resize((target_w, int(panel.height * scale)), Image.Resampling.LANCZOS)
    gap = 60
    combined = Image.new("RGB", (target_w, original.height + panel.height + gap), "white")
    combined.paste(original, (0, 0))
    combined.paste(panel, (0, original.height + gap))
    combined.save(FIG_DIR / "Figure_6_single_cell_validation.png")


def main() -> None:
    panel_c = make_pseudobulk_panel()
    merge_figure6(panel_c)
    print(panel_c)
    print(FIG_DIR / "Figure_6_single_cell_validation.png")


if __name__ == "__main__":
    main()
