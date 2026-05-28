#!/usr/bin/env python
"""Create a miRNA evidence-grading workflow figure for the manuscript."""

from __future__ import annotations

import os
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parents[1]
OUT_DIR = PROJECT_DIR / "results" / "mirna_evidence_workflow"
OUT_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(PROJECT_DIR / "results" / "_mpl_config"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd


def box(ax, x, y, w, h, text, fc="#F6F8FA", ec="#4A5568"):
    rect = patches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.012", fc=fc, ec=ec, lw=1)
    ax.add_patch(rect)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=8)
    return {"left": x, "right": x + w, "bottom": y, "top": y + h, "cx": x + w / 2, "cy": y + h / 2}


def arrow(ax, a, b):
    ax.annotate("", xy=b, xytext=a, arrowprops=dict(arrowstyle="->", lw=1.0, color="#4A5568", shrinkA=5, shrinkB=5))


def main() -> None:
    robust_qc = pd.read_csv(PROJECT_DIR / "results/robust_candidates/robust_candidate_qc.csv")
    axis_qc = pd.read_csv(PROJECT_DIR / "results/mirna_mrna_axes/mirna_mrna_axis_qc.csv").iloc[0]
    graded = pd.read_csv(PROJECT_DIR / "results/submission_enhancements/mirna_mrna_axes_evidence_graded.csv")

    robust_mirna = int(robust_qc.loc[robust_qc["data_type"].eq("miRNA"), "robust_strict_candidates"].iloc[0])
    robust_mrna = int(robust_qc.loc[robust_qc["data_type"].eq("bulk mRNA"), "robust_strict_candidates"].iloc[0])
    exact = int((graded["match_type"] == "exact").sum())
    arm = int((graded["match_type"] == "arm_agnostic").sum())
    main_text = int((graded["recommended_manuscript_role"] == "main_text_prioritized_axis").sum())

    fig, ax = plt.subplots(figsize=(9.2, 4.4))
    ax.axis("off")
    ax.text(0.5, 0.92, "miRNA-mRNA evidence grading workflow", ha="center", va="center", fontsize=13, fontweight="bold")

    b1 = box(ax, 0.04, 0.60, 0.18, 0.18, f"Robust candidates\n{robust_mirna} miRNAs\n{robust_mrna} mRNAs")
    b2 = box(ax, 0.29, 0.60, 0.18, 0.18, "miRTarBase 2025\ncandidate target\nmatching")
    b3 = box(ax, 0.54, 0.60, 0.18, 0.18, f"Opposite-direction\nfilter\n{int(axis_qc['negative_direction_axes'])} inverse axes")
    b4 = box(ax, 0.78, 0.60, 0.18, 0.18, f"Evidence grading\n{exact} exact\n{arm} arm-agnostic")
    b5 = box(ax, 0.29, 0.23, 0.24, 0.18, f"Main text\n{main_text} exact mature-miRNA axes\nhsa-miR-375 targets", fc="#EAF4F0")
    b6 = box(ax, 0.62, 0.23, 0.24, 0.18, "Supplementary\narm-agnostic axes\nrequire arm-specific validation", fc="#F8EFE6")

    for left, right in [(b1, b2), (b2, b3), (b3, b4)]:
        arrow(ax, (left["right"], left["cy"]), (right["left"], right["cy"]))
    arrow(ax, (b4["cx"], b4["bottom"]), (b5["cx"], b5["top"]))
    arrow(ax, (b4["cx"], b4["bottom"]), (b6["cx"], b6["top"]))

    ax.text(
        0.5,
        0.08,
        "Strict grading intentionally narrows the main-text miRNA result: exact mature-miRNA axes are prioritized, while arm-agnostic axes remain hypothesis-generating.",
        ha="center",
        va="center",
        fontsize=8,
    )
    fig.tight_layout()
    fig.savefig(OUT_DIR / "mirna_evidence_grading_workflow.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_DIR / "mirna_evidence_grading_workflow.pdf", bbox_inches="tight")
    print(OUT_DIR)


if __name__ == "__main__":
    main()
