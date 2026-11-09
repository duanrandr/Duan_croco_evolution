"""Shared builder for supplementary full-treatment bar plots."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from ncc_style import COLORS, TRAIT_LABELS, clean_axis, mm_to_in, save_figure, set_ncc_style
from utils import compact_letter_display, get_trait_summary, group_order_all, prefix_label, prefix_order_all, write_source_data


def build_supp_bar(
    excel_path: str | Path,
    out_dir: str | Path,
    trait: str,
    stem: str,
    panel_title: str | None = None,
    dpi: int = 600,
) -> list[Path]:
    set_ncc_style()
    summary = pd.read_excel(excel_path, sheet_name="Summary_All")
    posthoc = pd.read_excel(excel_path, sheet_name="PostHoc")
    source_dir = Path(out_dir).parent / "source_data"

    prefixes = prefix_order_all()
    groups = group_order_all()
    trait_df = get_trait_summary(summary, trait, groups)
    letters = compact_letter_display(posthoc, "All_Data", trait, groups)

    x = np.arange(len(prefixes), dtype=float)
    width = 0.34
    offsets = {"HFe": -width / 2, "LFe": width / 2}

    fig, ax = plt.subplots(1, 1, figsize=(mm_to_in(180), mm_to_in(82)))
    source_rows = []

    for i, prefix in enumerate(prefixes):
        for fe in ["HFe", "LFe"]:
            group = f"{prefix}{fe}"
            row = trait_df.loc[trait_df["Treatment_Group"].astype(str).eq(group)]
            if row.empty:
                continue
            mean = float(row["Mean"].iloc[0])
            sd = float(row["SD"].iloc[0])
            n = int(row["n"].iloc[0])
            xpos = x[i] + offsets[fe]

            hatch = None
            if str(prefix).startswith("T0"):
                hatch = "///"
            elif "-" in str(prefix):
                hatch = "..."

            ax.bar(
                xpos,
                mean,
                width=width,
                color=COLORS[fe],
                edgecolor="black",
                linewidth=0.35,
                hatch=hatch,
                zorder=2,
            )
            ax.errorbar(
                xpos,
                mean,
                yerr=sd,
                fmt="none",
                ecolor="black",
                elinewidth=0.45,
                capsize=1.4,
                capthick=0.45,
                zorder=3,
            )
            letter = letters.get(group, "")
            if letter:
                ax.text(
                    xpos,
                    mean + sd + 0.030 * max(trait_df["Mean"].max(), 1e-9),
                    letter,
                    ha="center",
                    va="bottom",
                    fontsize=4.3 if len(letter) > 3 else 4.8,
                    rotation=0,
                )
            source_rows.append(
                {
                    "figure": stem,
                    "trait": trait,
                    "treatment_group": group,
                    "mean": mean,
                    "sd": sd,
                    "n": n,
                    "letters": letter,
                }
            )

    for divider in [3.5, 7.5]:
        ax.axvline(divider, color=COLORS["light_grey"], lw=0.55, ls="--", zorder=1)
    ax.set_ylabel(TRAIT_LABELS[trait])
    ax.set_xticks(x)
    ax.set_xticklabels([prefix_label(p) for p in prefixes], rotation=45, ha="right")
    ax.set_ylim(bottom=0)
    if panel_title:
        ax.set_title(panel_title, pad=3)
    clean_axis(ax)

    handles = [
        Patch(facecolor=COLORS["HFe"], edgecolor="black", linewidth=0.35, label="HFe"),
        Patch(facecolor=COLORS["LFe"], edgecolor="black", linewidth=0.35, label="LFe"),
        Patch(facecolor="white", edgecolor="black", linewidth=0.35, hatch="///", label="T0 ancestral"),
        Patch(facecolor="white", edgecolor="black", linewidth=0.35, hatch="...", label="short-term transfer"),
    ]
    fig.legend(handles=handles, frameon=False, loc="upper center", ncol=4, bbox_to_anchor=(0.55, 1.02), columnspacing=0.9, handlelength=1.2)
    fig.subplots_adjust(left=0.07, right=0.995, bottom=0.26, top=0.84)
    write_source_data(pd.DataFrame(source_rows), source_dir, f"{stem}_source_data")
    return save_figure(fig, out_dir, stem, dpi=dpi)
