"""Generate Fig. 2: ancestral (T0) and long-term selected physiological responses."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.patches import Patch

from ncc_style import COLORS, TRAIT_LABELS, TRAIT_ORDER, clean_axis, mm_to_in, panel_label, save_figure, set_ncc_style
from utils import compact_letter_display, get_trait_summary, prefix_label, prefix_order_longterm, summary_group_order_longterm, write_source_data


def build_figure(excel_path: str | Path, out_dir: str | Path, dpi: int = 600) -> list[Path]:
    set_ncc_style()
    summary = pd.read_excel(excel_path, sheet_name="Summary_All")
    posthoc = pd.read_excel(excel_path, sheet_name="PostHoc")
    source_dir = Path(out_dir).parent / "source_data"

    prefixes = prefix_order_longterm()
    group_order = summary_group_order_longterm()
    x = np.arange(len(prefixes), dtype=float)
    width = 0.34
    offsets = {"HFe": -width / 2, "LFe": width / 2}

    fig, axes = plt.subplots(1, 3, figsize=(mm_to_in(180), mm_to_in(72)), sharex=True)

    source_rows = []
    panel_labels = ["a", "b", "c"]
    for ax, trait, label in zip(axes, TRAIT_ORDER, panel_labels):
        trait_df = get_trait_summary(summary, trait, group_order)
        # In the T0_and_Selected post-hoc sheet, selected groups are named
        # Selected_23HFe rather than 23HFe. Map them back to the plotting names.
        posthoc_group_map = {
            g: (g if g.startswith("T0_") else f"Selected_{g}") for g in group_order
        }
        raw_letters = compact_letter_display(
            posthoc, "T0_and_Selected", trait, [posthoc_group_map[g] for g in group_order]
        )
        letters = {g: raw_letters.get(posthoc_group_map[g], "") for g in group_order}

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
                hatch = "///" if str(prefix).startswith("T0") else None
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
                    capsize=1.5,
                    capthick=0.45,
                    zorder=3,
                )
                letter = letters.get(group, "")
                if letter:
                    ax.text(
                        xpos,
                        mean + sd + 0.035 * max(trait_df["Mean"].max(), 1e-9),
                        letter,
                        ha="center",
                        va="bottom",
                        fontsize=4.8 if len(letter) > 3 else 5.2,
                    )
                source_rows.append(
                    {
                        "figure": "Fig2",
                        "panel": label,
                        "trait": trait,
                        "treatment_group": group,
                        "mean": mean,
                        "sd": sd,
                        "n": n,
                        "letters": letter,
                    }
                )

        ax.set_ylabel(TRAIT_LABELS[trait])
        ax.set_xticks(x)
        ax.set_xticklabels([prefix_label(p) for p in prefixes], rotation=45, ha="right")
        ax.axvline(1.5, color=COLORS["light_grey"], lw=0.55, ls="--", zorder=1)
        ax.axvline(3.5, color=COLORS["light_grey"], lw=0.55, ls="--", zorder=1)
        ax.set_ylim(bottom=0)
        clean_axis(ax)
        panel_label(ax, label)

    handles = [
        Patch(facecolor=COLORS["HFe"], edgecolor="black", linewidth=0.35, label="HFe"),
        Patch(facecolor=COLORS["LFe"], edgecolor="black", linewidth=0.35, label="LFe"),
        Patch(facecolor="white", edgecolor="black", linewidth=0.35, hatch="///", label="T0 ancestral"),
    ]
    fig.legend(handles=handles, frameon=False, loc="upper center", ncol=3, bbox_to_anchor=(0.52, 1.02), handlelength=1.2, columnspacing=1.0)

    fig.subplots_adjust(left=0.065, right=0.99, bottom=0.24, top=0.86, wspace=0.34)
    write_source_data(pd.DataFrame(source_rows), source_dir, "Fig2_source_data")
    return save_figure(fig, out_dir, "Fig2_longterm_responses", dpi=dpi)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", default="data/Duan_Supp_Table1_052226.xlsx")
    parser.add_argument("--out", default="outputs/figures")
    parser.add_argument("--dpi", default=600, type=int)
    args = parser.parse_args()
    paths = build_figure(args.data, args.out, dpi=args.dpi)
    for p in paths:
        print(p)
