"""Generate Fig. 3: short-term reciprocal transfer percent-change responses."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ncc_style import COLORS, TRAIT_SHORT_LABELS, TRAIT_ORDER, clean_axis, mm_to_in, panel_label, role_legend, save_figure, set_ncc_style
from utils import p_to_stars, write_source_data

TRANSFER_ORDER = [
    ("Warming", 23, 28),
    ("Warming", 23, 32),
    ("Warming", 28, 32),
    ("Cooling", 32, 28),
    ("Cooling", 32, 23),
    ("Cooling", 28, 23),
]


def build_figure(excel_path: str | Path, out_dir: str | Path, dpi: int = 600) -> list[Path]:
    set_ncc_style()
    pc = pd.read_excel(excel_path, sheet_name="Percent_Change")
    source_dir = Path(out_dir).parent / "source_data"

    categories: list[dict] = []
    for direction, src, dst in TRANSFER_ORDER:
        for fe in ["HFe", "LFe"]:
            categories.append({"direction": direction, "source": src, "destination": dst, "fe": fe})

    x = np.arange(len(categories), dtype=float)
    width = 0.34
    offsets = {"Source": -width / 2, "Destination": width / 2}
    role_order = ["Source", "Destination"]

    fig, axes = plt.subplots(3, 1, figsize=(mm_to_in(180), mm_to_in(118)), sharex=True)
    panel_labels = ["a", "b", "c"]
    source_rows = []

    for ax, trait, label in zip(axes, TRAIT_ORDER, panel_labels):
        sub = pc.loc[pc["Trait"].astype(str).eq(trait)].copy()
        y_values = []
        for i, cat in enumerate(categories):
            for role in role_order:
                row = sub.loc[
                    sub["Direction"].astype(str).eq(cat["direction"])
                    & sub["Source_Temp"].eq(cat["source"])
                    & sub["Destination_Temp"].eq(cat["destination"])
                    & sub["Fe"].astype(str).eq(cat["fe"])
                    & sub["Reference_Role"].astype(str).eq(role)
                ]
                if row.empty:
                    continue
                pct = float(row["Percent_Change"].iloc[0])
                p = float(row["Tukey_Adjusted_P"].iloc[0])
                sig = str(row["Tukey_Significant"].iloc[0]) == "Yes"
                xpos = x[i] + offsets[role]
                ax.bar(
                    xpos,
                    pct,
                    width=width,
                    color=COLORS[role],
                    edgecolor="black",
                    linewidth=0.35,
                    zorder=2,
                )
                star = p_to_stars(p) if sig else ""
                if star:
                    # Place stars just outside the bar, above positive bars and below negative bars.
                    sign = 1 if pct >= 0 else -1
                    ax.text(
                        xpos,
                        pct + sign * 4,
                        star,
                        ha="center",
                        va="bottom" if pct >= 0 else "top",
                        fontsize=5.2,
                    )
                y_values.append(pct)
                source_rows.append(
                    {
                        "figure": "Fig3",
                        "panel": label,
                        "trait": trait,
                        "direction": cat["direction"],
                        "transfer": f"{cat['source']}->{cat['destination']}",
                        "fe": cat["fe"],
                        "reference_role": role,
                        "percent_change": pct,
                        "tukey_adjusted_p": p,
                        "significant": sig,
                    }
                )

        ax.axhline(0, color="black", lw=0.55, zorder=1)
        # Divider between warming and cooling categories.
        ax.axvline(5.5, color=COLORS["light_grey"], lw=0.55, ls="--", zorder=1)
        ax.set_ylabel(f"{TRAIT_SHORT_LABELS[trait]}\n% change")
        clean_axis(ax)
        panel_label(ax, label, x=-0.07, y=1.02)
        # Leave room for stars on both sides.
        if y_values:
            ymax = max(y_values)
            ymin = min(y_values)
            span = max(ymax - ymin, 1)
            ax.set_ylim(ymin - 0.12 * span - 8, ymax + 0.12 * span + 8)

    labels = [f"{c['source']}→{c['destination']}\n{c['fe']}" for c in categories]
    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels(labels, rotation=45, ha="right")
    axes[0].text(0.24, 1.04, "Warming transfers", transform=axes[0].transAxes, ha="center", va="bottom", fontsize=6)
    axes[0].text(0.76, 1.04, "Cooling transfers", transform=axes[0].transAxes, ha="center", va="bottom", fontsize=6)
    role_legend(axes[0], loc="upper right")

    fig.subplots_adjust(left=0.085, right=0.995, bottom=0.19, top=0.93, hspace=0.35)
    write_source_data(pd.DataFrame(source_rows), source_dir, "Fig3_source_data")
    return save_figure(fig, out_dir, "Fig3_percent_change", dpi=dpi)


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
