"""Generate Supplementary Fig. 2: heatmap of standardized physiological responses."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ncc_style import mm_to_in, save_figure, set_ncc_style
from utils import group_order_all, to_group_label, write_source_data

TRAITS = ["Growth_rates", "Nfix_PON", "Cfix_POC"]
TRAIT_DISPLAY = ["Growth\nrate", "N₂ fix.\n/PON", "C fix.\n/POC"]


def build_figure(excel_path: str | Path, out_dir: str | Path, dpi: int = 600) -> list[Path]:
    set_ncc_style()
    summary = pd.read_excel(excel_path, sheet_name="Summary_All")
    source_dir = Path(out_dir).parent / "source_data"

    matrix = (
        summary.loc[summary["Trait"].isin(TRAITS)]
        .pivot_table(index="Treatment_Group", columns="Trait", values="Mean", aggfunc="mean")
        .reindex(group_order_all())
        .dropna(how="any")
    )
    matrix = matrix[TRAITS]
    z = (matrix - matrix.mean(axis=0)) / matrix.std(axis=0, ddof=1)

    fig, ax = plt.subplots(1, 1, figsize=(mm_to_in(92), mm_to_in(120)))
    vmax = np.nanmax(np.abs(z.values))
    im = ax.imshow(z.values, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)

    ax.set_xticks(np.arange(len(TRAITS)))
    ax.set_xticklabels(TRAIT_DISPLAY)
    ax.set_yticks(np.arange(len(z.index)))
    ax.set_yticklabels([to_group_label(g) for g in z.index])
    ax.tick_params(axis="both", length=0)

    # White grid lines between cells.
    ax.set_xticks(np.arange(-0.5, len(TRAITS), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(z.index), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.5)

    for i in range(z.shape[0]):
        for j in range(z.shape[1]):
            val = z.iat[i, j]
            color = "white" if abs(val) > 1.4 else "black"
            ax.text(j, i, f"{val:.1f}", ha="center", va="center", fontsize=5.2, color=color)

    cbar = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03)
    cbar.set_label("z score")
    cbar.ax.tick_params(width=0.55, length=2.5)
    for spine in ax.spines.values():
        spine.set_visible(False)

    fig.subplots_adjust(left=0.38, right=0.91, bottom=0.06, top=0.99)
    source = z.reset_index().rename(columns={"index": "Treatment_Group"})
    write_source_data(source, source_dir, "SuppFig2_zscore_source_data")
    return save_figure(fig, out_dir, "SuppFig2_standardized_heatmap", dpi=dpi)


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
