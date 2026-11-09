"""Generate Supplementary Fig. 1: hierarchical clustering of treatment means."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist

from ncc_style import COLORS, mm_to_in, save_figure, set_ncc_style
from utils import group_order_all, to_group_label, write_source_data

TRAITS = ["Growth_rates", "Nfix_PON", "Cfix_POC"]


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
    z = (matrix - matrix.mean(axis=0)) / matrix.std(axis=0, ddof=1)

    dist = pdist(z.values, metric="euclidean")
    link = linkage(dist, method="average")

    fig, ax = plt.subplots(1, 1, figsize=(mm_to_in(145), mm_to_in(120)))
    dendrogram(
        link,
        labels=[to_group_label(g) for g in z.index],
        orientation="right",
        leaf_font_size=5.6,
        color_threshold=0,
        above_threshold_color=COLORS["dark_grey"],
        ax=ax,
    )
    ax.set_xlabel("Clustering distance")
    ax.tick_params(axis="y", length=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.55)
    ax.tick_params(axis="x", width=0.55, length=2.5)
    fig.subplots_adjust(left=0.25, right=0.99, bottom=0.08, top=0.98)

    source = z.reset_index().rename(columns={"index": "Treatment_Group"})
    write_source_data(source, source_dir, "SuppFig1_zscore_source_data")
    return save_figure(fig, out_dir, "SuppFig1_hierarchical_clustering", dpi=dpi)


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
