"""Generate Supplementary Fig. 2: heatmap of standardized physiological responses."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ncc_style import mm_to_in, save_figure, set_ncc_style
from utils import write_source_data

TRAITS = ["Growth_rates", "Nfix_PON", "Cfix_POC"]
TRAIT_DISPLAY = ["Growth\nrate", "N$_2$ fix.\n/PON", "C fix.\n/POC"]


def build_figure(out_dir: str | Path, dpi: int = 600, vector_only: bool = False) -> list[Path]:
    # Step 1: Load the supplied treatment means and verify standardized values.
    set_ncc_style()
    plt.rcParams.update({"font.family": "Arial", "font.sans-serif": ["Arial"], "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold"})
    source_dir = Path(out_dir) / "source_data"
    metadata = pd.read_csv(source_dir / "Treatment_mean_matrix.csv")
    matrix = metadata.set_index("Treatment_Group")[TRAITS]
    assert matrix.shape == (25, 3)
    assert [g for g in matrix.index if g.startswith("T0_")] == ["T0_28HFe"]
    assert len([g for g in matrix.index if g.startswith("IA_")]) == 6
    assert np.isfinite(matrix.values).all()
    z = (matrix - matrix.mean(axis=0)) / matrix.std(axis=0, ddof=1)
    expected_z = pd.read_csv(source_dir / "Standardized_matrix.csv").set_index("Treatment_Group")[TRAITS]
    assert np.max(np.abs(z.values-expected_z.loc[z.index].values)) < 1e-12

    # Step 2: Plot the standardized physiological profiles.
    fig, ax = plt.subplots(1, 1, figsize=(mm_to_in(92), mm_to_in(150)))
    vmax = np.nanmax(np.abs(z.values))
    im = ax.imshow(z.values, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)

    ax.set_xticks(np.arange(len(TRAITS)))
    ax.set_xticklabels(TRAIT_DISPLAY)
    ax.set_yticks(np.arange(len(z.index)))
    ax.set_yticklabels(metadata["Heatmap_Label"].tolist())
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
    # Step 3: Save the exact plotted matrix and figure files.
    source = z.reset_index().rename(columns={"index": "Treatment_Group"})
    write_source_data(source, source_dir, "SuppFig2_zscore_source_data")
    return save_figure(fig, out_dir, "SuppFig2_standardized_heatmap", dpi=dpi, formats=("pdf", "svg") if vector_only else ("pdf", "svg", "tif", "png"))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default=str(Path(__file__).resolve().parents[2] / "results/reproduced/multivariate"))
    parser.add_argument("--vector-only", action="store_true")
    parser.add_argument("--dpi", default=600, type=int)
    args = parser.parse_args()
    paths = build_figure(args.out, dpi=args.dpi, vector_only=args.vector_only)
    for p in paths:
        print(p.name)
