"""Run all Crocosphaera experimental-data figure scripts."""
from __future__ import annotations

from pathlib import Path

from fig02_longterm import build_figure as build_fig2
from fig03_percent_change import build_figure as build_fig3
from supp_fig1_cluster import build_figure as build_supp1
from supp_fig2_heatmap import build_figure as build_supp2
from supp_fullbars import build_supp_bar


def run_all(data: str | Path, out: str | Path, dpi: int = 600) -> list[Path]:
    out = Path(out)
    outputs: list[Path] = []
    outputs += build_fig2(data, out, dpi=dpi)
    outputs += build_fig3(data, out, dpi=dpi)
    outputs += build_supp1(data, out, dpi=dpi)
    outputs += build_supp2(data, out, dpi=dpi)
    outputs += build_supp_bar(data, out, trait="Growth_rates", stem="SuppFig3_growth_responses", dpi=dpi)
    outputs += build_supp_bar(data, out, trait="Nfix_PON", stem="SuppFig4_nfix_responses", dpi=dpi)
    outputs += build_supp_bar(data, out, trait="Cfix_POC", stem="SuppFig5_cfix_responses", dpi=dpi)
    return outputs


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", default="data/Duan_Supp_Table1_052226.xlsx")
    parser.add_argument("--out", default="outputs/figures")
    parser.add_argument("--dpi", default=600, type=int)
    args = parser.parse_args()
    paths = run_all(args.data, args.out, dpi=args.dpi)
    for p in paths:
        print(p)
