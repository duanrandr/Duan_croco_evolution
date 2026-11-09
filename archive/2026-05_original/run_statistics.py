"""Recompute ANOVA and Tukey-Kramer statistics from raw replicate data.

This is a lightweight wrapper around ``generate_supp_table1.py``. It writes the
main statistical outputs as CSV files for quick inspection, while
``generate_supp_table1.py`` writes the complete multi-sheet Supplementary Table 1
workbook.
"""
from __future__ import annotations

from pathlib import Path

from generate_supp_table1 import load_raw_data, run_anova, run_tukey


def main(data: str | Path, out: str | Path) -> list[Path]:
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    raw = load_raw_data(data)
    anova = run_anova(raw)
    tukey = run_tukey(raw, anova)
    paths = [out / "anova_recomputed.csv", out / "tukey_recomputed.csv"]
    anova.to_csv(paths[0], index=False)
    tukey.to_csv(paths[1], index=False)
    return paths


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", default="data/raw_croco_physiology.csv", help="Raw replicate-level CSV, or an Excel workbook with a Raw_All_Data sheet.")
    parser.add_argument("--out", default="outputs/stats")
    args = parser.parse_args()
    for p in main(args.data, args.out):
        print(p)
