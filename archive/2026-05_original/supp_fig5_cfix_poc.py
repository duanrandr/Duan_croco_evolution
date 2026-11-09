"""Generate Supplementary Fig. 5: C-specific carbon fixation ."""
from __future__ import annotations

from supp_fullbars import build_supp_bar

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", default="data/Duan_Supp_Table1_052226.xlsx")
    parser.add_argument("--out", default="outputs/figures")
    parser.add_argument("--dpi", default=600, type=int)
    args = parser.parse_args()
    paths = build_supp_bar(args.data, args.out, trait="Cfix_POC", stem="SuppFig5_cfix_responses", dpi=args.dpi)
    for p in paths:
        print(p)
