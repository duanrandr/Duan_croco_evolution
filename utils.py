"""Shared data-loading and statistics utilities for Crocosphaera figures."""
from __future__ import annotations

from itertools import combinations
from pathlib import Path
from string import ascii_lowercase
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

TRAITS = ["Growth_rates", "Nfix_PON", "Cfix_POC"]
FE_ORDER = ["HFe", "LFe"]


def read_tables(excel_path: str | Path) -> dict[str, pd.DataFrame]:
    """Load all manuscript data tables from the Excel workbook."""
    excel_path = Path(excel_path)
    required = [
        "Raw_All_Data",
        "Summary_All",
        "Percent_Change",
        "ANOVA",
        "PostHoc",
        "PostHoc_Significant",
        "Key_Comparisons",
    ]
    tables = {sheet: pd.read_excel(excel_path, sheet_name=sheet) for sheet in required}
    return tables


def to_group_label(group: str) -> str:
    """Convert internal treatment group names to compact figure labels."""
    label = str(group)
    label = label.replace("T0_", "T0 ")
    label = label.replace("-", "→")
    label = label.replace("HFe", " HFe")
    label = label.replace("LFe", " LFe")
    return label


def prefix_label(prefix: str) -> str:
    label = str(prefix).replace("T0_", "T0 ").replace("-", "→")
    return label


def summary_group_order_longterm() -> list[str]:
    """Treatment order used in Fig. 2."""
    prefixes = ["T0_23", "23", "T0_28", "28", "T0_32", "32"]
    return [f"{prefix}{fe}" for prefix in prefixes for fe in FE_ORDER]


def prefix_order_longterm() -> list[str]:
    return ["T0_23", "23", "T0_28", "28", "T0_32", "32"]


def group_order_all() -> list[str]:
    """Treatment order used for supplementary full-treatment bar plots."""
    prefixes = [
        "T0_23",
        "23",
        "23-28",
        "23-32",
        "T0_28",
        "28",
        "28-23",
        "28-32",
        "T0_32",
        "32",
        "32-23",
        "32-28",
    ]
    return [f"{prefix}{fe}" for prefix in prefixes for fe in FE_ORDER]


def prefix_order_all() -> list[str]:
    return [
        "T0_23",
        "23",
        "23-28",
        "23-32",
        "T0_28",
        "28",
        "28-23",
        "28-32",
        "T0_32",
        "32",
        "32-23",
        "32-28",
    ]


def get_trait_summary(summary: pd.DataFrame, trait: str, groups: Sequence[str]) -> pd.DataFrame:
    """Return summary rows for a trait in a specified treatment order."""
    sub = summary.loc[summary["Trait"].eq(trait)].copy()
    sub["Treatment_Group"] = sub["Treatment_Group"].astype(str)
    order = pd.Categorical(sub["Treatment_Group"], categories=list(groups), ordered=True)
    sub = sub.assign(_order=order).dropna(subset=["_order"]).sort_values("_order")
    return sub.drop(columns="_order")


def format_p(p: float) -> str:
    if pd.isna(p):
        return "NA"
    if p < 0.001:
        return f"{p:.1e}"
    return f"{p:.3f}"


def p_to_stars(p: float) -> str:
    if pd.isna(p):
        return ""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return ""


def _letters(n: int) -> list[str]:
    """Generate a, b, ..., z, aa, ab, ..."""
    labels: list[str] = []
    base = len(ascii_lowercase)
    i = 0
    while len(labels) < n:
        x = i
        s = ""
        while True:
            s = ascii_lowercase[x % base] + s
            x = x // base - 1
            if x < 0:
                break
        labels.append(s)
        i += 1
    return labels


def compact_letter_display(
    posthoc: pd.DataFrame,
    analysis_set: str,
    trait: str,
    group_order: Sequence[str],
    alpha: float = 0.05,
) -> dict[str, str]:
    """Create a compact, conservative letter display from Tukey comparisons.

    Groups sharing the same letter are not significantly different at `alpha`.
    The algorithm prioritizes legibility in journal figures: each group receives
    one letter, and no significant pair is allowed to share a letter. Exact
    adjusted P values are retained in the Source Data CSV files.
    """
    requested_groups = [str(g) for g in group_order]
    if not requested_groups:
        return {}

    ph = posthoc.loc[
        posthoc["Analysis_Set"].astype(str).eq(analysis_set)
        & posthoc["Trait"].astype(str).eq(trait)
    ].copy()
    if ph.empty:
        return {g: "" for g in requested_groups}

    # Use means to sort groups from high to low response so the letters are
    # stable and intuitive. Fall back to the requested plotting order.
    means: dict[str, float] = {}
    for _, r in ph.iterrows():
        means[str(r["Group_1"])] = float(r["Mean_1"])
        means[str(r["Group_2"])] = float(r["Mean_2"])

    groups = [g for g in requested_groups if g in means]
    groups = sorted(groups, key=lambda g: (-means.get(g, 0), requested_groups.index(g)))

    sig: dict[frozenset[str], bool] = {}
    for _, r in ph.iterrows():
        g1, g2 = str(r["Group_1"]), str(r["Group_2"])
        reject = bool(r["Reject_H0"]) if not pd.isna(r["Reject_H0"]) else bool(r["Adjusted_P"] < alpha)
        sig[frozenset([g1, g2])] = reject

    def significant(g1: str, g2: str) -> bool:
        if g1 == g2:
            return False
        return sig.get(frozenset([g1, g2]), False)

    letter_sets: list[set[str]] = []
    assignments: dict[str, str] = {}
    letters = _letters(len(groups) + 5)

    for g in groups:
        assigned = False
        for idx, s in enumerate(letter_sets):
            if all(not significant(g, h) for h in s):
                s.add(g)
                assignments[g] = letters[idx]
                assigned = True
                break
        if not assigned:
            letter_sets.append({g})
            assignments[g] = letters[len(letter_sets) - 1]

    return {g: assignments.get(g, "") for g in requested_groups}

def write_source_data(df: pd.DataFrame, out_dir: str | Path, stem: str) -> Path:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{stem}.csv"
    df.to_csv(path, index=False)
    return path
