"""Generate Supplementary Table 1 from raw Crocosphaera replicate data.

This script rebuilds the multi-sheet Supplementary Table 1 workbook from the
raw replicate-level physiology table only. It generates:

- README
- Raw_All_Data
- Summary_All
- Percent_Change
- ANOVA
- PostHoc
- PostHoc_Significant
- Key_Comparisons

Input can be either:
1. a CSV file containing the raw replicate table, or
2. an Excel workbook containing a sheet named ``Raw_All_Data``.

Example
-------
python scripts/generate_supp_table1.py \
    --raw data/raw_croco_physiology.csv \
    --out outputs/tables/Supplementary_Table_1_generated.xlsx
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import studentized_range
from statsmodels.formula.api import ols

TRAITS = ["Growth_rates", "Nfix_PON", "Cfix_POC"]
FE_LEVELS = ["HFe", "LFe"]
ALPHA = 0.05

RAW_COLUMNS = [
    "Replicate_ID",
    "Treatment_Class",
    "Selection_Status",
    "Thermal_Treatment",
    "Culture_Treatment",
    "Treatment_Group",
    "Fe",
    "Replicate",
    "Source_Temp",
    "Endpoint_Temp",
    "Selection_Temp",
    "Direction",
    "Growth_rates",
    "Nfix_PON",
    "Cfix_POC",
]

TRANSFER_ORDER = [
    (23, 28),
    (23, 32),
    (28, 23),
    (28, 32),
    (32, 23),
    (32, 28),
]

LEVELS = {
    "All_Data": ["T0_23", "T0_28", "T0_32", "23", "28", "32", "23-28", "23-32", "28-23", "28-32", "32-23", "32-28"],
    "T0_and_Selected": ["T0_23", "T0_28", "T0_32", "Selected_23", "Selected_28", "Selected_32"],
    "Selected_Only": ["23", "28", "32"],
    "T0_Only": ["T0_23", "T0_28", "T0_32"],
    "Selected_and_Switching": ["23", "28", "32", "23-28", "23-32", "28-23", "28-32", "32-23", "32-28"],
    "Switching_Only": ["23-28", "23-32", "28-23", "28-32", "32-23", "32-28"],
    "T0_vs_Selected_23C": ["T0", "Selected"],
    "T0_vs_Selected_28C": ["T0", "Selected"],
    "T0_vs_Selected_32C": ["T0", "Selected"],
}


@dataclass(frozen=True)
class AnalysisSet:
    name: str
    factor_a: str
    data: pd.DataFrame
    description: str
    levels_a: list[str]


def _clean_factor(series: pd.Series) -> pd.Series:
    """Convert categorical labels to stable strings while preserving missing values."""
    out = series.copy()
    out = out.astype(str).str.replace(r"\.0$", "", regex=True)
    out = out.replace({"nan": np.nan, "None": np.nan, "NaN": np.nan, "<NA>": np.nan})
    return out


def load_raw_data(path: str | Path) -> pd.DataFrame:
    """Load raw replicate data from CSV or from the Raw_All_Data Excel sheet."""
    path = Path(path)
    if path.suffix.lower() in {".csv", ".txt"}:
        raw = pd.read_csv(path)
    elif path.suffix.lower() in {".xlsx", ".xlsm", ".xls"}:
        raw = pd.read_excel(path, sheet_name="Raw_All_Data")
    else:
        raise ValueError(f"Unsupported raw-data file type: {path.suffix}")

    missing = [col for col in RAW_COLUMNS if col not in raw.columns]
    if missing:
        raise ValueError(f"Raw data are missing required columns: {missing}")

    raw = raw[RAW_COLUMNS].copy()
    for col in ["Treatment_Class", "Selection_Status", "Thermal_Treatment", "Culture_Treatment", "Treatment_Group", "Fe", "Direction"]:
        raw[col] = _clean_factor(raw[col])
    for col in ["Replicate", "Source_Temp", "Endpoint_Temp", "Selection_Temp", *TRAITS]:
        raw[col] = pd.to_numeric(raw[col], errors="coerce")
    return raw


def export_raw_csv_from_workbook(workbook: str | Path, out_csv: str | Path) -> Path:
    """Extract only the Raw_All_Data sheet from a workbook into a CSV file."""
    raw = load_raw_data(workbook)
    out_csv = Path(out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    raw.to_csv(out_csv, index=False)
    return out_csv


def build_readme(raw: pd.DataFrame) -> pd.DataFrame:
    n_records = int(raw.shape[0])
    n_groups = int(raw["Treatment_Group"].nunique())
    return pd.DataFrame(
        [
            ("Raw records", f"{n_records} replicate rows; treatment groups = {n_groups}."),
            ("Traits analyzed", ", ".join(TRAITS)),
            ("ANOVA", "All two-way ANOVA outputs are combined in one sheet. Models use Factor_A × Fe, with the analysis set defining Factor_A."),
            ("PostHoc", "All Tukey-Kramer HSD pairwise comparisons are combined in one sheet. Use Analysis_Set and Trait filters to view specific tests."),
            ("PostHoc_Significant", "Filtered view of PostHoc_Combined where adjusted P < 0.05."),
            ("Percent_Change", "Switching-culture percent changes relative to source and destination selected references. Positive values are salmon; negative values are blue."),
            ("Key_Comparisons", "Focused direct comparisons: switching vs source/destination selected references, and selected vs T0 ancestral lines."),
            ("ANOVA method", "Two-way ANOVA via partial F-tests: main effects from additive model comparisons and interaction from additive vs full model; missing values excluded."),
            ("Post-hoc method", "Tukey-Kramer HSD using pooled within-group MSE for unequal n; all reported P values are adjusted."),
        ],
        columns=["Item", "Description"],
    )


def summary_treatment_order() -> list[str]:
    """Treatment order used in Supplementary Table 1 summary sheets."""
    prefixes = [
        "T0_23", "T0_28", "T0_32",
        "23", "28", "32",
        "23-28", "23-32", "28-23", "28-32", "32-23", "32-28",
    ]
    return [f"{prefix}{fe}" for prefix in prefixes for fe in FE_LEVELS]


def build_summary(raw: pd.DataFrame) -> pd.DataFrame:
    """Summarize raw replicates as mean, SD and SEM for each treatment and trait."""
    metadata_cols = [
        "Treatment_Group",
        "Treatment_Class",
        "Selection_Status",
        "Thermal_Treatment",
        "Culture_Treatment",
        "Fe",
        "Source_Temp",
        "Endpoint_Temp",
        "Selection_Temp",
        "Direction",
    ]
    records = []
    available = set(raw["Treatment_Group"].astype(str))
    group_order = [g for g in summary_treatment_order() if g in available]
    # Preserve any unexpected groups after the predefined order.
    group_order += [g for g in pd.unique(raw["Treatment_Group"].astype(str)) if g not in group_order]
    for group in group_order:
        g = raw.loc[raw["Treatment_Group"].astype(str).eq(str(group))].copy()
        meta = {col: g[col].iloc[0] for col in metadata_cols}
        for trait in TRAITS:
            values = pd.to_numeric(g[trait], errors="coerce").dropna()
            n = int(values.shape[0])
            mean = float(values.mean()) if n else np.nan
            sd = float(values.std(ddof=1)) if n > 1 else np.nan
            sem = float(sd / np.sqrt(n)) if n > 1 else np.nan
            records.append({**meta, "Trait": trait, "n": n, "Mean": mean, "SD": sd, "SEM": sem})
    return pd.DataFrame.from_records(records)


def get_analysis_sets(raw: pd.DataFrame) -> list[AnalysisSet]:
    df = raw.copy()
    for col in ["Thermal_Treatment", "Culture_Treatment", "Selection_Status", "Fe"]:
        df[col] = _clean_factor(df[col])

    sets = [
        AnalysisSet(
            "All_Data",
            "Thermal_Treatment",
            df.copy(),
            "All data: T0 + long-term selected + short-term switching",
            LEVELS["All_Data"],
        ),
        AnalysisSet(
            "T0_and_Selected",
            "Culture_Treatment",
            df.loc[df["Selection_Status"].isin(["T0", "Selected"])].copy(),
            "T0 ancestral + long-term selected cultures",
            LEVELS["T0_and_Selected"],
        ),
        AnalysisSet(
            "Selected_Only",
            "Thermal_Treatment",
            df.loc[df["Selection_Status"].eq("Selected")].copy(),
            "Long-term selected cultures only",
            LEVELS["Selected_Only"],
        ),
        AnalysisSet(
            "T0_Only",
            "Thermal_Treatment",
            df.loc[df["Selection_Status"].eq("T0")].copy(),
            "T0 ancestral cultures only",
            LEVELS["T0_Only"],
        ),
        AnalysisSet(
            "Selected_and_Switching",
            "Thermal_Treatment",
            df.loc[df["Selection_Status"].isin(["Selected", "Switching"])].copy(),
            "Long-term selected + short-term switching cultures",
            LEVELS["Selected_and_Switching"],
        ),
        AnalysisSet(
            "Switching_Only",
            "Thermal_Treatment",
            df.loc[df["Selection_Status"].eq("Switching")].copy(),
            "Short-term switching cultures only",
            LEVELS["Switching_Only"],
        ),
    ]
    for temp in [23, 28, 32]:
        name = f"T0_vs_Selected_{temp}C"
        sets.append(
            AnalysisSet(
                name,
                "Selection_Status",
                df.loc[
                    df["Endpoint_Temp"].eq(temp)
                    & df["Selection_Status"].isin(["T0", "Selected"])
                ].copy(),
                f"T0 ancestral vs long-term selected at endpoint {temp} °C",
                LEVELS[name],
            )
        )
    return sets


def _ordered_levels_present(series: pd.Series, preferred: Iterable[str]) -> list[str]:
    present = set(series.dropna().astype(str))
    ordered = [x for x in preferred if x in present]
    extras = [x for x in pd.unique(series.dropna().astype(str)) if x not in ordered]
    return ordered + extras


def run_anova(raw: pd.DataFrame) -> pd.DataFrame:
    records = []
    for aset in get_analysis_sets(raw):
        for trait in TRAITS:
            d = aset.data.dropna(subset=[trait, aset.factor_a, "Fe"]).copy()
            if d.empty or d[aset.factor_a].nunique() < 2 or d["Fe"].nunique() < 2:
                continue
            levels_a = _ordered_levels_present(d[aset.factor_a], aset.levels_a)
            d[aset.factor_a] = pd.Categorical(d[aset.factor_a], categories=levels_a, ordered=True)
            d["Fe"] = pd.Categorical(d["Fe"], categories=FE_LEVELS, ordered=True)
            formula = f"{trait} ~ C({aset.factor_a}) * C(Fe)"
            model = ols(formula, data=d).fit()
            table = sm.stats.anova_lm(model, typ=2)
            residual = table.loc["Residual"]
            residual_ss = float(residual["sum_sq"])
            residual_df = int(residual["df"])
            residual_ms = float(residual_ss / residual_df)
            for effect_name, row in table.drop(index="Residual").iterrows():
                effect = effect_name.replace(f"C({aset.factor_a})", aset.factor_a).replace("C(Fe)", "Fe")
                ss = float(row["sum_sq"])
                df_effect = int(row["df"])
                records.append(
                    {
                        "Analysis_Set": aset.name,
                        "Dataset_Description": aset.description,
                        "Trait": trait,
                        "Model": f"Response ~ {aset.factor_a} * Fe",
                        "Factor_A": aset.factor_a,
                        "Factor_B": "Fe",
                        "Effect": effect,
                        "SS": ss,
                        "df": df_effect,
                        "MS": ss / df_effect,
                        "F": float(row["F"]),
                        "p_value": float(row["PR(>F)"]),
                        "Residual_SS": residual_ss,
                        "Residual_df": residual_df,
                        "Residual_MS": residual_ms,
                        "n": int(d.shape[0]),
                        "Levels_A": ", ".join(levels_a),
                        "Levels_B": ", ".join(FE_LEVELS),
                        "Type": "Type II partial F-test; interaction tested from additive to full; missing values excluded",
                    }
                )
    return pd.DataFrame.from_records(records)


def _anova_lookup(anova: pd.DataFrame) -> dict[tuple[str, str], dict[str, float]]:
    out: dict[tuple[str, str], dict[str, float]] = {}
    for (aset, trait), sub in anova.groupby(["Analysis_Set", "Trait"], sort=False):
        first = sub.iloc[0]
        out[(str(aset), str(trait))] = {
            "MSE": float(first["Residual_MS"]),
            "Residual_df": int(first["Residual_df"]),
        }
    return out


def run_tukey(raw: pd.DataFrame, anova: pd.DataFrame, alpha: float = ALPHA) -> pd.DataFrame:
    """Run Tukey-Kramer pairwise tests using each ANOVA model's residual MSE.

    P values and confidence intervals are calculated from the studentized range
    distribution. The implementation vectorizes the distribution call within
    each analysis set to keep the script reasonably fast while preserving exact
    adjusted P values.
    """
    anova_info = _anova_lookup(anova)
    records = []
    for aset in get_analysis_sets(raw):
        for trait in TRAITS:
            d = aset.data.dropna(subset=[trait, aset.factor_a, "Fe"]).copy()
            if d.empty or d[aset.factor_a].nunique() < 2 or d["Fe"].nunique() < 2:
                continue
            levels_a = _ordered_levels_present(d[aset.factor_a], aset.levels_a)
            if aset.factor_a == "Selection_Status":
                group_order = [f"{level}_{fe}" for level in levels_a for fe in FE_LEVELS]
                d["Group"] = _clean_factor(d[aset.factor_a]) + "_" + _clean_factor(d["Fe"])
            else:
                group_order = [f"{level}{fe}" for level in levels_a for fe in FE_LEVELS]
                d["Group"] = _clean_factor(d[aset.factor_a]) + _clean_factor(d["Fe"])
            present = set(d["Group"])
            group_order = [g for g in group_order if g in present]
            if len(group_order) < 2:
                continue

            info = anova_info[(aset.name, trait)]
            mse = float(info["MSE"])
            residual_df = int(info["Residual_df"])
            number_of_groups = len(group_order)
            q_crit = float(studentized_range.ppf(1 - alpha, number_of_groups, residual_df))
            grouped = d.groupby("Group", sort=False)[trait]
            means = grouped.mean().to_dict()
            ns = grouped.count().to_dict()

            pair_data = []
            q_values = []
            for i, group_1 in enumerate(group_order):
                for group_2 in group_order[i + 1 :]:
                    n1 = int(ns[group_1])
                    n2 = int(ns[group_2])
                    mean_1 = float(means[group_1])
                    mean_2 = float(means[group_2])
                    mean_diff = mean_2 - mean_1
                    se = float(np.sqrt(mse * 0.5 * ((1 / n1) + (1 / n2))))
                    q_stat = float(abs(mean_diff) / se) if se > 0 else np.nan
                    pair_data.append((group_1, group_2, mean_1, mean_2, mean_diff, se, q_stat, n1, n2))
                    q_values.append(q_stat)

            p_values = studentized_range.sf(np.asarray(q_values, dtype=float), number_of_groups, residual_df)
            for (group_1, group_2, mean_1, mean_2, mean_diff, se, q_stat, n1, n2), p_adj in zip(pair_data, p_values):
                p_adj = float(p_adj)
                ci_lower = float(mean_diff - q_crit * se)
                ci_upper = float(mean_diff + q_crit * se)
                records.append(
                    {
                        "Analysis_Set": aset.name,
                        "Dataset_Description": aset.description,
                        "Trait": trait,
                        "Group_Field": f"{aset.factor_a}+Fe",
                        "Group_1": group_1,
                        "Group_2": group_2,
                        "Mean_1": mean_1,
                        "Mean_2": mean_2,
                        "Mean_Difference_group2_minus_group1": mean_diff,
                        "SE": se,
                        "q_statistic": q_stat,
                        "Adjusted_P": p_adj,
                        "CI_Lower": ci_lower,
                        "CI_Upper": ci_upper,
                        "Reject_H0": bool(p_adj < alpha) if np.isfinite(p_adj) else False,
                        "Alpha": alpha,
                        "n_1": n1,
                        "n_2": n2,
                        "MSE": mse,
                        "Residual_df": residual_df,
                        "Number_of_Groups": number_of_groups,
                        "Method": "Tukey-Kramer HSD using pooled within-group MSE; adjusted p-values",
                    }
                )
    return pd.DataFrame.from_records(records)

def _lookup_mean(summary: pd.DataFrame, group: str, trait: str) -> float:
    row = summary.loc[
        summary["Treatment_Group"].astype(str).eq(str(group))
        & summary["Trait"].astype(str).eq(str(trait))
    ]
    if row.empty:
        raise KeyError(f"No summary mean for group={group}, trait={trait}")
    return float(row["Mean"].iloc[0])


def _lookup_tukey(posthoc: pd.DataFrame, analysis_set: str, trait: str, group_1: str, group_2: str) -> tuple[float, str]:
    sub = posthoc.loc[
        posthoc["Analysis_Set"].astype(str).eq(str(analysis_set))
        & posthoc["Trait"].astype(str).eq(str(trait))
    ]
    match = sub.loc[
        ((sub["Group_1"].astype(str).eq(str(group_1))) & (sub["Group_2"].astype(str).eq(str(group_2))))
        | ((sub["Group_1"].astype(str).eq(str(group_2))) & (sub["Group_2"].astype(str).eq(str(group_1))))
    ]
    if match.empty:
        raise KeyError(f"No Tukey result for {analysis_set}, {trait}, {group_1} vs {group_2}")
    row = match.iloc[0]
    p = float(row["Adjusted_P"])
    return p, "Yes" if bool(row["Reject_H0"]) else "No"


def build_percent_change(summary: pd.DataFrame, posthoc: pd.DataFrame) -> pd.DataFrame:
    records = []
    for source_temp, destination_temp in TRANSFER_ORDER:
        direction = "Warming" if destination_temp > source_temp else "Cooling"
        for fe in FE_LEVELS:
            transfer_group = f"{source_temp}-{destination_temp}{fe}"
            for role, ref_temp in [("Source", source_temp), ("Destination", destination_temp)]:
                reference_group = f"{ref_temp}{fe}"
                for trait in TRAITS:
                    transfer_mean = _lookup_mean(summary, transfer_group, trait)
                    reference_mean = _lookup_mean(summary, reference_group, trait)
                    pct = ((transfer_mean - reference_mean) / reference_mean) * 100
                    p_adj, sig = _lookup_tukey(posthoc, "Selected_and_Switching", trait, reference_group, transfer_group)
                    records.append(
                        {
                            "Transfer_Group": transfer_group,
                            "Direction": direction,
                            "Fe": fe,
                            "Source_Temp": source_temp,
                            "Destination_Temp": destination_temp,
                            "Reference_Role": role,
                            "Reference_Group": reference_group,
                            "Trait": trait,
                            "Transfer_Mean": transfer_mean,
                            "Reference_Mean": reference_mean,
                            "Percent_Change": pct,
                            "Tukey_Analysis_Set": "Selected_and_Switching",
                            "Tukey_Adjusted_P": p_adj,
                            "Tukey_Significant": sig,
                        }
                    )
    return pd.DataFrame.from_records(records)


def build_key_comparisons(summary: pd.DataFrame, percent_change: pd.DataFrame, posthoc: pd.DataFrame) -> pd.DataFrame:
    records = []
    for _, row in percent_change.iterrows():
        records.append(
            {
                "Comparison_Type": "Switching vs source/destination selected line",
                "Analysis_Set": "Selected_and_Switching",
                "Trait": row["Trait"],
                "Group_1": row["Reference_Group"],
                "Group_2": row["Transfer_Group"],
                "Mean_1": row["Reference_Mean"],
                "Mean_2": row["Transfer_Mean"],
                "Percent_Change_group2_vs_group1": row["Percent_Change"],
                "Adjusted_P": row["Tukey_Adjusted_P"],
                "Significant": row["Tukey_Significant"],
                "Reference_Role": row["Reference_Role"],
                "Note": "Percent change = (switching - reference) / reference × 100",
            }
        )

    for temp in [23, 28, 32]:
        analysis_set = f"T0_vs_Selected_{temp}C"
        for fe in FE_LEVELS:
            t0_group = f"T0_{temp}{fe}"
            selected_group = f"{temp}{fe}"
            for trait in TRAITS:
                mean_1 = _lookup_mean(summary, t0_group, trait)
                mean_2 = _lookup_mean(summary, selected_group, trait)
                pct = ((mean_2 - mean_1) / mean_1) * 100
                # In this analysis set, Tukey group labels are Selection_Status + "_" + Fe.
                p_adj, sig = _lookup_tukey(posthoc, analysis_set, trait, f"T0_{fe}", f"Selected_{fe}")
                records.append(
                    {
                        "Comparison_Type": "Selected vs T0 ancestral line",
                        "Analysis_Set": analysis_set,
                        "Trait": trait,
                        "Group_1": t0_group,
                        "Group_2": selected_group,
                        "Mean_1": mean_1,
                        "Mean_2": mean_2,
                        "Percent_Change_group2_vs_group1": pct,
                        "Adjusted_P": p_adj,
                        "Significant": sig,
                        "Reference_Role": "T0",
                        "Note": "Percent change = (selected - T0) / T0 × 100",
                    }
                )
    return pd.DataFrame.from_records(records)


def build_supplementary_table(raw: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Build all Supplementary Table 1 sheets from raw replicate data."""
    summary = build_summary(raw)
    anova = run_anova(raw)
    posthoc = run_tukey(raw, anova)
    posthoc_sig = posthoc.loc[posthoc["Adjusted_P"] < ALPHA].copy()
    percent_change = build_percent_change(summary, posthoc)
    key = build_key_comparisons(summary, percent_change, posthoc)
    return {
        "README": build_readme(raw),
        "Raw_All_Data": raw.copy(),
        "Summary_All": summary,
        "Percent_Change": percent_change,
        "ANOVA": anova,
        "PostHoc": posthoc,
        "PostHoc_Significant": posthoc_sig,
        "Key_Comparisons": key,
    }


def write_workbook(tables: dict[str, pd.DataFrame], out_path: str | Path) -> Path:
    """Write all generated tables to one Excel workbook.

    The workbook is intentionally kept simple and reproducible: one dataframe per
    sheet, frozen header rows, and readable column widths. Formatting is minimal
    so that the file is generated quickly and remains easy to inspect in Excel.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        for sheet_name, df in tables.items():
            df.to_excel(writer, sheet_name=sheet_name, index=False)
            ws = writer.book[sheet_name]
            ws.freeze_panes = "A2"
            for col_idx, col_name in enumerate(df.columns, start=1):
                sample = [str(col_name)] + ["" if v is None or pd.isna(v) else str(v) for v in df[col_name].head(50)]
                width = min(max(max(len(x) for x in sample) + 2, 10), 42)
                ws.column_dimensions[ws.cell(row=1, column=col_idx).column_letter].width = width
    return out_path

def write_sheet_csvs(tables: dict[str, pd.DataFrame], out_dir: str | Path) -> list[Path]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for sheet_name, df in tables.items():
        path = out_dir / f"{sheet_name}.csv"
        df.to_csv(path, index=False)
        paths.append(path)
    return paths


def main(raw_path: str | Path, out_path: str | Path, csv_dir: str | Path | None = None) -> list[Path]:
    raw = load_raw_data(raw_path)
    tables = build_supplementary_table(raw)
    paths = [write_workbook(tables, out_path)]
    if csv_dir is not None:
        paths.extend(write_sheet_csvs(tables, csv_dir))
    return paths


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw", default="data/raw_croco_physiology.csv", help="Raw replicate-level CSV, or an Excel workbook with a Raw_All_Data sheet.")
    parser.add_argument("--out", default="outputs/tables/Supplementary_Table_1_generated.xlsx", help="Output Excel workbook path.")
    parser.add_argument("--csv-dir", default="outputs/tables/supplementary_table_1_csv", help="Optional folder for one CSV per sheet. Use an empty string to skip CSV export.")
    args = parser.parse_args()
    csv_dir = args.csv_dir if args.csv_dir else None
    for path in main(args.raw, args.out, csv_dir=csv_dir):
        print(path)
