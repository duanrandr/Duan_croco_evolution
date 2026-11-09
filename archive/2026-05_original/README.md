# Crocosphaera thermal selection figures and Supplementary Table 1

This repository contains Python code for regenerating the experimental-data figures, statistical outputs, and Supplementary Table 1 for the Crocosphaera thermal selection manuscript.

## What this code generates

### Figures

The script `scripts/run_all_figures.py` generates:

- `Fig2_longterm_responses`: ancestral (T0) and long-term selected physiological responses.
- `Fig3_percent_change`: short-term reciprocal transfer percent-change responses.
- `SuppFig1_hierarchical_clustering`: clustering of treatment-level physiological profiles.
- `SuppFig2_standardized_heatmap`: standardized z-score heatmap.
- `SuppFig3_growth_responses`: full-treatment growth-rate bar plot.
- `SuppFig4_nfix_responses`: full-treatment N-specific N2-fixation bar plot.
- `SuppFig5_cfix_responses`: full-treatment C-specific C-fixation bar plot.

### Supplementary Table 1

The script `scripts/generate_supp_table1.py` rebuilds the complete multi-sheet Supplementary Table 1 workbook from the raw replicate-level data in `data/raw_croco_physiology.csv`. The generated workbook contains:

- `README`
- `Raw_All_Data`
- `Summary_All`
- `Percent_Change`
- `ANOVA`
- `PostHoc`
- `PostHoc_Significant`
- `Key_Comparisons`

The same script can also export each sheet as a CSV file for easier version control and review.


## Data

`data/raw_croco_physiology.csv` is the raw replicate-level physiology table used to regenerate Supplementary Table 1. 
Model simulation data/code should be cited separately in the manuscript's Data/Code Availability statements.
