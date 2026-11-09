# Crocosphaera Thermal Experimental Evolution

**Code and data analysis for investigating long-term thermal selection and short-term physiological responses of the marine nitrogen-fixing cyanobacterium *Crocosphaera watsonii* under contrasting iron availability.**

## Overview

This repository accompanies the manuscript:

**Duan, R., Wrightson, L., Yang, N., et al.**  
*Long and short term thermal responses of iron limited ocean nitrogen fixers produce divergent biogeographic projections.*

The study investigates how long-term thermal selection and short-term temperature acclimation influence the physiological performance of *Crocosphaera watsonii* under iron-replete and iron-limited conditions, and how these timescale-dependent responses affect projections of marine nitrogen fixation.

The experimental work combines two years of laboratory selection, reciprocal temperature-transfer experiments, physiological measurements, and Earth system modeling using PISCES-QUOTA.

## Experimental Design

*Crocosphaera watsonii* WH0005 was cultured under a factorial combination of:

- **Temperatures:** 23°C, 28°C, and 32°C
- **Iron concentrations:** 5 nM (Fe-limited) and 250 nM (Fe-replete)
- **Replicates:** Three biological cultures per temperature-iron treatment
- **Long-term selection:** Approximately two years

The experimental dataset distinguishes four stages:

1. **Ancestral (T0):** The original stock culture sampled at 28°C under Fe-replete conditions at the beginning of the experiment.
2. **Initial acclimation:** Cultures characterized after the initial three-month period under each temperature–iron treatment.
3. **Long-term selection (Tfinal):** Cultures maintained under their respective temperature–iron conditions for approximately two years.
4. **Reciprocal temperature transfers:** Long-term-selected cultures transferred to alternative temperatures while maintaining their original iron treatments, followed by 14 days of acclimation at the target temperature.

Physiological responses were characterized using specific growth rates, nitrogen fixation normalized to particulate organic nitrogen (PON), and carbon fixation normalized to particulate organic carbon (POC).

## Repository Structure

| Directory | Description |
|---|---|
| `scripts/` | Current R and Python scripts for physiological data visualization and multivariate analyses |
| `archive/2026-05_original/` | Archived scripts from an earlier manuscript and figure-development stage |

The archived scripts are retained for historical reference and should not be assumed to reproduce the current manuscript figures.

## Analysis and Figure Scripts

### Main Figures

| Script | Description |
|---|---|
| `build_figure2.R` | Physiological responses of ancestral, initial-acclimation, and long-term-selected cultures under contrasting temperature and iron conditions |
| `build_figure3.R` | Growth, nitrogen fixation, and carbon fixation responses following reciprocal temperature transfers, displayed at matched measurement temperatures |

### Supplementary Figures

| Script | Description |
|---|---|
| `prepare_multivariate.py` | Constructs treatment-level physiological profiles, standardizes traits, and performs hierarchical clustering |
| `build_supp1_cluster.R` | Generates the hierarchical-clustering dendrogram of physiological response profiles (Supplementary Fig. 1) |
| `supp_fig2_heatmap.py` | Generates the standardized physiological-response heatmap (Supplementary Fig. 2) |
| `build_supplementary_3_5.R` | Generates full-treatment growth, nitrogen-fixation, and carbon-fixation plots (Supplementary Figs. 3–5) |
| `build_supplementary6.R` | Visualizes percentage changes relative to long-term-selected cultures at the original selection temperature (Supplementary Fig. 6) |
| `build_T0_letters.R` | Generates the growth-rate comparison between the ancestral T0 baseline and selected reciprocal-transfer treatments at 28°C under Fe-replete conditions (Supplementary Fig. 7) |

### Supporting Scripts

`vector_export.R` provides functions for exporting vector graphics in PDF and SVG formats.

## Statistical Analyses

The physiological analyses include comparisons across temperature treatments, iron regimes, experimental stages, and reciprocal-transfer histories.

The study employs analysis of variance (ANOVA) and Tukey–Kramer post hoc comparisons, with adjustments for multiple testing.

The manuscript's Supplementary Table 1 documents the underlying physiological measurements, treatment summaries, percentage changes, statistical tests, and pairwise comparisons.

Experimental replicates are treated as independent biological cultures. Missing measurements are not imputed in the multivariate preprocessing script.

## Software and Dependencies

The current visualization and analysis scripts use both **R** and **Python**.

Main R dependencies include:

- `ggplot2`
- `patchwork`
- `svglite`
- `ragg`

Python dependencies include:

- `numpy`
- `pandas`
- `scipy`
- `matplotlib`

The scripts also reference supporting plotting utilities and source-data files. Package versions and complete environment specifications should be recorded alongside the finalized reproducible analysis release.

## Data Availability and Reproducibility

The physiological dataset, including raw measurements, summary statistics, and statistical comparisons, is documented in Supplementary Table 1 of the associated manuscript.

The current repository contains the updated figure-generation and analysis scripts, but does not yet include all data inputs and supporting utilities required for complete end-to-end execution.

Individual scripts use specific input file paths and output directories. These paths and dependencies must be verified when assembling a fully reproducible analysis environment.

## Earth System Modeling

The Earth system modeling component uses a modified version of PISCES-QUOTA to examine the consequences of short- and long-timescale thermal responses for future marine nitrogen fixation.

The corresponding model implementation is available separately:

**PISCES-QUOTA-P6Z-EVO:**  
https://github.com/lewiswrightson/PISCES-QUOTA-P6Z-EVO

Supporting model datasets are available through Zenodo:

**Model data DOI:**  
https://doi.org/10.5281/zenodo.20395566

## Citation

If you use these analyses or scripts, please cite the associated manuscript once its publication details are available.

For the Earth system modeling framework, see:

Wrightson, L., Yang, N., Mahaffey, C., Hutchins, D. A., & Tagliabue, A. (2022). Integrating the impact of global change on the niche and physiology of marine nitrogen-fixing cyanobacteria. *Global Change Biology*, 28, 7078–7093. https://doi.org/10.1111/gcb.16399

## Contact

**Ran Duan**  
University of Southern California  
GitHub: https://github.com/duanrandr

For questions regarding the analyses or code, please open an issue in this repository.
