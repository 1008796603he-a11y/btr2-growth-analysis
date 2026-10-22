# Btr2 growth-rate analysis (Supplementary File S1)

R script for the statistics and boxplot of specific growth rates of the nine constructs in the thesis.

## Files
- `Supplementary_File_S1_final_growth_analysis.R` – analysis script
- `comb2_growth_rates.csv` – per-well specific growth rates (141 wells; columns `Construct`, `Run`, `Well`, `mu` in h^-1). Input of Section B of the script.
- `LICENSE` – MIT

## Run
```r
install.packages(c("data.table", "ggplot2", "multcompView"))
source("Supplementary_File_S1_final_growth_analysis.R")   # run from this directory
```
Writes `comb2_growth_rate_summary.csv` (mean, SD, n, Tukey letters per construct) and `growth_rate_boxplot_final.png`.
Tested with R 4.5.3, data.table 1.17.8, ggplot2 4.0.3, multcompView 0.1.12.

## What the script does
- Per well, ln(OD600) (not blank-subtracted) is regressed on time over the eligible window with the highest R^2 (at least 12 consecutive points, slope > 0, R^2 >= 0.95). The resulting value is a window-fit specific growth rate, not a maximum growth rate.
- One-way ANOVA (`aov(mu ~ Construct)`), then Tukey HSD and compact letter display. Run is not a model term.
- No normality, variance-homogeneity or residual diagnostics are performed in this script.

## Data composition
Two independent runs (27 Aug and 19 Sep) were combined. 9 constructs x 8 wells x 2 runs = 144 wells, minus 3 non-growing wells from the 27 Aug plate that were removed manually, gives the 141 wells in the table. The Aug 15 B2P3 clone was excluded (suspected wrong plasmid).

## Scope
Section B reproduces the reported statistics and letters from the rate table. Section A (rate estimation from raw plate-reader data) documents the window rule; the raw OD600 files are not included here.

## Construct labels
Data keys in the table differ from labels on the figure: `CO7` is shown as "2.5 kb insert control", `Btr1` as "5 kb insert control", and `B2P3-AGA` as "B2P3-T>A".
