# SMART-MC and MSCOR: Reproducibility Repository

**Paper Title:**  
SMART-MC: Characterizing the Dynamics of Multiple Sclerosis Therapy Transitions Using a Covariate-Based Markov Model

This repository provides the necessary code and documentation for reproducing the results in the article listed above. It implements the SMART-MC model for estimating covariate-driven treatment transitions in Multiple Sclerosis, and introduces MSCOR, a novel global optimization routine (over multiple spherically constrained parameter space) used in this context.

---

## 📄 Paper Description

The manuscript introduces:
- **SMART-MC**, a penalized covariate-driven Markov model to estimate dynamic treatment transitions in MS
- **MSCOR**, a spherically constrained optimization routine for high-dimensional non-convex estimation
- Applications to real-world EHR data and extensive simulation studies
- Figures, tables, and diagnostics for estimation performance and phenotypic treatment heterogeneity

---

## 🧮 Figures and Tables in the Paper

| Output        | Description                                | Script Path                            |
|---------------|--------------------------------------------|----------------------------------------|
| **Figure 1**  | Sankey diagrams of treatment sequences     | `Real Data Analysis/exploratory_analysis.R` |
| **Figures 5–6a** | Estimated transition probabilities (SMART-MC) | `Real Data Analysis/SMART_MC_Var_effect_plot.R` |
| **Figure 6b** | Odds ratios for across-treatment transitions | `SMART_MC_ODDS_ratio_calculation.m` → `SMART_MC_Odds_ratio_plot.R` |
| **Table 1**   | MSCOR benchmark results                    | `MSCOR Benchmark/MSCOR_Benchmark_comparison.m` → `MSCOR_post_evaluation.m` |
| **Tables S2–S5** | Simulation results                       | `Simulation Study/` (see details in `SMART_MC_Reproducibility_and_DEMO_instructions.pdf`)        |
| **Table S6**, **Figure S2** | SMART-MC estimated coefficients | `Real Data Analysis/SMART_MC_Real_data.m` |

---

## 📦 Requirements

- **R version**: 4.3.1  
- **MATLAB**: R2022a or later  
- **Required R packages**:
  - `ggplot2`, `dplyr`, `tidyr`, `readr`, `patchwork`, `RColorBrewer`, `ggalluvial`, `reshape2`, `purrr`

To install all required R packages:

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "readr", "RColorBrewer", 
                   "patchwork", "ggalluvial", "reshape2", "purrr"))
