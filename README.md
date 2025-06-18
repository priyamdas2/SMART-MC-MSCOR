# SMART-MC and MSCOR: Modeling Multiple Sclerosis Therapy Transitions

**Paper Title:**  
SMART-MC: Characterizing the Dynamics of Multiple Sclerosis Therapy Transitions Using a Covariate-Based Markovian Framework

## Overview

This repository contains all code and data needed to reproduce the key results and figures from the manuscript accepted for publication in the *Journal of the American Statistical Association (JASA)*. The paper develops SMART-MC, a covariate-based Markov model for estimating treatment transition probabilities in Multiple Sclerosis (MS), and proposes MSCOR, a novel global optimization routine tailored for sparse, non-convex estimation.

The analysis includes both simulated data and real-world DMT (Disease Modifying Therapy) sequences from the CLIMB cohort and other sources within the Mass General Brigham health system.

---

## Figures and Scripts

| Figure / Table | Description | Script |
|----------------|-------------|--------|
| Figure 1       | Sankey diagrams of treatment sequences | `R/Empirical_Sankey_Plot.R` |
| Figure 2       | Heatmaps of empirical transition probabilities | `R/Empirical_Heatmap.R` |
| Figure 3       | SMART-MC estimated transition curves by covariates | `R/SMART_MC_TransitionCurves.R` |
| Figure 4       | Initial treatment probabilities and odds ratios | `R/SMART_MC_Initial_OR.R` |
| Table S6       | Estimated coefficients from SMART-MC | `R/Estimate_SMART_MC_Coefficients.R` |
| Simulation Results | Simulation setup and performance evaluation | `Simulations/Run_Simulation_Results.R` |

Each script is documented with in-line comments explaining its functionality. Scripts automatically save figures to the `Plots/` folder.

---

## Requirements

- **R version**: 4.3.1  
- **R packages**:
  - `ggplot2 (>= 3.4.2)`
  - `dplyr (>= 1.1.2)`
  - `tidyr (>= 1.3.0)`
  - `readr (>= 2.1.4)`
  - `RColorBrewer (>= 1.1-3)`
  - `patchwork (>= 1.1.2)`
  - `ggalluvial (>= 0.12.5)`
  - `reshape2 (>= 1.4.4)`
  - `purrr (>= 1.0.1)`

You can install all required packages via:

```R
install.packages(c("ggplot2", "dplyr", "tidyr", "readr", "RColorBrewer", "patchwork", "ggalluvial", "reshape2", "purrr"))
