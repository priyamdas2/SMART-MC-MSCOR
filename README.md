# SMART-MC and MSCOR: Reproducibility Repository

**Paper Title:**  
SMART-MC: Characterizing the Dynamics of Multiple Sclerosis Therapy Transitions Using a Covariate-Based Markov Model 

This repository provides the necessary code and documentation for reproducing the results presented in the article above. It implements the **SMART-MC** model for estimating covariate-driven treatment transitions in Multiple Sclerosis, and introduces **MSCOR**, a novel global optimization routine designed for optimization over multiple spherically constrained parameter spaces.

---

## 📄 Paper Description

The manuscript introduces:

- **SMART-MC**, a covariate-driven Markov model for estimating dynamic treatment transition probabilities in MS as a function of patient covariates.
- **MSCOR**, a spherically constrained optimization routine for globally optimizing black-box functions over parameter spaces constrained to collections of unit spheres.
- A benchmark study demonstrating the superiority of **MSCOR** over existing global optimization techniques such as Genetic Algorithms and Simulated Annealing.
- A comprehensive simulation study of **SMART-MC**, powered by **MSCOR**, evaluating estimation performance under various scenarios.
- An application of **SMART-MC** to real-world data to identify the influence of key clinical and demographic factors on both within- and across-treatment transition probabilities for MS-DMTs.

---

## 🧮 Figures and Tables in the Paper

To reproduce the tables and figures presented in the paper, please refer to `SMART_MC_Reproducibility_and_DEMO_instructions.pdf`. A brief overview is provided below.


| Output        | Description                                | Script Path                            |
|---------------|--------------------------------------------|----------------------------------------|
| **Figure 1**  | Sankey diagrams of treatment sequences     | `Real Data Analysis/exploratory_analysis.R` |
| **Figures 5,6a** | Estimated transition probabilities (SMART-MC) | `Real Data Analysis/SMART_MC_Var_effect_plot.R` |
| **Figure 6b** | Odds ratios for across-treatment transitions | `SMART_MC_ODDS_ratio_calculation.m` → `SMART_MC_Odds_ratio_plot.R` |
| **Table 1, S1**, **Figure S1**   | MSCOR benchmark results                    | `MSCOR Benchmark/MSCOR_Benchmark_comparison.m` → `MSCOR_post_evaluation.m` |
| **Tables S2–S5** | Simulation results                       | `Simulation Study/` (see details in `SMART_MC_Reproducibility_and_DEMO_instructions.pdf`)        |
| **Table S6**, **Figure S2** | SMART-MC estimated coefficients; simulated treatment trajectory | `Real Data Analysis/SMART_MC_Real_data.m` |

---

## 📦 Requirements

- **R version**: 4.3.1  
- **MATLAB**: R2022a or later  
- **Required R packages**:
  - `ggplot2`, `dplyr`, `tidyr`, `readr`, `patchwork`, `RColorBrewer`, `ggalluvial`, `reshape2`, `purrr`

