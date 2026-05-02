# Absorbing Measurement Artifacts in Multi-Group SEM: Simulation Evidence and an Empirical Application of the Stabilizer Variable Test

[![Simulation Data](https://img.shields.io/badge/Simulation%20Data-Figshare-blue)](https://doi.org/10.6084/m9.figshare.31990182)
[![Theory Paper](https://img.shields.io/badge/Theory%20Paper-Mathematics%2014(6)%2C%201064-green)](https://doi.org/10.3390/math14061064)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository contains the analysis code for the empirical companion paper to [Yilmaz & Çene (2026)](https://doi.org/10.3390/math14061064). While the theory paper established the mathematical foundations and regression-based Monte Carlo validation of the Stabilizer Variable Test (SVT), this companion paper applies the SVT within a confirmatory factor analysis / structural equation modeling (CFA/SEM) framework using empirical data.

**Key finding:** Biological rhythm regularity (BRIAN) functions as a Type B (directional alignment) stabilizer variable for the chronotype–physical activity path, reducing cross-group parameter heterogeneity by 82.3% (CV: 175.8% → 31.1%) across 10 demographic moderators.

## Study Design

| Component | Details |
|---|---|
| Sample | N = 1,729 Turkish adults (optimized from N = 3,072 via AVE trajectory) |
| Predictor | Chronotype (MEQ, 10-item latent factor, 5 error covariances) |
| Outcome | Physical activity (YJIPAQ, Yeo-Johnson transformed IPAQ) |
| Stabilizer | Biological Rhythms (BRIAN total score) |
| Moderators | 10 demographic variables (age, gender, education, marital status, income, BMI, smoking, alcohol, chronic disease, medical nutrition therapy) |
| Estimator | MLR (lavaan); FIML used exclusively for missing data robustness condition |

## Validation Framework

The empirical application follows a five-stage validation protocol:

1. **CFA-Based Monte Carlo Simulation** (S6): 270,000 SVT applications across 540 conditions validating Type I error (mean = 0.046) and power (mean = 0.787) under latent variable estimation
2. **Primary SVT Analysis** (S7): BRIAN classified as stabilizer (Δℓ = 1.629, d = 0.981, p_boot = 0.0009, p_binom = 0.011)
3. **Comparative Negative Controls** (S8): PSQI and KIDMED tested as alternative candidates; dual-criterion alone insufficient for discrimination → C1/C2 verification required
4. **Condition Verification** (S9, S13): C1 dose-response confirmed (r = +0.405); C2 approximate (β = 0.12, within Phase 2D safe zone)
5. **Sensitivity Analyses** (S10–S14): AVE trajectory (2,061 steps), empirical error rates (Type I = 4.3%, power = 78.8%), item-level LOO (10/10 robust), robustness (6/6 dimensions robust)

## Cumulative Validation

| Source | Replications | Framework |
|---|---|---|
| Theory paper (Phases 0–3) | 949,100 | Regression-based |
| CFA simulation (S6) | 270,000 | CFA/SEM |
| AVE trajectory (S10) | 2,061 | Empirical CFA/SEM |
| Type I/Power (S11) | 2,000 | Empirical permutation/bootstrap |
| **Total** | **1,223,161** | **Both frameworks** |

## Repository Structure
├── 00_CFA_simulation.R              # CFA-based Monte Carlo (270K SVT runs, ~374 hrs)
├── 01_SVT_Empirical.R               # Main SVT analysis
├── 02_SVT_negative_controls.R       # Comparative stabilizer analysis
├── 03_C2_structural_independence.R   # C2 verification (OLS, SEM, TOST, per-group)
├── 04_AVE_trajectory.R              # Full SVT at every AVE step (~6 hrs)
├── 05_TypeI_TypeII_error.R          # Permutation/bootstrap error analysis (~4 hrs)
├── 06_item_level_analysis.R         # Leave-one-out item analysis
├── 07_MI_spectrum_analysis.R        # MI severity–stabilization relationship
└── 08_robustness.R                  # Comprehensive robustness (6 dimensions, ~3.6 hrs)

## Data Availability

The empirical datasets contain sensitive health information (sleep quality, BMI, chronic disease status, alcohol use) and are available from the corresponding author upon reasonable request, subject to ethics committee approval. The CFA-based Monte Carlo simulation results (270,000 SVT applications, fully synthetic data) are deposited on [Figshare](https://doi.org/10.6084/m9.figshare.31990182). To reproduce the empirical analyses, obtain the datasets from the corresponding author, place them in a `Datasets/` directory, and run the scripts in order (`00` through `08`).

## Key Results

### SVT Decision

| Metric | Value |
|---|---|
| Decision | **Stabilizer** |
| Mechanism | Type B (Directional Alignment) |
| Weighted mean Δℓ | 1.629 [95% CI: 0.614, 2.691] |
| Cohen's d | 0.981 |
| p (bootstrap) | 0.0009 |
| p (binomial) | 0.011 |
| Orientation Share | 0.900 |
| CV reduction | 175.8% → 31.1% (82.3% reduction) |
| Positive moderators | 9/10 |
| MI satisfied | 9/10 |

### Error Rates (Empirical)

| | Dual-Criterion | Bradley |
|---|---|---|
| Type I | 4.3% | PASS (within [2.5%, 7.5%]) |
| Power | 78.8% | Marginal (target: 80%) |

### Cross-Validation with Simulation

| Metric | Simulation (S6, M10) | Empirical (S11) |
|---|---|---|
| Type I | 0.046 | 0.043 |
| Power | 0.787 | 0.788 |

## Computational Requirements

| Analysis | Time | Cores | SVT Runs |
|---|---|---|---|
| CFA Simulation (00) | 373.9 hrs | 14 | 270,000 |
| AVE Trajectory (04) | 6.0 hrs | 14 | 2,061 |
| Type I/Power (05) | 4.1 hrs | 14 | 2,000 |
| Robustness (08) | 3.6 hrs | 14 | ~3,000 |
| **Total** | **~388 hrs** | | **~277,000** |

## Software

- **R** 4.5.2
- **lavaan** (CFA, SEM, multi-group MI testing)
- **boot** (bootstrap inference)
- **foreach / doParallel** (parallel computation)
- **readxl, dplyr, stringr** (data management)
- **ggplot2, ggrepel, gridExtra** (visualization)

## Reproducibility

1. Clone this repository
2. Obtain datasets from the corresponding author and place in `Datasets/`
3. Run scripts in order: `00` through `08`
4. Each script saves results to `RDS/` and `Outputs/` directories

**Note:** The CFA simulation (`00_CFA_simulation.R`) requires approximately 374 hours on 14 cores. Pre-computed simulation results are available on [Figshare](https://doi.org/10.6084/m9.figshare.31990182).

## Citation

```bibtex
@article{stabilizer_practice,
  title={Absorbing Measurement Artifacts in Multi-Group {SEM}: Simulation Evidence and an Empirical Application of the Stabilizer Variable Test},
  author=ANONYMIZED,
  journal={*****},
  year={2026},
  note={Under review}
}
```

### Theory Paper

```bibtex
@article{yilmaz2026stabilizer,
  title={Stabilizer Variables for Measurement Invariance--Induced Heterogeneity: Identification Theory and Testing in Multi-Group Models},
  author={Yilmaz, Salim and {\c{C}}ene, Erhan},
  journal={Mathematics},
  volume={14},
  number={6},
  pages={1064},
  year={2026},
  doi={10.3390/math14061064}
}
```

## Related

- **Theory Paper:** [Yilmaz & Çene (2026), *Mathematics*, 14(6), 1064](https://doi.org/10.3390/math14061064)
- **Theory Paper Repository:** [GitHub](https://github.com/sy142/stabilizer-variable-simulations)
- **Theory Paper Data:** [Figshare](https://doi.org/10.6084/m9.figshare.30731633)
- **Simulation Data (this study):** [Figshare](https://doi.org/10.6084/m9.figshare.31990182)

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

## Contact

**Salim Yilmaz** — Department of Healthcare Management, Acibadem Mehmet Ali Aydinlar University, Istanbul, Turkiye
