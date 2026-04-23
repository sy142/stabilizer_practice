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
  author={Yilmaz, Salim and {\c{C}}ene, Erhan},
  journal={Multivariate Behavioral Research},
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
