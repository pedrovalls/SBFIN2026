# Monte Carlo Valuation Experiment – Replication Package

This repository contains the code required to replicate the Monte Carlo results presented in the paper:

> "Stochastic Valuation with Regime Dynamics"

## 📦 Contents
# Monte Carlo Valuation Experiment – Replication Package

This repository contains the code required to replicate the Monte Carlo results presented in the paper:

> "Stochastic Valuation with Regime Dynamics"

## 📦 Contents

The script:

- Simulates cash-flow paths under a 3-state Markov chain
- Computes true firm value (finite-horizon present value)
- Estimates valuation using:
  - Traditional DCF
  - Regime-adaptive (RS) method
- Computes:
  - Bias
  - MAE
  - RMSE
  - MAPE
  - P(APE > 10%)
- Performs paired t-tests on quadratic loss
- Computes bootstrap confidence intervals
- Generates:
  - Histograms
  - Boxplots
  - Empirical CDF (3-panel figure)
- Exports LaTeX-ready tables

---

## 🔧 Software Requirements

- R version 4.2 or higher
- Base R packages only

Optional:
- ggplot2 (if alternative visualization is desired)

---

## ▶️ How to Run

1. Clone the repository: https://github.com/pedrovalls/SBFIN2026

2. Open `simulacao_com_tudo.R` in R or RStudio.

3. Run the script entirely.

The script will:

- Print the results table in the console
- Generate PDF figures in the working directory
- Produce LaTeX-ready table output

---

## 🎯 Simulation Design

Three environments are implemented:

1. High persistence
2. High volatility
3. High regime switching (low persistence)

Each scenario runs 10,000 Monte Carlo replications.

---

## 🔁 Reproducibility

The script fixes the random seed:
set.seed(42)

Results should replicate exactly across systems (minor floating-point differences may occur).

---

## 📊 Output

Generated files include:

- `figure3_cdf_3panels.pdf`
- Scenario-specific PDFs
- LaTeX table output printed in console

---


