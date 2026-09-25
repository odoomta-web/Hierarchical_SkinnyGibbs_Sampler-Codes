# Hierarchical Skinny Gibbs Sampler for Logistic Regression with Pólya–Gamma Latent Variables

Code and simulation studies accompanying:

> Hierarchical skinny Gibbs sampler in logistic regression using Pólya–Gamma latent variables.
> *Statistics and Its Interface*, 19(2), pp. 179–196.

## Overview

This repository implements two Gibbs samplers for Bayesian variable selection in
high-dimensional logistic regression under a hierarchical spike-and-slab prior,
using Pólya–Gamma data augmentation:

- **H-SGPG** — Hierarchical **Skinny** Gibbs sampler (scalable approximation).
- **H-EGPG** — Hierarchical **Exact** Gibbs sampler (full conditional updates).

The simulation study compares both samplers against the *t*-approximation Skinny
Gibbs sampler of Narisetty, Shen & He (2019, JASA), and three frequentist penalized
methods: lasso, MCP, and SCAD.

## Repository structure

```
R/
├── simulate_data.R          Functions generating high-dimensional logistic data
├── generate_datasets.R      Driver script: generates & saves simulation replicates
├── gibbs_skinny_pg.R        hsgpg_func() — Hierarchical Skinny Gibbs w/ Pólya-Gamma
├── gibbs_exact_pg.R         hegpg_func() — Hierarchical Exact Gibbs w/ Pólya-Gamma
├── evaluation_metrics.R     evaluation() — Sensitivity, Specificity, MCC, MSPE
└── run_model_comparison.R   Driver script: fits all methods & aggregates results
data/       (generated .RData datasets land here; not tracked in git)
results/    (aggregated result tables land here; not tracked in git)
```

## Requirements

R packages (all from CRAN except `skinnybasad`):

```r
install.packages(c("MASS", "invgamma", "BayesLogit", "mvnfast", "Matrix",
                    "truncnorm", "mvtnorm", "glmnet", "ncvreg"))
```

`skinnybasad` is **not on CRAN**. Install it from the supplementary materials of
Narisetty, Shen & He (2019), *"Skinny Gibbs: A Consistent and Scalable Gibbs
Sampler for Model Selection,"* JASA:
<https://www.tandfonline.com/doi/suppl/10.1080/01621459.2018.1482754?scroll=top>

## Reproducing the simulation study

1. Create output folders (or edit paths in the scripts to point elsewhere):
   ```r
   dir.create("data", showWarnings = FALSE)
   dir.create("results", showWarnings = FALSE)
   ```
2. Generate simulated datasets:
   ```r
   source("R/generate_datasets.R")
   ```
   This writes one `.RData` file per replicate under
   `data/<scenario>/dflist_<i>.RData`. The default scenario in the script is
   `p = 250`, 8 active covariates, strong signal, both independent and
   correlated designs (see comments in `generate_datasets.R` for how to change
   `p`, the number of active covariates, signal strength, correlation
   structure, and number of replicates).
3. Fit and compare all methods:
   ```r
   source("R/run_model_comparison.R")
   ```
   This loads the simulated datasets, fits H-SGPG, H-EGPG, *t*-approximation
   Skinny/Exact Gibbs, lasso, MCP, and SCAD, evaluates each against the known
   active-covariate set, and saves an aggregated results table to
   `results/result_p_<p>_<p_act>.RData`.


