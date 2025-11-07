# autodiff_phylosdm
Probabilistic Machine Learning for Phylogenetic Species Distribution Models

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Last Commit](https://img.shields.io/github/last-commit/shubhi124081/phyloSDMs2)](https://github.com/shubhi124081/phyloSDMs2/commits/main)
[![Stars](https://img.shields.io/github/stars/shubhi124081/phylo-sdms?style=social)](https://github.com/shubhi124081/phyloSDMs2/stargazers)
[![Last Commit](https://img.shields.io/github/last-commit/shubhi124081/autodiff_phylosdm)](https://github.com/shubhi124081/autodiff_phylosdm/commits/main) [![Stars](https://img.shields.io/github/stars/shubhi124081/autodiff_phylosdm?style=social)](https://github.com/shubhi124081/autodiff_phylosdm/stargazers)

## Overview

This repository implements a probabilistic framework for phylogenetic species distribution models (SDMs) using Template Model Builder (TMB) for high-performance automatic differentiation and scalable likelihood optimization. It integrates spatial, phylogenetic, and observation-effort information to jointly learn environmental niches and phylogenetic covariance across species.

## Key features

- Poisson Log-Gaussian Cox process for occurrence data
- Phylogenetic Gaussian Processes to borrow strength across related species
- TMB backend: C++ automatic differentiation and Laplace approximation
- Flexible offsets for observation effort and expert range maps
- Posterior mean intensity and relative occurrence maps
- HPC-ready: parallel prediction and memory-safe tiling

## Citation

If you use this code, please cite:

```bibtex
@misc{sharma2025phyloSDM,
    title = {Probabilistic Machine Learning for Phylogenetic Species Distribution Models},
    author = {Sharma, Shubhi},
    year = {2025},
    note = {Manuscript in preparation}
}
```

## Model summary

| Term                      | Description                                                                                    |
|---------------------------|------------------------------------------------------------------------------------------------|
| X_i β_s                   | Environmental predictors with species-specific coefficients                                   |
| f_s                       | Phylogenetic Gaussian Process random effect                                                   |
| offset_i                  | Observation effort or expert range prior                                                      |
| Σ_phylo                   | Covariance from phylogenetic distances (exponential or squared-exponential kernels)           |

The model uses Laplace approximation (via TMB) to integrate over random effects for computational efficiency.

## Implementation highlights

- C++ autodiff backend via TMB
- Parallel raster predictions and tiling for large maps
- Configurable phylogenetic kernels (exponential or squared-exponential)
- Supports simulated, temporal, and expert-informed experiments

## Repository structure

phylo-sdms/
├── scripts/
│   ├── 000-phyloGenie_functions.R
│   ├── 010-generate_data.R
│   ├── 020-tmb_model.cpp
│   ├── 030-fit_model.R
│   ├── 033-spatial_prediction.R
│   ├── 040-evaluation.R
│   └── 050-visualization.R
├── data/
├── raw_data/
├── res/
├── analysis/
└── README.md

## Installation

Install required R packages and compile the TMB model:

```r
install.packages(c("TMB", "terra", "Matrix", "ggplot2", "parallel", "data.table"))
TMB::compile("scripts/020-tmb_model.cpp")
dyn.load(dynlib("scripts/020-tmb_model"))
```
