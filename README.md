# autodiff_phylosdm
Probabilistic Machine Learning for Phylogenetic Species Distribution Models

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Last Commit](https://img.shields.io/github/last-commit/shubhi124081/autodiff_phylosdm)](https://github.com/shubhi124081/autodiff_phylosdm/commits/main) 
---

### ⚠️ Project Status: Under Active Development

This repository is **under active development** and not yet production-ready.  
The current version represents a **transition phase from a Stan-based Bayesian framework to a Template Model Builder (TMB)** implementation.  

- Core algorithms and probabilistic formulations are stable, but **interfaces and functions may change**.  
- Expect **frequent updates**, restructuring, and improvements in computational efficiency.  
- Model behavior and outputs should be treated as **experimental** until the TMB workflow is fully validated against the original Stan models.  

Researchers and collaborators are welcome to explore the framework, but please **use caution for applied or large-scale analyses** until the codebase reaches a stable release.

---


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

```text
phylo-sdms/
├── scripts/
│   ├── 000-phyloGenie_functions.R
│   ├── 000-config.yaml
│   ├── 011-gen_data.R
│   ├── 010-workflow.sh
│   ├── 020-lgcp_background.cpp
│   ├── 021-buildjobs.R
│   ├── 023-main.R
│   ├── 024-runjobs.R
│   ├── 031-cond_pred.R
│   ├── 031-recursive_cond_pred.R
│   ├── 032-eval.R
│   ├── 032-recursive_eval.R
│   ├── 032-eval.R
│   ├── 034-soft_clip.R
│   ├── 034-spatial_pred.R
│   └── 035-thresholds.R
├── data/
├── raw_data/
├── res/
├── analysis/
└── README.md
```

## Installation

Install required R packages and compile the TMB model:

```r
install.packages(c("TMB", "terra", "Matrix", "ggplot2", "parallel", "data.table"))
TMB::compile("scripts/021-lgcp_background.cpp")
dyn.load(dynlib("scripts/021-lgcp_background.cpp"))
```
