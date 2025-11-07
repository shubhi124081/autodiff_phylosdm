# autodiff_phylosdm
# Probabilistic Machine Learning for Phylogenetic Species Distribution Models

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Last Commit](https://img.shields.io/github/last-commit/shubhi124081/phylo-sdms)](https://github.com/shubhi124081/phylo-sdms/commits/main)
[![Stars](https://img.shields.io/github/stars/shubhi124081/phylo-sdms?style=social)](https://github.com/shubhi124081/phylo-sdms/stargazers)

---

## 🧠 Overview

This repository implements a **probabilistic machine learning framework** for **phylogenetic species distribution modeling (SDM)** using **Template Model Builder (TMB)** for high-performance automatic differentiation and scalable likelihood optimization.  

The framework integrates **spatial**, **phylogenetic**, and **observation-effort** information into a unified generative model that estimates environmental responses and their covariance across species. It bridges ecology, statistics, and AI—allowing **joint learning of environmental niches and phylogenetic structure** at scale.

---

## ✨ Key Features

- **Probabilistic foundation:** Poisson Log-Gaussian Cox process for occurrence data  
- **Phylogenetic Gaussian Processes:** Borrow strength across related species via covariance kernels  
- **TMB backend:** C++ automatic differentiation for scalable gradient-based optimization  
- **Flexible offsets:** Incorporates expert range maps and observation effort as additive log-offsets  
- **Interpretability:** Produces posterior mean intensity and relative occurrence maps  
- **HPC-ready:** Efficiently fits thousands of models across clusters  

---

## 📄 Citation

If you use this code, please cite:

```bibtex
@misc{sharma2025phyloSDM,
  title = {Probabilistic Machine Learning for Phylogenetic Species Distribution Models},
  author = {Sharma, Shubhi},
  year = {2025},
  note = {Manuscript in preparation}
}

# Model Summary 
| Term                      | Description                                                                                    |
| ------------------------- | ---------------------------------------------------------------------------------------------- |
| ( X_i \beta_s )           | Environmental predictors with species-specific coefficients                                    |
| ( f_s )                   | Phylogenetic Gaussian Process random effect                                                    |
| ( \text{offset}_i )       | Observation effort or expert range prior                                                       |
| ( \Sigma_{\text{phylo}} ) | Covariance derived from phylogenetic distance using exponential or squared-exponential kernels |

The model uses Laplace approximation within TMB to efficiently integrate over random effects, combining Bayesian flexibility with computational scalability.

🚀 Implementation Highlights

C++ autodiff backend via TMB

Parallel raster predictions and memory-safe tiling

Configurable phylogenetic kernels (exponential or squared-exponential)

Supports artificial, temporal, and expert-informed experiments

Designed for large-scale ecological applications