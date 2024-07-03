# Robust Linear Mixed Models using Hierarchical Gamma Divergence

This repository provides R code implementing robust linear mixed models, as proposed by the following paper.

Sugasawa, S., Hui, K. F. C. and Welsh, A. H. (2024). Robust Linear Mixed Models using Hierarchical Gamma-Divergence. *(arXiv)* 

The repository includes the following files.

- `RLMM-HGD.R` : Implementation of the proposed robust linear mixed models (RLMM) 
- `Example-Sim.R`: Example of fitting RLMM to simulated data
- `Example-AIDS.R`: Example of fitting RLMM to a multi-center AIDS cohort study
- `Simulation-Estimation.R`: Monte Carlo simulation study for parameter estimation (presented in Section 5)
- `Simulation-Inference.R`: Monte Carlo simulation study for confidence intervals (presented in Section 5)
- `LMM-DPD-function.R`: Implementation of robust linear mixed models using marginal density power divergence 
