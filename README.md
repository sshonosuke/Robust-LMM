# Robust Linear Mixed Models using Hierarchical Gamma Divergence

This repository provides R code implementing robust linear mixed models, as proposed by the following paper.

Sugasawa, S., Hui, K. F. C. and Welsh, A. H. (2024). Robust Linear Mixed Models using Hierarchical Gamma-Divergence. *(arXiv:2407.01883)* 

The main files are as follows: 
 
- `RLMM_HGD.R` : Implementation of the proposed robust linear mixed models with hierarchical gamma divergence 
- `Example.R`: Oneshot example of fitting RLMM to simulated data

This respository includes the following files to replicate the results in the paper:

- `Section-5-MSE.R`: Script to replicate the results of mean squared errors in Section 5
- `Section-5-MSE-Figure.R`: Script to make Figure 1 in Section 5 and Figure S1
- `Section-5-CI.R`: Script to replicate the results of credible intervals in Section 5
- `Section-5-CI-Figure.R`: Script to make Figure 2 in Section 5
- `Section-6.R`: Script for application to AIDS Cohort Study in Section 6
- `Section-S5-MSE.R`: Script to replicate the results of mean squared errors in Section S5
- `Section-S5-MSE-Figure.R`: Script to make Figure S3 in Section S5
- `Section-S6.R`: Script for application to AIDS Cohort Study with transformed response variables in Section S6
- `LMM-DPD-function.R`: Implementation of robust linear mixed models using marginal density power divergence 
