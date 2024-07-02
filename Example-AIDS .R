rm(list=ls())

## load packages and functions
library(catdata)
library(robustlmm)
source("RLMM-function.R")

## load dataset
data(aids)
Y <- aids$cd4/100
aids[,2:7] <- apply(aids[,2:7], 2, scale)
X <- as.matrix(cbind(aids[,2:7], aids$time^2, aids$time^3, aids$age^2, aids$age^3, aids$cesd^2, aids$cesd^3))
X <- cbind(1, X)
Z <- cbind(1, X[,2])
p <- dim(X)[2]

ID <- aids$person
m <- length(unique(ID))
ni <- table(ID)
ID <- rep(1:m, ni)
N <- sum(ni)


## Estimation 
# ML (non-robust method)
fit_ML <- lmer(Y~X[,-1]+(Z[,-1]|ID), REML=F)


# HGD (proposed)
gam.set <- seq(0, 0.2, by=0.01)
select <- LMM_HGD_select(Y, X, Z, ID, gam.set=gam.set, print=T)
gam.opt <- max(select$gam1, select$gam2)

fit_HGD <- LMM_HGD(Y, X, Z, ID, gam=gam.opt)

set.seed(1)
boot_HGD <- LMM_HGD_boot(Y, X, Z, ID, gam=gam.opt, B=1000)


# 'robustlmm' package
fit_rml <- rlmer(Y~X[,-1]+(Z[,-1]|ID))


