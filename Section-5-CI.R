###------------------------------------------------------###
###           R code for simulation study of             ###
###         confidence intervals in Section 5            ###
###------------------------------------------------------###
rm(list=ls())

## packages
library(MASS)
library(progress)
library(robustlmm)
source("RLMM_HGD.R")

set.seed(1)
R <- 200    # number of replications 

## Scenarios
s <- 7   # 1~9
Om <- cbind(sort(rep(c(0,0.05,0.1),3)), rep(c(0,0.05,0.1),3))
om1 <- Om[s,1]
om2 <- Om[s,2]

aa <- 9    # outlier location 


## Settings
m <- 50    # 50 or 100
nn <- c(10, 15, 20, 25, 30)   # cluster-size
meth <- c("aHGD", "HGD", "ML", "RML")
L <- length(meth)

## True parameters
V <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
Sig <- 1.5
Beta <- c(0.5, 0.3, 0.5, 0.8)

## Cluster size
ni <- rep(nn, rep(m/5, 5))
N <- sum(ni)
ID <- rep(1:m, ni)

## Covariate 
mat <- 0.4*matrix(1,3,3)+0.6*diag(3)
XX <- mvrnorm(N, rep(0, 3), mat)
x1 <- XX[,1]
x2 <- XX[,2]
x3 <- XX[,3]
X <- cbind(1, XX)
Z <- cbind(1, x2)

## Objects to store results 
CP <- array(NA, c(R, length(Beta), L))
AL <- array(NA, c(R, length(Beta), L))
IS <- array(NA, c(R, length(Beta), L))
dimnames(CP)[[3]] <- dimnames(AL)[[3]] <- dimnames(IS)[[3]] <- meth
simY <- matrix(NA, N, R)


## Monte Carlo Replications 
for(r in 1:R){
  # data generation
  U <- mvrnorm(m, c(0, 0), V)
  ch1 <- matrix(rbinom(2*m, 1, om1), m, 2)
  U[ch1==1] <- U[ch1==1] + aa
  ch2 <- rbinom(N, 1, prob=om2*2/(1+exp(-3*x1)) )
  Ep <- (1-ch2)*rnorm(N, 0, Sig) + ch2*rnorm(N, aa, 1)
  Y <- as.vector(X%*%Beta + U[ID,1]+ Z[,2]*U[ID,2] + Ep)
  simY[,r] <- Y
  
  # Hierarchical gamma-divergence with selected tuning parameter (aHGD) 
  select <- RLMM_HGD_select(Y, X, Z, ID)
  hgam <- select$gam
  rfit <- RLMM_HGD(Y, X, Z, ID, gam=hgam)
  Est <- rfit$Beta
  boot <- RLMM_HGD_boot(Y, X, Z, ID, gam=hgam, B=100)
  CI1 <- apply(boot$Boot.Beta, 2, quantile, prob=c(0.025, 0.975))
  
  # Hierarchical gamma-divergence with fixed tuning parameter (HGD) 
  rfit <- RLMM_HGD(Y, X, Z, ID, gam=0.5)
  Est <- rfit$Beta
  boot <- RLMM_HGD_boot(Y, X, Z, ID, gam=0.5, B=100)
  SE <- apply(boot$Boot.Beta, 2, sd)
  CI2 <- apply(boot$Boot.Beta, 2, quantile, prob=c(0.025, 0.975))

  # standard maximum likelihood (ML) using "lme4" package
  B <- 100
  Boot_beta_ML <- matrix(NA, B, dim(X)[2])
  pb <- progress_bar$new(total=B)   # progress bar 
  for(b in 1:B){
    BW <- as.vector(rdirichlet(1,rep(1,m)))*m
    fit_boot <- lmer(Y~X[,-1]+(Z[,-1]|ID), REML=F, weights=BW[ID])
    Boot_beta_ML[b,] <- summary(fit_boot)$coefficients[,1]
    pb$tick()  # progress
  }
  CI3 <- apply(Boot_beta_ML, 2, quantile, prob=c(0.025, 0.975))
  
  # robust maximum likelihood (RML) using "robustlmm" package
  Boot_beta_RML <- matrix(NA, B, dim(X)[2])
  pb <- progress_bar$new(total=B)   # progress bar 
  for(b in 1:B){
    BW <- as.vector(rdirichlet(1,rep(1,m)))*m
    try( rfit_boot <- rlmer(Y~X[,-1]+(Z[,-1]|ID), weights=BW[ID], rel.tol=0.1) )
    Boot_beta_RML[b,] <- summary(rfit_boot)$coefficients[,1]
    pb$tick()  # progress
  }
  CI4 <- apply(Boot_beta_RML, 2, quantile, prob=c(0.025, 0.975))

  # evaluation 
  lCI <- cbind(CI1[1,], CI2[1,], CI3[1,], CI4[1,]) 
  uCI <- cbind(CI1[2,], CI2[2,], CI3[2,], CI4[2,]) 
  CP[r,,] <- ifelse(lCI<Beta, 1, 0) * ifelse(uCI>Beta, 1, 0)
  AL[r,,] <- uCI - lCI
  IS[r,,] <- AL[r,,] + 40*ifelse(Beta>uCI, 1, 0)*(Beta-uCI) + 40*ifelse(Beta<lCI, 1, 0)*(lCI-Beta)
  
  print(r)
  print( apply(IS[r,,], 2, mean) )
}


save(list=ls(), file=paste0("sim-CI-result(s=",s,",m=",m,").RData"))


