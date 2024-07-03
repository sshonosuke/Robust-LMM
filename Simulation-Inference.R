###------------------------------------------------------###
###           Monte Carlo simulation for                 ###
###        evaluating confidence intervals               ###
###------------------------------------------------------###
rm(list=ls())

# load packages 
library(MASS)
library(robustlmm)
source("RLMM_HGD.R")

# Preparation 
s <- 9    # Scenarios: 1-9
m <- 50   # number of clusters 

set.seed(1)
R <- 200    # number of replications 


## Settings
# contamination ratio
Om <- cbind(sort(rep(c(0,0.05,0.1),3)), rep(c(0,0.05,0.1),3))
om1 <- Om[s,1]
om2 <- Om[s,2]
aa <- 9    # outlier location 

# cluster and within-cluster sample size
m <- 50    # number of clusters (50 or 100)
nn <- c(10, 15, 20, 25, 30)   # cluster-size
meth <- c("ML", "HGD1", "HGD1-p", "HGD2", "HGD2-p", "RML")
L <- length(meth)

# true parameters 
V <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
Sig <- 1.5
Beta <- c(0.5, 0.3, 0.5, 0.8)

# cluster ID 
ni <- rep(nn, rep(m/5, 5))
N <- sum(ni)
ID <- rep(1:m, ni)

# covariate matrix
mat <- 0.4*matrix(1,3,3)+0.6*diag(3)
XX <- mvrnorm(N, rep(0, 3), mat)
x1 <- XX[,1]
x2 <- XX[,2]
x3 <- XX[,3]
X <- cbind(1, XX)
Z <- cbind(1, x2)

# Objects for storing estimation results 
CP <- array(NA, c(R, length(Beta), L))
AL <- array(NA, c(R, length(Beta), L))
IS <- array(NA, c(R, length(Beta), L))
dimnames(CP)[[3]] <- dimnames(AL)[[3]] <- dimnames(IS)[[3]] <- meth
simY <- matrix(NA, N, R)


## Monte Carlo replications 
for(r in 1:R){
  # data generation
  U <- mvrnorm(m, c(0, 0), V)
  ch1 <- matrix(rbinom(2*m, 1, om1), m, 2)
  U[ch1==1] <- U[ch1==1] + aa
  ch2 <- rbinom(N, 1, prob=om2*2/(1+exp(-3*x1)) )
  Ep <- (1-ch2)*rnorm(N, 0, Sig) + ch2*rnorm(N, aa, 1)
  Y <- as.vector(X%*%Beta + U[ID,1]+ Z[,2]*U[ID,2] + Ep)
  simY[,r] <- Y
  
  # standard LMM
  fit <- lmer(Y~X[,-1]+(Z[,-1]|ID), REML=F)
  Est <-  summary(fit)$coefficients[,1]
  SE <- summary(fit)$coefficients[,2]
  CI1 <- cbind(Est-1.96*SE, Est+1.96*SE)
  
  # proposed method 1
  select <- RLMM_HGD_select(Y, X, Z, ID)
  gam_opt <- select$gam
  rfit <- RLMM_HGD(Y, X, Z, ID, gam=gam_opt)
  Est <- rfit$Beta
  boot <- RLMM_HGD_boot(Y, X, Z, ID, gam=gam_opt, B=100)
  SE <- apply(boot$Boot.Beta, 2, sd)
  CI2 <- cbind(Est-1.96*SE, Est+1.96*SE)
  CI2.p <- apply(boot$Boot.Beta, 2, quantile, prob=c(0.025, 0.975))
  
  # proposed method 2
  rfit <- RLMM_HGD(Y, X, Z, ID, gam=0.5)
  Est <- rfit$Beta
  boot <- RLMM_HGD_boot(Y, X, Z, ID, gam=0.5, B=100)
  SE <- apply(boot$Boot.Beta, 2, sd)
  CI3 <- cbind(Est-1.96*SE, Est+1.96*SE)
  CI3.p <- apply(boot$Boot.Beta, 2, quantile, prob=c(0.025, 0.975))

  # 'robustlmm' package
  rfit <- fit
  try(rfit <- rlmer(Y~X[,-1]+(Z[,-1]|ID)))
  Est <- summary(rfit)$coefficients[,1]
  SE <- summary(rfit)$coefficients[,2]
  CI4 <- cbind(Est-1.96*SE, Est+1.96*SE)
  
  # Summary
  lCI <- cbind(CI1[,1], CI2[,1], CI2.p[1,], CI3[,1], CI3.p[1,], CI4[,1]) 
  uCI <- cbind(CI1[,2], CI2[,2], CI2.p[2,], CI3[,2], CI3.p[2,], CI4[,2]) 
  CP[r,,] <- ifelse(lCI<Beta, 1, 0) * ifelse(uCI>Beta, 1, 0)
  AL[r,,] <- uCI - lCI
  IS[r,,] <- AL[r,,] + 40*ifelse(Beta>uCI, 1, 0)*(Beta-uCI) + 40*ifelse(Beta<lCI, 1, 0)*(lCI-Beta)
  
  # Displays the current number of replications
  print(r)
}


# Save results
save(list=ls(), file=paste0("sim-CI-result(s=", s, ",m=", m, ").RData"))


