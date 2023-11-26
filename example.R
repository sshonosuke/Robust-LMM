rm(list=ls())
set.seed(1)

# load packages 
source("RLMM_HGD.R")
library(MASS)

# Scenarios
s <- 9   # 1~9
Om <- cbind(sort(rep(c(0,0.05,0.1),3)), rep(c(0,0.05,0.1),3))
om1 <- Om[s,1]
om2 <- Om[s,2]
aa <- 9    # outlier location 

## Settings
# sample size
m <- 50    # 50 or 100
nn <- c(10, 15, 20, 25, 30)   # cluster-size

ni <- rep(nn, rep(m/5, 5))
N <- sum(ni)
ID <- rep(1:m, ni)

# true parameters
V <- matrix(c(1, 0.3, 0.3, 1), 2, 2)   
Sig <- 1.5
Beta <- c(0.5, 0.3, 0.5, 0.8)

# covariates
mat <- 0.4*matrix(1,3,3)+0.6*diag(3)
XX <- mvrnorm(N, rep(0, 3), mat)
x1 <- XX[,1]
x2 <- XX[,2]
x3 <- XX[,3]
X <- cbind(1, XX)
Z <- cbind(1, x2)


## data generation
U <- mvrnorm(m, c(0, 0), V)
ch1 <- matrix(rbinom(2*m, 1, om1), m, 2)
U[ch1==1] <- U[ch1==1] + aa
ch2 <- rbinom(N, 1, prob=om2*2/(1+exp(-3*x1)) )
Ep <- (1-ch2)*rnorm(N, 0, Sig) + ch2*rnorm(N, aa, 1)
Y <- as.vector(X%*%Beta + U[ID,1]+ Z[,2]*U[ID,2] + Ep)


## standard LMM
fit_ML <- lmer(Y~X[,-1]+(Z[,-1]|ID), REML=F)
Beta_ML <-  summary(fit_ML)$coefficients[,1]
Sig_ML <- summary(fit_ML)$sigma    # error standard deviation
R_ML <- as.matrix( summary(fit_ML)$varcor$ID )    # covariance matrix of random effects
RE_ML <- cbind(coef(fit_ML)$ID[,2]-Beta_ML[1], coef(fit_ML)$ID[,1])   # estimated random effects


## Hierarchical gamma-divergence 
select <- RLMM_HGD_select(Y, X, Z, ID)
(gam_opt <- select$gam)
fit_HGD <- RLMM_HGD(Y, X, Z, ID, gam=gam_opt)
Beta_HGD <- fit_HGD$Beta
Sig_HGD <- fit_HGD$Sig
R_HGD <- fit_HGD$R
RE_HGD <- matrix(fit_HGD$RE, m, 2)


## Estimates 
cbind(Beta, Beta_ML, Beta_HGD)   # coefficients for fixed effects
c(Sig, Sig_ML, Sig_HGD)    # error standard deviation 
list(V, R_ML, R_HGD)      # covariance matrix of random effects


## MSE of random effects
mean((RE_ML[,1]-U[,1])^2)
mean((RE_HGD[,1]-U[,1])^2)
mean((RE_ML[2]-U[,2])^2)
mean((RE_HGD[,2]-U[,2])^2)


## plot (estimated random effects)
par(mfcol=c(1,2))
# random intercept 
plot(U[,1], RE_ML[,1], col=4, main="random intercept")
points(U[,1], RE_HGD[,1], col=2, pch=20)
legend("topleft", c("ML", "HGD"), col=c(4,2), pch=c(1,20))
abline(0, 1)
# random slope
plot(U[,2], RE_ML[,2], col=4, main="random slope")
points(U[,2], RE_HGD[,2], col=2, pch=20)
legend("topleft", c("ML", "HGD"), col=c(4,2), pch=c(1,20))
abline(0, 1)





## confidence intervals 
SE_ML <- summary(fit_ML)$coefficients[,2]
CI_ML <- cbind(Beta_ML-qnorm(0.975)*Beta_ML, Beta_ML+qnorm(0.975)*Beta_ML)

HGD_Boot <- RLMM_HGD_boot(Y, X, Z, ID, gam=gam_opt, B=100)  # Bootstrap of HGD
CI_HGD <- t( apply(HGD_Boot$Boot.Beta, 2, quantile, prob=c(0.025, 0.975)) )

CI_ML
CI_HGD
