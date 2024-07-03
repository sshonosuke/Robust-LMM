###------------------------------------------------------###
###           Monte Carlo simulation for                 ###
###          evaluating point estimation of              ###
###        model parameters and random effects           ###
###------------------------------------------------------###
## Remark: 'heavy' package should be directly downloaded from https://cran.r-project.org/web/packages/heavy/index.html
rm(list=ls())


# load packages 
library(MASS)
library(robustlmm)
library(heavy)   
source("RLMM_HGD.R")
source("LMM-DPD-function.R")


# Preparation 
s <- 9    # Scenarios: 1-9
m <- 50   # number of clusters 

set.seed(1)
R <- 500    # number of replications 


## Settings
# contamination ratio
Om <- cbind(sort(rep(c(0,0.05,0.1),3)), rep(c(0,0.05,0.1),3))
om1 <- Om[s,1]
om2 <- Om[s,2]
aa <- 7     # outlier location 

# cluster and within-cluster sample size
nn <- c(10, 15, 20, 25, 30)   # cluster-size
meth <- c("ML", "aHGD", "HGD", "RLMM", "HT", "DPD")
L <- length(meth)

# true parameters    
V <- matrix(c(1, 0.3, 0.3, 1), 2, 2)    # variance-covariance matrix of random effects
Sig <- 1.5    # error standard deviation 
Beta <- c(0.5, 0.3, 0.5, 0.8)     # regression coefficient 

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


## Objects for storing estimation results 
Gam <- c()
Beta.est <- array(NA, c(R, 4, L))
V.est <- array(NA, c(R, 3, L))
Sig.est <- matrix(NA, R, L)
MSE.RE1 <- matrix(NA, R, L)
MSE.RE2 <- matrix(NA, R, L)
CPT <- matrix(NA, R, L)
dimnames(Beta.est)[[3]] <- dimnames(V.est)[[3]]<- dimnames(Sig.est)[[2]] <- meth
dimnames(MSE.RE1)[[2]] <- dimnames(MSE.RE1)[[2]] <- dimnames(CPT)[[2]] <- meth


## Monte Carlo replications 
for(r in 1:R){
  # data generation
  U <- mvrnorm(m, c(0, 0), V)
  ch1 <- matrix(rbinom(2*m, 1, om1), m, 2)
  U[ch1==1] <- U[ch1==1] + aa
  ch2 <- rbinom(N, 1, prob=om2*2/(1+exp(-3*x1)) )
  Ep <- (1-ch2)*rnorm(N, 0, Sig) + ch2*rnorm(N, aa, 1)
  Y <- as.vector(X%*%Beta + U[ID,1]+ Z[,2]*U[ID,2] + Ep)
  
  # standard LMM
  tt <- proc.time()
  fit <- lmer(Y~X[,-1]+(Z[,-1]|ID), REML=F)
  CPT[r, 1] <- ( proc.time() - tt )[3]
  
  hBeta <-  summary(fit)$coefficients[,1]
  Beta.est[r,,1] <- hBeta
  hRE <- cbind(coef(fit)$ID[,2]-hBeta[1], coef(fit)$ID[,1]) 
  MSE.RE1[r, 1] <- mean((hRE[,1]-U[,1])^2)
  MSE.RE2[r, 1] <- mean((hRE[,2]-U[,2])^2)
  Sig.est[r,1] <- summary(fit)$sigma
  hV <- as.matrix( summary(fit)$varcor$ID )
  V.est[r,,1] <- hV[c(1,2,4)]
  
  # proposed method 1 (estimated gamma)
  tt <- proc.time()
  select <- RLMM_HGD_select(Y, X, Z, ID)
  gam_opt <- select$gam
  fit1 <- RLMM_HGD(Y, X, Z, ID, gam=gam_opt)
  CPT[r, 2] <- ( proc.time() - tt )[3]
  
  Gam[r] <- gam_opt
  Beta.est[r,,2] <- fit1$Beta
  hRE <- matrix(fit1$RE, m, 2)
  MSE.RE1[r, 2] <- mean((hRE[,1]-U[,1])^2)
  MSE.RE2[r, 2] <- mean((hRE[,2]-U[,2])^2)
  Sig.est[r, 2] <- fit1$Sig
  V.est[r,,2] <- fit1$R[c(1,2,4)]
  
  # proposed method 2 (fixed gamma)
  tt <- proc.time()
  fit2 <- RLMM_HGD(Y, X, Z, ID, gam=0.5)
  CPT[r, 3] <- ( proc.time() - tt )[3]
  
  Beta.est[r,,3] <- fit2$Beta
  hRE <- matrix(fit2$RE, m, 2)
  MSE.RE1[r, 3] <- mean((hRE[,1]-U[,1])^2)
  MSE.RE2[r, 3] <- mean((hRE[,2]-U[,2])^2)
  Sig.est[r, 3] <- fit2$Sig
  V.est[r,,3] <- fit2$R[c(1,2,4)]
  
  # 'robustlmm' package
  rfit <- fit
  tt <- proc.time()
  try( rfit <- rlmer(Y~X[,-1]+(Z[,-1]|ID)) )
  CPT[r, 4] <- ( proc.time() - tt )[3]

  qq <- 2
  hBeta <- summary(rfit)$coefficients[,1]
  v <- as.matrix( coef(rfit)$ID[,1:qq]-hBeta[qq:1] )
  v.mat <- as.matrix( cbind(v[,qq], v[,-qq]) )
  hRE <- as.vector(v.mat)
  hSig <- summary(rfit)$sigma
  hV <- as.matrix( summary(rfit)$varcor$ID )
  Beta.est[r,,4] <- hBeta
  hRE <- matrix(hRE, m, 2)
  MSE.RE1[r, 4] <- mean((hRE[,1]-U[,1])^2)
  MSE.RE2[r, 4] <- mean((hRE[,2]-U[,2])^2)
  Sig.est[r, 4] <- hSig
  V.est[r,,4] <- hV[c(1,2,4)]
  
  # marginal DPD (density power divergence)
  tt <- proc.time()
  fit2 <- LMM.DPD(Y, X, Z, ID, gam=0.5, maxitr=100)
  CPT[r, 5] <- ( proc.time() - tt )[3]
  
  Beta.est[r,,5] <- fit2$Beta
  hRE <- matrix(fit2$RE, m, 2)
  MSE.RE1[r, 5] <- mean((hRE[,1]-U[,1])^2)
  MSE.RE2[r, 5] <- mean((hRE[,2]-U[,2])^2)
  Sig.est[r, 5] <- fit2$Sig
  V.est[r,,5] <- fit2$R[c(1,2,4)]
  
  # heavy-tailed distribution 
  tt <- proc.time()
  hfit <- heavyLme(Y~x1+x2+x3, random=~1+x2, groups=~ID)
  CPT[r, 6] <- ( proc.time() - tt )[3]
  
  hBeta <- coef(hfit)
  Beta.est[r,,6] <- hBeta
  hRE <- hfit$ranef
  MSE.RE1[r, 6] <- mean((hRE[,1]-U[,1])^2)
  MSE.RE2[r, 6] <- mean((hRE[,2]-U[,2])^2)
  
  # print
  print(r)
}



## save 
save(list=ls(), file=paste0("sim-res", s, "-m", m, ".RData"))
