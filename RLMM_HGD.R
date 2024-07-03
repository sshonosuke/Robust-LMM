###----------------------------------------------------------###
###   This script includes the following three functions:    ###
###    - RLMM_HGD: Fitting robust linear mixed models        ###
###    - RLMM_HGD_select: Tuning parameter selection         ###
###    - RLMM_HGD_boot: Bootstrap for inference              ###
###----------------------------------------------------------###

# packages 
library(mvtnorm)
library(lme4)
library(SparseM)
library(MCMCpack)


###----------------------------------------------------------###
###              Robust linear mixed models                  ###
###          via hierarchical gamma-divergence               ###
###----------------------------------------------------------###
## INPUT 
# Y: response variable (vector)
# X: covariate matrix for fixed effects
# Z: covariate matrix for random effects
# ID: group labels (vector)
# gam: tuning parameter of gamma-divergence (from 0 to 1)
# bw (optional): weight for observations 
# maxitr (optional): maximum number of iterations for MM-algorithm
# init.fit (optional): initial fit of the linear mixed models

## OUTPUT
# Beta: coefficient of fixed effects (vector)
# Sig: standard deviation of the error term 
# R: covariance matrix of random effects
# RE: random effects
# Mu: mean vector of response variable
# itr: number of iterations at convergence

RLMM_HGD <- function(Y, X, Z, ID, gam=0.5, bw=NULL, maxitr=500, init.fit=NULL){
  ## preparation 
  qq <- dim(Z)[2]    # dimension of random effects
  m <- length(unique(ID))
  ni <- table(ID)
  ID <- rep(1:m, ni)
  N <- sum(ni)
  tr <- function(mat){ sum(diag(mat)) }
  
  ## initial fit (standard maximum likelihood)
  if(is.null(NULL)){
    fit <- lmer(Y~X[,-1]+(Z[,-1]|ID), REML=F)
    Beta <-  summary(fit)$coefficients[,1]
    V <- as.matrix( coef(fit)$ID[,1:qq]-Beta[qq:1] )
    V.mat <- as.matrix( cbind(V[,qq], V[,-qq]) )
    V <- as.vector(V.mat)
    Sig <- summary(fit)$sigma
    R <- as.matrix( summary(fit)$varcor$ID )[1:qq, 1:qq] 
    R <- R + 0.01*diag(qq)     # regularization to avoid numerical error
    init.fit <- list(Beta=Beta, RE=V, V=R, Sig=Sig)
  }else{
    Beta <- init.fit$Beta
    V <- init.fit$RE
    R <- init.fit$V
    Sig <- init.fit$Sig
  }
  
  ##  design matrix for random effects
  ZZ <- matrix(0, N, m*qq)  
  for(i in 1:m){ 
    for(k in 1:qq){
      ZZ[ID==i, m*(k-1)+i] <- Z[ID==i, k]
    }
  }
  ZZ <- as(ZZ, "sparseMatrix")
  
  ##  weight (used for cluster bootsrap)
  if(is.null(bw)){ bw <- rep(1, m) }
  bbw <- bw[ID]
  
  ##   MM algorithm   ##
  for(k in 1:maxitr){
    Beta0 <- Beta
    V0 <- V
    
    # weight
    mu <- as.vector(X%*%Beta+ZZ%*%V)
    val <- dnorm(Y, mu, Sig)^(gam)
    ww <- val/mean(val)*bbw
    
    # update beta
    resid <- as.vector(Y-ZZ%*%V)
    Beta <- as.vector( solve(t(ww*X)%*%X)%*%t(X)%*%(ww*resid) )
    
    # update sigma 
    resid <- as.vector(Y-X%*%Beta-ZZ%*%V)
    IS <- list()
    for(i in 1:m){
      sZ <- Z[ID==i,]
      if( is.null(dim(sZ)) ){ sZ <- matrix(sZ, 1, qq) }
      IS[[i]] <- solve( Sig^2*diag(ni[i]) + sZ%*%R%*%t(sZ) )
    }
    N2 <- Sig^2*sum(unlist(lapply(IS, tr))) - N*gam/(1+gam)
    Sig <- sqrt( sum(ww*resid^2)/N2 )
    
    # update random effects
    V.mat <- matrix(V, m, qq)
    log.uu <- gam*dmvnorm(V.mat, rep(0,qq), R, log=T)
    log.uu <- log.uu - max(log.uu)
    val <- exp(log.uu)
    uu <- val/mean(val)*bw
    resid <- as.vector(Y-X%*%Beta)
    V.mat <- matrix(NA, m, qq)
    for(i in 1:m){
      sZ <- Z[ID==i,]
      if( is.null(dim(sZ)) ){ sZ <- matrix(sZ, 1, qq) }
      V.mat[i,] <- as.vector( solve(t(ww[ID==i]*sZ)%*%sZ + uu[i]*Sig^2*solve(R))%*%t(sZ)%*%(ww[ID==i]*resid[ID==i]) )
    }
    V <- as.vector(V.mat)
    
    # update R
    mat <- 0
    for(i in 1:m){
      sZ <- Z[ID==i,]
      if( is.null(dim(sZ)) ){ sZ <- matrix(sZ, 1, qq) }
      IS <- solve( Sig^2*diag(ni[i]) + sZ%*%R%*%t(sZ) )
      mat <- mat + t(sZ)%*%IS%*%sZ
    }
    R <- (1+gam)*(t(uu*V.mat)%*%V.mat - R%*%mat%*%R + m*R)/m
    R <- R + 0.01*diag(qq)     # regularization to avoid numerical error
  
    # checking convergence
    dd <- sum(abs(Beta - Beta0)) / sum(abs(Beta0)+0.0001)
    if( dd<10^(-6) ){ break }
  }
  Mu <- as.vector(X%*%Beta+ZZ%*%V)
  
  ## Result
  Res <- list(Beta=Beta, Sig=Sig, R=R, RE=V, Mu=Mu, itr=k)
  return(Res)
}





###----------------------------------------------------------###
###        Tuning parameter selection by H-score             ###
###----------------------------------------------------------###
## INPUT
# Y: response variable (vector)
# X: covariate matrix for fixed effects
# Z: covariate matrix for random effects
# ID: group labels (vector)
# gam.set: set of candidate values for gamma
# maxitr (optional): maximum number of iterations for MM-algorithm
# print: If 'T', H-scores under each candidate value is shown

## OUTPUT
# gam: optimal values of gamma 
# HS1: set of H-scores for conditional distributions of observations
# HS2: set of H-scores for marginal distribuiton of random effects 

RLMM_HGD_select <- function(Y, X, Z, ID, gam.set=NULL, maxitr=50, print=F){
  ## preparation
  qq <- dim(Z)[2]
  if(is.null(gam.set)){
    gam.set <- seq(0, 0.5, by=0.05)
  }
  L <- length(gam.set)
  
  ## H-score
  HS1 <- c()
  HS2 <- c()
  for(j in 1:L){
    gam <- gam.set[j]
    fit <- RLMM_HGD(Y=Y, X=X, Z=Z, ID=ID, gam=gam, maxitr=maxitr)
    Beta <- fit$Beta
    Sig <- fit$Sig
    V <- fit$R
    RE <- matrix(fit$RE, m, qq)
    Mu <- fit$Mu
    # H1
    ww <- dnorm(Y, Mu, Sig)^(gam)
    C <- ( (2*pi)^(-gam/2)*Sig^(-gam)*(1+gam)^(-1/2) )^(-gam/(1+gam))
    res <- Y - Mu
    T1 <- C*ww*gam*res^2/Sig^4 - C*ww/Sig^2
    T2 <- C^2*ww^2*res^2/Sig^4
    HS1[j] <- sum(2*T1) + sum(T2)
    # H2
    Om <- solve(V)
    dens <- (2*pi)^(-qq/2)*det(Om)^(1/2)*exp(-0.5*diag(RE%*%Om%*%t(RE)))
    ww <- dens^(gam)
    C <- ( (2*pi)^(-qq*gam/2)*det(Om)^(gam/2)*(1+gam)^(-qq/2) )^(-gam/(1+gam))
    T1 <- C*ww*gam*apply((Om%*%t(RE))^2, 2, sum) - C*ww*sum(diag(Om)) 
    T2 <- C^2*ww^2*apply((Om%*%t(RE))^2, 2, sum) 
    HS2[j] <- sum(2*T1) + sum(T2)
    if(print){ print(c(HS1[j], HS2[j])) }
  }
  
  ## Result 
  gam.opt1 <- gam.set[which.min(HS1)]
  gam.opt2 <- gam.set[which.min(HS2)]
  gam.opt <- max(gam.opt1, gam.opt2)
  Res <- list(gam=gam.opt, HS1=HS1, HS2=HS2)
}





###----------------------------------------------------------###
###       Cluster weighted bootstrap for inference           ###
###----------------------------------------------------------###
## INPUT
# Y: response variable (vector)
# X: covariate matrix for fixed effects
# Z: covariate matrix for random effects
# ID: group labels (vector)
# B: number of bootstrap 
# gam: tuning parameter of gamma-divergence (from 0 to 1)

## OUTPUT
# Beta: bootstrap samples of coefficient of fixed effects
# Sig: bootstrap samples of standard deviation of the error term 
# R: bootstrap samples of covariance matrix of random effects

RLMM_HGD_boot <- function(Y, X, Z, ID, B=200, gam=0.5, print=T){
  m <- length(unique(ID))
  p <- dim(X)[2]
  qq <- dim(Z)[2]
  Boot.Beta <- matrix(NA, B, p)
  Boot.R <- array(NA, c(B, qq, qq))
  Boot.Sig <- c()
  for(b in 1:B){
    BW <- as.vector(rdirichlet(1,rep(1,m)))*m   # Dirichlet weight
    bfit <- RLMM_HGD(Y, X, Z, ID, gam=gam, bw=BW)
    Boot.Beta[b,] <- bfit$Beta
    Boot.Sig[b] <- bfit$Sig
    Boot.R[b,,] <- bfit$R
    if(round(10*b/B)==(10*b/B)){ 
      print( paste0("Bootstrap: ", round(100*b/B), "% done") ) 
    }
  }
  Res <- list(Boot.Beta=Boot.Beta, Boot.R=Boot.R, Boot.Sig=Boot.Sig)
  return(Res)
}



