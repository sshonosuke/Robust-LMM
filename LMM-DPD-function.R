###----------------------------------------------------------###
###       This script includes LMM.DPD for fitting           ###
###       robust linear mixed models via marginal            ###
###              density power divergence                    ###        
###----------------------------------------------------------###

## INPUT 
# Y: response variable (vector)
# X: covariate matrix for fixed effects
# Z: covariate matrix for random effects
# ID: group labels (vector)
# gam: tuning parameter of density power divergence (from 0 to 1)
# maxitr (optional): maximum number of iterations for MM-algorithm
# init.fit (optional): initial fit of the linear mixed models

## OUTPUT
# Beta: coefficient of fixed effects (vector)
# Sig: standard deviation of the error term 
# R: covariance matrix of random effects
# RE: random effects
# Mu: mean vector of response variable
# weight: density power weight for each observation 
# itr: number of iterations at convergence

LMM.DPD <- function(Y, X, Z, ID, gam=0.25, maxitr=500, init=NULL){
  ## preparation 
  qq <- dim(Z)[2]    # dimension of random effects
  m <- length(unique(ID))
  ni <- table(ID)
  ID <- rep(1:m, ni)
  N <- sum(ni)
  tr <- function(mat){ sum(diag(mat)) }
  
  ##  design matrix for random effects
  ZZ <- matrix(0, N, m*qq)  
  for(i in 1:m){ 
    for(k in 1:qq){
      ZZ[ID==i, m*(k-1)+i] <- Z[ID==i, k]
    }
  }
  ZZ <- as(ZZ, "sparseMatrix")
  
  ## initial fit (standard maximum likelihood)
  if(is.null(NULL)){
    fit <- lmer(Y~X[,-1]+(Z[,-1]|ID), REML=F)
    Beta <-  summary(fit)$coefficients[,1]
    V <- as.matrix( coef(fit)$ID[,1:qq]-Beta[qq:1] )
    V.mat <- as.matrix( cbind(V[,qq], V[,-qq]) )
    V <- as.vector(V.mat)
    Sig <- summary(fit)$sigma
    R <- as.matrix( summary(fit)$varcor$ID )[1:qq, 1:qq] + 0.01*diag(qq)
    init.fit <- list(Beta=Beta, RE=V, V=R, Sig=Sig)
  }else{
    Beta <- init.fit$Beta
    V <- init.fit$RE
    R <- init.fit$V
    Sig <- init.fit$Sig
  }
  
  ## objective function 
  Q <- function(beta, sig, R){
    val <- c()
    mu <- as.vector( X%*%beta )
    for(i in 1:m){
      Z_i <- Z[ID==i,]
      Sigma_i <- Z_i%*%R%*%t(Z_i) + sig^2*diag(ni[i])
      dens <- dmvnorm(Y[ID==i], mean=Mu[ID==i], sigma=Sigma_i)^gam
      val[i] <- det(Sigma_i)^(-gam/2)*(2*pi)^(-0.5*gam*ni[i])*(1+gam)^(-0.5*ni[i]) - (1+1/gam)*dens
    }
    return(mean(val))
  }
  
  # iteration 
  for(k in 1:maxitr){
    Beta0 <- Beta
    
    # weight 
    Mu <- as.vector( X%*%Beta )
    ww <- c()
    Sigma_i <- list()
    IS_i <- list()
    for(i in 1:m){
      Z_i <- Z[ID==i,]
      Sigma_i[[i]] <- Z_i%*%R%*%t(Z_i) + Sig^2*diag(ni[i])
      IS_i[[i]] <- solve(Sigma_i[[i]])
      ww[i] <- dmvnorm(Y[ID==i], mean=Mu[ID==i], sigma=Sigma_i[[i]])^gam
    }
    
    # update beta
    denom <- 0
    num <- 0
    for(i in 1:m){
      X_i <- X[ID==i,]
      denom <- denom + ww[i]*t(X_i)%*%IS_i[[i]]%*%X_i
      num <- num + ww[i]*t(X_i)%*%IS_i[[i]]%*%Y[ID==i]
    }
    Beta <- as.vector( solve(denom)%*%num )
    
    # update sigma 
    Q_sigma <- function(x){ Q(beta=Beta, sig=x, R=R) }
    Sig <- optim(fn=Q_sigma, par=Sig, lower=0.0001, upper=10, method="L-BFGS-B")$par
    
    # update R (only when qq=2)
    Q_R <- function(x){ 
      v1 <- x[1]
      v2 <- x[2]
      rho <- x[3]
      Rx <- matrix(c(v1, rep(rho*sqrt(v1*v2), 2), v2), 2, 2)
      return( Q(beta=Beta, sig=Sig, R=Rx) )
    }
    R_par <- optim(fn=Q_R, par=c(R[1], R[4], R[2]/sqrt(R[1]*R[4])), lower=c(0.01, 0.01, -0.99), 
                   upper=c(10, 10, 0.99), method="L-BFGS-B")$par
    R <- matrix(c(R_par[1],  rep(R_par[3]*sqrt(R_par[1]*R_par[2]), 2), R_par[2]), 2, 2)
    
    # difference
    dd <- sum(abs(Beta - Beta0)) / sum(abs(Beta0)+10^(-5))
    if( dd<10^(-6) ){ break }
  }
  
  # random effects
  RE <- matrix(NA, m, qq)
  resid <- Y - Mu
  for(i in 1:m){
    Z_i <- Z[ID==i,]
    Sigma_i <- Z_i%*%R%*%t(Z_i) + Sig^2*diag(ni[i])
    IS_i <- solve(Sigma_i)
    RE[i,] <- as.vector( R%*%t(Z_i)%*%IS_i%*%resid[ID==i] )
  }
  Mu <- as.vector(X%*%Beta+ZZ%*%V)
  
  ## Result
  Res <- list(Beta=Beta, Sig=Sig, R=R, RE=RE, Mu=Mu, weight=ww, itr=k)
  return(Res)
}



