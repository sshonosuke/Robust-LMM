###------------------------------------------------------###
###             R code for application to                ###
###          AIDS Cohort Study in Section S6             ###
###           (Transformed response values)              ###
###------------------------------------------------------###
rm(list=ls())

## Packages
library(catdata)
library(tidyverse)
library(lme4)
library(robustlmm)
library(colorspace)
source("RLMM_HGD.R")


## transformation functions
trans <- "log"   # Three options: "log", "square_root", "cube_root"

## load dataset
data(aids)
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

## transformation of response values
Y_origin <- aids$cd4/100
if(trans=="log"){ Y <- log(Y_origin) }
if(trans=="square_root"){ Y <- sqrt(Y_origin) }
if(trans=="cube_root"){ Y <- (Y_origin)^(1/3) }




## Estimation 
# ML (non-robust method)
fit_ML <- lmer(Y~X[,-1]+(Z[,-1]|ID), REML=F)

Sig_ml <- summary(fit_ML)$sigma
Beta_ml <- summary(fit_ML)$coefficients[,1]
Pred <- predict(fit_ML)
Resid <- (Y-Pred)/Sig_ml
R_ml <- summary(fit_ML)$varcor[[1]]
U_ml <- t(matrix(fit_ML@u, 2, m))

# HGD (proposed)
gam.set <- seq(0, 0.5, by=0.02)
select <- RLMM_HGD_select(Y, X, Z, ID, gam.set=gam.set, print=T)
(gam.opt <- select$gam)

fit_HGD <- RLMM_HGD(Y, X, Z, ID, gam=gam.opt)

# Robust LMM using 'robustlmm' package
fit_rml <- rlmer(Y~X[,-1]+(Z[,-1]|ID))


## Random effects
U_rml <- t(matrix(fit_rml@b.s, 2, m))
U_HGD <- cbind(fit_HGD$RE[1:m], fit_HGD$RE[-(1:m)])




##-----------------
## Figure S4 (QQ-plot)
##-----------------
all_results <- bind_rows(
  data.frame(criterion="Residual", value=Resid),
  data.frame(criterion="Random intercept", value=U_ml[,1]/sd(U_ml[,1])),
  data.frame(criterion="Random slope", value=U_ml[,2]/sd(U_ml[,2]))) %>% 
  mutate(criterion=fct_inorder(criterion))

ggplot(all_results, aes(sample=value)) +
  geom_qq() +
  geom_abline(intercept=0, slope=1, linetype=2) +
  facet_wrap(. ~ criterion, scales="free", nrow=1) +
  theme_bw()

ggsave(paste0("Figure-S4-", trans,".pdf"), height=4, width=10)



##-----------------
## Figure S5 (Scatter plots of predicted random effects)
##-----------------
all_results <- bind_rows(
  data.frame(ML=U_ml[,1], HGD=U_HGD[,1], RML=U_rml[,1], effect="Intercept") %>% pivot_longer(HGD:RML, names_to="method"),
  data.frame(ML=U_ml[,2], HGD=U_HGD[,2], RML=U_rml[,2], effect="Slope") %>% pivot_longer(HGD:RML, names_to="method"))

ggplot(all_results, aes(x=ML, y=value, group=method, color=method, shape=method)) + 
  geom_point(size=1) +
  scale_color_discrete_qualitative(palette="Set 2", nmax=6, order=c(1,4)) +
  geom_abline(intercept=0, slope=1, linetype=2) +
  facet_wrap(. ~ effect, nrow=1, scales="free") +
  labs(x="ML", y="Robust method", group="Method", color="Method", shape="Method") +
  theme_bw() +
  theme(legend.position="bottom") +
  guides(colour=guide_legend(nrow=1))

ggsave(paste0("Figure-S5-", trans,".pdf"), height=4, width=10)



