###------------------------------------------------------###
###            R code for application to                 ###
###          AIDS Cohort Study in Section 6              ###
###------------------------------------------------------###
rm(list=ls())

## Packages
library(catdata)
library(tidyverse)
library(lme4)
library(robustlmm)
library(colorspace)
source("RLMM_HGD.R")


## Dataset
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
select <- RLMM_HGD_select(Y, X, Z, ID, gam.set=gam.set, print=T)
gam.opt <- select$gam

fit_HGD <- RLMM_HGD(Y, X, Z, ID, gam=gam.opt)

set.seed(1)
boot_HGD <- RLMM_HGD_boot(Y, X, Z, ID, gam=gam.opt, B=100)

# Robust LMM using 'robustlmm' package
fit_rml <- rlmer(Y~X[,-1]+(Z[,-1]|ID))



# HGD
CI_HGD <- t(apply(boot_HGD$Boot.Beta, 2, quantile, prob=c(0.025, 0.975)))
Beta_HGD <- fit_HGD$Beta
U_HGD <- cbind(fit_HGD$RE[1:m], fit_HGD$RE[-(1:m)])


# ML
Sig_ml <- summary(fit_ML)$sigma
Beta_ml <- summary(fit_ML)$coefficients[,1]
Pred <- predict(fit_ML)
Resid <- (Y-Pred)/Sig_ml
R_ml <- summary(fit_ML)$varcor[[1]]
U_ml <- t(matrix(fit_ML@u, 2, m))

zz <- qnorm(0.975)
SE_ml <- summary(fit_ML)$coefficients[,2]
CI_ml <- cbind(Beta_ml-zz*SE_ml, Beta_ml+zz*SE_ml)


# RML
Beta_rml <- summary(fit_rml)$coefficients[,1]
SE_rml <- summary(fit_rml)$coefficients[,2]
CI_rml <- cbind(Beta_rml-zz*SE_rml, Beta_rml+zz*SE_rml)
Sig_rml <- summary(fit_rml)$sigma
R_rml <- summary(fit_rml)$varcor$ID
U_rml <- t(matrix(fit_rml@b.s, 2, m))


##-----------------
## Table 2 (Regression coefficients)
##-----------------
CI_ml <- round(CI_ml, 2)
CI_ml <- paste0("(",CI_ml[,1],", ",CI_ml[,2],")")
CI_rml <- round(CI_rml, 2)
CI_rml <- paste0("(",CI_rml[,1],", ",CI_rml[,2],")")
CI_HGD <- round(CI_HGD, 2)
CI_HGD <- paste0("(",CI_HGD[,1],", ",CI_HGD[,2],")")
result <- cbind(round(Beta_ml,2), CI_ml, round(Beta_rml,2), CI_rml, round(Beta_HGD,2), CI_HGD)
result <- result[c(1,3,4,5,2,8,9,6,12,13,7,10,11),]

write.csv(result, file="Table-2.csv")



##-----------------
## Table 2 (Variance parameters)
##-----------------
mat <- rbind(c(Sig_ml^2, as.vector(R_ml)[c(1,2,4)]),
             c(Sig_rml^2, as.vector(R_rml)[c(1,2,4)]),
             c(fit_HGD$Sig^2, as.vector(fit_HGD$R)[c(1,2,4)]))

write.csv(mat, file="Table-3.csv")


##-----------------
## Figure 3 (QQ-plot)
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

ggsave("Figure-3.pdf", height=4, width=10)



##-----------------
## Figure 4 (non-linear effects)
##-----------------
ran <- range(X[,2])
xx <- seq(ran[1], ran[2], length=200)
XX <- cbind(xx, xx^2, xx^3)
all_results <- data.frame(variable="Time",
                          covariate=xx,
                          HGD=XX%*%Beta_HGD[c(2,8,9)], 
                          ML=XX%*%Beta_ml[c(2,8,9)], 
                          RML=XX%*%Beta_rml[c(2,8,9)]) %>% 
  pivot_longer(HGD:RML, names_to="method")

ran <- range(X[,7])
xx <- seq(ran[1], ran[2], length=200)
XX <- cbind(xx, xx^2, xx^3)
all_results <- bind_rows(
  all_results,
  data.frame(variable="Age",
             covariate=xx, 
             HGD=XX%*%Beta_HGD[c(7,10,11)], 
             ML=XX%*%Beta_ml[c(7,10,11)], 
             RML=XX%*%Beta_rml[c(7,10,11)]) %>% 
    pivot_longer(HGD:RML, names_to="method"))

ran <- range(X[,6])
xx <- seq(ran[1], ran[2], length=200)
XX <- cbind(xx, xx^2, xx^3)
all_results <- bind_rows(
  all_results,
  data.frame(variable="Cesd",
             covariate=xx, 
             HGD=XX%*%Beta_HGD[c(6,12,13)], 
             ML=XX%*%Beta_ml[c(6,12,13)], 
             RML=XX%*%Beta_rml[c(6,12,13)]) %>% 
    pivot_longer(HGD:RML, names_to="method"))

ggplot(all_results %>% 
         mutate(variable=fct_inorder(variable)) %>% 
         mutate(method=fct_relevel(method, c("HGD", "ML", "RML"))), 
       aes(x=covariate, y=value, color=method, group=method)) +
  geom_line(size=1) +
  scale_color_discrete_qualitative(palette="Set 2", nmax=6, order=c(1,3,4)) +
  facet_wrap(. ~ variable, nrow=1, scales="free") +
  labs(x="Covariate", y="Estimated conditional effect", group="Method", color="Method", shape="Method") +
  theme_bw() +
  theme(legend.position="bottom") +
  guides(colour=guide_legend(nrow=1))

ggsave("Figure-4.pdf", width=10, height=5)



##-----------------
## Figure 5 (predicted random effects)
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

ggsave("Figure-5.pdf", width=8, height=5)




