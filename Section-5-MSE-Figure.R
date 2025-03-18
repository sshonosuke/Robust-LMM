###----------------------------------------------------------###
###     R code for reproducing Figure 1 (and Figure S1)      ###
###              in Sugasawa et al. (2024)                   ###
###----------------------------------------------------------###
## Remark: This code requires multiple outputs (.RData) 
##           obtained by running "Section-5-MSE.R".

rm(list=ls())
library(tidyverse)
library(colorspace)
library(ggplot2)

R <- 500   # number of Monte Carlo replications

# true parameter
V <- matrix(c(1, 0.3, 0.3, 1), 2, 2)[c(1,2,4)]
Sig <- 1.5
Beta <- c(0.5, 0.3, 0.5, 0.8)

s.set <- c(1, 4, 7, 2, 5, 8, 3, 6, 9)
meth <- c("ML", "aHGD", "HGD", "RML", "mDPD", "HT")
M <- length(meth)

MSE <- array(NA, c(9, 6, M))
sMSE <- array(NA, c(9, 3, M))
RMSE <- matrix(NA, 9, M)
Lab <- c(paste0("Beta",0:3), "V", "Sigma")
dimnames(MSE)[[1]] <- dimnames(RMSE)[[1]] <- dimnames(sMSE)[[1]] <- paste0("S",1:9)
dimnames(MSE)[[2]] <- Lab
dimnames(MSE)[[3]] <- dimnames(RMSE)[[2]] <- dimnames(sMSE)[[3]] <- meth
gam.est <- c()


# m=100
m <- 100  
for(ss in 1:9){
  load(paste0("sim-res", s.set[ss],"-m", m,".RData"))
  sub <- is.na(apply(MSE.RE1, 1, sum))==F
  RMSE[ss,] <-  apply(MSE.RE1[sub,], 2, mean) + apply(MSE.RE2[sub,], 2, mean) 
  for(j in 1:4){
    MSE[ss,j,] <- apply((Beta.est[sub,j,]-Beta[j])^2, 2, mean)
  }
  sMSE[ss,1,] <- apply(MSE[ss,1:4,], 2, sum)
  MSE[ss,6,] <- sMSE[ss,2,] <- apply((Sig.est[sub,]-Sig)^2, 2, mean)
  for(l in 1:M){
    MSE[ss, 5, l] <- mean( apply((t(V.est[sub,,l])-V[c(1,2,4)])^2, 2, sum) )
  }
  sMSE[ss,3,] <- MSE[ss,5,]
  gam.est[ss] <- mean(Gam)
}

sMSE_100 <- sMSE
RMSE_100 <- RMSE
write.csv(t(gam.est), file=paste0("gam(m=",m,").csv"))

# m=50
m <- 50 
for(ss in 1:9){
  load(paste0("sim-res", s.set[ss],"-m", m,".RData"))
  sub <- is.na(apply(MSE.RE1, 1, sum))==F
  RMSE[ss,] <-  apply(MSE.RE1[sub,], 2, mean) + apply(MSE.RE2[sub,], 2, mean) 
  for(j in 1:4){
    MSE[ss,j,] <- apply((Beta.est[sub,j,]-Beta[j])^2, 2, mean)
  }
  sMSE[ss,1,] <- apply(MSE[ss,1:4,], 2, sum)
  MSE[ss,6,] <- sMSE[ss,2,] <- apply((Sig.est[sub,]-Sig)^2, 2, mean)
  for(l in 1:M){
    MSE[ss, 5, l] <- mean( apply((t(V.est[sub,,l])-V[c(1,2,4)])^2, 2, sum) )
  }
  sMSE[ss,3,] <- MSE[ss,5,]
  gam.est[ss] <- mean(Gam)
}

sMSE_50 <- sMSE
RMSE_50 <- RMSE
write.csv(t(gam.est), file=paste0("gam(m=",m,").csv"))



# Objects for ggplot 
dimnames(sMSE_50)[[2]] <- dimnames(sMSE_100)[[2]] <- c("Beta", "Sigma", "R")

all_results <- bind_rows(
  cbind(as.data.frame.table(sMSE_50), m="m=50"),
  cbind(as.data.frame.table(sMSE_100), m="m=100"))
colnames(all_results) <- c("scenario", "parameter", "method", "value", "m")

all_ranef_results <- bind_rows(
  data.frame(as.data.frame.table(RMSE_50), m="m=50", parameter="random_effects"),
  data.frame(as.data.frame.table(RMSE_100), m="m=100", parameter="random_effects"))
colnames(all_ranef_results) <- c("scenario", "method", "value", "m", "parameter")
all_ranef_results <- all_ranef_results %>% 
  relocate(colnames(all_results))


all_results <- bind_rows(all_results, all_ranef_results) 
rm(all_ranef_results)
rm(list=ls(pattern="RMSE")) 
rm(list=ls(pattern="sMSE")) 

all_results <- all_results %>% 
  mutate(m=fct_inorder(m)) %>%
  mutate(parameter=fct_inorder(parameter)) %>%
  mutate(parameter=fct_recode(parameter,
                                "Regression coefficients"="Beta",
                                "Random effects covariance"="R",
                                "Random effects"="random_effects",
                                "Error variance"="Sigma")) %>% 
  mutate(isHGD=factor(method %in% c("HGD", "aHGD"), labels=c("No", "Yes"))) %>% 
  mutate(method=fct_relevel(method, c("aHGD", "HGD", "ML", "RML", "mDPD", "HT")))


## Settings of plot
meth <- c("ML", "aHGD", "HGD", "RML", "mDPD", "HT")
M <- length(meth)
Ind <- c(1, 4, 5, 6, 2, 3)
col <- c(1, 1, 3, 4, 2, 2)
lty <- c(1, 2, 1, 1, 1, 2)
lwd <- c(1, 1, 1, 1, 1, 1)
scale <- 1.4


##-----------------
##  Figure 1 (Beta, random effect, random effect variance)
##-----------------

ggplot(data=all_results %>% dplyr::filter(parameter != "Error variance"), aes(x=scenario, y=value, group=method, color=method, shape=method)) +
  geom_point(size=2.5) +
  geom_line() +
  scale_color_discrete_qualitative(palette="Set 2") +
  scale_shape_manual(labels=levels(all_results$method), values=c(16,16,17,17,17,17)) +
  scale_y_log10() +
  facet_wrap(parameter ~ m, nrow=3, scales="free_y") +
  labs(x="Scenario", y="MSE", group="Method", color="Method", shape="Method") +
  theme_bw() +
  theme(legend.position="bottom") +
  guides(colour=guide_legend(nrow=2))

ggsave("Figure-1.pdf", width=8, height=10)



##-----------------
##  Figure S1 (error variance)
##-----------------
ggplot(data=all_results %>% dplyr::filter(parameter == "Error variance"), aes(x=scenario, y=value, group=method, color=method, shape=method)) +
  geom_point(size=2.5) +
  geom_line() +
  scale_color_discrete_qualitative(palette="Set 2") +
  scale_shape_manual(labels=levels(all_results$method), values=c(16,16,17,17,17,17)) +
  scale_y_log10() +
  facet_wrap(. ~ m, nrow=1, scales="free_y") +
  labs(x="Scenario", y="MSE (error variance)", group="Method", color="Method", shape="Method") +
  theme_bw() +
  theme(legend.position="bottom") +
  guides(colour=guide_legend(nrow=1))

ggsave("Figure-S1.pdf", width=8, height=3.8)
