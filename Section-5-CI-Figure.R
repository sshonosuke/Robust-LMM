###----------------------------------------------------------###
###           R code for reproducing Figure 2                ###
###              in Sugasawa et al. (2024)                   ###
###----------------------------------------------------------###
## Remark: This code requires multiple outputs (.RData) 
##           obtained by running "Section-5-CI.R".

rm(list=ls())
library(tidyverse)
library(colorspace)
library(ggplot2)

## Summary of results 
s.set <- c(1, 4, 7, 2, 5, 8, 3, 6, 9)
meth <- c("aHGD", "HGD", "ML", "RML")
L <- length(meth)

IS1 <- IS2 <- CP1 <- CP2 <- array(NA, c(9, 4, L))
dimnames(CP1)[[3]] <- dimnames(CP2)[[3]] <- meth
dimnames(IS1)[[3]] <- dimnames(IS2)[[3]] <- meth

for(ss in 1:9){
  load(paste0("sim-CI-result(s=",s.set[ss],",m=50).RData"))
  IS1[ss,,] <- apply(IS, c(2,3), mean)
  CP1[ss,,] <- apply(CP, c(2,3), mean)
  
  load(paste0("sim-CI-result(s=",s.set[ss],",m=100).RData"))
  IS2[ss,,] <- apply(IS, c(2,3), mean)
  CP2[ss,,] <- apply(CP, c(2,3), mean)
}


# Objects for ggplot 
all_results1 <- bind_rows(
  cbind(CP1[,1,] %>% as.data.frame, scenario=paste0("S", 1:9), criterion="Coverage probability", m="m=50"),
  cbind(CP2[,1,] %>% as.data.frame, scenario=paste0("S", 1:9), criterion="Coverage probability", m="m=100"),
  cbind(IS1[,1,] %>% as.data.frame, scenario=paste0("S", 1:9), criterion="Interval score", m="m=50"),
  cbind(IS2[,1,] %>% as.data.frame, scenario=paste0("S", 1:9), criterion="Interval score", m="m=100")) %>% 
  pivot_longer(aHGD:RML, names_to="method") %>% 
  mutate(criterion=fct_inorder(criterion)) %>% 
  mutate(m=fct_inorder(m)) %>% 
  mutate(method=fct_relevel(method, c("aHGD", "HGD", "ML", "RML"))) 

all_results2 <- bind_rows(
  cbind(CP1[,2,] %>% as.data.frame, scenario=paste0("S", 1:9), criterion="Coverage probability", m="m=50"),
  cbind(CP2[,2,] %>% as.data.frame, scenario=paste0("S", 1:9), criterion="Coverage probability", m="m=100"),
  cbind(IS1[,2,] %>% as.data.frame, scenario=paste0("S", 1:9), criterion="Interval score", m="m=50"),
  cbind(IS2[,2,] %>% as.data.frame, scenario=paste0("S", 1:9), criterion="Interval score", m="m=100")) %>% 
  pivot_longer(aHGD:RML, names_to="method") %>% 
  mutate(criterion=fct_inorder(criterion)) %>% 
  mutate(m=fct_inorder(m)) %>% 
  mutate(method=fct_relevel(method, c("aHGD", "HGD", "ML", "RML"))) 

all_results3 <- bind_rows(
  cbind(CP1[,3,] %>% as.data.frame, scenario=paste0("S", 1:9), criterion="Coverage probability", m="m=50"),
  cbind(CP2[,3,] %>% as.data.frame, scenario=paste0("S", 1:9), criterion="Coverage probability", m="m=100"),
  cbind(IS1[,3,] %>% as.data.frame, scenario=paste0("S", 1:9), criterion="Interval score", m="m=50"),
  cbind(IS2[,3,] %>% as.data.frame, scenario=paste0("S", 1:9), criterion="Interval score", m="m=100")) %>% 
  pivot_longer(aHGD:RML, names_to="method") %>% 
  mutate(criterion=fct_inorder(criterion)) %>% 
  mutate(m=fct_inorder(m)) %>% 
  mutate(method=fct_relevel(method, c("aHGD", "HGD", "ML", "RML"))) 

all_results4 <- bind_rows(
  cbind(CP1[,4,] %>% as.data.frame, scenario=paste0("S", 1:9), criterion="Coverage probability", m="m=50"),
  cbind(CP2[,4,] %>% as.data.frame, scenario=paste0("S", 1:9), criterion="Coverage probability", m="m=100"),
  cbind(IS1[,4,] %>% as.data.frame, scenario=paste0("S", 1:9), criterion="Interval score", m="m=50"),
  cbind(IS2[,4,] %>% as.data.frame, scenario=paste0("S", 1:9), criterion="Interval score", m="m=100")) %>% 
  pivot_longer(aHGD:RML, names_to="method") %>% 
  mutate(criterion=fct_inorder(criterion)) %>% 
  mutate(m=fct_inorder(m)) %>% 
  mutate(method=fct_relevel(method, c("aHGD", "HGD", "ML", "RML"))) 

all_results <- bind_rows(
  all_results1 %>% mutate(dataset="Beta 1"),
  all_results2 %>% mutate(dataset="Beta 2"),
  all_results3 %>% mutate(dataset="Beta 3"),
  all_results4 %>% mutate(dataset="Beta 4")
)



##-----------------
##  Figure 2: Coverage probability and interval score (4 by 4 panels)
##-----------------
ggplot(all_results, aes(x=scenario, y=value, group=method, color=method, shape=method)) +
  geom_point(size=2.5) +
  geom_line() +
  scale_color_discrete_qualitative(palette="Set 2", nmax=6) +
  scale_shape_manual(labels=levels(all_results$method), values=c(16,16,17,17)) +
  facet_wrap(vars(dataset, criterion, m), scales="free_y") + 
  labs(x="Scenario", y="Value", color="Method", shape="Method") +
  theme_bw() +
  theme(legend.position="bottom") +
  guides(colour=guide_legend(nrow=1))

ggsave("Figure-2.pdf", width=13, height=13)




