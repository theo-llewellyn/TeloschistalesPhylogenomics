#compare MCMCTree estimates to Reltime estimates

setwd("/Users/tbl19/The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/PHYLOGENOMICS/TIME_CALIBRATION/")

library(ggplot2)
library(dplyr)
library(magrittr)
library(cowplot)
library(ape)
library(ggtree)
library(tibble)
library(tidyr)
library(colorspace)
library(MCMCtreeR)
library(caper)
library(treeio)
library(stringr)
library(phytools)
library(ggpubr)
library(ggrepel)

#read in trees with 95% confidence intervals
MCMCtree_phy <- readMCMCtree("MCMCTREE/LICHENS_MCMCTree_84T/FigTree.tre")
ggtree(MCMCtree_phy$apePhy) + geom_text(aes(label=node), hjust=-.3)

Reltime_OLS_phy <- read.mega("RELTIME/Telos85T_4parts_OLS_method.nexus")
OLS_data <- na.omit(get.data(Reltime_OLS_phy))
Reltime_OLS_topology <- read.tree("RELTIME/Telos85T_4parts_OLS_method.tre")
ggtree(Reltime_OLS_topology) + geom_text(aes(label=node), hjust=-.3)

Reltime_branchlength_phy <- read.mega("RELTIME/Telos85T_4parts_branchlength_method.nexus")
branchlength_data <- na.omit(get.data(Reltime_branchlength_phy))
Reltime_BL_topology <- read.tree("RELTIME/Telos85T_4parts_branchlength_method.tre")
ggtree(Reltime_BL_topology) + geom_text(aes(label=node), hjust=-.3)


df <- MCMCtree_phy$nodeAges*100
df <- cbind(df, MCMCTree_diff=NA,OLS_mean=NA,OLS_lower=NA,OLS_upper=NA,OLS_diff=NA,BL_mean=NA,BL_lower=NA,BL_upper=NA,BLS_diff=NA)

#match node numbers in Reltime and MCMCTree
matched_nodes_OLS <- as.data.frame(matchNodes(MCMCtree_phy$apePhy,Reltime_OLS_topology))
matched_nodes_BL <- as.data.frame(matchNodes(MCMCtree_phy$apePhy,Reltime_BL_topology))

for(i in 4:nrow(df)){
  #get difference between CI
  df[i,4] <- df[i,3] - df[i,2]
  #get node number for that clade in the RelTime trees
  nodenumber_OLS <- filter(matched_nodes_OLS, tr1 == rownames(df)[i])[[2]]
  nodenumber_bl <- filter(matched_nodes_BL, tr1 == rownames(df)[i])[[2]]
  #save mean OLS div time for that node
  df[i,5] <- as.numeric(filter(OLS_data, node == nodenumber_OLS)[3])
  #save OLS confidence intervals and difference
  CI_OLS <- filter(OLS_data, node == nodenumber_OLS)[4][[1]]
  df[i,6] <- CI_OLS[[1]][2]
  df[i,7] <- CI_OLS[[1]][1]
  df[i,8] <- df[i,7] - df[i,6]
  #save mean bl div time for that node
  df[i,9] <- as.numeric(filter(branchlength_data, node == nodenumber_bl)[2])
  #save OLS confidence intervals and difference
  CI_BL <- filter(branchlength_data, node == nodenumber_bl)[3][[1]]
  df[i,10] <- CI_BL[[1]][2]
  df[i,11] <- CI_BL[[1]][1]
  df[i,12] <- df[i,11] - df[i,10]
}
#add node numbers
df <- as.data.frame(cbind(df, MCMCTree_nodenumber = rownames(df)))
#convert to numeric
df2 <- as.data.frame(sapply(df, as.numeric ))

#convert to long format
df_long <- pivot_longer(na.omit(df2)[,c(4,8,12,13)], !MCMCTree_nodenumber)


### STATS
###

#MCMCTree vs RelTime BL
#correlation to see if mean estimates themselves
mean_corr_BL <- cor.test(df2$mean, df2$BL_mean, 
                         method = c("pearson"),
                         conf.level = 0.95)
#paired t-test to see if MCMCTree and ReltimeOLS mean estimates are different
mean_diff_BL <- t.test(df2$mean, df2$BL_mean, paired = TRUE, alternative = "two.sided")

#paired t-test to see if MCMCTree and ReltimeBL error ranges are different
t.test(df2$BLS_diff, df2$MCMCTree_diff, paired = TRUE, alternative = "greater")
#just correlation test to see if they correlate
cor.test(df2$MCMCTree_diff, df2$BLS_diff, 
         method = c("pearson"),
         conf.level = 0.95)

mean_BL_plot <- ggplot(df2, aes(x = mean, y = BL_mean)) +
  geom_point() +
  theme_bw() +
  geom_abline(slope=1, intercept = 0) +
  ylab("RelTime BL estimate (Million years)") +
  xlab("MCMCTree estimate (Million years)") +
  xlim(0,110) +
  ylim(0,110) +
  annotate("text", x = 5, y = 100,  hjust = 0, 
           label = paste("R=",round(mean_corr_BL$estimate,2),", p<0.0001",sep = "")) +
  annotate("text", x = 5, y = 95, hjust = 0, 
           label = paste("mean difference=",round(mean_diff_BL$estimate,2),", p<0.0001",sep = ""))


#MCMCTree vs RelTime OLS
#correlation to see if mean estimates themselves
mean_corr_OLS <- cor.test(df2$mean, df2$OLS_mean, 
                          method = c("pearson"),
                          conf.level = 0.95)
#paired t-test to see if MCMCTree and ReltimeOLS mean estimates are different
mean_diff_OLS <- t.test(df2$mean, df2$OLS_mean, paired = TRUE, alternative = "two.sided")
#just correlation test to see if range sizes correlate
cor.test(df2$MCMCTree_diff, df2$OLS_diff, 
         method = c("pearson"),
         conf.level = 0.95)
#paired t-test to see if MCMCTree and ReltimeOLS error ranges are different
t.test(df2$OLS_diff, df2$MCMCTree_diff, paired = TRUE, alternative = "greater")

mean_OLS_plot <- ggplot(df2, aes(x = mean, y = OLS_mean)) +
  geom_point()  +
  theme_bw() +
  geom_abline(slope=1, intercept = 0) +
  ylab("RelTime OLS estimate (Million years)") +
  xlab("MCMCTree estimate (Million years)") +
  xlim(0,110) +
  ylim(0,110) +
  annotate("text", x = 5, y = 100,  hjust = 0, 
           label = paste("R=",round(mean_corr_OLS$estimate,2),", p<0.0001",sep = "")) +
  annotate("text", x = 5, y = 95, hjust = 0, 
           label = paste("mean difference=",round(mean_diff_OLS$estimate,2),", p=",round(mean_diff_OLS$p.value,4),sep = ""))


#boxplots of ranges
node_age_boxplot <- ggplot(data = df_long, aes(x = name, y = value , col = name)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  scale_colour_manual(values=c("black", "#E7B800","grey")) +
  geom_point() +
  geom_line(aes(group = MCMCTree_nodenumber), color="grey", alpha = .5) +
  theme_bw() +
  #geom_text_repel(aes(label = MCMCTree_nodenumber), size = 2, col = "black") +
  theme(aspect.ratio = 1, legend.position = "none") +
  ylab("95% HPD interval size (Million years)") +
  xlab("Method") +
  scale_x_discrete(breaks=c("BLS_diff","MCMCTree_diff","OLS_diff"),
                   labels=c("RelTime BL", "MCMCTree","RelTime OLS"))

mean_plots <- plot_grid(mean_BL_plot, mean_OLS_plot, ncol = 2, labels = c("(a)","(b)"))

png("/Users/tbl19/The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/FIGURES/RelTime_MCMCTreee_nodeages_comparison.png",  res = 600, width = 210, height = 210, unit = "mm")
plot_grid(mean_plots, node_age_boxplot, ncol = 1,labels = c("","(c)"))
dev.off()
