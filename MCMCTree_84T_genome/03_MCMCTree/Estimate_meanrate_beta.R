#----------------#
# LOAD LIBRARIES #
#----------------#
library( ape )
library( rstudioapi )

#source function that estimates beta from gene trees
source("~/The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/04_R_SCRIPTS/myRpackages.R")
setwd("~/The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/PHYLOGENOMICS/TIME_CALIBRATION/MCMCTREE/00_Estimate_beta/Telos84T_genetrees/")

#function needs concatenated file of rooted gene trees and the min and max age of node in Mya divided by 100 i.e. 137 My = 1.37
estimate_beta("rooted_trees.tree",1.37,2.20)
