library(ape)
library(tidyverse)
library(phytools)
library(caper)

setwd("~")
tree_84T     <- read.tree("The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/PHYLOGENOMICS/TIME_CALIBRATION/MCMCTREE/LICHENS_MCMCTree_553T/00_Esters1a_calibration/FigTree_84sp_nodelabels_renamed.tree" )
#if using raxml tree
tree_553T <- read.tree("The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/PHYLOGENOMICS/MARKER_GENES/TREES/Telos553T_muscle_raxml_constrain.raxml.support.rooted.tree")
#if IQTree
tree_553T <- read.tree("The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/PHYLOGENOMICS/MARKER_GENES/TREES/Telos553T_genomeconstrain_v3.treefile.rooted.tree")
#remove branch lengths and support
tree_553T$edge.length<-NULL
tree_553T$node.label<-rep("",length(tree_553T$node.label))
#remove dashes from tip labels
tree_553T$tip.label <- gsub(pattern = "'", replacement = "", tree_553T$tip.label)

#save unlabelled tree for BASEML
#get phylip header
phylip_header <- paste( length(tree_553T$tip.label), "  1", sep = "" )
#write tree file
write( x = phylip_header, file = "The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/PHYLOGENOMICS/TIME_CALIBRATION/MCMCTREE/Telos553T_rooted_baseml.tree" )
write.tree(tree_553T, file = "The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/PHYLOGENOMICS/TIME_CALIBRATION/MCMCTREE/Telos553T_rooted_baseml.tree", append = TRUE )
#if IQTree
write( x = phylip_header, file = "The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/PHYLOGENOMICS/TIME_CALIBRATION/MCMCTREE/LICHENS_MCMCTree_553T/00_Esters1a_calibration/Telos553T_v3_IQTree_rooted_baseml.tree" )
write.tree(tree_553T, file = "The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/PHYLOGENOMICS/TIME_CALIBRATION/MCMCTREE/LICHENS_MCMCTree_553T/00_Esters1a_calibration/Telos553T_v3_IQTree_rooted_baseml.tree", append = TRUE )


SkewTs <- read_tsv("The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/PHYLOGENOMICS/TIME_CALIBRATION/MCMCTREE/LICHENS_MCMCTree_553T/00_Esters1a_calibration/ST.fitted.dists.G2.40_rownames.tsv")

#for each node in genome tree
for (i in 1:nrow(SkewTs)){
  #get node number
  nodenumber <- as.numeric(gsub(pattern = "t_n", replacement = "", x = SkewTs$node[i]))
  #get taxa descending from that node
  descendents <- clade.members(nodenumber, tree_84T, tip.labels = TRUE)
  #remove the 1_ format so that they match the 553T tree
  descendents <- gsub(pattern = "\\d+?_", replacement = "", descendents)
  #get node label in 553T tree
  nodename_553T <- getMRCA(tree_553T, descendents)
  #get SkewT calibration of that node
  SkewTvalue <- SkewTs[[i,6]]
  #replace parentheses and commas with asterisks otherwise ape cant interpret as node
  SkewTvalue <- gsub(pattern = "\\(", replacement = "*", x = SkewTvalue)
  SkewTvalue <- gsub(pattern = "\\)", replacement = "£", x = SkewTvalue)
  SkewTvalue <- gsub(pattern = ",", replacement = "^", x = SkewTvalue)
  #add that SkewT as a node label in the marker gene tree
  tree_553T$node.label[nodename_553T-length(tree_553T$tip.label)] <- paste("'",SkewTvalue,"'",sep="")
}


#write tree file
write( x = phylip_header, file = "The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/PHYLOGENOMICS/TIME_CALIBRATION/MCMCTREE/LICHENS_MCMCTree_553T/00_Esters1a_calibration/Telos553T_v3_MCMCtree_calib_symbols.tree" )
write.tree( phy = tree_553T, file = "The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/PHYLOGENOMICS/TIME_CALIBRATION/MCMCTREE/LICHENS_MCMCTree_553T/00_Esters1a_calibration/Telos553T_v3_MCMCtree_calib_symbols.tree", append = TRUE )

#we have to read it back in again to replace the node label special characters with the original symbols, this way we circumvent apes reformatting problem
tree_553T_txt <- readLines(con = "The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/PHYLOGENOMICS/TIME_CALIBRATION/MCMCTREE/LICHENS_MCMCTree_553T/00_Esters1a_calibration/Telos553T_v3_MCMCtree_calib_symbols.tree")
tree_553T_txt <- gsub(replacement = "\\(", pattern = "\\*", x = tree_553T_txt)
tree_553T_txt <- gsub(replacement = "\\)", pattern = "£", x = tree_553T_txt)
tree_553T_txt <- gsub(replacement = ",", pattern = "\\^", x = tree_553T_txt)
write(x = tree_553T_txt, file = "The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/PHYLOGENOMICS/TIME_CALIBRATION/MCMCTREE/LICHENS_MCMCTree_553T/00_Esters1a_calibration/Telos553T_v3_MCMCtree_calib.tree")
