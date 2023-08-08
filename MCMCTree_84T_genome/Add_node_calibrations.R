#Add node calibrations from Gaya et al 2015 to whole genome tree 85T
#add node calibrations
setwd("~/The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/PHYLOGENOMICS/TIME_CALIBRATION/")

#read in calibrations file
calibrations <- read.table( file = "Calibrations_lichens_maintree.txt",
                            stringsAsFactors = F, sep = "|", blank.lines.skip = T )
colnames( calibrations ) <- c( "name", "MCMCtree calib" )

#read in tree
tree.file <- list.files(pattern = "Telos85T_4parts_partitions12.raxml.support_rooted_labelled_maintree.tre" )
trees.calibs <- readLines(tree.file)

#replace node names with calibrations
for( j in 1:dim( calibrations )[1] ){
  trees.calibs <- gsub( pattern = paste0("\\[",calibrations[j,1],"\\]"),
                        x = trees.calibs,
                        replacement = paste( "'", calibrations[j,2], "'", sep = "" ) )
}

#save as new file
write(x = trees.calibs, file = "Telos85T_4parts_partitions12.raxml.support_rooted_labelled_MCMCtree_calib_1a.tree")
