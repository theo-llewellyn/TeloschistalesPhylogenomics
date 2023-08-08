library(ape)

#read in full species tree
lichen_tree <- read.tree("Telos85T_4parts_partitions12.raxml.rooted.tre")

for(gene in list.files()){
  #change to that directory
  setwd(gene)
  #read in alignment
  alignment <- read.dna(dir(pattern='partitions12_Telos.aln')[1],format = "sequential")
  #pull taxa names
  taxa <- rownames(alignment)
  #keep only tips in alignment
  subsetted_tree <- keep.tip(lichen_tree, taxa)
  #save subsetted tree
  write.tree(subsetted_tree,"Telos85T_4parts_partitions12.raxml.rooted.subset.tre")
  setwd("../")
}
