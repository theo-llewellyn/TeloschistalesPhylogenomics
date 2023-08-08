#Subset species tree for partial alignments for PAML

library(ape)


#for MCMCTree
setwd("~/The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/SELECTION_ANALYSIS/SPECIES_TREES/Leca117T/")

#read in full species tree
lichen_tree <- read.tree("Leca117T_concat_OF_tAl_rooted.treefile")

setwd("../../ALIGNMENTS/Leca117T/")

alignmenthsp90 <- read.dna('Hsp90_Leca117T_muscle5_msa_codon_trimmed_renamed.fa',format = "fasta")
#pull taxa names
taxa <- rownames(alignmenthsp90)
#keep only tips in alignment
subsetted_tree <- keep.tip(lichen_tree, taxa)
#unroot tree
unrooted_tree <- unroot(subsetted_tree)
#save subsetted tree
write.tree(unrooted_tree,"../../SPECIES_TREES/Leca117T/SpeciesTree_rooted_renamed_Hsp90.tre")
