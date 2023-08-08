library(tidyr)
library(dplyr)
library(readr)
library(tibble)
library(janitor)

setwd("~/The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/ORTHOFINDER/Results_Leca125T/")

#load in data
Orthogroups <- read_tsv("Orthogroups.GeneCount.tsv")
head(Orthogroups)

#filter for rows where all columns are either 0 or 1 and then where the total number of genes is at least half the number of taxa
#ncol(Orthogroups) - 2 gives the number of taxa as each taxon is a column but there's two extra columns with orthogroup and total which need to be removed
#can change the 0.75 to higher proportion to reduce amount of missing data
Orthogroups %>% filter(across(ends_with(c("renamed","proteins")), ~ .x < 2), Total >= (ncol(Orthogroups) -2) * 0.5) -> single_copy_orthologues_50
nrow(single_copy_orthologues_50)

#filtering just for taxa in the Teloschistales and Caliciales and removing genomes with >5% duplicated BUSCOs and >70% Complete single copy BUSCOs
Teloschistales_dup_comp <- Orthogroups[,c(1,6,8:10,12:17,19:30,39,41:53,55:57,62:66,73,77:84,86:95,102:115,119:120,122:126)]
#same as above but filtering for genes present in at least 10 genomes and the orthogroup must contain Xanthoria parietina and T. chrysopthalmus 125 to be able to calculate pairwise distances
Teloschistales_dup_comp %>% filter(across(ends_with(c("renamed","proteins")), ~ .x < 2)) %>%
  filter(Xanthoria_parietina_renamed ==1 & Teloschistes_chrysophthalmus_125.proteins == 1) %>%
  adorn_totals(c("row", "col")) %>%
  filter(Total >= 10) -> teloschistales_dup_comp_orthologues_10sp
nrow(teloschistales_dup_comp_orthologues_10sp)

#save orthogroup column as a text file
write.table(teloschistales_dup_comp_orthologues_10sp$Orthogroup,"Orthologues_Telos85T_dup_comp_10sp.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
