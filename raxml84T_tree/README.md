## 4. RAxML 83 taxon genome-scale phylogenetic tree
The following steps produce a genome-scale species tree for the 83 taxon Teloschistales dataset

`cd raxml83T_tree`
1. `Rscript Orthogroups_10sp.R` uses the `Orthogroups.GeneCount.tsv` file from OrthoFinder to filter genomes with <5% duplicated BUSCOs and >70% Complete single copy BUSCOs and then extracts a list of single copy orthologues present in at least 10 taxa
2. `qsub pull_Telos_seqs.sh` pull orthogroups for Teloschistales genomes and copy to new directory
3. `qsub pull_CDS.sh` pull CDS sequences
4. `qsub check_markers.sh` check none of the marker genes are in the orthogroups
5. `qsub muscle5.sh` align proteins with MUSCLE5 [Muscle5]([https://mafft.cbrc.jp/alignment/software/](https://github.com/rcedgar/muscle))
6. `qsub pal2nal.sh` convert to nucleotides using PAL2NAL
7. `qsub partition_alignment.sh` split alignment into two partitions (1st and 2nd codons) and (3rd codons)
8. `qsub raxml_genetrees.sh` make gene trees
9. Remove genes with outlier long branches and order from slow to fast as in Alvarez-Carretero et al 2022. Submit the following scripts: `00_Filtering_genes.sh`, `00_Get_filtered_genes_in_dir.R`, `01.1_Get_data_for_baseml.sh`, `01.1_Get_data_for_baseml_partitioned.sh`, `01.2_Get_mouse_human_aln.sh`, `01_Analysis_filtered_genes.R`, `03_Concatenate_genes_separated_for_partition.sh`, `fasta-phylip-partitions.sh`
10. `raxml84T_speciestree.sh` run raxml on the concatenated alignment with 4 partitions of slow to fast evolving genes using GTRGAMMA model
