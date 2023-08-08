# *Teloschistales* Metagenomics
Bioinformatic scripts/code for Chapter 3: Phylogenomics and molecular systematics of the _Teloschistaceae_

Scripts were run on the Imperial College London High Performance Computer or Royal Botanic Gardens Kew HPC. These HPCs use the PBS queueing system, therefore core/RAM/runtimes in .sh scripts are specified in PBS format. All scripts are written for a single genome file (replacing the word ACCESSION for the name of the sequence) but can be converted into array scripts to handle multiple genomes.

The following is an overview of the main steps in the pipeline. For detailed instructions for each step see the README files in each subdirectory

## 1. Metagenome Assembly
`cd assembly`
1. Quality assessment of Illumina reads with fastqc and trimmomatic
2. Metagenome assembly with megahit
3. Metagenome assessment with quast and BUSCO

## 2. Mycobiont read filtering
The following steps filter the metagenome to retrieve only the contigs belonging to the Lecanoromycete mycobiont.
`cd mycobiont_filtering`
1. BlobTools (round 1)
2. Binning with CONCOCT
3. Bin Refinement with BinArena
4. BlobTools (round 2)
5. Mycobiont assembly cleaning with redundans

## 3. Annotation and orthology inference
### 3.1 Repeat Masking
`cd annotation_orthology/repeat_masking`
1. `qsub funnanotate_sort.sh` sorts contigs using [funnanotate](https://github.com/nextgenusfs/funannotate) pipeline
2. `qsub repeatmodeler.sh` generates repeat content library using [RepeatModeler](https://www.repeatmasker.org/RepeatModeler/)
3. `qsub repeatmasker.sh` uses custom repeat library to softmask genome using [RepeatMasker](https://www.repeatmasker.org/RepeatMasker/)

### 3.2 Gene Prediction
`cd annotation_orthology/gene_prediction`
1. `qsub funnanotate_masked.sh` predicts genes from sorted, masked genome using funnanotate. This uses ESTs and proteins from Xanthoria parietina, Cladonia grayii and Usnea florida as evidence downloaded from JGI mycocosm.
2. `qsub rename_downloaded_proteins.sh` renames the protein headers from assemblies downloaded from NCBI/JGI so that they are unique between proteomes

### 3.3 Orthology inference
`cd annotation_orthology/orthology`
1. `cp *_CONCOCT_genes_masked/predict_results/*proteins.fa formatted_proteomes_117T` copies all predicted proteomes to a new directory called `formatted_proteomes_117T`
2. `qsub orthofinder.sh` runs orthology inference using [OrthoFinder](https://github.com/davidemms/OrthoFinder)

## 4. RAxML 83 taxon genome-scale phylogenetic tree
The first step is to infer a genome-scale species tree for the 83 taxon Teloschistales dataset
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


