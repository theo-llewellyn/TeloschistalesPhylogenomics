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
`cd annotation_orthology`
1. Repeat Modelling and masking
2. Gene prediction with funannotate
3. Orthology inference with Orthfinder

## 4. RAxML 83 taxon genome-scale phylogenetic tree
`cd raxml84T_tree`
1. Orthogroup filtering
2. Protein alignment
3. Convert AA to nucleotide alignment
4. Partition codons
5. Make gene trees
6. Filter gene trees
7. Make species tree

