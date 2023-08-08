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

### 2.1 BlobTools (round 1)
Uses a DIAMOND blast of the contigs against the UniRef90 database which can be downloaded [here](https://ftp.expasy.org/databases/uniprot/current_release/uniref/uniref90/uniref90.fasta.gz) and a BLASTn against all Lecanoromycetes genomes in NCBI and three JGI Mycocosm genomes ([Xanthoria parietina](https://mycocosm.jgi.doe.gov/Xanpa2/Xanpa2.home.html), [Cladonia grayii](https://mycocosm.jgi.doe.gov/Clagr3/Clagr3.home.html) and [Usnea florida](https://mycocosm.jgi.doe.gov/Usnflo1/Usnflo1.home.html))  
`cd mycobiont_filtering/BlobTools`
1. `qsub diamond.sh` [DIAMOND](https://github.com/bbuchfink/diamond) blast against UniRef90
2. `qsub blastn.sh` [BLASTn](https://blast.ncbi.nlm.nih.gov/Blast.cgi) against Lecanoromycete database
3. `qsub bbmap.sh` calculates read coverage using [BBTools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) bbmap function
4. `qsub taxlist_diamond.sh` link Uniref taxids to taxa
5. `qsub BlobTools.sh` uses [BlobTools](https://github.com/DRL/blobtools) to visualise coverage, GC-content and blast results of contigs. Requires taxonomy file to taxify the output of blast and DIAMOND searches. NCBI taxlist can be downloaded [here](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&ved=2ahUKEwjmq-eb9qT3AhXcQkEAHZ4TCCMQFnoECAMQAQ&url=https%3A%2F%2Fftp.ncbi.nih.gov%2Fpub%2Ftaxonomy%2Faccession2taxid%2Fnucl_wgs.accession2taxid.gz&usg=AOvVaw2Oeb-8gVxs3HaGSJh4Ck4L) and supplemented with the taxids for JGI genomes by adding rows at the bottom.

### 2.2 CONCOCT
`cd mycobiont_filtering/CONCOCT`
1. `qsub bwa_gatk.sh` align reads to contigs using [BWA-mem](https://github.com/lh3/bwa) and convert to .bam with [GATK](https://gatk.broadinstitute.org/hc/en-us)
2. `qsub concoct.sh` bins metagenome contigs into MAGs using [CONCOCT](https://github.com/BinPro/CONCOCT)
3. `qsub make_cov_gc.sh` makes a coverage and gc_content file from the bbmap output to be used in the following step

### 2.3 BinArena
Visualises the results of Blobtools and CONCOCT in order to merge bins belonging to mycobiont
`cd mycobiont_filtering/BinArena`
1. `qsub make_binarena_tsv.sh` extract gc and coverage from Blobtools output and concoct bins and saves into tsv
2. Open this file in BinArena using the settings x=gc , y=cov to 4rt, colour=tax
3. `qsub seqkit_array.sh` pull contigs using headers
4. `qsub filter_BlobTools.sh` filter blobtools input for filtered concoct contigs

### 2.4 Blobtools (round 2)
`cd mycobiont_filtering/BlobTools_round2`
1. `qsub BlobTools_round2.sh` rerun BlobTools
2. `qsub remove_contam.sh` remove final contaminants
3. `qsub pull_reads.sh` pull reads

### 2.4 Mycobiont assembly cleaning
`cd mycobiont_filtering/cleaning`
1. `qsub redundans.sh` uses [redundans](https://github.com/lpryszcz/redundans) to remove redundant contigs, scaffold and close gaps

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


