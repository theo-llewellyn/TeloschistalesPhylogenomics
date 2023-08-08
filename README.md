# *Teloschistales* Metagenomics
Bioinformatic scripts/code for Chapter 3: Phylogenomics and molecular systematics of the _Teloschistaceae_

Scripts were run on the Imperial College London High Performance Computer or Royal Botanic Gardens Kew HPC. These HPCs use the PBS queueing system, therefore core/RAM/runtimes in .sh scripts are specified in PBS format. All scripts are written for a single genome file (replacing the word ACCESSION for the name of the sequence) but can be converted into array scripts to handle multiple genomes.

## 1. Metagenome Assembly
### 1.1 Quality assessment of Illumina reads
Uses fastq.gz paired end Illumina raw reads. Read trimming requires `TruSeq3-PE-2.fa` for TruSeq Nano Library prep and `NexteraPE-PE.fa` for Nextera XT Library prep. Both .fa adpater files are included in trimmmotatic v0.36 within the adapters directory.  
`cd assembly/QC`
1. `qsub fastqc.sh` assesses raw read quality using [FastQC](https://github.com/s-andrews/FastQC)
2. `qsub trimmomatic.sh` trims low-quality bases and adapters with [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
3. `qsub fastqc_trimmed.sh` assesses read quality post-trimming

### 1.2 Metagenome assembly
`cd assembly/metagenome_assembly`
1. `qsub megahit.sh` metagenome assembly using [MEGAHIT](https://github.com/voutcn/megahit)

### 1.3 Metagenome assessment
`cd assembly/assessment`
1. `qsub quast.sh` assembly contiguity using [QUAST](https://github.com/ablab/quast)
2. `qsub busco.sh` assembly completeness using [BUSCO](https://busco.ezlab.org/) and the Ascomycota dataset

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

## 4. Phylogenomics
`cd phylogenomcis`
1. `Rscript Orthogroups_50percent.R` uses the `Orthogroups.GeneCount.tsv` file from OrthoFinder to extract a list of single copy orthologues present in at least 50% of taxa
2. `qsub extract_50_orthogroups.sh` pull 50% orthogroups and copy to new directory
3. `qsub mafft_trimAL_loop.sh` uses [MAFFT](https://mafft.cbrc.jp/alignment/software/) to align each orthogroup and [TrimAl](http://trimal.cgenomics.org/) to remove ambiguous regions, script also removes trailing information on protein headers so that only the species name remains. This is needed in order for tree building tools to recognise which sequences belong to the same genome
5. `qsub iqtree.sh` produces a concatenated maximum likelihood tree from all orthgroups alignments and also individual orthogroup 'gene trees' for each orthogroup separately using [IQ-Tree](https://github.com/iqtree/iqtree2)

## 5. Selection analysis
`cd selection_analysis`
The following scripts pull the genes of interest from the genomes and perform selection analysis to test to positive selection in the _Teloschistales_. The scripts are for the Hsp90 gene but are the same for all three genes

1. `qsub pull_gene.sh` Pull gene of interest if know it already or do a BlastP against predicted proteins using a query
2. `qsub muscle5.sh` Align with MUSCLE5
3. `qsub pull_CDS.sh` pull the CDS regions for those sequences
4. `qsub pal2nal.sh` run PAL2NAL using the aligned AA and unaligned CDSs to get aligned CDSs for CODEML
5. Manually inspect the alignment for strange things and replace the terminal gaps for ?. Check if there are any sequences that were split into two and merge if they are identical at overlapping regions. Save as GENE_Leca118T_muscle5_msa_codon_trimmed_renamed.fa
6. `qsub format_species_tree.sh` root with pxrr and remove branch lengths and support values
7. `qsub edit_protein_headers.sh` edit fasta headers in gene alignment so they just contain name of taxon
8. `qsub Drop_tips_PAML.R` drop tips from species tree that arent in gene alignments
9. Label Teloschistaceae branch with '#1'
10. `qsub codeml.sh` Selection analysis with PAML. We will test null and alternative models for both branch and branchsite models using Teloschistaceae as the foreground branch. Therefore we have four analysis per gene. A template directory (TEMPLATE) is included which can be used for all genes replacing the name of the alignment and species tree in the .ctl files. All four models can then be run. Remember to make a new copy of TEMPLATE for each gene to avoid overwriting outputs.
11. likelihood ratio test to compare the null and the alternative models: M0 vs branch & branch-site null vs. branch. Calculate t-test statistic: 2 x (log likelihood branch model - log likelihood M0), df = 1, chi2: `paml4.9j/bin/chi2 1 LRT`

