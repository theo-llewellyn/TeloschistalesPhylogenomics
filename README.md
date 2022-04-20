# *Teloschistales* Metagenomics
Bioinformatic scripts/code for Llewellyn et al. (2022) Metagenomics shines light on the evolution og 'sunscreen' pigments in the *Teloschistales* (lichen-forming Ascomycota)

All scripts (except .R scripts) were run on the Imperial College London High Performance Computer. This HPC uses the PBS queueing system, therefore core/RAM/runtimes in .sh scripts are specified in PBS format. All scripts are written in the format of a single genome (replacing the word ACCESSION for the name of the sequence) but can be converted into array scripts to handle multiple genomes.

## 1. Metagenome Assembly
### 1.1 Quality assessment of Illumina reads
Uses fastq.gz paired end Illumina raw reads. Read trimming requires `TruSeq3-PE-2.fa` for TruSeq Nano Library prep and `NexteraPE-PE.fa` for Nextera XT Library prep. Both .fa adpater files are included in trimmmotatic v0.36 within the adapters directory.
`cd assembly/QC`
1. `qsub fastqc.sh` assesses raw read quality using FastQC
2. `qsub trimmomatic.sh` trims low-quality bases and adapters with Trimmomatic
3. `qsub fastqc_trimmed.sh` assesses read quality post-trimming

### 1.2 Metagenome assembly
`cd assembly/metagenome_assembly`
1. `qsub megahit.sh` metagenome assembly using MEGAHIT
2. `qsub metaspades.sh` metagenome assembly using MetaSPAdes

### 1.3 Metagenome assessment
`cd assembly/assessment`
1. `qsub quast.sh` assembly contiguity using QUAST
2. `qsub busco.sh` assembly completeness using BUSCO and the Ascomycota dataset

## 2. Mycobiont read filtering
The following steps filter the metagenome to retrieve only the contigs belonging to the Lecanoromycete mycobiont.

### 2.1 BlobTools (round 1)
Uses a DIAMOND blast of the contigs against the UniRef90 database which can be downloaded here (INSERT LINK) and a BLASTn against all Lecanoromycetes genomes in NCBI and three JGI Mycocosm genomes (Xanthoria parietina, Cladonia grayii and Usnea florida)
`cd mycobiont_filtering/BlobTools`
1. `qsub diamond.sh` DIAMOND blast against UniRef90
2. `qsub blastn.sh` BLASTn against Lecanoromycete database
3. `qsub bbamp.sh` calculates read coverage using BBTools bbmap function
4. `qsub make_covfile.sh` converts bbmap output to a 'covfile' for use in BlobTools
5. `qsub BlobTools.sh` uses BlobTools to visualise coverage, GC-content and blast results of contigs. Requires TAXONOMY file to taxify the output of blast and DIAMOND searches

### 2.2 CONCOCT
`cd mycobiont_filtering/CONCOCT`
1. `qsub concoct.sh` bins metagenome contigs into MAGs
2. `qsub make_cov_gc.sh` makes a coverage and gc_content file from the bbmap output to be used in the following step
3. `Rscript concoct_mags_plot.r` visualises CONCOCT binning and BlobTools blasts to identify Ascomycota bins. Requires the 'clustering_merged.csv' file in the concoct_output
4. `cat bin1.fa bin2.fa bin3.fa > Lecanoromycete_MAG.fa` merge potential mycobiont bins into a single file

### 2.3 Blobtools (round 2)
The `Lecanoromycete_MAG.fa` can then be run through the steps in 2.1 using the exact same scripts but replacing the metagenome assembly for `Lecanoromycete_MAG.fa` taking care to change the output file names so as not to overwrite the first round of BlobTools. The results can then be used to remove any remaining non-mycbiont reads as follows
`cd mycobiont_filtering/BlobTools_round2`
1. `qsub remove_contam.sh` extracts contig headers of contigs with a top blast of Ascomycota or 'no hit'. Requires the `.bestsum.table.txt` file from BlobTools
2. `qsub seqkit_bbmap.sh` extract contigs based on headers file from previous step using SeqKit and then pulls reads which map to those contigs using bbmap

### 2.4 Mycobiont assembly cleaning
`cd mycobiont_filtering/cleaning`
1. `qsub redundans.sh` uses redundans to remove redundant contifs, scaffold and close gaps

## 3. Annotation and orthology inference
### 3.1 Repeat Masking
`cd annotation_orthology/repeat_masking`
1. `qsub repeatmodeler.sh` generates repeat content library using RepeatModeler
2. `qsub funnanotate_sort.sh` sorts contigs using funnanotate pipeline
3. `qsub repeatmasker.sh` uses custom repeat library to softmask genome using RepeatMasker

### 3.2 Gene Prediction
`cd annotation_orthology/gene_prediction`
1. `qsub funnanotate_masked.sh` predicts genes from sorted, masked genome using funnanotate. This uses ESTs and proteins from Xanthoria parietina, Cladonia grayii and Usnea florida as evidence downloaded from JGI mycocosm.

### 3.3 Orthology inference
`cd annotation_orthology/orthology`
1. `cp *_CONCOCT_genes_masked/predict_results/*proteins.fa formatted_proteomes_45T` copies all predicted proteomes to a new directory called `ormatted_proteomes_45T`
2. `qsub orthofinder.sh` runs orthology inference using OrthoFinder

## 4. Phylogenomics
`cd phylogenomcis`
1. `Rscript Orthogroups_75percent.R` uses the `Orthogroups.GeneCount.tsv` file from OrthoFinder to extract a list of single copy orthologues present in at least 75% of taxa
2. `qsub extract_75_orthogroups.sh` pull 75% orthogroups and copy to new directory
3. `qsub mafft_trimAL_loop.sh` uses (MAFFT)[] to align each orthogroup and (TrimAL)[] to remove ambiguous regions
4. `qsub edit_protein_headers.sh` removes trailing information on protein headers so that only the species name remains. This is needed in order for tree building tools to recognise which sequences belong to the same genome
5. `qsub iqtree.sh` produces a concatenated maximum likelihood tree from all orthgroups alignments and also individual orthogroup 'gene trees' for each orthogroup separately using (IQTree)[]
6. `qsub iqtree_gfc.sh` calculates gene- and site-concordance factors using IQTree
7. `qsub astral.sh` produces coalescent-based species tree using (ASTRAL)[]

## 5. Secondary metabolite biosynthetic gene cluster (BGC) analysis
`cd SMGC_analysis`
1. `qsub antismash.sh` predicts BGCs from funnanotate predictions using (antismash)[]
2. `qsub edit_antismash_gbks.sh` edits the antismash .gbk files so that each one includes the name of the genome it came from. This prevent identical named files from being removed in the next stage
3. `cp *.gbk antismash_Leca45T_gbks` copies renamed .gbks from all genomes into a new directory
4. `qsub bigscape.sh` identifies BGC families using (BiG-SCAPE)[]
5. `Rscript BGC_cluster_analysis.r` Mantel tests and PCoA of BGC data in R

## 6. Anthraquinone BGC analysis
### 6.1 Functional annotation
Once anthraquinone BGCFs have been identified from BiG-SCAPE output the orthofinder orthogroup to which those PKSs belong can be further characterised.
`cd anthraquinones/functional_annotation`
1. `sed 's/*//g' orthogroup.fa > file_noasterisk.fa` removes any asterisks from orthgroup file as these can't be read by Interpro
2. `qsub interproscan.sh` functional annotation of orthogroup using InterProScan
3. `qsub extract_PANTHER.sh` extract PANTHER annotations from interpro output

### 6.2 PT domain analysis
The following steps combine the PT domains from our anthraquinone PKSs with the Pfam seed alignment for that domain which can be downloaded [here](https://pfam.xfam.org/family/PF14765#tabview=tab3). The orthogroup of interest may also contain non-anthraqunone PKSs, therefore we need to make a note of which ones are putative anthraquinone PKSs by hovering over the relevant PKSs in the BiG-SCAPE PKSI html file. The names of the putative anthraquinone PKSs need to be saved into a file called `OG12_anthraquinone_headers.txt`

`cd anthraquinones/hmmalign`
1. `qsub hmm_convert.sh` converts the Pfam seed alignment to an HMM file
2. `qsub pull_anthraquinone_pks.sh` extract just anthraquinone PKSs from orthogroup of interest
3. `qsub blastp.sh` extracts PT domains from PKSs using aptA Aspergillus PT domain as a query
4. `qsub convert_2_fasta.sh` extract PT coordinates from blastp and convert to a fast file [bedtools]()
5. The above two steps are repeated for the Liu et al. (2015) sequences and then combined into a single file of PT domains with our putative anthraquinone PT domains
6. `qsub hmmalign.sh` uses (this)[https://github.com/reubwn/scripts/blob/master/hmmsearch-easy.pl] script from Reuben Nowell to align the PT domains to those in the Pfam seed alignment `PF14765_seed.txt`
7. `qsub iqtree_hmmalign.sh` uses IQTree to build a ML tree from the clustal alignment
8. `qsub hmmalign_taxify.sh` uses (this)[https://github.com/reubwn/scripts/blob/master/taxify_uniprot_treefile.pl] script from Reuben Nowell to taxify the output. Requires the Uniref90 taxlist used for BlobTools.
