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

#
