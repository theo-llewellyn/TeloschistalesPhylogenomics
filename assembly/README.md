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
