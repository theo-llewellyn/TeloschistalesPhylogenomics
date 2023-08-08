#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=1:mem=1gb
#PBS -J 1-7

module load blast+

#get name of contigs file
GENE=$(head -n $PBS_ARRAY_INDEX /rds/general/user/tbl19/home/GenePull/Genes.txt|tail -n 1)

cd $PBS_O_WORKDIR

#converts the genomes to a searchable blast database made up of multiple files named -out
cat *genome.fa > Leca118T_genomes.fa

makeblastdb \
-in Leca118T_genomes.fa \
-out Leca118T \
-parse_seqids \
-dbtype prot

blastn \
-num_threads 8 \
-query /rds/general/user/tbl19/home/GenePull/${GENE}.fasta \
-db Leca118T_genomes \
-evalue 0.05 \
-word_size 11 \
-gapopen 5 \
-gapextend 2 \
-max_target_seqs 1000 \
-max_hsps 1 \
-outfmt '6 sseqid sseq' \
-out /rds/general/project/theollewellynproject/ephemeral/${GENE}.vs.Leca108T_genomes.blastn.out

#turn output into multifasta files of all the hits
cat  /rds/general/project/theollewellynproject/ephemeral/${GENE}.vs.Leca108T_genomes.blastn.out | awk '{print">"$1"\n"$2}' >  /rds/general/project/theollewellynproject/ephemeral/Leca108T_${GENE}.fa
