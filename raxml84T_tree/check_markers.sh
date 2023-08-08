#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=1:mem=1gb
#PBS -J 1-7

module load blast+

#checks whether any of the marker genes ITS, mitSSU, nucLSU, RPB1, RPB2, RET1 and MCMC7 are in genomic data. Uses LIQ125TCHR as example
cd $PBS_O_WORKDIR

GENE=$(head -n $PBS_ARRAY_INDEX Genes.txt|tail -n 1)

blastn \
-num_threads 8 \
-query /rds/general/project/theollewellynproject/ephemeral/${GENE}_query.fa \
-db /rds/general/project/theollewellynproject/ephemeral/LIQ125TCHR_cds \
-evalue 0.05 \
-word_size 11 \
-gapopen 5 \
-gapextend 2 \
-max_target_seqs 100 \
-max_hsps 1 \
-outfmt '6 sseqid sseq' \
-out /rds/general/project/theollewellynproject/ephemeral/${GENE}.vs.LIQ125TCHR_cds.blastn.out

HEADER=$(cut -f 1 ${GENE}.vs.LIQ125TCHR_cds.blastn.out)

ORTHGROUGROUP=$(grep $HEADER /rds/general/project/theollewellynproject/ephemeral/Single_Copy_Orthologues_10sp_CDS/* | grep 'OG.*fa' -o)

echo Orthogroup is $ORTHOGROUP
