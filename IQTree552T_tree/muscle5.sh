#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=8:mem=10gb
#PBS -J 1-7

module load anaconda3/personal
source activate muscle5-env

#get name of contigs file
GENE=$(head -n $PBS_ARRAY_INDEX /rds/general/user/tbl19/home/Genes1.txt|tail -n 1)

cd /rds/general/user/tbl19/home/Telos553T_multifastas

/rds/general/user/tbl19/home/software/muscle5/muscle5.1.linux_intel64 \
 -align Telos553T_genomes85_${GENE}.fa \
 -output /rds/general/user/tbl19/home/Telos553T_msa/Telos553T_${GENE}_muscle5_msa.fa \
 -threads 8
