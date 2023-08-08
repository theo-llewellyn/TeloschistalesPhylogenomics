#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=16:mem=124gb
#PBS -J 1-83

#get name of contigs file
ACCESSION=$(head -n $PBS_ARRAY_INDEX /rds/general/user/tbl19/home/macrogen_2022_genomes.txt|tail -n 1)

module load anaconda3/personal
source activate seqkit-env

cd $PBS_O_WORKDIR

#opens the assembly, searches for particular contigs and puts into a new file
cat /rds/general/project/theollewellynproject/live/megahit_metagenome_assembly/megahit_${ACCESSION}/final.contigs.fa | /rds/general/user/tbl19/home/anaconda3/envs/seqkit-env/bin/seqkit grep -f /rds/general/project/theollewellynproject/live/filtered_Asco_contigs/${ACCESSION}_concoct_merged_headers.txt > /rds/general/project/theollewellynproject/live/filtered_Asco_contigs/${ACCESSION}_concoct_merged_contigs.fa

conda deactivate
