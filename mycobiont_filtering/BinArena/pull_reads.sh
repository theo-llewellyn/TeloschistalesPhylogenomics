#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=16:mem=124gb

#get name of contigs file
ACCESSION=ACCESSION
module load anaconda3/personal
source activate seqkit-env

cd $PBS_O_WORKDIR

#opens the assembly, searches for particular contigs and puts into a new file
cat /rds/general/project/theollewellynproject/live/filtered_Asco_contigs/${ACCESSION}_concoct_merged_contigs.fa | /rds/general/user/tbl19/home/anaconda3/envs/seqkit-env/bin/seqkit grep -f /rds/general/project/theollewellynproject/live/filtered_Asco_contigs/${ACCESSION}_concoct_Asco_headers.txt > /rds/general/project/theollewellynproject/live/filtered_Asco_contigs/${ACCESSION}_concoct_Asco_contigs.fa

conda deactivate

#pull reads mapping to those contigs with bbmap
module load bbmap

bbmap.sh -Xmx62g \
nodisk=t \
ref=/rds/general/project/theollewellynproject/live/filtered_Asco_contigs/${ACCESSION}_concoct_Asco_contigs.fa \
in1=/rds/general/project/theollewellynproject/live/data/Trimmed_reads/${ACCESSION}_Trimmed_1P.fq.gz \
in2=/rds/general/project/theollewellynproject/live/data/Trimmed_reads/${ACCESSION}_Trimmed_2P.fq.gz \
threads=16 \
outm=/rds/general/project/theollewellynproject/ephemeral/filtered_Asco_reads/${ACCESSION}_concoct_Asco_R#.fq.gz
