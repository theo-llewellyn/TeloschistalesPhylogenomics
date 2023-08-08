#PBS -l walltime=00:20:00
#PBS -l select=1:ncpus=1:mem=1gb

cd $PBS_O_WORKDIR

/rds/general/user/tbl19/home/software/pal2nal.v14.pal2nal.pl \
 /rds/general/project/theollewellynproject/live/DNA_Repair_Analysis/Hsp90_Leca117T_muscle5_msa.fa \
 /rds/general/project/theollewellynproject/live/DNA_Repair_Analysis/Hsp90_CDS_Leca117T.fa \
 -output fasta
