#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=8:mem=10gb

module load anaconda3/personal
source activate muscle5-env
/rds/general/user/tbl19/home/software/muscle5/muscle5.1.linux_intel64 \
 -align /rds/general/project/theollewellynproject/live/DNA_Repair_Analysis/Hsp90_Leca118T.fa \
 -output /rds/general/project/theollewellynproject/live/DNA_Repair_Analysis/Hsp90_Leca118T_muscle5_msa.fa \
 -threads 8
