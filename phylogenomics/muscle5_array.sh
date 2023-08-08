#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=8:mem=10gb
#PBS -J 1-715

module load anaconda3/personal
source activate muscle5-env

cd $PBS_O_WORKDIR

#get name of OG file
OG=$(head -n $PBS_ARRAY_INDEX /rds/general/project/theollewellynproject/live/DNA_Repair_Analysis/formatted_proteomes_Leca118T/OrthoFinder/Results_Results_Leca118T/Orthologues_Leca118T_50p.txt|tail -n 1)

/rds/general/user/tbl19/home/software/muscle5/muscle5.1.linux_intel64 \
 -align /rds/general/project/theollewellynproject/live/DNA_Repair_Analysis/formatted_proteomes_Leca118T/OrthoFinder/Results_Results_Leca118T/Single_Copy_Orthologues_50percent/${OG}.fa \
 -output /rds/general/user/tbl19/home/Leca118T_msa/${OG}_muscle5_msa.fa \
 -threads 8
