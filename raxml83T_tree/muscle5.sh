#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=8:mem=10gb
#PBS -J 1-1524

module load anaconda3/personal
source activate muscle5-env

#get name of OG file
OG=$(head -n $PBS_ARRAY_INDEX /rds/general/user/tbl19/home/OGs_Telos85T_10sp.txt|tail -n 1)
PREFIX=${OG%%.fa}

/rds/general/user/tbl19/home/software/muscle5/muscle5.1.linux_intel64 \
 -align /rds/general/project/theollewellynproject/live/OrthoFinder/Results_Leca125T/Single_Copy_Orthologues_10sp/${OG} \
 -output /rds/general/user/tbl19/home/Telos85T_10sp_msa/${PREFIX}_muscle5_msa.fa \
 -threads 8
