#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=32:mem=124gb

module load anaconda3/personal
source activate iqtree-env 


mkdir /rds/general/project/theollewellynproject/live/IQTree/IQTree_Leca117T_50p_tAL
cd /rds/general/project/theollewellynproject/live/IQTree/IQTree_Leca117T_50p_tAL

#make concatenated species tree
iqtree2 -p /rds/general/user/tbl19/home/Leca117T_msa_tAl \
 --prefix concat_OF_tAl \
 -B 1000 \
 -T 32

conda deactivate
