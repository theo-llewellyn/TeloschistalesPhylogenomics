#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=8:mem=100gb

module load anaconda3/personal
source activate iqtree-env 

cd /rds/general/project/theollewellynproject/live/IQTree/IQTree_Telos552T

#untrimmed
iqtree2 \
 -s ~/Telos553T_msa/Telos552T_muscle_concatenated.fa \
 -p ~/Telos553T_msa/Telos552T_muscle_partitions.txt \
 -m GTR+I+G \
 -g /rds/general/project/theollewellynproject/live/raxml-ng/Telos85T_4parts_partitions12/*.support_renamed \
 --prefix Telos552T_genomeconstrain \
 -B 1000 \
 -T 8 \
 --undo \
 -nm 2000

#trimmed
 iqtree2 \
 -s ~/Telos553T_msa/Telos552T_muscle_tAl_concatenated.fa \
 -p ~/Telos553T_msa/Telos552T_muscle_tAl_partitions.txt \
 -m GTR+I+G \
 -g /rds/general/project/theollewellynproject/live/raxml-ng/Telos85T_4parts_partitions12/*.support_renamed \
 --prefix Telos552T_genomeconstrain_tAl \
 -B 1000 \
 -T 8 \
 --undo \
 -nm 2000
