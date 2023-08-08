#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=32:mem=124gb

module load anaconda3/personal
source activate OrthoFinder-env

cd $PBS_O_WORKDIR

/rds/general/user/tbl19/home/software/OrthoFinder/orthofinder \
 -f /rds/general/project/theollewellynproject/live/DNA_Repair_Analysis/formatted_proteomes_Leca117T \
 -t 32 \
 -n Leca117T

conda deactivate
