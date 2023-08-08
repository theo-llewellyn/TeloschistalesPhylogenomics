#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=8:mem=100gb

module load anaconda3/personal
source activate modeltest-ng-env

cd /rds/general/user/tbl19/home/modeltest-ng/Telos553T_markers
cp /rds/general/user/tbl19/home/Telos553T_msa/Telos553T_muscle*_concatenated.fa .
cp /rds/general/user/tbl19/home/Telos553T_msa/Telos553T_muscle*_partitions.txt .

~/software/modeltest-ng-static \
 -d nt \
 -i Telos553T_muscle_concatenated.fa \
 -q Telos553T_muscle_partitions.txt \
 -p 6 \
 -T raxml

~/software/modeltest-ng-static \
 -d nt \
 -i Telos553T_muscle_tAl_concatenated.fa \
 -q Telos553T_muscle_tAl_partitions.txt \
 -p 6 \
 -T raxml
