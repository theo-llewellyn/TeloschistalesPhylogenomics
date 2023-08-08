#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=99:mem=100gb

mkdir /rds/general/project/theollewellynproject/live/raxml-ng/Telos85T_4parts_partitions12
cd /rds/general/project/theollewellynproject/live/raxml-ng/Telos85T_4parts_partitions12

#parse trimmed alignment
~/software/raxml-ng --parse \
         --msa ~/Telos85T_10sp_msa/000_alignments/Telos85T_4parts_partitions12.fa \
         --model ~/Telos85T_10sp_msa/000_alignments/Telos85T_4parts_partitions12.txt \
         --prefix Telos85T_4parts_partitions12

#ML tree search trimmed
~/software/raxml-ng --all \
         --msa Telos85T_4parts_partitions12.raxml.rba \
         --prefix Telos85T_4parts_partitions12 \
         --seed 2 \
         --threads 99 \
         --bs-trees 100
#BS convergence
~/software/raxml-ng --bsconverge \
         --bs-trees Telos85T_4parts_partitions12.raxml.bootstraps \
         --prefix Telos85T_4parts_partitions12_convergence_test \
         --seed 2
