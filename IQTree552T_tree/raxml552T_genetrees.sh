#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=8:mem=1gb
#PBS -J 1-7

mkdir cd /rds/general/project/theollewellynproject/live/raxml-ng/Telos552T/genetrees
cd /rds/general/project/theollewellynproject/live/raxml-ng/Telos552T/genetrees
GENE=$(head -n $PBS_ARRAY_INDEX GENES.txt|tail -n 1)
mkdir ${GENE}
cd ${GENE}


#ML tree search with bootstrapping untrimmed and RAxML genome tree as constraint
~/software/raxml-ng --all \
         --msa ~/Telos552T_msa/Telos552T_${GENE}_muscle5_msa.fa \
	 --model GTR+I+G4 \
         --prefix Telos552T_${GENE} \
         --seed 2 \
         --threads 8 \
         --bs-trees autoMRE
