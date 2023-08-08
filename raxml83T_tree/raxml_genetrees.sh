#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=8:mem=50gb
#PBS -J 1-1524

#get name of OG file
OG=$(head -n $PBS_ARRAY_INDEX /rds/general/user/tbl19/home/OGs_Telos85T_10sp.txt|tail -n 1)
PREFIX=${OG%%.fa}

#make a directory for that gene
mkdir /rds/general/project/theollewellynproject/live/raxml-ng/Telos85T_10sp_genetrees/${PREFIX}
cd /rds/general/project/theollewellynproject/live/raxml-ng/Telos85T_10sp_genetrees/${PREFIX}

#run raxml tree search for that gene using best model, untrimmed
/rds/general/user/tbl19/home/software/raxml-ng --search \
	--msa /rds/general/user/tbl19/home/Telos85T_10sp_msa/${PREFIX}/${PREFIX}_muscle5_msa_nucl.fa \
	--prefix ${PREFIX} \
	--model /rds/general/user/tbl19/home/Telos85T_10sp_msa/${PREFIX}/partition.txt \
  --seed 2 \
  --threads 8
