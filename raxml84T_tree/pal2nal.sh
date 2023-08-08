#PBS -l walltime=00:20:00
#PBS -l select=1:ncpus=1:mem=1gb
#PBS -J 1-1524

#get name of OG file
OG=$(head -n $PBS_ARRAY_INDEX /rds/general/user/tbl19/home/OGs_Telos85T_10sp.txt|tail -n 1)
PREFIX=${OG%%.fa}

mkdir /rds/general/user/tbl19/home/Telos85T_10sp_msa/$PREFIX
cd /rds/general/user/tbl19/home/Telos85T_10sp_msa/$PREFIX

/rds/general/user/tbl19/home/software/pal2nal.v14.pal2nal.pl \
 /rds/general/user/tbl19/home/Telos85T_10sp_msa/${PREFIX}_muscle5_msa.fa \
 /rds/general/project/theollewellynproject/ephemeral/Single_Copy_Orthologues_10sp_CDS/${PREFIX}_CDS.fa \
 -output fasta
