#PBS -l walltime=00:20:00
#PBS -l select=1:ncpus=1:mem=1gb
#PBS -J 1-1524

#get name of OG file
OG=$(head -n $PBS_ARRAY_INDEX /rds/general/user/tbl19/home/OGs_Telos85T_10sp.txt|tail -n 1)
PREFIX=${OG%%.fa}

cd /rds/general/user/tbl19/home/Telos85T_10sp_msa/$PREFIX

#get alignment length
~/software/AMAS/amas/AMAS.py summary -f fasta -d dna -i ${PREFIX}_muscle5_nucl.fa
length=$(cut -f 3 summary.txt | tail -n 1)

#make partition file
echo pos12 = 1-$length'\'3, 2-$length'\'3 > partition.txt
echo pos3 = 3-$length'\'3 >> partition.txt
#Add gene models

~/software/AMAS/amas/AMAS.py split -f fasta -d dna -i ${PREFIX}_muscle5_nucl.fa -l partition.txt
sed -i 's/^/GTR+G, /' partition.txt
