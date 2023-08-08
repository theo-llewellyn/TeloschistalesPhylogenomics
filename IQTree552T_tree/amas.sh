#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=8:mem=240gb

cd ~/Telos552T_msa

~/software/AMAS/amas/AMAS.py concat \
 -p Telos552T_muscle_partitions.txt \
 -t Telos552T_muscle_concatenated.fa \
 -i Telos552T_*_muscle5_msa.fa \
 -f fasta -d dna -c 8

~/software/AMAS/amas/AMAS.py concat \
 -p Telos552T_muscle_tAl_partitions.txt \
 -t Telos552T_muscle_tAl_concatenated.fa \
 -i Telos552T_*_muscle5_msa_tAl.fa \
 -f fasta -d dna -c 8

#add DNA, to the beginning to be read by modeltest-ng
sed -i 's/^/DNA, /' Telos552T_muscle_partitions.txt
sed -i 's/^/DNA, /' Telos552T_muscle_tAl_partitions.txt
