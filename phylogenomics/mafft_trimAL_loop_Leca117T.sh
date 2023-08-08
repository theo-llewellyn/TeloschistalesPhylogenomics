#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=8:mem=96gb

module load mafft/7.271 trimal

cd /rds/general/project/theollewellynproject/live/OrthoFinder/Results_Leca117T/Single_Copy_Orthologues_50percent
mkdir /rds/general/user/tbl19/home/Leca117T_msa

for i in *.fa
do
ORTHOGROUP=${i%%.fa}
mafft --auto $i > /rds/general/user/tbl19/home/Leca117T_msa/${ORTHOGROUP}_msa.fa
done

cd /rds/general/user/tbl19/home/Leca117T_msa

for i in *.fa
do 
ORTHGROUP=${i%%.fa}
trimal -in $i -out ${ORTHGROUP}_tAl.fa -fasta -automated1
done

for i in *.fa; do sed -i 's/_[0-9]*-T1//g;s/_prot[0-9]*[^:]//g' $i; done
