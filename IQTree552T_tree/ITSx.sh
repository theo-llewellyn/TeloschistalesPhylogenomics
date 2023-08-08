#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=16:mem=8gb

#get name of contigs file
ACCESSION=ACCESSION
module load hmmer/3.1

mkdir /rds/general/project/theollewellynproject/live/ITSx/${ACCESSION}
cd /rds/general/project/theollewellynproject/live/ITSx/${ACCESSION}

/rds/general/user/tbl19/home/software/ITSx_1.1.3/ITSx \
 -i /rds/general/project/theollewellynproject/live/redundans/redundans_${ACCESSION}_concoct/scaffolds_filled_sorted.fa \
 -o /rds/general/project/theollewellynproject/live/ITSx/${ACCESSION}_concoct_ITSx.out \
 -t fungi \
 --cpu 16

#move files to said directory, not need if already changed to that directory
#mv /rds/general/project/theollewellynproject/live/ITSx/${ACCESSION}_* /rds/general/project/theollewellynproject/live/ITSx/${ACCESSION}

#add species accession to the fasta headers to make them identifiable
sed -i "s/scaffold/${ACCESSION}_scaffold/g" ${ACCESSION}*.fasta
