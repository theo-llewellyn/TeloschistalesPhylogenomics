#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=16:mem=124gb
module load anaconda3/personal
source activate blobtools

#get name of contigs file
ACCESSION=NAME_OF_SEQUENCE

cd $PBS_O_WORKDIR

mkdir /rds/general/project/theollewellynproject/live/BlobTools/${ACCESSION}_concoct_blobtools.out

#create the blobDB with input, -t hits files for uniprot and lichens
/rds/general/user/tbl19/home/blobtools/blobtools create \
 -i /rds/general/project/theollewellynproject/live/filtered_Asco_contigs/${ACCESSION}_concoct_merged_contigs.fa \
 -c /rds/general/project/theollewellynproject/live/BlobTools/bbmap_alignments/${ACCESSION}_concoct.sam.cov \
 -t /rds/general/project/theollewellynproject/live/BlobTools/diamond_output/${ACCESSION}_concoct.vs.uniref90.diamond.taxified.out \
 -t /rds/general/project/theollewellynproject/live/BlobTools/blast_output/${ACCESSION}_concoct.vs.Lecanoromycetes_genomes.blastn.taxified.out \
 -o ${ACCESSION}_concoct_blobtools.out/${ACCESSION}_concoct_blobplot

#create a view of the blobDB, at family rank, default is phylum
/rds/general/user/tbl19/home/blobtools/blobtools view \
 -i ${ACCESSION}_concoct_blobtools.out/${ACCESSION}_concoct_blobplot.blobDB.json \
 -r phylum \
 -o ${ACCESSION}_concoct_blobtools.out/

#create the blobplot colouring by family and plotting unwanted ones first
/rds/general/user/tbl19/home/blobtools/blobtools plot \
 -i ${ACCESSION}_concoct_blobtools.out/${ACCESSION}_concoct_blobplot.blobDB.json \
 -o ${ACCESSION}_concoct_blobtools.out/ \
 -r phylum \
 -m \
 --sort_first no-hit,other 

conda deactivate
