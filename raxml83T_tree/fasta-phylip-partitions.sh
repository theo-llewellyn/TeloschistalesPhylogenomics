
#This script was written by Sandra Alvarez-Carretero and adapted for the Teloschistales dataset by Theo Llewellyn

# Open a terminal within the 01_SeqBayes_S1/00_Gene_filtering/000_alignments
# directory and run the following code:

# Loop over all the "partX" directories, move there, and run the pipeline 
# within each directory to obtain individual concatenated alignments for  
# each of the partitions
PATH=/rds/general/user/tbl19/home/software/fasta-phylip-partitions/src/:$PATH
for i in `seq 1 4` 
do 

cd part$i 

# 1. Create `species_names.txt`. In order to create a file listing the 72 species 
#    for which most of the molecular data can be found, we are going to use the 
#    alignment we generated for one of the genes within the `filtered_genes_step2_all72sp` 
#    directory. Therefore, we can run the following command which will extract 
#    the first column of this file - the one with the taxa names - and then output it 
#    in the main directory from which you should run this command (`000_alignments`) 
#    in a file called `species_names.txt`, which will be later used by the pipeline:

sed -n '2,85p' ../../filtered_genes_step2_all72sp/1/partitions12_*.aln | awk '{print $1}' | sed 's/_prot.*//' | sed 's/_[0-9].*-T1//' > species_names.txt

# 2. Now, run the pipeline. If you wanted more information about how to install this pipeline 
#    and the arguments used,
#    please read the tutorial [here](https://github.com/sabifo4/fasta-phylip-partitions/blob/main/README.md)
#    NOTE: If you are running this code in a directory that you have synched to Dropbox or another 
#    cloud storage and you have issues, just move the folder out of the synched directory and run the 
#    command below again.

Run_tasks.sh . Telos partN

# 3. Move back to main directory for next iteration to proceed 
cd ..

done

# Open a terminal within the 01_SeqBayes_S1/00_Gene_filtering/000_alignments
# directory and run the following code:
wd=$( pwd )
for i in `seq 1 4`
do

cd part$i/phylip_format/02_concatenated_alignments 
cat *aln >> $wd"/alignment_4parts.aln" 
cd $wd 

done  

#or with AMAS
~/software/AMAS/amas/AMAS.py concat \
 -p Telos85T_4parts_partitions12.txt \
 -t Telos85T_4parts_partitions12.fa \
 -i part*/phylip_format/02_concatenated_alignments/*aln \
 -f phylip -d dna

sed -i 's/^/GTR+I+G4, /' Telos85T_4parts_partitions12.txt
