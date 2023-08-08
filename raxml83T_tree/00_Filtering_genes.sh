#!/bin/bash

# This script was written by Sandra Alvarez-Carretero and adapted for 
# the Teloschistales dataset by Theo Llewellyn

# Start counter
count=0

# Move to main dir, where this script saved in `00_Gene_filtering/scripts`
# is run, i.e., `00_Gene_filtering`
curr_dir=$( pwd )
cd $curr_dir
cd ../../src
src=$( pwd )
cd $curr_dir

mkdir out_data
# Create log file
printf "Gene,Num_seqs,TCHR_presence,XANPA_presence,Num_codons_hum,Length_seq\n" > out_data/genes.csv

# Loop over genes (15,904 genes) and filter
for i in *_Telos
do

# Start general counter
count=$(( count + 1 ))

# If the cds.aln.fasta exists...
if [ -f $i/*muscle5_msa_nucl.fa ]
then

# Get gene name
gene=$( echo $i | sed 's/..*\///')
# Check num sequences
num_seq=$( grep '>' $i/*muscle5_msa_nucl.fa | wc -l )
# Check Teloschistes chrysopthalmus is there 
TCHR=$( grep '>LIQ125TCHR2' $i/*muscle5_msa_nucl.fa | wc -l )
# Check Xanthoria parietina is there 
XANPA=$( grep '>Xanthoria_parietina' $i/*muscle5_msa_nucl.fa | wc -l )
# Check length of codons (1 codon = 3 bps)
# The perl scripts gets the fasta sequences in one line 
# so it is easier to count the amount of nucleotides,
# then this file is removed.
src/one_line_fasta.pl $i/*muscle5_msa_nucl.fa
num_nucs=$( wc -L $i/*one_line.fa | awk '{print $1}' )
num_codons=$( grep -A1 'LIQ125TCHR2' $i/*one_line.fa | sed 's/\-//g' | sed -n '2,2p' | wc -L | awk '{print $1/3}' )
rm $i/*_one_line.fa

else 
gene=$( echo NULL )
num_seq=$( echo NULL )
TCHR=$( echo NULL )
XANPA=$( echo NULL )
num_nucs=$( echo NULL )

fi 

echo There are $count genes visited to generate summary statistic
printf $gene","$num_seq","$TCHR","$XANPA","$num_codons","$num_nucs"\n" >> out_data/genes.csv
printf $gene","$num_seq","$TCHR","$XANPA","$num_codons","$num_nucs"\n"

done
