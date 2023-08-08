#!/bin/bash
#This script was written by Sandra Alvarez-Carretero and adapted for the Teloschistales dataset by Theo Llewellyn

# 1. Start counter
count=0

# 2. Move to main dir, where this script saved in `00_Gene_filtering/scripts`
#    is run, i.e., `00_Gene_filtering`
curr_dir=$( pwd )
cd $curr_dir
cd src
src=$( pwd )
cd $curr_dir

# 3. Loop over filtered genes and prepare them for baseml 
#    analyses
mkdir baseml
mkdir out_logs
for i in *_Telos
do

	# 3.0. Start general counter & create dir
	count=$(( count + 1 ))
	mkdir -p baseml/$count

	# 3.1 Get gene name  
	gene=$( echo $i | sed 's/..*\///')
	# 3.2 Check num sequences
	num_seq=$( grep '>' $i/*muscle5_msa_nucl.fa | wc -l )

	echo Parsing gene $count":" $gene 
	echo Parsing gene $count":" $gene >> out_logs/log_02_baseml_getdataformat.txt

	# 3.3 Get one line sequences
	$src/one_line_fasta.pl $i/*muscle5_msa_nucl.fa
	mv $i/*one_line.fa baseml/$count
	# 3.4. Get sequence next to header 
	$src/00_get_seq_next_to_header.pl baseml/$count/*one_line*
	# 3.5. Get alignment in PHYLIP format 
	$src/01_concatenate_genes.pl baseml/$count/*_tab.aln
	# 3.6. Remove unnecessary files 
	rm baseml/$count/*one_line*

	# 3.7. Get tree in PHYLIP format
	printf $num_seq" 1\n" > baseml/$count/$gene".tree"
	cat /rds/general/project/theollewellynproject/live/raxml-ng/Telos85T_10sp_genetrees/$i/*bestTree >> baseml/$count/$gene".tree"

done

echo There are $count genes visited to generate summary statistic
echo There are $count genes visited to generate summary statistic >> out_logs/log_02_baseml_getdataformat.txt
