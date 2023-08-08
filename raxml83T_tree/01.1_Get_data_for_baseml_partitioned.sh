#!/bin/bash

#This script was written by Sandra Alvarez-Carretero and adapted for the Teloschistales dataset by Theo Llewellyn

# 1. Get path to main dir,`00_Gene_filtering/scripts`,
#    and then get the path to `src`. Move to baseml
curr_dir=$( pwd )
cd $curr_dir
cd src
src=$( pwd )
cd $curr_dir/baseml

# 2. Run `src/partition_alignment` for each of the genes 
#    in the `baseml` directory
for i in *
do
cd $i
printf "Parsing "$i" ...\n\n"
# The usage of this perl script is the following:
#   <path_to_script>/partition_alignments.pl <alignment_file> <lines_to_skip> <separator>
$src/partition_alignments.pl OG*[0-9].aln 2 "\s{6}"
 for filename in partitions12_*
 do
 cp $filename ${filename%.*}_Telos.aln
 done
cd ..
done
