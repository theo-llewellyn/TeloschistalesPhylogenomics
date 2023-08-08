#!/bin/bash

#This script was written by Sandra Alvarez-Carretero and adapted for the Teloschistales dataset by Theo Llewellyn

# Go through each alignment, extract mouse and human sequences 
# and create file to be read in R  
for i in baseml/*
do

# 1. Get gene name 
gene_name=$( ls $i/*tree | sed 's/..*\///' | sed 's/\..*//' )
echo Parsing $gene_name
echo Parsing $gene_name >> out_logs/log_03_get_XanPa_TelChr.txt

#> # 2. Create the mouse_human file to be later read in R 
grep 'LIQ125TCHR2' $i/partitions12_*aln > $i/$gene_name"_TCHR_XANPA.aln"
grep 'Xanthoria_parietina' $i/partitions12_*.aln >> $i/$gene_name"_TCHR_XANPA.aln"
sed -i 's/      /\t/g' $i/$gene_name"_TCHR_XANPA.aln"

done
