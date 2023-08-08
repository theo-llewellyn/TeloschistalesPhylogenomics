#!/bin/bash

filename='macrogen_2022_genomes.txt'

while read p; do 
    echo $p
    perl -lane 'if(($F[5]=~"Ascomycota")|($F[5]=~"no-hit")){print $F[0]}' ${p}_concoct_blobtools.out/${p}_concoct_blobplot.blobDB.bestsum.table.txt > /rds/general/project/theollewellynproject/live/filtered_Asco_contigs/${p}_concoct_Asco_headers.txt
done < $filename
