#!/bin/bash

#get a list of the filtered genes 
ls ~/Telos85T_10sp_msa/000_alignments/part* | grep -o -P '(?<=_)OG.*Telos(?=.fasta)' > OGs_Telos85T_4parts_partitions12_1254genes.txt

#keep only accession number
cd rds/general/project/theollewellynproject/live/raxml-ng/Telos85T_10sp_genetrees/

#for loop goes through each gene, edits the taxon names and saves in a new file
LINES=$(cat ~/OGs_Telos85T_4parts_partitions12_1254genes.txt)
for LINE in $LINES
 do
  #get the fasta header and any sequence until the next fasta then remove the final line and change header to just species name. This is saved to a separate file for each genome
  sed 's/_prot[0-9]*[^:]//g' ${LINE}/*bestTree | sed 's/_[0-9]*-T1//g' > ${LINE}/${LINE}.raxml_edited.bestTree
done

#make astral directory
mkdir /rds/general/project/theollewellynproject/live/ASTRAL/Telos85T_10sp_1254genes

#concatenate all into a single multi-tree file and then remove the redundant files
cat */*_edited.bestTree > /rds/general/project/theollewellynproject/live/ASTRAL/Telos85T_10sp/TTelos85T_10sp_1254genes.muscle.raxml.bestTree.allgenetrees

cd /rds/general/project/theollewellynproject/live/ASTRAL/Telos85T_10sp_1254genes

java -jar ~/software/ASTRAL/Astral/astral.5.7.8.jar \
 -i Telos85T_10sp_1254genes.muscle.raxml.bestTree.allgenetrees \
 -t 3 \
 -o Telos85T_10sp_1254genes.muscle.ASTRAL.tre \
 2>ASTRAL_out.log
