#!/bin/bash

#concatenate all gene trees but removing the newick header with number of taxa and trees
awk FNR!=1 filtered_genes_step2/*/*.tree > genetrees.tree 

#keep only gene trees with the Caliciales outgroup
grep -E '(LIQ199HEOB|LIQ203DILE) > genetrees_root.tree

#count number genes
NGENES=$(wc -l genetrees_root.tree | awk '{print $1}')

for i in {1..$NGENES}
do
#save the taxa names of the outgroup in each gene
grep -Eo '(LIQ199HEOB.*?T1|LIQ203DILE.*?T1)' <(head -n $i genetrees_root.tree | tail -n 1) > outgroups_${i}.txt
#reroot the tree using that outgroup file
pxrr -f outgroups_${i}.txt -t <(head -n $i genetrees_root.tree | tail -n 1) > rooted_tree_${i}.tree
done

rm outgroups_*.txt
cat rooted_tree*.tree > rooted_trees.tree
