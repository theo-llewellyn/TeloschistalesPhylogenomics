#!/bin/bash

for i in *.support
do
grep -Eo '(LIQ199HEOB|LIQ203DILE)' $i > outgroups_${i}.txt
pxrr -f outgroups_${i}.txt -t $i > ${i}.rooted
done

rm outgroups_*.txt
cat *.rooted > rooted_trees.tree
