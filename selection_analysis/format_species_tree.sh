#!/bin/bash
#remove branch lengths and support from iqtree species tree and reroot

sed 's/:[0-9]*\.[0-9]*//g' concat_OF_tAl.treefile | sed 's/100//g;s/96,/,/g;s/71)/)/g;s/78)/)/g' > tree.tre
pxrr -t tree.tre -g LIQ171ARSP,LIQ211ARSP,LIQ180TREL,LIQ214TRSU > Leca117T_concat_OF_tAl_rooted.treefile
