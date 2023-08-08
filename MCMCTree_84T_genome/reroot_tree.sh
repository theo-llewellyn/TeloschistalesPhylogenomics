#!/bin/bash

#reroot tree and remove branch lengths and support
sed 's/:[0-9]*\.[0-9]*//g' Telos85T_4parts_partitions12.raxml.support | sed 's/100//g' | sed 's/99,/,/g' > tree.tre
pxrr -t tree.tre -g LIQ203DILE,LIQ199HEOB > Telos85T_4parts_partitions12.raxml.rooted.tre
