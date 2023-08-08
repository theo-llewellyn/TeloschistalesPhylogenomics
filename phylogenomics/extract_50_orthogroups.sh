#!/bin/bash

mkdir Single_Copy_Orthologues_50percent

cat Orthologues_Leca117T_50p.txt | while read line; do cp Orthogroup_Sequences/${line}.fa Single_Copy_Orthologues_50percent; done

