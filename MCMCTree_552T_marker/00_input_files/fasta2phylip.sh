#!/bin/bash

sed -i 's/DNA, //' Telos553T_muscle_partitions_v2.txt 
~/software/AMAS/amas/AMAS.py split -l Telos553T_muscle_partitions_v2.txt -u phylip -f fasta -d dna -i Telos553T_muscle_concatenated_v2.fa
