#!/bin/bash

sed -i 's/DNA, //' Telos553T_muscle_partitions_v2.txt 
~/software/AMAS/amas/AMAS.py split -l Telos553T_muscle_partitions_v2.txt -u phylip -f fasta -d dna -i Telos553T_muscle_concatenated_v2.fa

cat Telos553T_muscle_concatenated_p1_Telos553T_ITS_muscle5_msa-out.phy Telos553T_muscle_concatenated_p2_Telos553T_LSU_muscle5_msa-out.phy Telos553T_muscle_concatenated_p3_Telos553T_MCM7_muscle5_msa-out.phy Telos553T_muscle_concatenated_p4_Telos553T_RET1_muscle5_msa-out.phy Telos553T_muscle_concatenated_p5_Telos553T_RPB1_muscle5_msa-out.phy Telos553T_muscle_concatenated_p6_Telos553T_RPB2_muscle5_msa-out.phy Telos553T_muscle_concatenated_p7_Telos553T_SSU_muscle5_msa-out.phy > Telos553T_7parts.aln
