#!/bin/bash
# Edit protein headers so they just contain the accession code

sed '/^>/s/_prot.*//' Hsp90_Leca118T_muscle5_msa_codon_trimmed.fa | sed '/^>.*/s/_[0-9].*-T1//' > Hsp90_Leca118T_muscle5_msa_codon_trimmed_renamed.fa
