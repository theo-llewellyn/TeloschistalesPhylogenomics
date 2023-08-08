      seqfile = /rds/general/project/theollewellynproject/live/DNA_Repair_Analysis/Leca117T/Hsp90_Leca117T_muscle5_msa_codon_trimmed_renamed.fa            * Path to the alignment file
     treefile = /rds/general/project/theollewellynproject/live/DNA_Repair_Analysis/Leca117T/SpeciesTree_rooted_renamed_Hsp90.tre           * Path to the tree file
      outfile = /rds/general/project/theollewellynproject/live/DNA_Repair_Analysis/PAML/Leca117T/Hsp90/out_branchsite.txt            * Path to the output file
   
        noisy = 3              * How much rubbish on the screen
      verbose = 1              * More or less detailed report

      seqtype = 1              * Data type
        ndata = 1           * Number of data sets or loci
        icode = 0              * Genetic code 
    cleandata = 0              * Remove sites with ambiguity data?
		
        model = 2         * Models for ω varying across lineages
	  NSsites = 2          * Models for ω varying across sites
    CodonFreq = 7        * Codon frequencies
	  estFreq = 0        * Use observed freqs or estimate freqs by ML
        clock = 0          * Clock model
    fix_omega = 0         * Estimate or fix omega
        omega = 0.5        * Initial or fixed omega
