## 5. Multilocus tree 552T 84T backbone

The following code allows you to pull 7 marker genes from a wole genome dataset, combine with marker gene datasets and then make a phylogeny constrained to whole genome tree. The steps to obtain concatenated 7-locus phylogeny and individual gene trees for the 553T dataset are as follows:

1. `qsub pull_marker_genes.sh` pull genes with gene pull
2. `qsub ITSx.sh` pull ITS
3. Blast them on ncbi to ID and check top hit is Teloschistales
4. concatenate markers pulled from genomes with bigger marker gene dataset with a simple cat command
5. `qsub muscle5.sh` align marker gene multifastas with muscle5
6. `qsub trimAL.sh` trim
7. Change terminal gaps (-) to missing (?) so raxml doesnt interpret them as indels. Do this manually in Aliview or other sequence viewers
8. `qsub amas.sh` concatenate both trimmed and untrimmed
9. `qsub modeltest-ng.sh` choose the best substitution model using modeltest-ng and the AIC resuklts
10. `rename_taxa.sh` edit names of whole genome trees so matches marker genes
11. `qsub IQTree_Telos552T_genomeconstrain.sh` run ML IQTRee of 552T marker gene dataset with 84T genome tree as a constraint. This is dont for both trimmed and untrimmed alignments
12. `qsub raxml552T_genetrees.sh` make gene trees for each marker
13. Use gene trees to look for unusual placement or branch lengths. If any errors spotted, remove sequences and repeat from step 5.
