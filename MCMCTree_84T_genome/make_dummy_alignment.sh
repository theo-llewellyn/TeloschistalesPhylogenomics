#!/bin/bash

#make dummy alignment
awk '{print $1}' ~/LICHENS_Hessian/alignments/Telos_concat_part1.aln | sed -r '/^\s*$/d' | sed "s/$/  AT/" | sed 's/84  AT/84 2/' > dummy_aln.aln

#make the directories for the 5 independent runs, copy the control file, alignment and tree with node priors made in step 2
#run within 00_MCMCtree_prior directory
for i in `seq 1 5`
do
mkdir mcmc$i
cp mcmctree_GBM_calibs_LICHENS.ctl mcmc$i
cp Telos85T_4parts_partitions12.raxml.support_rooted_labelled_MCMCtree_calib.tree mcmc$i
cp dummy_aln.aln mcmc$i
done
