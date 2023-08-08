#!/bin/bash

#make dummy alignment
awk '{print $1}' Telos553T_muscle_concatenated_p1_Telos553T_ITS_muscle5_msa-out.phy | sed -r '/^\s*$/d' | sed "s/$/  AT/" | sed 's/553  AT/553 2/' > dummy_aln.aln

#make the directories for the 5 independent runs, copy the control file, alignment and tree with node priors made in step 2
#run within 00_MCMCtree_prior directory
for i in `seq 1 5`
do
mkdir mcmc$i
cp mcmctree_GBM_calibs_LICHENS.ctl mcmc$i
cp Telos553T_MCMCtree_calib.tree mcmc$i
cp dummy_aln.aln mcmc$i
done
