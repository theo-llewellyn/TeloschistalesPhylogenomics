Once we have a RAxML 85T tree we can begin the dating analyses.

Steps again follow Alvarez-Carretero 2022.

## 1. Root RAxML 85T tree using Caliciales as outgroup using pxrr
`qsub reroot_tree.sh`

## 2. Calibrate nodes 
`cd 00_calibrate_nodes`
For this we are using Gaya et al 2015 (https://www.pnas.org/doi/full/10.1073/pnas.1507072112) estimates for Teloschistales and Prieto and Wedin 2017 (https://link.springer.com/article/10.1007/s13225-016-0372-y) for the node connecting Teloschistales and Caliciales.

IMPORTANT!! Make sure to convert ages into 100Myr units, i.e. 237Mya -> 2.37 in the files as MCMCTree can estimate better values for the priors when it is in this format.

### 2.1 Label nodes on the tree with [NODENAME]
Do manually
### 2.2 Make file with calibrations in the format [NODENAME]|'B(MIN,MAX)'

### 2.3 Replace node names with calibrations
`qsub Add_node_calibrations.R`

We can then manually add the newick header with no. taxa and no. trees i.e. 84 1. Specifying the node calibrations with this format in MCMCtree will tell the programme these are soft bounds with a default probability of 0.025 that the bounds will be violated.

## 3. Estimate marginal likelihood of relaxed clock model
`cd 01_clock_model_testing`

We will use mcmc3r to test the fit of three clock models (strict clock, ILN and GBM) in order to assess which is the best for our data for the proper dating with mcmctree. To do this we follow this tutorial
https://dosreislab.github.io/2017/10/24/marginal-likelihood-mcmc3r.html

### 3.1 Generate list of 32 beta points
We generate a list of 32 beta points to estimate the marginal likelihood under the strict clock (clk), the independent log normal (ILN) and generalised brownian motion (GBM) models and make separate directories for each
`qsub Generate_betapoints.R`

### 3.2 Repeat for genes shared by all taxa
We will do this for the 4 genes shared by all taxa (can just copy the control files created in step 3.1 and use same beta points). We will make a directory with clk, iln and gbm for each gene, within which we will have the 32 subdirectories for the beta points. Therefore we have to replace the alignment name with the name of the gene alignment. We also have to remove the protein suffixes from the alignment files as they were copied directly from the muscle alignments before we renamed them.
`qsub rename_alignments.sh`

In total we have 4 genes x 32 beta points x 3 clock models = 384 MCMCTree jobs

### 3.2 Run MCMCTree using those 32 b values.
We will run 4 array jobs for each of the four genes
`qsub MCMCTree_betapoints.sh` script shown is for one gene, repeat for other four

### 3.3 Calculate marginal likelihood and standard error
We then parse the MCMCTree output in R and calculate the log marginal likelihood and standard error for each model. We repeat this for each of the four models and then calculate the Bayes factors and posterior probabilities to identify the best model. We can do this in a loop:
`qsub Calculate_marginal_likelihood.R`

This will print how many times each clock model was the best fitting model.

## 4. Calculate gradient and Hessian to use in approximate likelihood calculation
`cd 02_hessian`

### 4.1 First we copy the 85T tree (no bl and labels) alignments for the 4 partitions (slow to fast) and the concatenated alignment to a new directory.

### 4.2 Make directories for the 4 partitions and then soft link the alignments and tree to directories
`qsub prepare_baseml_directories.sh`

We then change into each of the four p0* directories and run `mcmctree *ctl` to generate the input files for baseml. When we see "X site patterns read, X sites" we can cancel the job. We remove unecessary files with
```
for i in *
do cd $i
rm rst* out.* 2base.t rub tmp*out
sed -i 's/method\ \=\ 0/method\ \=\ 1/' tmp*ctl
cd ..
done
```
###  4.3 run baseml
`qsub baseml.sh`

### 4.4 Concatenate rst2 files
```
for i in `seq 1 4`
do
cat p0$i/rst2 >> in.BV.01".p1-4"
printf "\n\n" >> in.BV.01".p1-4"
done
```


## 5. MCMCTreee
`cd 03_MCMCTree`

Here we are estimating divergence times using the approximate likelihood method. We will use the GBM clock model (autocorrelated, i.e. closer branch are more likely to have similar rates) as clades within the Teloschistales are likely to be evolving more or less the same. Its a smallish clade (~1000sp) and they are all lichens so no huge lifestyle shifts that could dramatically affect evolution rates. We will also specify a uniform node age priors by specifying birth death parameters of 1, 1 and 0.1. For the substitution rate priors we will use a shape parameter (alpha) of 2 and estimate the scaling parameter (beta). This can be done as follows:

Using the tree height of the gene trees and the divergence time of the root. We need to convert root age from My to 100Myr (i.e. 137Mya will become 1.37) as it leads to better values for rate priors in MCMCtree. We need rooted gene trees with the Caliciales outgroup. We can get this with the following code:
`qsub root_genetrees.sh`

We can then use the rooted tree file in the following script `Estimate_meanrate_beta.R` which first looks for outlier gene trees where a single branch is >60% of total tree length as we did to make the species tree. We calculate the relative branch lengths by dividing largest branch by tree length. We order them from slow to fast and we take the average evolutionary rate. If alpha = 2 we can then calculate beta by doing alpha/mean rate
B = 2/0.8840532
  = 2.262307


### 5.1 Prior
We first run MCMCTree without data (just sampling from prior) to check MCMCTree is actually sampling from the priors we provided. We will run 5 independent runs on a dummy alignment as we're not using the data anyway so it will speed up runtime.

`qsub make_dummy_alignment.sh`

Now we run mcmctree 5 times `MCMCTree_prior_array.sh` and combine results using `Combine_MCMC.sh` script to generate a single tracer file. We run MCMCtree again on the combined data to get a tree with branch lengths (needs to be run using paml4.9j otherwise the print -1 option doesnt work).

### 5.2 Posterior
Run with the data. We will run 16 independent chains using the GBM and ILN models
`qsub MCMCTree_posterior_array.sh`

### 5.3 Compare prior and posterior
Plot prior and posterior estimates to see how different chains performed and how post compares to prior in R using `00_Check_MCMCs_MCMCtreeR_prior_post.R`

### 5.4 Calculate ESS 
with script `Calculate_ESS.R`

### 5.5 Concatenate chains
Concatenate the chains into a single tracer file in the same way as with the prior
`qsub Combine_posterior_chains.sh`

Can check this in Tracer to have a look.

We also can create a summary of the chains using `mcmc_tracer.txt`, `dummy_aln.aln` and `Telos85T_4parts_partitions12.raxml.support_rooted_labelled_MCMCtree_calib_1a.tree` and running with MCMCTree and the -1 option.


## 6. Fit Skew-T distribution to mean posterior estimates.  
`cd 04_SkewT`

Run `00_Fit_skewT.R` script to estimate skew-T distrib. Can then check that the skew-T distribution are good estimates of the posterior from previous MCMCtree run

### 6.1 Add Distributions onto tree
Run `01_Add_STcalibs_to_tree.R`.
Needs three input files:
- tree file with node number labels (can be found in out.txt from posterior mcmctree
- node calibrations file used for mcmctree with an additional column showing which node numbers they refer to
- SkewT distributions from previous step

### 6.2 Run MCMCtree again using just the priors
We do this to see if MCMCtrees effective priors overlap well with the skew-T priors we specified. As in step 5.1 we will run it without any data. We will run two independent chains.

### 6.3 Compare effective priors with skewT
When this is finished we compare the effective priors estimated by mcmctree against the skewT distributions and the samples from the posterior using `02_Eval_skewT.R`. This will show us if the skewT distribution priors are effectively capturing the posterior estimates from the whole genome data which will then be used as calibration points in the big multilocus phylogeny dating

### 6.4 Compare with RelTime estimates
I also estimated dates using the rapid relative divergence tool RelTime. Both OLS and branch-length methods in MEGA. I provided the same secondary calibrations as for MCMCTree with uniform probability and hard bounds as RelTime does not allow for soft-bounds. For OLS I used the Maximum Composite Likelihood substitution model + G with a gamma parameter of 2 as RelTime OLS does not allow for HKY+G5. For the branch length method I set max rate ratio to 100. Once finished I saved as nexus, converted to newick and then compared the results in R. I compared results to MCMCTree estimates with `RelTime_MCMCTree_comparison.R`
