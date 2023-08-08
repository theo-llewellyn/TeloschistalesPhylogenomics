Once we have a RAxML 85T tree we can begin the dating analyses.

Steps again follow Alvarez-Carretero 2022.

## 1. Root RAxML 85T tree using Caliciales as outgroup using pxrr

```
#remove branch lengths and support
 sed 's/:[0-9]*\.[0-9]*//g' Telos85T_4parts_partitions12.raxml.support | sed 's/100//g' | sed 's/99,/,/g' > tree.tre
pxrr -t tree.tre -g LIQ203DILE,LIQ199HEOB > Telos85T_4parts_partitions12.raxml.rooted.tre
```
## 2. Calibrate nodes 
For this we are using Gaya et al 2015 (https://www.pnas.org/doi/full/10.1073/pnas.1507072112) estimates for Teloschistales and Prieto and Wedin 2017 (https://link.springer.com/article/10.1007/s13225-016-0372-y) for the node connecting Teloschistales and Caliciales.

IMPORTANT!! Make sure to convert ages into 100Myr units, i.e. 237Mya -> 2.37 in the files as MCMCTree can estimate better values for the priors when it is in this format.

### 2.1 Label nodes on the tree with [NODENAME]
Do manually
### 2.2 Make file with calibrations in the format [NODENAME]|'B(MIN,MAX)'

### 2.3 Replace node names with calibrations
We can use the R script `Add_node_calibrations.R`

```
setwd("~/The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/PHYLOGENOMICS/TIME_CALIBRATION/")

#read in calibrations file
calibrations <- read.table( file = "Calibrations_lichens_maintree.txt",
                            stringsAsFactors = F, sep = "|", blank.lines.skip = T )
colnames( calibrations ) <- c( "name", "MCMCtree calib" )

#read in tree
tree.file <- list.files(pattern = "Telos85T_4parts_partitions12.raxml.support_rooted_labelled_maintree.tre" )
trees.calibs <- readLines(tree.file)

#replace node names with calibrations
for( j in 1:dim( calibrations )[1] ){
  trees.calibs <- gsub( pattern = paste0("\\[",calibrations[j,1],"\\]"),
                        x = trees.calibs,
                        replacement = paste( "'", calibrations[j,2], "'", sep = "" ) )
}

#save as new file
write(x = trees.calibs, file = "Telos85T_4parts_partitions12.raxml.support_rooted_labelled_maintree_MCMCtree_calib.tree")
```
We can then manually add the newick header with no. taxa and no. trees i.e. 84 1. Specifying the node calibrations with this format in MCMCtree will tell the programme these are soft bounds with a default probability of 0.025 that the bounds will be violated.

## 3. Estimate marginal likelihood of relaxed clock model

We will use mcmc3r to test the fit of three clock models (strict clock, ILN and GBM) in order to assess which is the best for our data for the proper dating with mcmctree. To do this we follow this tutorial
https://dosreislab.github.io/2017/10/24/marginal-likelihood-mcmc3r.html

### 3.1 Generate list of 32 beta points
We generate a list of 32 beta points to estimate the marginal likelihood under the strict clock (clk), the independent log normal (ILN) and generalised brownian motion (GBM) models and make separate directories for each
```
##R-script
b = mcmc3r::make.beta(n=8, a=5, method="step-stones")
mcmc3r::make.bfctlf(b, ctlf="mcmctree.ctl", betaf="beta.txt")

dir.create("../../iln")
setwd("../../iln")
#then copy the clk .ctl file to iln and change clock model to 2
mcmc3r::make.bfctlf(b, ctlf="mcmctree.ctl", betaf="beta.txt")
dir.create("../gbm")
setwd("../gbm")
#then copy the clk .ctl file to iln and change clock model to 3
mcmc3r::make.bfctlf(b, ctlf="mcmctree.ctl", betaf="beta.txt")

setwd("../clk"); clk <- mcmc3r::stepping.stones()
setwd("../iln"); iln <- mcmc3r::stepping.stones()
setwd("../gbm"); gbm <- mcmc3r::stepping.stones()

models <- c("CLK","ILN","GBM")
logML <- c(clk$logml,iln$logml,gbm$logml)
SE <- c(clk$se,iln$se,gbm$se)
data.frame(models,logML,SE)
```
### 3.2 Repeat for genes shared by all taxa
We will do this for the 4 genes shared by all taxa (can just copy the control files created in step 3.1 and use same beta points). We will make a directory with clk, iln and gbm for each gene, within which we will have the 32 subdirectories for the beta points. Therefore we have to replace the alignment name with the name of the gene alignment. 
```
cd OG0003067
sed -i 's/\*\.aln/partitions12_OG0003067_Telos.aln/' */*ctl
cd ../OG0001822
sed -i 's/\*\.aln/partitions12_OG0001822_Telos.aln/' */*ctl
cd ../OG0002170
sed -i 's/\*\.aln/partitions12_OG0002170_Telos.aln/' */*ctl
cd ../OG0002971
sed -i 's/\*\.aln/partitions12_OG0002971_Telos.aln/' */*ctl
```
We also have to remove the protein suffixes from the alignment files as they were copied directly from the muscle alignments before we renamed them.
```
sed -i -e 's/_prot[0-9]*[^:]//g' */partitions12_*_Telos.aln
sed -i -e 's/_[0-9].*-T1//' */partitions12_*_Telos.aln
```

In total we have 4 genes x 32 beta points x 3 clock models = 384 MCMCTree jobs

### 3.2 Run MCMCTree using those 32 b values.
We will run 4 array jobs for each of the four genes
```
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=100gb
#PBS -J 1-32

cd /rds/general/user/tbl19/home/MCMC3r/clk/OG0001822/$PBS_ARRAY_INDEX

PATH=/rds/general/user/tbl19/home/software/paml-4.10.5/bin:$PATH

mcmctree > /dev/null

```
### 3.3 Calculate marginal likelihood and standard error
We then parse the MCMCTree output in R and calculate the log marginal likelihood and standard error for each model. We repeat this for each of the four models and then calculate the Bayes factors and posterior probabilities to identify the best model. We can do this in a loop:
```
setwd("~/The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/PHYLOGENOMICS/TIME_CALIBRATION/MCMCTREE/MCMC3r/")
#set up blank dataframe to store results
BF_df <- data.frame(gene_name=character(),models=character(),logML=numeric(),SE=numeric(), bf=numeric(),logbf=numeric(),pr=numeric(),prior=numeric(),pr.ci.2.5.=numeric(),pr.ci.97.5.=numeric())

#set up counters to count number of times each model was the best model
CLK_count <- 0
ILN_count <- 0
GBM_count <- 0

#for each gene tested
for(gene in list.files(pattern="^OG*")){
  #change to that directory
  setwd(gene)
  #pull the log-marginal likelihoods for each model
  setwd("clk"); clk <- mcmc3r::stepping.stones()
  setwd("../iln"); iln <- mcmc3r::stepping.stones()
  setwd("../gbm"); gbm <- mcmc3r::stepping.stones()
  #names of model
  models <- c("CLK","ILN","GBM")
  #name of gene
  gene_name <- rep(gene,3)
  #gets log marginal likelihoods and standard errors for each
  logML <- c(clk$logml,iln$logml,gbm$logml)
  SE <- c(clk$se,iln$se,gbm$se)
  #adds logML and SE to Bayes factors, log BF, Posterior probabilities, priors and confidence intervals for Posterior probs.
  tmp <- data.frame(gene_name,models,logML,SE, mcmc3r::bayes.factors(clk, iln, gbm))
  
  #order probabilities low to high and get the name of the gene with the highest
  best_model <- tmp[order(tmp$pr)]$models[3]
  #if a particular model is the best it will add 1 to the count for that model
  if(best_model == "CLK"){
    CLK_count <- CLK_count + 1
  }
  if(best_model == "ILN"){
    ILN_count <- ILN_count + 1
  }
  if(best_model == "GBM"){
    GBM_count <- GBM_count + 1
  }
  #adds the results from that gene to dataframe
  BF_df <- rbind(BF_df,tmp)
  setwd("../..")
}

paste("CLK is best",CLK_count,"times")
paste("ILN is best",ILN_count,"times")
paste("GBM is best",GBM_count,"times")
```
This will print how many times each clock model was the best fitting model. If this works we can then repeat the model testing for all the genes.

### 3.4 Repeat for all genes
To repeat model testing for all genes we have to use subtrees to make sure the taxa in the tree match those in the alignment. To do that we can use the drop.tip function of the ape package.


First we prepare the directories
```
cp -r ~/Telos85T_10sp_msa/filtered_genes_step2 ~/MCMC3r
cd ~/MCMC3r/filtered_genes_step2

#rename taxon headers
for i in *
do
cd $i
sed -e 's/_prot[0-9]*[^:]//g' partitions12_*_Telos.aln > partitions12_Telos.aln
sed -i -e 's/_[0-9].*-T1//' partitions12_Telos.aln
cd ..
done
```
We then make subtrees for each gene using drop.tip()
```
library(ape)

#read in full species tree
lichen_tree <- read.tree("Telos85T_4parts_partitions12.raxml.rooted.tre")

for(gene in list.files()){
  #change to that directory
  setwd(gene)
  #read in alignment
  alignment <- read.dna(dir(pattern='partitions12_Telos.aln')[1],format = "sequential")
  #pull taxa names
  taxa <- rownames(alignment)
  #keep only tips in alignment
  subsetted_tree <- keep.tip(lichen_tree, taxa)
  #save subsetted tree
  write.tree(subsetted_tree,"Telos85T_4parts_partitions12.raxml.rooted.subset.tre")
  setwd("../")
}
```
We need to reformat the root node label as drop.tip changes the brackets to dashes and also add the number of taxa and trees to make it a newick tree
```
for i in *
do
cd $i
#reformat root age label
sed -i 's/B-0.999-1.001-/B(0.999,1.001)/' *.tre
#add root age label to taxa that dont have the caliciales outgroup
sed -i "s/;/\'B(0.999,1.001)\';/" *.tre
#add newick header to tree
cat <(head -n 1 *.tree) Telos85T_4parts_partitions12.raxml.rooted.subset.tre > Telos85T_4parts_partitions12.raxml.rooted.subset.formatted.tre
rm Telos85T_4parts_partitions12.raxml.rooted.subset.tre
cd ..
done
```
We can then use the same template directory structure from the 4 genes for all genes but changing the alignment and tree names in the control files. We can reformat .ctl as follows:
```
sed -i 's/partitions12_.*/partitions12_Telos.aln/' */*/mcmctree.ctl
sed -i 's/partitions12_.*/partitions12_Telos.aln/' */mcmctree.ctl
sed -i 's/\.\.\/Telos85T_.*/Telos85T_4parts_partitions12.raxml.rooted.subset.formatted.tre/' */*/mcmctree.ctl
sed -i 's/\.\.\/Telos85T_.*/Telos85T_4parts_partitions12.raxml.rooted.subset.formatted.tre/' */mcmctree.ctl 
```
Now we can run MCMCTree again on all genes, so 1254 genes x 32 beta points x 3 clock models = 120,384 MCMCTree jobs. To allow this to run on the HPC we will split it into 32 arrays for each beta point, each consisting of 1254 jobs and looping through the three models

```
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=10gb
#PBS -J 1-1254


PATH=/rds/general/user/tbl19/home/software/paml-4.10.5/bin:$PATH

LINES=$(cat /rds/general/user/tbl19/home/MCMC3r/all_genes_scripts/clockmodels.txt)

#for each of the three models
for LINE in $LINES
do
    cd /rds/general/user/tbl19/home/MCMC3r/filtered_genes_step2/${PBS_ARRAY_INDEX}/${LINE}/1
    mcmctree > /dev/null
done

```
We can then copy this script for the other 31 beta points
```
for i in {2..32}
do
cp MCMCTree_clockmodeltest_array_1.sh MCMCTree_clockmodeltest_array_${i}.sh
sed -i "s/\/1/\/${i}/" MCMCTree_clockmodeltest_array_${i}.sh
done
```
We then repeat step 3.3 for all genes to see which model is best across all genes.



## 4. Calculate gradient and Hessian to use in approximate likelihood calculation

### 4.1 First we copy the 85T tree (no bl and labels) alignments for the 4 partitions (slow to fast) and the concatenated alignment to a new directory.

### 4.2 Make directories for the 4 partitions and then soft link the alignments and tree to directories
```
cd ~/LICHENS_Hessian/
wd=$( pwd )

for i in `seq 1 4`
do ln -s $wd/alignments/Telos_concat_part$i".aln" p0$i
ln -s $wd/trees/*.tre p0$i
cp $wd/control_file/mcmctree_LICHENS.ctl baseml_method1/p0$i
cd baseml_method1/p0$i
tree_name=$( ls *tre )
aln_name=$( ls *aln )
sed -i -e 's/ALN/'${aln_name}'/' *ctl 
sed -i -e 's/TREE/'${tree_name}'/' *ctl
cd $wd; done
```
We then change into each of the four p0* directories and run `mcmctree *ctl` to generate the input files for baseml. When we see "X site patterns read, X sites" we can cancel the job. We remove unecessary files with
```
for i in *
do cd $i
rm rst* out.* 2base.t rub tmp*out
sed -i 's/method\ \=\ 0/method\ \=\ 1/' tmp*ctl
cd ..
done
```
4.3 run baseml
```
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=5gb
#PBS -J 1-4

cd /rds/general/user/tbl19/home/LICHENS_Hessian/baseml_method1/p0${PBS_ARRAY_INDEX}

PATH=/rds/general/user/tbl19/home/software/paml-4.10.5/bin:$PATH

baseml tmp*ctl
```
### 4.4 Concatenate rst2 files
```
for i in `seq 1 4`
do
cat p0$i/rst2 >> in.BV.01".p1-4"
printf "\n\n" >> in.BV.01".p1-4"
done
```


## 5. MCMCTreee

Here we are estimating divergence times using the approximate likelihood method. We will use the GBM clock model (autocorrelated, i.e. closer branch are more likely to have similar rates) as clades within the Teloschistales are likely to be evolving more or less the same. Its a smallish clade (~1000sp) and they are all lichens so no huge lifestyle shifts that could dramatically affect evolution rates. We will also specify a uniform node age priors by specifying birth death parameters of 1, 1 and 0.1. For the substitution rate priors we will use a shape parameter (alpha) of 2 and estimate the scaling parameter (beta). This can be done two ways:

Method 1: using the equation: (alpha x root time)/mean tip to root distance. Mean tip to root can be estimated in R:
```
newtree <- read.tree("/Users/tbl19/Library/CloudStorage/OneDrive-SharedLibraries-TheRoyalBotanicGardens,Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/PHYLOGENOMICS/TREE_INFERENCE/RAxML_Telos85T_4parts/Telos85T_4parts_partitions12.raxml.support_renamed.rooted.tre")
mean(distRoot(newtree))
```
This give 0.50 therefore beta = (2x171)/0.5
beta = 684

Method 2: Using the tree height of the gene trees and the divergence time of the root. We need to convert root age from My to 100Myr (i.e. 137Mya will become 1.37) as it leads to better values for rate priors in MCMCtree. We need rooted gene trees with the Caliciales outgroup. We can get this with the following code:
```
#concatenate all gene trees but removing the newick header with number of taxa and trees
awk FNR!=1 filtered_genes_step2/*/*.tree > genetrees.tree 

#keep only gene trees with the Caliciales outgroup
grep -E '(LIQ199HEOB|LIQ203DILE) > genetrees_root.tree

#count number genes
NGENES=$(wc -l genetrees_root.tree | awk '{print $1}')

for i in {1..$NGENES}
do
#save the taxa names of the outgroup in each gene
grep -Eo '(LIQ199HEOB.*?T1|LIQ203DILE.*?T1)' <(head -n $i genetrees_root.tree | tail -n 1) > outgroups_${i}.txt
#reroot the tree using that outgroup file
pxrr -f outgroups_${i}.txt -t <(head -n $i genetrees_root.tree | tail -n 1) > rooted_tree_${i}.tree
done

rm outgroups_*.txt
cat rooted_tree*.tree > rooted_trees.tree
```
We can then use the rooted tree file in the following script `Estimate_meanrate_beta.R` which first looks for outlier gene trees where a single branch is >60% of total tree length as we did to make the species tree. We calculate the relative branch lengths by dividing largest branch by tree length. We order them from slow to fast and we take the average evolutionary rate. If alpha = 2 we can then calculate beta by doing alpha/mean rate
B = 2/0.8840532
  = 2.262307


### 5.1 Prior
run without data (just sampling from prior) to check MCMCTree is actually sampling from the priors we provided. We will run 5 independent runs on a dummy alignment as we're not using the data anyway so it will speed up runtime.

```
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
```
run mcmctree 5 times and combine results using Sandras Combine_MCMC.sh script to generate a single tracer file. We run MCMCtree again on the combined data to get a tree with branch lengths (needs to be run using paml4.9j otherwise the print -1 option doesnt work).

### 5.2 Posterior
Run with the data. We will run 16 independent chains using the GBM and ILN models
```
for i in `seq 1 16`
do
mkdir run$i
mkdir run$i/mcmctree_GBM
mkdir run$i/mcmctree_ILN
cp mcmctree_GBM.ctl run$i/mcmctree_GBM
cp mcmctree_ILN.ctl run$i/mcmctree_ILN
cp *tree run$i/mcmctree_GBM
cp *tree run$i/mcmctree_ILN
cp in.BV* run$i/mcmctree_GBM/in.BV
cp in.BV* run$i/mcmctree_ILN/in.BV
done
```
### 5.3 Compare prior and posterior
Plot prior and posterior estimates to see how different chains performed and how post compares to prior in R using 00_Check_MCMCs_MCMCtreeR_prior_post.R

### 5.4 Calculate ESS 
with script Calculate_ESS.R

### 5.5 Concatenate chains
Concatenate the chains into a single tracer file in the same way as with the prior
```
for i in {1..16}
 do
 printf "Parsing run"$i"/mcmctree_GBM ... ... \n"
 end=$( wc -l run$i/mcmctree_GBM/mcmc.txt | sed 's/ ..*//' )
 if [[ $i -eq 1 ]]
  then begin=1
  else begin=2
 fi
 sed -n ''${begin}','${end}'p' run$i/mcmctree_GBM/mcmc.txt >> mcmc_files/mcmc_tracer.txt
 sed -n '1,'${end}'p' run$i/mcmctree_GBM/mcmc.txt >>  run$i/mcmctree_GBM/mcmc_clean.txt
done 
```
Can check this is Tracer to have a look.

We also can create a summary of the chains using `mcmc_tracer.txt`, `dummy_aln.aln' and `Telos85T_4parts_partitions12.raxml.support_rooted_labelled_MCMCtree_calib_1a.tree` and running with MCMCTree and the -1 option.

## 6. Fit Skew-T distribution to mean posterior estimates.  
Run `00_Fit_skewT.R` script to estimate skew-T distrib. Can then check that the skew-T distribution are good estimates of the posterior from previous MCMCtree run

### 6.1 Add Distributions onto tree with `01_Add_STcalibs_to_tree.R`.
Needs three input files:
- tree file with node number labels (can be found in out.txt from posterior mcmctree
- node calibrations file used for mcmctree with an additional column showing which node numbers they refer to
- SkewT distributions from previous step

### 6.2 Run MCMCtree again using just the priors
We do this to see if MCMCtrees effective priors overlap well with the skew-T priors we specified. As in step 5.1 we will run it without any data. We will run two independent chains.

### 6.3 Compare effective priors with skewT
When this is finished we compare the effective priors estimated by mcmctree against the skewT distributions and the samples from the posterior using `02_Eval_skewT.R`. This will show us if the skewT distribution priors are effectively capturing the posterior estimates from the whole genome data which will then be used as calibration points in the big multilocus phylogeny dating

### 6.4 Compare with RelTime estimates
I also estimated dates using the rapid relative divergence tool RelTime. Both OLS and branch-length methods in MEGA. I provided the same secondary calibrations as for MCMCTree with uniform probability and hard bounds as RelTime does not allow for soft-bounds. For OLS I used the Maximum Composite Likelihood substitution model + G with a gamma parameter of 2 as RelTime OLS does not allow for HKY+G5. For the branch length method I set max rate ratio to 100. Once finished I saved as nexus, converted to newick and then compared the results in R.
