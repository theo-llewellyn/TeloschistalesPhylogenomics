## 1. Prepare Input files
`cd 00_input_files`

### 1.1 Calibration file
We will prepare the calibration file of dates from the whole genome skewT distributions. This takes the format CLADE|'ST(NUM,NUM,NUM,NUM)'. We can get this from the ST.fitted.dists.G2.40.tsv file. We can use this to label the nodes in the 553T tree using an R script. First we need to edit the table to give it another column header 'node' so R can read the rownames as a new column. Then we use the mrca function of ape to find the MRCAs of all nodes in the 84T tree, find same node in the 553T and then add the skew-t calibration at that node (see step 1.3).

### 1.2 Alignment
The alignment is the 7 gene partitioned alignment we used to make the Telos553T tree. However it needs to be in phylip format for PAML. We can split the concatenated fasta alignment into 7 phylip alignments using AMAS and then cat them again to make a partitioned one.
`qsub fasta2phylip.sh`

### 1.3 Tree
The tree is the ML topology from raxml of the Telos553T partitioned dataset using the Telos85T genome tree as a constraint. We need to add the skewT node calibrations from the Telos85T posterior samples. To do this we first need to edit the tip names in 84T tree to match the 553T

`qsub edit_tip_names.sh`

We also need to root the 553T using the Caliciales
```
pxrr -g LIQ203DILE,LIQ199HEOB -t MARKER_GENES/TREES/Telos553T_muscle_raxml_constrain.raxml.support > MARKER_GENES/TREES/Telos553T_muscle_raxml_constrain.raxml.support.rooted.tree
```
Then we can use the R script `Calibrations_Telos553T.R` which takes the ST distribution from the `ST.fitted.dists.G2.40_rownames.tsv` file, finds which node that corresponds to in the 553T tree and adds it as a node label. This will also generate the unlabelled tree `Telos553T_rooted_baseml.tree` which is needed for the baseml step.


## 2. Baseml to estimate Hessian and Gradient
`cd 01_baseml`
### 2.1 File prep
This is done pretty much the same process as with the genome tree. We make 7 directories for the 7 genes, each with their own control file.

p01 - ITS

p02 - LSU

p03 - MCM7

p04 - RET1

p05 - RPB1

p06 - RPB2

p07 - SSU

We make a blank control file and then replace the alignment name for each gene
`qsub prepare_directories.sh`

We then run mcmctree within each directory `mcmctree *ctl` to generate the input files for baseml. We may need to recompile mcmctree as the max number of taxa is set to 500. This can be done be just editing line 25 in the src mcmctree.c file and then recompiling. When we see "X site patterns read, X sites" we can cancel the job. We remove unecessary files with
```
for i in p*
do cd $i
rm rst* out.* 2base.t rub tmp*out
sed -i 's/method\ \=\ 0/method\ \=\ 1/' tmp*ctl
cd ..
done
```

### 2.2 Run BASEML

`qsub baseml.sh`

Then we concatenate the rst files for the genes
```
for i in `seq 1 7`
do
cat p0$i/rst2 >> in.BV.01".p1-7"
printf "\n\n" >> in.BV.01".p1-7"
done
```

## 3. MCMCTree for divergence times
`cd 02_MCMCTree`

We estimate the average substitution rate the same as with the genome trees. However we have less file wrangling to do as the gene trees are already in the right format.
`qsub root_genetrees.sh`

### 3.1 Prior
Same as genome dataset
`qsub make_dummy_alignment.sh`

We run 5 chains using a dummy alignment and skewT calibrated 553T tree using the `mcmctree_GBM_calibs_LICHENS.ctl` for GBM model.

`qsub MCMCTree_prior_553T.sh`

As with the 84T tree we combine results of 5 chains using `Combine_MCMC.sh` script to generate a single tracer file. We run MCMCtree again on the combined data to get a tree with branch lengths (needs to be run using paml4.9j otherwise the print -1 option doesnt work).

### 3.2 Posterior
We run it using both ILN and GBM models. First we prepare directories with `prepare_directories.sh` and the submit jobs with `MCMCTree_posterior_553T_array_ILN.sh` with `mcmctree_ILN.ctl ` for the ILN model and `MCMCTree_posterior_553T_array.sh` with `mcmctree.ctl` for the GBM model.

### 3.3 Check convergence and post. vs prior

Once these have run we can combine them using `Combine_posterior_chains.sh` script. If we want to check ESS in Tracer we need to run a few extra commands so that the generation time is consecutive `prepare_tracer_files.sh`

To generate the labelled tree we use mcmctree with the mcmc_tracer_gens.txt, dummy alignment and ML tree with the -1 option.
