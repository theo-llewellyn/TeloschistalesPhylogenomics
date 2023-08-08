## 1. Prepare Input files

### 1.1 Calibration file
We will prepare the calibration file of dates from the whole genome skewT distributions. This takes the format CLADE|'ST(NUM,NUM,NUM,NUM)'. We can get this from the ST.fitted.dists.G2.40.tsv file. We can use this to label the nodes in the 553T tree using an R script. First we need to edit the table to give it another column header 'node' so R can read the rownames as a new column. Then we use the mrca function of ape to find the MRCAs of all nodes in the 84T tree, find same node in the 553T and then add the skew-t calibration at that node (see step 1.3).

### 1.2 Alignment
The alignment is the 7 gene partitioned alignment we used to make the Telos553T tree. However it needs to be in phylip format for PAML. We can split the concatenated fasta alignment into 7 phylip alignments using AMAS
```
sed -i 's/DNA, //' Telos553T_muscle_partitions_v2.txt 
~/software/AMAS/amas/AMAS.py split -l Telos553T_muscle_partitions_v2.txt -u phylip -f fasta -d dna -i Telos553T_muscle_concatenated_v2.fa
```
We can just cat them to make the partitioned one.
```
cat Telos553T_muscle_concatenated_p1_Telos553T_ITS_muscle5_msa-out.phy Telos553T_muscle_concatenated_p2_Telos553T_LSU_muscle5_msa-out.phy Telos553T_muscle_concatenated_p3_Telos553T_MCM7_muscle5_msa-out.phy Telos553T_muscle_concatenated_p4_Telos553T_RET1_muscle5_msa-out.phy Telos553T_muscle_concatenated_p5_Telos553T_RPB1_muscle5_msa-out.phy Telos553T_muscle_concatenated_p6_Telos553T_RPB2_muscle5_msa-out.phy Telos553T_muscle_concatenated_p7_Telos553T_SSU_muscle5_msa-out.phy > Telos553T_7parts.aln
```

### 1.3 Tree
The tree is the ML topology from raxml of the Telos553T partitioned dataset using the Telos85T genome tree as a constraint. We need to add the skewT node calibrations from the Telos85T posterior samples. To do this we first need to edit the tip names in 84T tree to match the 553T
```
sed 's/L6R42/LIQ80XSP/g;s/LQ337/LIQ153CAOA/g;s/LQ341/LIQ74CAUR/g;s/LQ352/LIQ145TFLA/g;s/L6R38/LIQ146XSP/g;s/L6R39/LIQ85CALI/g;s/LQ348/LIQ92SELA_2/g;s/LQ345/LIQ72TVIL/g;s/L6R36/LIQ73XASTE-2/g;s/LQ339/LIQ75XAME-2/g;s/LQ351/LIQ93LETR/g;s/LQ338/LIQ69CCAR/g;s/L6R37/LIQ70TSP/g;s/L6R41/LIQ106LELE/g;s/LQ346/LIQ82CAATT/g;s/LQ349/LIQ109XAAU_2/g;s/L6R35/LIQ143CAAG/g;s/LQ343/LIQ76CEHR/g;s/LQ350/LIQ78TCHR-2/g;s/L6R40/LIQ81XMZF-2/g;s/LQ342/LIQ94LETR-2/g;s/LQ344/LIQ71TLAC/g;s/Xanthoria_parietina/Xanpa2/g;s/LQ340/LIQ79DIDIA/g;s/LQ347/LIQ84UMVE/g' FigTree_84sp_nodelabels.tree > FigTree_84sp_nodelabels_renamed.tree
```
We also need to root the 553T using the Caliciales
```
pxrr -g LIQ203DILE,LIQ199HEOB -t MARKER_GENES/TREES/Telos553T_muscle_raxml_constrain.raxml.support > MARKER_GENES/TREES/Telos553T_muscle_raxml_constrain.raxml.support.rooted.tree
```
Then we can use the R script `Calibrations_Telos553T.R` which takes the ST distribution from the `ST.fitted.dists.G2.40_rownames.tsv` file, finds which node that corresponds to in the 553T tree and adds it as a node label. This will also generate the unlabelled tree `Telos553T_rooted_baseml.tree` which is needed for the baseml step.


## 2. Baseml to estimate Hessian and Gradient
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
```
cd ~/LICHENS_553T_Hessian/
wd=$( pwd )

for i in p*
do
cp mcmctree_LICHENS.ctl $i
#cp Telos553T_rooted_baseml.tree $i #if using raxml tree
cp Telos553T_IQTree_rooted_baseml.tree $i # if using IQTree
cd $i
tree_name=$( ls *tree )
aln_name=$( ls *phy )
sed -i -e 's/ALN/'${aln_name}'/' *ctl 
sed -i -e 's/TREE/'${tree_name}'/' *ctl
cd ..
done
```
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

```
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=5gb
#PBS -J 1-7

cd /rds/general/user/tbl19/home/LICHENS_553T_Hessian/p0${PBS_ARRAY_INDEX}

PATH=/rds/general/user/tbl19/home/software/paml-4.10.5/bin:$PATH

baseml tmp*ctl
```
Then we concatenate the rst files for the genes
```
for i in `seq 1 7`
do
cat p0$i/rst2 >> in.BV.01".p1-7"
printf "\n\n" >> in.BV.01".p1-7"
done
```

## 3. MCMCTree for divergence times
We estimate the average substitution rate the same as with the genome trees. However we have less file wrangling to do as the gene trees are already in the right format.
```
for i in *.support
do
grep -Eo '(LIQ199HEOB|LIQ203DILE)' $i > outgroups_${i}.txt
pxrr -f outgroups_${i}.txt -t $i > ${i}.rooted
done

rm outgroups_*.txt
cat *.rooted > rooted_trees.tree
```

### 3.1 Prior
Same as genome dataset
```
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
```
We run 5 chains using a dummy alignment and skewT calibrated 553T tree

### 3.2 Posterior
We run it using both ILN and GBM models
```
for i in `seq 1 32`
do
mkdir run$i
mkdir run$i/mcmctree_GBM
mkdir run$i/mcmctree_ILN
cp mcmctree.ctl run$i/mcmctree_GBM
cp mcmctree_ILN.ctl run$i/mcmctree_ILN
cp in.BV* run$i/mcmctree_GBM/in.BV
cp in.BV* run$i/mcmctree_ILN/in.BV
done
```
### 3.3 Check convergence and post. vs prior

Once these have run we can combine them using `Combine_MCMC.sh` script. If we want to check ESS in Tracer we need to run a few extra commands so that the generation time is consecutive.
```
seq -f '%.7g' 0 500 3459500 > gens.txt
paste -d' ' gens.txt mcmc_tracer.txt > mcmc_tracer_gens.txt
gcut --complement -f 2 mcmc_tracer_gens.txt > mcmc_tracer_gens.txt1
mv mcmc_tracer_gens.txt1 > mcmc_tracer_gens.txt
```
To generate the labelled tree we use mcmctree with the mcmc_tracer_gens.txt, dummy alignment and ML tree with the -1 option.
