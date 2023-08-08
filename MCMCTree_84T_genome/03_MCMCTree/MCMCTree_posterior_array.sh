#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=5gb
#PBS -J 1-16

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

cd 01_MCMCtree_posterior/run${PBS_ARRAY_INDEX}/mcmctree_GBM

PATH=/rds/general/user/tbl19/home/software/paml-4.10.5/bin:$PATH

mcmctree *ctl
