#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=100gb
#PBS -J 1-32

cd /rds/general/user/tbl19/home/MCMC3r/clk/OG0001822/$PBS_ARRAY_INDEX

PATH=/rds/general/user/tbl19/home/software/paml-4.10.5/bin:$PATH

mcmctree > /dev/null
