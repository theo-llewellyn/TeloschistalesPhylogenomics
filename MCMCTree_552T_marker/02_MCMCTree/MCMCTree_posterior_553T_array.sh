#!/bin/bash
#SBATCH -c 1
#SBATCH -p long
#SBATCH -J MCMCTree_post
#SBATCH -t 30-00:00:00
#SBATCH -o /data/projects/gaya_lab/Teloschistaceae/LICHENS_MCMCTree_553T/00_Esters1a_calibration/01_MCMCtree_posterior_v2/MCMCTree_posterior.out
#SBATCH -e /data/projects/gaya_lab/Teloschistaceae/LICHENS_MCMCTree_553T/00_Esters1a_calibration/01_MCMCtree_posterior_v2/MCMCTree_posterior.err
#SBATCH -a 1-32

module load paml

cd /data/projects/gaya_lab/Teloschistaceae/LICHENS_MCMCTree_553T/00_Esters1a_calibration/01_MCMCtree_posterior_v2/run${SLURM_ARRAY_TASK_ID}/mcmctree_GBM

/home/tl23kg/software/paml4.9j/bin/mcmctree *.ctl
