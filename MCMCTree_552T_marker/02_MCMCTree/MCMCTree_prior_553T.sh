#!/bin/bash
#SBATCH -c 1
#SBATCH -p long
#SBATCH -J MCMCTree_prior
#SBATCH -t 30-00:00:00
#SBATCH -o /data/projects/gaya_lab/Teloschistaceae/LICHENS_MCMCTree_553T/00_Esters1a_calibration/00_MCMCtree_prior_v2/MCMCTree_prior_553T.out
#SBATCH	-e /data/projects/gaya_lab/Teloschistaceae/LICHENS_MCMCTree_553T/00_Esters1a_calibration/00_MCMCtree_prior_v2/MCMCTree_prior_553T.err
#SBATCH -a 1-5

cd /data/projects/gaya_lab/Teloschistaceae/LICHENS_MCMCTree_553T/00_Esters1a_calibration/00_MCMCtree_prior_v2/mcmc${SLURM_ARRAY_TASK_ID}

/home/tl23kg/software/paml4.9j/bin/mcmctree *.ctl 
