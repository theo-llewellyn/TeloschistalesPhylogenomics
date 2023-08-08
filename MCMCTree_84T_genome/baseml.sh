#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=5gb
#PBS -J 1-4

cd /rds/general/user/tbl19/home/LICHENS_Hessian/baseml_method1/p0${PBS_ARRAY_INDEX}

PATH=/rds/general/user/tbl19/home/software/paml-4.10.5/bin:$PATH

baseml tmp*ctl
