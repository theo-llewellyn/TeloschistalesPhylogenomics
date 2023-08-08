#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=96gb
#PBS -J 1-4

model=$(head -n $PBS_ARRAY_INDEX /rds/general/project/theollewellynproject/live/DNA_Repair_Analysis/PAML/Leca117T/models.txt | tail -n 1)

cd /rds/general/project/theollewellynproject/live/DNA_Repair_Analysis/PAML/Leca117T/Hsp90/${model}

~/software/paml4.9j/bin/codeml *.ctl 
