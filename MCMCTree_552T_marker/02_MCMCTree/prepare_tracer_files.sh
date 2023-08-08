#!/bin/bash

cd mcmc_files
#count number of combined generations
lines=$(wc -l mcmc_tracer.txt | cut -f 1 -d ' ')
gens=$((($lines - 1)*500))

#make a column with generation number every 500
seq -f '%.7g' 0 500 $gens > gens.txt
sed -i '1 s/0/1/' gens.txt
sed -i '1 i\Gen' gens.txt
#add it as a column to first file
paste -d'\t' gens.txt mcmc_tracer.txt > mcmc_tracer_gens.txt
cut --complement -f 2 mcmc_tracer_gens.txt > mcmc_tracer_gens.txt1
mv mcmc_tracer_gens.txt1 mcmc_tracer_gens.txt
sed -i '$ d' mcmc_tracer_gens.txt

lines=$(wc -l mcmc_ILN_tracer.txt | cut -f 1 -d ' ')
gens=$((($lines - 1)*500))

seq -f '%.7g' 0 500 $gens > gens.txt
sed -i '1 s/0/1/' gens.txt
sed -i '1 i\Gen' gens.txt
paste -d'\t' gens.txt mcmc_ILN_tracer.txt > mcmc_ILN_tracer_gens.txt
cut --complement -f 2 mcmc_ILN_tracer_gens.txt > mcmc_tracer_gens.txt1
mv mcmc_tracer_gens.txt1 mcmc_ILN_tracer_gens.txt
sed -i '$ d' mcmc_ILN_tracer_gens.txt
rm gens.txt
