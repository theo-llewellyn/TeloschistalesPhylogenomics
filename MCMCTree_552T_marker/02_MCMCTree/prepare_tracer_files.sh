#!/bin/bash

cd mcmc_files
lines=$(wc -l mcmc_tracer.txt | cut -f 1 -d ' ')
gens=$((($lines - 1)*500))

seq -f '%.7g' 0 500 $gens > gens.txt
paste -d' ' gens.txt mcmc_tracer.txt > mcmc_tracer_gens.txt
cut --complement -f 2 mcmc_tracer_gens.txt > mcmc_tracer_gens.txt1
mv mcmc_tracer_gens.txt1 mcmc_tracer_gens.txt

lines=$(wc -l mcmc_tracer_ILN.txt | cut -f 1 -d ' ')
gens=$((($lines - 1)*500))

seq -f '%.7g' 0 500 $gens > gens.txt
paste -d' ' gens.txt mcmc_ILN_tracer.txt > mcmc_ILN_tracer_gens.txt
cut --complement -f 2 mcmc_ILN_tracer_gens.txt > mcmc_tracer_gens.txt1
mv mcmc_tracer_gens.txt1 mcmc_ILN_tracer_gens.txt

rm gens.txt
