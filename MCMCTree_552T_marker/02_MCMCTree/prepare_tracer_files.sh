#!/bin/bash

seq -f '%.7g' 0 500 3459500 > gens.txt
paste -d' ' gens.txt mcmc_tracer.txt > mcmc_tracer_gens.txt
gcut --complement -f 2 mcmc_tracer_gens.txt > mcmc_tracer_gens.txt1
mv mcmc_tracer_gens.txt1 > mcmc_tracer_gens.txt
