#!/bin/bash

for i in {1..32}
 do
 printf "Parsing run"$i"/mcmctree_GBM ... ... \n"
 end=$( wc -l run$i/mcmctree_GBM/mcmc.txt | sed 's/ ..*//' )
 if [[ $i -eq 1 ]]
  then begin=1
  else begin=2
 fi
 sed -n ''${begin}','${end}'p' run$i/mcmctree_GBM/mcmc.txt >> mcmc_files/mcmc_tracer.txt
 sed -n '1,'${end}'p' run$i/mcmctree_GBM/mcmc.txt >>  run$i/mcmctree_GBM/mcmc_clean.txt
done 

for i in {1..32}
 do
 printf "Parsing run"$i"/mcmctree_ILN ... ... \n"
 end=$( wc -l run$i/mcmctree_ILN/mcmc.txt | sed 's/ ..*//' )
 if [[ $i -eq 1 ]]
  then begin=1
  else begin=2
 fi
 sed -n ''${begin}','${end}'p' run$i/mcmctree_ILN/mcmc.txt >> mcmc_files/mcmc_ILN_tracer.txt
 sed -n '1,'${end}'p' run$i/mcmctree_ILN/mcmc.txt >>  run$i/mcmctree_ILN/mcmc_clean.txt
done 
