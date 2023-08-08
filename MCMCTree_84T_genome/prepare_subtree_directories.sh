#!/bin/bash

cp -r ~/Telos85T_10sp_msa/filtered_genes_step2 ~/MCMC3r
cd ~/MCMC3r/filtered_genes_step2

#rename taxon headers
for i in *
do
cd $i
sed -e 's/_prot[0-9]*[^:]//g' partitions12_*_Telos.aln > partitions12_Telos.aln
sed -i -e 's/_[0-9].*-T1//' partitions12_Telos.aln
cd ..
done
