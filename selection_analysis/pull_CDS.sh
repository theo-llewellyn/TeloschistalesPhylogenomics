#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=1:mem=1gb

#rename ncbi and jgi protein headers so they are unique and contain the name of the taxon
for i in *fa; do PREFIX=${i%%.cds*}; sed -i "s/FUN/${PREFIX}/" $i; done
cat Usnflo1_GeneCatalog_CDS_20160419.fasta | awk '/^>/{print ">Usnea_florida_prot" ++i; next}{print}' > Usnflo1_GeneCatalog_CDS_20160419.renamed.fa
cat Xanpa2_GeneCatalog_CDS_20140928.fasta | awk '/^>/{print ">Xanthoria_parietina_prot" ++i; next}{print}' > Xanpa2_GeneCatalog_CDS_20140928.renamed.fa
cat Clagr3_GeneCatalog_CDS_20111121.fasta | awk '/^>/{print ">Cladonia_grayi_prot" ++i; next}{print}' > Clagr3_GeneCatalog_CDS_20111121.renamed.fa

#copy all CDS to new directory and cat into a single file
cat *cds* > Leca118T_CDS_all.fa

#take the headers of the guidance filtered alignment
cat CDS/Leca118T_CDS_all.fa | seqkit grep -f Hsp90_headers.txt > Hsp90_CDS_Leca118T.fa
