#rename ncbi and jgi protein headers so they are unique and contain the name of the taxon
for i in *fa; do PREFIX=${i%%.proteins*}; sed -i "s/FUN/${PREFIX}/" $i; done
cat Usnflo1_GeneCatalog_proteins_20160419.aa.fasta | awk '/^>/{print ">Usnea_florida_prot" ++i; next}{print}' > Usnflo1_GeneCatalog_proteins_20160419.aa.renamed.fa
cat Xanpa2_GeneCatalog_proteins_20140928.aa.fasta | awk '/^>/{print ">Xanthoria_parietina_prot" ++i; next}{print}' > Xanpa2_GeneCatalog_proteins_20140928.aa.renamed.fa
cat Clagr3_GeneCatalog_proteins_20111121.aa.fasta | awk '/^>/{print ">Cladonia_grayi_prot" ++i; next}{print}' > Clagr3_GeneCatalog_proteins_20111121.aa.renamed.fa

#concatenate all proteins into a single file
cat *.fa > Leca117T_proteins.fa

#make into blastdb
makeblastdb \
-in Leca117T_proteins.fa \
-out Leca117caT \
-parse_seqids \
-dbtype prot

#blastp to pull homologues
QUERY=Hsp90_Anidulans.fasta

blastp \
-num_threads 8 \
-query queries/${QUERY} \
-db blastdb/Leca117T \
-outfmt "6 sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore sseq" \
-evalue 1 \
-out Hsp90.vs.Leca117T.blastp.out

#filter to get matches >35% pident and 40% cov
cat Hsp90.vs.Leca117T.blastp.out | awk '$2>=35 && $3>=40' > Hsp90.vs.Leca117T.blastp.out.filtered

# extract fastas
cut -f 1 Hsp90.vs.Leca117T.blastp.out.filtered | uniq > Hsp90_headers.txt
cat blastdb/Leca117T_proteins.fa | seqkit grep -f Hsp90_headers.txt > Hsp90_Leca117T.fa
