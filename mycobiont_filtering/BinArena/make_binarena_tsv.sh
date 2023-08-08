#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=1:mem=1gb

#get name of contigs file
ACCESSION=ACCESSION

cd $PBS_O_WORKDIR

#filter BlobTools summary for gc, coverage, read length and taxonomy
cut -f 1,2,3,5,6 BlobTools/${ACCESSION}_blobtools.out/*_blobplot.blobDB.bestsum.table.txt | awk 'NR>11' | awk '{if($2>500)print $0}' | sort > ${ACCESSION}_gc_cov2.tsv

#format the concoct clustering file
sed '1d' concoct/${ACCESSION}/concoct_output/clustering_merged.csv | tr ',' '\t' | sort > ${ACCESSION}_clustering_merged.tsv

#join tables and convert to tsv again
join ${ACCESSION}_clustering_merged.tsv ${ACCESSION}_gc_cov2.tsv | tr ' ' '\t' > ${ACCESSION}_gc_cov2_bins.tsv

#make bin column into text for binarena and edit format of particular taxonomy that has a space in
awk 'BEGIN{FS=OFS="\t"}{$2="bin_"$2}1' ${ACCESSION}_gc_cov2_bins.tsv | sed 's/Candidatus\t/Candidatus_/' > ${ACCESSION}_gc_cov2_bins_binarena.tsv

#add column headers
echo -e "contig_id\tmetagenome_bin\tlength\tGC\tcov\tTaxonomy" | cat - ${ACCESSION}_gc_cov2_bins_binarena.tsv > ${ACCESSION}_gc_cov2_bins_binarena_final.tsv

rm ${ACCESSION}_gc_cov2.tsv ${ACCESSION}_clustering_merged.tsv ${ACCESSION}_gc_cov2_bins.tsv ${ACCESSION}_gc_cov2_bins_binarena.tsv
mv ${ACCESSION}_gc_cov2_bins_binarena_final.tsv binarena
