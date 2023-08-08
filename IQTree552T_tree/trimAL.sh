#bin/bash

for i in *.fa; do GENE=${i%%.fa}; trimal -in $i -out ${GENE}_tAl.fa -fasta -automated1; done
