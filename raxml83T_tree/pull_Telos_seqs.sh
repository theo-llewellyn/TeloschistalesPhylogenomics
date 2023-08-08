mkdir Single_Copy_Orthologues_50percent
#search for orthogroup headers and save into new directory
cat Orthologues_Telos85T_dup_comp_10sp.txt | while read line; do cp Orthogroup_Sequences/${line}.fa Single_Copy_Orthologues_10sp; done

FILENAME="Telos85T.txt"

#loop through each orthologue file
for i in *.fa
 #extract just the species name from filename
 do PREFIX=${i%%.fa}
   #read lines of Telos85T
    LINES=$(cat $FILENAME)
     #loop through Telos. taxa fasta headers
     for LINE in $LINES
      #for each taxon pull all sequences in said orthogroup and save to new file
      do sed -n -e "/$LINE/,/>/p" $i | sed '1b;/>.*/d' >> ${PREFIX}_Telos.fa
     done
 done
 
find . -type f ! -name '*Telos.fa' -delete
