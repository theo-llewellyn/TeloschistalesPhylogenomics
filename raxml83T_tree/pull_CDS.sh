#PBS -l walltime=00:20:00
#PBS -l select=1:ncpus=1:mem=1gb
#PBS -J 1-1524

#get name of OG file
OG=$(head -n $PBS_ARRAY_INDEX /rds/general/user/tbl19/home/OGs_Telos85T_10sp.txt|tail -n 1)
PREFIX=${OG%%.fa}

mkdir Single_Copy_Orthologues_10sp_CDS
cd Single_Copy_Orthologues_10sp_CDS

#pull fasta headers from protein orthologue file
grep '>' /rds/general/project/theollewellynproject/live/OrthoFinder/Results_Leca125T/Single_Copy_Orthologues_10sp/${OG} > ${PREFIX}_headers.fa
#loop through protein headers
LINES=$(cat ${PREFIX}_headers.fa)
for LINE in $LINES
 do
  #get the fasta header and any sequence until the next fasta then remove the final line and change header to just species name. This is saved to a separate file for each genome
  sed -n -e "/$LINE/,/>/p" /rds/general/project/theollewellynproject/ephemeral/Leca125T_CDS_all.fa | sed '1b;/>.*/d' | sed "s/>.*/${LINE}/" >> ${PREFIX}_CDS.fa
done

rm ${PREFIX}_headers.fa
