#!/bin/bash


#This script was written by Sandra Alvarez-Carretero and adapted for the Teloschistales dataset by Theo Llewellyn

# 1. Set counters (1254/4 = 313.5 genes per partition)
i=0
jinc=314
j=314
count=0
count_per_parts=0
count_del_dirs=0

# 2. Move to main dir `00_Gene_filtering`, from where 
#    this script is run
curr_dir=$( pwd )
cd $curr_dir

# 3. Create `000_alignments` if not there 

if [ ! -d 000_alignments ] 
then 
mkdir 000_alignments 
fi 

# 4. Combine all the different fasta alignment in a unique
#    fasta file for each partition
for partition in `seq 1 4`
do 
	# Create folder for each partition 
	mkdir -p 000_alignments/part$partition
	
	# If it is not the last partition
	if [ $partition != 4 ]
	then
		## 210622 - modified to fit new data 								
		#while [ $i -lt $j ]
		while [ $count_per_parts -lt $jinc ]
		do
			i=$(( i + 1 ))
			## 210622 - new while loop to account for dirs removed 
			##          in step 3
			## If this directory was removed, increase i again 
			## until next directory is found
			while [[ ! -d filtered_genes_step2/$i  ]]
			do
			count_del_dirs=$(( count_del_dirs + 1 ))
			printf "Have to increase i plus 1 as dir "$i" does not exist\n"
			printf "Have to increase i plus 1 as dir "$i" does not exist\n"  >> out_logs/log_12_file_generate_alns.txt
			i=$(( i + 1 ))
			done
			count=$(( count + 1 ))
			count_per_parts=$(( count_per_parts + 1 ))
			name_gene=$( echo `ls filtered_genes_step2/${i}/*aln` | sed 's/..*\/partitions12\_//' | sed 's/\.aln//'  )
			# Convert to fasta-phylip-partitions needed input format (fasta format)
			cp filtered_genes_step2/${i}/*aln 000_alignments/part$partition/gene$count"_"$name_gene.fasta
			sed -i -e '1,1d' 000_alignments/part$partition/gene$count"_"$name_gene.fasta          # Remove phylip header
			sed -i -e 's/^/\>/' 000_alignments/part$partition/gene$count"_"$name_gene.fasta       # Add ">" at the beginning of each sequence 
			sed -i -e 's/\s\{6\}/\n/' 000_alignments/part$partition/gene$count"_"$name_gene.fasta # Replace spaces between taxa names and seqs with new line
			echo Gene $count visited
			echo Gene $count visited >> out_logs/log_12_file_generate_alns.txt
		done 
		
		# Update counters for genes
		printf "Total of genes for partition "$partition": "$count_per_parts"\n"
		printf "Total of genes for partition "$partition": "$count_per_parts"\n" >> out_logs/log_12_file_generate_alns.txt
		## 210622 - This is not true anymore. Now count 
		## does not equal to i as there are missing dirs !
		#i=$( echo $count )
		#j=$(( i + jinc ))
		j=$(( count + jinc ))
		count_per_parts=$( echo 0 )
	
	# Otherwise, just reach until the end of genes available 
	# in the directory
	else
		tot=$( ls filtered_genes_step2 | wc -l )
		k=$(( tot ))
		## 210622 - Not needed anymore
		# Set i as the last gene !
		#i=$( echo $count )
		#while [ $i -lt $k ]
		while [ $count -lt $k ]
		do
			i=$(( i + 1 ))
			## 210622 - new new while loop to account for dirs removed 
			##          in step 3			
			## If this directory was removed, increase i again 
			## until next directory is found
			while [[ ! -d filtered_genes_step2/$i  ]]
			do
			count_del_dirs=$(( count_del_dirs + 1 ))
			printf "Have to increase i plus 1 as dir "$i" does not exist\n"
			printf "Have to increase i plus 1 as dir "$i" does not exist\n"  >> out_logs/log_12_file_generate_alns.txt
			i=$(( i + 1 ))
			done
			count=$(( count + 1 ))
			name_gene=$( echo `ls filtered_genes_step2/${i}/*aln` | sed 's/..*\/partitions12\_//' | sed 's/\.aln//'  )
			cp filtered_genes_step2/${i}/*aln 000_alignments/part$partition/gene$count"_"$name_gene.fasta
			# Convert to fasta-phylip-partitions needed input format (fasta format)
			sed -i -e '1,1d' 000_alignments/part$partition/gene$count"_"$name_gene.fasta          # Remove phylip header
			sed -i -e 's/^/\>/' 000_alignments/part$partition/gene$count"_"$name_gene.fasta       # Add ">" at the beginning of each sequence 
			sed -i -e 's/\s\{6\}/\n/' 000_alignments/part$partition/gene$count"_"$name_gene.fasta # Replace spaces between taxa names and seqs with new line
			sed -i -e 's/_prot[0-9]*[^:]//g' 000_alignments/part$partition/gene$count"_"$name_gene.fasta #removes protein number from species names so they are same across genes
			sed -i -e '/^>.*/s/_[0-9].*-T1//' 000_alignments/part$partition/gene$count"_"$name_gene.fasta
			echo Gene $count visited
			echo Gene $count visited >> out_logs/log_12_file_generate_alns.txt
			count_per_parts=$(( count_per_parts + 1 ))
		done 
		printf "Total of genes for partition "$partition": "$count_per_parts"\n"
		printf "Total of genes for partition "$partition": "$count_per_parts"\n" >> out_logs/log_12_file_generate_alns.txt
	
	fi
	
	printf "Partition "${partition}" finishes !\n"

done 

printf "Total of dirs that had to be skipped: "$count_del_dirs"\n"
printf "Total of dirs that had to be skipped: "$count_del_dirs"\n" >> out_logs/log_12_file_generate_alns.txt
