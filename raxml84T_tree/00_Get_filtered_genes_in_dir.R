# This script was written by Sandra Alvarez-Carretero and adapted for the Teloschistales dataset by Theo Llewellyn

#----------------------#
# 0. CLEAN ENVIRONMENT #
#----------------------#
rm( list = ls( ) )

#--------------------------#
# 1. SET WORKING DIRECTORY #
#--------------------------#
library( rstudioapi ) 
# Get the path to current open R script and find main dir "00_Gene_filtering"
path_to_file <- getActiveDocumentContext()$path
script_wd <- paste( dirname( path_to_file ), "/", sep = "" )
wd <- gsub( pattern = "/scripts", replacement = "", x = script_wd )
setwd( wd )

#---------------------------------------------#
# 2. READ FILE WITH SUMM STATISTICS FOR GENES #
#---------------------------------------------#
# Read file with summary genes 
sum_genes <- read.csv( file = "out_data/genes.csv" )

#--------------------------------------#
# 3. APPLY RESTRICTIONS TO FILTER DATA #
#--------------------------------------#
# Find which genes do not have TCHR sequence
no_TCHR_ind  <- which( sum_genes$TCHR_presence == 0 )
length( no_TCHR_ind )   # 0 

# Find which genes do not have XANPA sequence 
no_XANPA_ind <- which( sum_genes$XANPA_presence == 0 )
length( no_XANPA_ind )   # 0


# Find which genes have less than 10 sequences 
less_10seq_ind  <- which( sum_genes$Num_seqs < 10 ) 
length( less_10seq_ind ) # 0

# Find out which genes do not meet at least one of the previous 
# restrictions 
restrictions <- unique( c( no_TCHR_ind, no_XANPA_ind, less_10seq_ind ) )
length( restrictions )   # 0

# Get filtered genes
filt_genes      <- sum_genes
dim( filt_genes )[1] # 1300 genes
name_filt_genes <- sum_genes[,1]

# >> CHECK | Check we have processed them alright 
all.equal( which( sum_genes$Gene %in% sum_genes$Gene[restrictions] ),
           sort( restrictions ) )
# [1] TRUE
# >> END CHECK

# Check if directory where to copy files exist. Otherwise, create it
if ( ! dir.exists( "filtered_genes/" ) ){
  dir.create( "filtered_genes/" )
}

# Copy filtered genes in a new dir "filtered_genes"
counter = 0 
for( i in seq( 1:length( name_filt_genes ) ) ){
  
  counter = counter + 1
  # Copy files
  file.copy( paste(as.character( name_filt_genes[i] ), "/", sep = "" ),
             paste( "filtered_genes/", sep = ""),
             recursive = TRUE )
  cat( "Gene ", counter, " is: ", as.character( name_filt_genes[i] ), "\n" )
  write( paste( "Gene ", counter, " is: ",
                as.character( name_filt_genes[i] ), "\n", sep = "" ),
         "out_logs/log_01_copy_filtered_genes.txt",
         append = TRUE, sep = "\n" )
}

#-------------------------------------------------------------#
# 4. APPLY EXTRA RESTRICTION: GENES PRESENT IN ALL 84 SPECIES #
#-------------------------------------------------------------#
# Find which genes are present in the 84 species
# and meet the requirements
ind_all_84  <- which( sum_genes$Num_seqs == 84 )
all84sp     <- sum_genes[c(ind_all_84),] # 4 genes !

# >> CHECK | Check all the genes are present in 84 species!
all.equal( length( which( all84sp$Num_seqs == 84 ) ), length( ind_all_84 ) )
#[1] TRUE 
# >> CHECK END

# Find which genes do not have human sequence
no_human_84sp_ind  <- which( all84sp$XANPA_presence == 0 ) 
length( no_human_84sp_ind ) # 0 

# Find which genes do not have mouse sequence
no_mouse_84sp_ind  <- which( all84sp$TCHR_presence == 0 ) 
length( no_mouse_84sp_ind ) # 0 

# Find which genes have under 100 codons 
less_cod_84sp_ind  <- which( all84sp$Num_codons_hum < 100 ) 
length( less_cod_84sp_ind ) # 0 

# Find which genes have less than 10 sequences 
less_10seq_84sp_ind  <- which( all84sp$Num_seqs < 10 )
length( less_10seq_84sp_ind ) #0 

# Get genes names
name_genes_84sp <- all84sp[,1]
length( name_genes_84sp ) # 4

# Check if directory where to copy files exist. Otherwise, create it
if ( ! dir.exists( "filtered_genes_all84sp/" ) ){
  dir.create( "filtered_genes_all84sp/" )
}

# Copy genes that are within the 84sp in a new dir "filtered_genes_all84sp"
counter = 0 
for( i in seq( 1:length( name_genes_84sp ) ) ){
  
  counter = counter + 1
  # Copy files
  file.copy( paste( as.character( name_genes_84sp[i] ), "/", sep = "" ),
             paste( "filtered_genes_all84sp/", sep = ""),
             recursive = TRUE )
  cat( "Gene ", counter, " is: ", as.character( name_genes_84sp[i] ), "\n" )
  write( paste( "Gene ", counter, " is: ",
                as.character( name_genes_84sp[i] ), "\n", sep = "" ),
         "out_logs/log_01_copy_filtered_all84sp_genes.txt",
         append = TRUE, sep = "\n" )
}
