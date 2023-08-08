
#This script was written by Sandra Alvarez-Carretero and adapted for the Teloschistales dataset by Theo Llewellyn

#----------------------#
# 0. CLEAN ENVIRONMENT #
#----------------------#
rm( list = ls( ) )

#-------------------#
# 1. LOAD LIBRARIES #
#-------------------#
library( ape )

#--------------------------#
# 2. SET WORKING DIRECTORY #
#--------------------------#
library( rstudioapi ) 
# Get the path to current open R script and find main dir "00_Gene_filtering"
path_to_file <- getActiveDocumentContext()$path
script_wd <- paste( dirname( path_to_file ), "/", sep = "" )
wd <- gsub( pattern = "/scripts", replacement = "", x = script_wd )
setwd( wd )

#-----------#
# 3. TASKS  #
#-----------#
## 3.1. Read trees from "filtered_genes" and save them in a list
trees.dir <- "/rds/general/project/theollewellynproject/live/raxml-ng/Telos85T_10sp_genetrees/"
num.genes <- length( list.files( trees.dir ) )

lichen.trees <- vector( mode = "list", num.genes )
names( lichen.trees ) <- list.files( trees.dir )

count <- 0
for ( i in 1:num.genes ){
  
  files <- list.files( paste( trees.dir,
                              list.files( trees.dir )[i],
                              sep = "" ) )
  check.files <-  grep( pattern = "bestTree", files )
  
  if( length( check.files ) != 0 ){
    
    count <- count + 1
    lichen.trees[[ i ]] <- ape::read.tree( file = paste( trees.dir,
                                                         list.files( trees.dir )[i],
                                                         "/",
                                                         list.files( trees.dir )[i],
                                                         ".raxml.bestTree",
                                                         sep = "" ) )
    # Ensure the name of the entry in the list is the same gene that has been added
    # to this entry
    names( lichen.trees )[i] <- list.files( trees.dir )[i]
  }
  
  cat( count, "\n" )
  
}

# Save RData objects with filtered genes
save( lichen.trees, file = "out_RData/lichen.trees.RData" )

## 3.2. Sum branch lengths for each tree to get the tree length 
blengths  <- as.data.frame( matrix( 0, nrow = length( lichen.trees ), ncol = 5 ) )
rownames( blengths ) <- names( lichen.trees )
colnames( blengths ) <- c( "Sum.blengths", "Largbr.rel2tree",
                           "Tip.name", "H-M.dist",
                           "Delete?" )

for( i in 1:length( lichen.trees ) ){

  blengths[i,1] <- sum( lichen.trees[[i]]$edge.length )

}


#-------------------------------------#
# 4. CALCULATE RELATIVE BRANCH LENGTH #
#-------------------------------------#
# Take each branch length (i.e., lichen.trees[[i]]$edge.length)
# and divide it into the corresponding tree length that was stored
# in column 1 for each gene (i.e., blenghts[i,1]). These will result 
# into a vector of relative branch lengths which, if multiplied
# by 100, you can use a useful statistic to detect possible
# alignment errors. If something went very
# very wrong (threshold set to relative branch lengths >= 60%),
# then the alignment makes no sense and the tree has
# a very incredible long branch for that species.
counter <- 0
for( i in 1:length( lichen.trees ) ){
  
  # Divide all branches into the tree length so we can obtain the relative 
  # branch lengths, i.e., r_{ij} = b{ij} / SUM_j (b_{ij}); as detailed in
  # dos Reis et al. 2012.
  # Then, take the largest relative branch length (in %) for each gene and store 
  # it in second  column in data.frame `blengths`.
  blengths[i,2] <- max( lichen.trees[[i]]$edge.length/blengths[i,1] )*100
  print( max( lichen.trees[[i]]$edge.length/blengths[i,1] )*100 )
  
  
  # If the largest relative branch length is equal to or larger than 60%...
  if ( round( blengths[i,2] ) >= 60 ){
    cat( "Gene ", rownames( blengths)[i], " - row", i,
         "- has a branch length longer than 60% of the total tree length\n" )
    write( paste( "Gene ", rownames( blengths)[i], " - row", i,
                  "- has a branch length longer than 60% of the total tree length: ",
                  blengths[i,2], "%\n", sep = ""),
           file = "out_logs/log_04_R_genes_blenghth_longer_60pc.txt",
           append = T )
    counter <- counter + 1
    blengths[i,5] <- c("Y")
    
    # 1. Find which position has largest branch length
    ind1   <- which( round( lichen.trees[[i]]$edge.length/blengths[i,1] ) >= 0.6 )
    # 2. In case more than one branch length was larger than 0.6 and stored as
    #   `ind1`, find the largest and save it as `ind1.1`
    ind1.1 <- max( c(lichen.trees[[i]]$edge.length/blengths[i,1])[ind1] )
    # 3. Find the position of this largest branch which index is `ind1.1`, just
    #    in case there were more than one saved in step 1
    ind1.2 <- which( lichen.trees[[i]]$edge.length/blengths[i,1] == ind1.1 )
    # 4. Find the node position of the individual which branch length is `ind1.2`
    ind2   <- lichen.trees[[i]]$edge[ind1.2,2]
    blengths[i,3] <- lichen.trees[[i]]$tip.label[ind2]
    
    # 5. Now we want to find the corresponding taxon name. But if
    #    ind2 is larger than the amount of taxa available in
    #    this gene tree, then increase ind1.2 by 1 until
    #    it finds the correct taxon name.
    if( ind2 >= length( lichen.trees[[i]]$tip.label ) ){
      while( ind2 >= length( lichen.trees[[i]]$tip.label ) ){
        ind1.2 <- ind1.2 + 1
        ind2   <- lichen.trees[[i]]$edge[ind1.2,2]
      }
      blengths[i,3] <- paste( lichen.trees[[i]]$tip.label[ind2], ",",
                              lichen.trees[[i]]$tip.label[ind2+1], sep = "" )
    }else{
      blengths[i,3] <- lichen.trees[[i]]$tip.label[ind2]
    }
  }else{
    blengths[i,3] <- "NULL"
    blengths[i,5] <- c("N")
  }
  
  # Append that the analysis has finished in the log file!
  if( i == length( lichen.trees ) ){
    cat("\nA total of ", counter, "genes had a branch length longer than 60% of
        the total tree length\n")
    write( paste( "\nA total of ", counter, "genes had a branch length longer than 60% of
                  the total tree length\n", sep = "" ),
           file = "out_logs/log_04_R_genes_blenghth_longer_60pc.txt",
           append = T )
  }
  
}

# Now, identify the 133 genes identified to have at least one branch length 
# longer than 60% of the total tree length.
# This will be later erased, but not now.
ind.delete <- which( blengths[,5] == "Y" )

# Save RData objects with filtered genes
save( blengths, file = "out_RData/mammals_summary_matrix_filtstep1.RData" )

#------------------------------------------------------#
# 5. Calculating pairwise distance from mouse to human #
#------------------------------------------------------#
## 5.1. Read *_XanPa_TelChr.aln files in "baseml" and save them in a list
baseml.dir <- c( "baseml/" )

XT_aln.aln <- vector( mode = "list", num.genes )
names( XT_aln.aln ) <- list.files( trees.dir )

count <- 0
for ( i in 1:num.genes ){
  
  files <- list.files( paste( baseml.dir, i, sep = "" ) )
  check.files <-  grep( pattern = "TCHR_XANPA", files )
  
  if( length( check.files ) != 0 ){
    
    count <- count + 1
    # Read the file
    aln <- read.table( file = paste( baseml.dir, i,
                                     "/", files[check.files], sep = "" ),
                       sep = "\t", stringsAsFactors = FALSE )
    # Get species names
    sp_names <- aln[,1]
    
    # Get length of sequences and construct the matrix that will be later
    # used to calculate the H-M distance
    num_seqs <- length( unlist( strsplit(aln[1,2], split = "" ) ) )
    aln_mat  <- matrix( 0, nrow = dim( aln )[1], ncol = num_seqs )
    rownames( aln_mat ) <- sp_names
    
    # Append to matrix
    for( j in seq( 1:dim(aln)[1] ) ){
      
      aln_mat[j,] <- unlist( strsplit(aln[j,2], split = "" ) )
      
    }
    
    # Put in the list
    XT_aln.aln[[ i ]] <- ape::as.DNAbin( aln_mat )
    
  }
  
  cat( count, "\n" )
  
}

## 5.2. Save RData objects with XT alignments as DNAbin format.
##      If you have already generated this file, you can run the 
##      second command and load the corresponding RData file
save( XT_aln.aln, file = "out_RData/XTalns_DNAbin.RData" )

## 5.3. Go through every matrix in the XT_aln.aln list 
##      and compute the XT distance
## We picked the TN93 model because this is the one that BASEML uses:
## Page 26 PAML doc:
## ** The pairwise sequence distances are included in the output as well, and also in
## ** a separate file called 2base.t. This is a lower-diagonal distance matrice, readable
## ** by the NEIGHBOR program in Felesenstein's PHYLIP package (Felsenstein 2005).
## ** For models JC69, K80, F81, F84, the appropriate distance formulas are used,
## ** while for more complex models, the TN93 formula is used. BASEML is mainly a maximum
## ** likelihood program, and the distance matrix is printed out for convenience and really
## ** has nothing to do with the later likelihood calculation.

# Use model TN93
count <- 0
for( i in 1:length( XT_aln.aln ) ){
  
  count <- count + 1
  XT.dist <- ape::dist.dna( XT_aln.aln[[ i ]], model = "TN93" )
  blengths[i,4] <- XT.dist
  
  # Safety check
  files <- list.files( paste( baseml.dir, i, sep = "" ) )
  check.files <- grep( pattern = "XanPa_TelChr", files )
  cat( "Gene as name in list: ", names( XT_aln.aln )[i], "\n",
       "Gene as name in order it was appended: ", files[check.files], "\n",
       "Gene as name in matrix blengths: ", rownames( blengths )[i], "\n" )
  cat( count, "\n" )
  
  write( paste( "Gene as name in list: ", names( XT_aln.aln )[i], "\n",
                "Gene as name in order it was appended: ", files[check.files], "\n",
                "Gene as name in matrix blengths: ", rownames( blengths )[i], "\n",
                sep = "" ),
         file = "out_logs/log_05_R_check_genes_included_XT_aln_TN93.txt",
         append = TRUE )
  
}

# Do the same with JC69 and raw models. So extend the data frame two more columns.
# This is why the object `blengths2` is created here.
blengths2       <- as.data.frame( matrix( 0, nrow = length( lichen.trees ), ncol = 6 ) )
blengths2[,1:4] <- blengths[,1:4]
rownames( blengths2 ) <- names( lichen.trees )
colnames( blengths2 ) <- c( "Sum.blengths", "Largbr.rel2tree", "Tip.name",
                            "H-M.dist.TN93", "H-M.dist.JC69", "H-M.dist.raw" )

count <- 0
for( i in 1:length( XT_aln.aln ) ){
  
  count <- count + 1
  XT.dist.JC69   <- ape::dist.dna( XT_aln.aln[[ i ]], model = "JC69" )
  blengths2[i,5] <- XT.dist.JC69
  XT.dist.raw    <- ape::dist.dna( XT_aln.aln[[ i ]], model = "raw" )
  blengths2[i,6] <- XT.dist.raw
  
  # Safety check
  files <- list.files( paste( baseml.dir, i, sep = "" ) )
  check.files <-  grep( pattern = "XanPa_TelChr", files )
  cat( "Gene as name in list: ", names( XT_aln.aln )[i], "\n",
       "Gene as name in order it was appended: ", files[check.files], "\n",
       "Gene as name in matrix blengths2: ", rownames( blengths2 )[i], "\n" )
  cat( count, "\n" )
  write( paste( "Gene as name in list: ", names( XT_aln.aln )[i], "\n",
                "Gene as name in order it was appended: ", files[check.files], "\n",
                "Gene as name in matrix blengths2: ", rownames( blengths2 )[i], "\n",
                sep = "" ),
         file = "out_logs/log_05_R_check_genes_included_XT_aln_JC69_raw.txt",
         append = TRUE )
  
}

# Save the object so you can avoid running the loop above
save( blengths2,
      file = "out_RData/mammals_summary_matrix.TN93.JC69.raw_filtstep2.RData" )

#------------------------------------------------------------#
# 6. Find those genes which the largest single branch is     #
#    larger than or equal to 60% relative to the tree length #
#------------------------------------------------------------#
# Find index of largest bl 
ind.larg.b  <- which( round(blengths2$Largbr.rel2tree) >= 60 )
# >> Check that we have 133 and that they are the same than found before
length( ind.larg.b ) # 14
setdiff( ind.larg.b, ind.delete ) # 0
# >> SUCCESS! 
filt.blengths <- blengths2[-c( ind.larg.b ),]
# >> Safety check: Double check that now all rows in col 3 are "NULL"
length( which( filt.blengths[,3] == "NULL" ) ) == dim( filt.blengths )[1]
# Double check and OK!
filt.df.blengths <- filt.blengths[,-3] # Remove 3rd column as double check done!
# >> SUCCESS!

# Which are the tips which bl is too large?
large.blengths <- blengths2[ind.larg.b,]
tips.longb <- factor( blengths2[ind.larg.b,3] ) #24 different
plot( tips.longb, las = 3 )

#-------------------------------------------#
# 7. Order filtered genes from slow to fast #
#-------------------------------------------#
# Create objects with ordered genes (fast- to slow-evolving genes) according to 
# the three models tested
ord.blengths.f2s.TN93 <- filt.df.blengths[order( round( filt.df.blengths$`H-M.dist.TN93`, 
                                                        15) ),]
ord.blengths.f2s.JC69 <- filt.df.blengths[order( round( filt.df.blengths$`H-M.dist.JC69`, 
                                                        15) ),]
ord.blengths.f2s.raw <- filt.df.blengths[order( round( filt.df.blengths$`H-M.dist.raw`, 
                                                       15) ),]

# Create objects with ordered genes (slow- to fast-evolving genes) according to 
# the three models tested
ord.blengths.s2f.TN93 <- filt.df.blengths[order( round( filt.df.blengths$`H-M.dist.TN93`, 
                                                        15 ),decreasing = TRUE ),]
ord.blengths.s2f.JC69 <- filt.df.blengths[order( round( filt.df.blengths$`H-M.dist.JC69`, 
                                                        15 ),decreasing = TRUE ),]
ord.blengths.s2f.raw <- filt.df.blengths[order( round( filt.df.blengths$`H-M.dist.raw`, 
                                                       15 ),decreasing = TRUE ),]

# >> CHECK | Which are NaN or Inf --> very different sequences.
#    Do this check in the three data frames previously created with ordered genes. 

# 2. Use object `ord.blengths.s2f.JC69`, which genes were ordered from slow- 
#    to fast-evolving according to the model `JC69`.
# Column H-M.dist.TN93 (index = 3)
which( ord.blengths.s2f.JC69[,3] == "Inf" )
# integer(0)
which( ord.blengths.s2f.JC69[,3] == "NaN" ) 
#[1] 1 15436
which( ord.blengths.s2f.JC69[,3] >= 0.75 )
#[1] 4 8
ord.blengths.s2f.JC69[c( 1, 4, 8, 15436 ),]
rownames( ord.blengths.s2f.JC69 )[ c( 1, 4, 8, 15436 )]
# [1] "OG0006811_Telos" "OG0006783_Telos" "OG0008730_Telos" "OG0007581_Telos"
# [5] "OG0008302_Telos" "OG0007175_Telos" "OG0007697_Telos" "OG0007635_Telos"
# [9] "OG0005828_Telos" "OG0007167_Telos" "OG0006156_Telos" "OG0005251_Telos"
# [13] "OG0005235_Telos" "OG0005781_Telos" "OG0007086_Telos" "OG0006403_Telos"
# [17] "OG0006105_Telos" "OG0007489_Telos" "OG0007088_Telos" "OG0007808_Telos"
# [21] "OG0006717_Telos" "OG0007390_Telos" "OG0006762_Telos" "OG0007753_Telos"
# [25] "OG0007797_Telos" "OG0006334_Telos" "OG0006414_Telos" "OG0006740_Telos"
# [29] "OG0004013_Telos" "OG0008076_Telos" "OG0008436_Telos" "OG0005859_Telos"

#> Column H-M.dist.JC69 (index = 4)
which( ord.blengths.s2f.JC69[,4] == "Inf" ) 
# integer(0)
which( ord.blengths.s2f.JC69[,4] == "NaN" ) 
# integer(0)
which( ord.blengths.s2f.JC69[,4] >= 0.75 )
# [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
# [26] 26 27 28 29 30
ord.blengths.s2f.JC69[c( 1:30),]
rownames( ord.blengths.s2f.JC69 )[ c( 1:30)]
# [1] "OG0006811_Telos" "OG0006783_Telos" "OG0008730_Telos" "OG0007581_Telos"
# [5] "OG0008302_Telos" "OG0007175_Telos" "OG0007697_Telos" "OG0007635_Telos"
# [9] "OG0005828_Telos" "OG0007167_Telos" "OG0005251_Telos" "OG0005235_Telos"
# [13] "OG0006105_Telos" "OG0006403_Telos" "OG0006156_Telos" "OG0005781_Telos"
# [17] "OG0007086_Telos" "OG0007088_Telos" "OG0007808_Telos" "OG0007753_Telos"
# [21] "OG0006762_Telos" "OG0006334_Telos" "OG0007797_Telos" "OG0007489_Telos"
# [25] "OG0006717_Telos" "OG0006740_Telos" "OG0007390_Telos" "OG0004013_Telos"
# [29] "OG0006414_Telos" "OG0008076_Telos"

#> Column H-M.dist.raw (index = 5)
which( ord.blengths.s2f.JC69[,5] == "NaN" ) 
# integer(0)
which( ord.blengths.s2f.JC69[,5] == "Inf" )
# integer(0)
which( ord.blengths.s2f.JC69[,5] >= 0.75 )
# integer(0)

# 3. Use object `ord.blengths.s2f.raw`, which genes were ordered from slow-
#     to fast-evolving according to the model `raw`.
#> Column H-M.dist.TN93 (index = 3)
which( ord.blengths.s2f.raw[,3] == "Inf" )
# integer(0)
which( ord.blengths.s2f.raw[,3] == "NaN" ) 
# integer(0)
which( ord.blengths.s2f.raw[,3] >= 0.75 )
# [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
# [26] 26 27 28 29 30 31 32
ord.blengths.s2f.raw[c( 1:32),]
rownames( ord.blengths.s2f.raw )[ c( 1:32)]
# [5] "OG0008302_Telos" "OG0007175_Telos" "OG0007697_Telos" "OG0007635_Telos"
# [9] "OG0005828_Telos" "OG0007167_Telos" "OG0005251_Telos" "OG0005235_Telos"
# [13] "OG0006105_Telos" "OG0006403_Telos" "OG0006156_Telos" "OG0005781_Telos"
# [17] "OG0007086_Telos" "OG0007088_Telos" "OG0007808_Telos" "OG0007753_Telos"
# [21] "OG0006762_Telos" "OG0006334_Telos" "OG0007797_Telos" "OG0007489_Telos"
# [25] "OG0006717_Telos" "OG0006740_Telos" "OG0007390_Telos" "OG0004013_Telos"
# [29] "OG0006414_Telos" "OG0008076_Telos" "OG0008436_Telos" "OG0005859_Telos"

#> Column H-M.dist.JC69 (index = 4)
which( ord.blengths.s2f.raw[,4] == "NaN" ) #column H-M.dist.JC69
# integer(0)
which( ord.blengths.s2f.raw[,4] == "Inf" ) #column H-M.dist.JC69
# integer(0)
which( ord.blengths.s2f.raw[,4] >= 0.75 )
# [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
# [26] 26 27 28 29 30
ord.blengths.s2f.raw[c( 1:30),]
rownames( ord.blengths.s2f.raw )[ c( 1:30)]
# [1] "OG0006811_Telos" "OG0006783_Telos" "OG0008730_Telos" "OG0007581_Telos"
# [5] "OG0008302_Telos" "OG0007175_Telos" "OG0007697_Telos" "OG0007635_Telos"
# [9] "OG0005828_Telos" "OG0007167_Telos" "OG0005251_Telos" "OG0005235_Telos"
# [13] "OG0006105_Telos" "OG0006403_Telos" "OG0006156_Telos" "OG0005781_Telos"
# [17] "OG0007086_Telos" "OG0007088_Telos" "OG0007808_Telos" "OG0007753_Telos"
# [21] "OG0006762_Telos" "OG0006334_Telos" "OG0007797_Telos" "OG0007489_Telos"
# [25] "OG0006717_Telos" "OG0006740_Telos" "OG0007390_Telos" "OG0004013_Telos"
# [29] "OG0006414_Telos" "OG0008076_Telos"

#> Column H-M.dist.raw (index = 5)
which( ord.blengths.s2f.raw[,5] == "Inf" ) #column H-M.dist.raw
# integer(0)
which( ord.blengths.s2f.raw[,5] == "NaN" ) #column H-M.dist.raw
# integer(0)
which( ord.blengths.s2f.raw[,5] >= 0.75 )
# integer(0)

# >> END CHECK

# See the genes for which at least one of the models had issues to compute the distance
# between Xanpa and TelChr sequences or this distance was equal or higher than 0.75.
ord.blengths.s2f.TN93[c( 1:32),]

#                Sum.blengths Largbr.rel2tree H-M.dist.TN93 H-M.dist.JC69
#OG0006811_Telos     7.619425        7.966165     2.8432477     2.3250692
#OG0006783_Telos     8.769181       38.296427     2.1679567     1.7045491
#OG0008730_Telos     6.417053       11.893980     1.4681349     1.3383650
#OG0007581_Telos     9.094273       11.347955     1.2432654     1.1911204
#OG0008302_Telos    10.960978       11.756168     1.2402994     1.1075088
#OG0007175_Telos    10.314096        9.565802     1.1344589     1.0997528
#OG0007697_Telos     6.016370       14.150426     1.1177640     1.0316993
#OG0007635_Telos     8.308688        8.478571     1.0003461     0.9913169
#OG0005828_Telos    16.609983       18.179441     0.9600876     0.9430205
#OG0007167_Telos    13.489911        7.974760     0.9213760     0.8889349
#OG0006156_Telos    10.480141        6.044117     0.8889015     0.8351535
#OG0005251_Telos    18.342696        5.017403     0.8766294     0.8577995
#OG0005235_Telos    16.323069       15.753177     0.8737800     0.8392140
#OG0005781_Telos    13.122461        6.590151     0.8593503     0.8315350
#OG0007086_Telos    11.081498        7.642730     0.8525848     0.8152635
#OG0006403_Telos    12.466850       10.142626     0.8468812     0.8352376
#OG0006105_Telos    13.607668        5.604009     0.8442645     0.8384524
#OG0007489_Telos     9.865341        7.657759     0.8283420     0.7826644
#OG0007088_Telos     7.277304        8.189104     0.8276181     0.8133207
#OG0007808_Telos     5.427359       10.886326     0.8238612     0.8030276
#OG0006717_Telos    11.525118        7.584426     0.8189955     0.7792651
#OG0007390_Telos     9.773490        8.917940     0.8179443     0.7739407
#OG0006762_Telos     9.003578        7.973075     0.8165011     0.8001477
#OG0007753_Telos     7.421411       11.203584     0.8152945     0.8013444
#OG0007797_Telos     6.998134        9.554476     0.8129096     0.7831883
#OG0006334_Telos    14.780844        6.840389     0.8078377     0.7832516
#OG0006414_Telos    16.228352        6.047589     0.7973883     0.7644768
#OG0006740_Telos     9.576985        5.028743     0.7927488     0.7755553
#OG0004013_Telos    18.771031        6.952010     0.7786462     0.7678583
#OG0008076_Telos     6.389503       19.262375     0.7741748     0.7644041
#OG0008436_Telos    11.622895        9.317137     0.7553027     0.7437301
#OG0005859_Telos    18.860961        8.987654     0.7552263     0.7354044
#                H-M.dist.raw3
#OG0006811_Telos    0.7162162
#OG0006783_Telos    0.6848473
#OG0008730_Telos    0.6240876
#OG0007581_Telos    0.5967742
#OG0008302_Telos    0.5787037
#OG0007175_Telos    0.5769231
#OG0007697_Telos    0.5604839
#OG0007635_Telos    0.5500000
#OG0005828_Telos    0.5366984
#OG0007167_Telos    0.5207469
#OG0006156_Telos    0.5037037
#OG0005251_Telos    0.5110294
#OG0005235_Telos    0.5050336
#OG0005781_Telos    0.5025126
#OG0007086_Telos    0.4970845
#OG0006403_Telos    0.5037313
#OG0006105_Telos    0.5047847
#OG0007489_Telos    0.4858491
#OG0007088_Telos    0.4964286
#OG0007808_Telos    0.4929245
#OG0006717_Telos    0.4846491
#OG0007390_Telos    0.4827586
#OG0006762_Telos    0.4919355
#OG0007753_Telos    0.4923469
#OG0007797_Telos    0.4860335
#OG0006334_Telos    0.4860558
#OG0006414_Telos    0.4793651
#OG0006740_Telos    0.4833333
#OG0004013_Telos    0.4805825
#OG0008076_Telos    0.4793388
#OG0008436_Telos    0.4717742
#OG0005859_Telos    0.4686684
# Delete genes
ord.blengths.s2f.TN93.filt <- ord.blengths.s2f.TN93[ -c( 1:32),]
dim( ord.blengths.s2f.TN93.filt ) # [1] 1254     5

# Write table with the matrix
write.table( ord.blengths.s2f.TN93.filt,
             file = "out_RData/filtered_ordered_matrix.csv",
             quote = F, sep = "," )

## 7.1. Plot filtered and not filtered matrix with relative and absolute branch lengths
# Output a pdf with the plot
pdf( file = "out_RData/check_relblVStreelength.pdf", paper = "a4" )
plot( ord.blengths.s2f.TN93.filt[,1], ord.blengths.s2f.TN93.filt[,2],
      main = "Tree length VS Relative branch lengths",
      xlab = "Tree length", ylab = "Relative branch lengths" )
dev.off()
pdf( file = "out_RData/check_logtreelengthVSlargestbl.pdf", paper = "a4" )
plot( log(ord.blengths.s2f.TN93.filt[,1]), ord.blengths.s2f.TN93.filt[,2],
      main = "Quality control check of relative branch lengths",
      xlab = "Total tree length (log-scale)", ylab = "Relative branch lengths (%)" )
dev.off()

# There are no outliers

# Final size
dim( ord.blengths.s2f.TN93.filt )[1] # 1254

save( ord.blengths.s2f.TN93.filt,
      file = "out_RData/mammals_summary_matrix.TN93_filtstep3.RData" )


#--------------------------------------------#
# 8. Get a summary statistics of the filters #
#--------------------------------------------#
# Start counter with starting number of genes. Then proceed with the sum stats.
start.num.genes <- 1300
##============================================================================##
## UPDATE -- We added an extra filter to account for those genes that had to  ##
## be removed because they are present in the second data set too. This is    ##
## the new extra column. We know that there are a total of 15,268 genes that  ##
## are not overlapping with the second data set (explained in the README.md   ##
## file in the GitHub repository as bash scripts were used for that purpose). ##
## Therefore, the number in the last column is manualy added to `sum.stats`   ##
## as you can see in line 564 to compute the summary statistics.              ##
##============================================================================##
#sum.stats <- as.data.frame( matrix( 0, nrow = 5, ncol = 4 ) )
#colnames(sum.stats)   <- c( "Raw", "First_filter", "Second_filter", "Third filter" )
sum.stats <- as.data.frame( matrix( 0, nrow = 5, ncol = 5 ) )
colnames(sum.stats)   <- c( "Raw", "First_filter", "Second_filter", "Third filter", "Fourth filter" )
rownames( sum.stats ) <- c( "Num_genes", "Lost_genes", "%_Loss",
                            "Cum_lost_genes", "%_Cum_loss" )
sum.stats[1,] <- c( start.num.genes, dim( blengths)[1], dim(filt.df.blengths)[1],
                    dim(ord.blengths.s2f.TN93.filt)[1], 1300 )
sum.stats[2,] <- c( sum.stats[1,1]-sum.stats[1,1],
                    sum.stats[1,1]-sum.stats[1,2],
                    sum.stats[1,2]-sum.stats[1,3],
                    sum.stats[1,3]-sum.stats[1,4],
                    sum.stats[1,4]-sum.stats[1,5])
sum.stats[3,] <- c( round( (100 - (sum.stats[1,1]*100/sum.stats[1,1])),3 ),
                    round( (100 - (sum.stats[1,2]*100/sum.stats[1,1])),3 ),
                    round( (100 - (sum.stats[1,3]*100/sum.stats[1,2])),3 ),
                    round( (100 - (sum.stats[1,4]*100/sum.stats[1,3])),3 ),
                    round( (100 - (sum.stats[1,5]*100/sum.stats[1,4])),3 ))
sum.stats[4,] <- c( start.num.genes-sum.stats[1,1],
                    start.num.genes-sum.stats[1,2],
                    start.num.genes-sum.stats[1,3],
                    start.num.genes-sum.stats[1,4],
                    start.num.genes-sum.stats[1,5])
sum.stats[5,] <- c( round( (100 - (sum.stats[1,1]*100/start.num.genes)),3 ),
                    round( (100 - (sum.stats[1,2]*100/start.num.genes)),3 ),
                    round( (100 - (sum.stats[1,3]*100/start.num.genes)),3 ),
                    round( (100 - (sum.stats[1,4]*100/start.num.genes)),3 ),
                    round( (100 - (sum.stats[1,5]*100/start.num.genes)),3 ))

# Write summary table
#write.table(file = "out_RData/filtering_sum_stats.csv", sum.stats, quote = F, sep = ",")
write.table(file = "out_RData/filtering_sum_stats_nonucgenes.csv", sum.stats, quote = F, sep = ",")

#----------------------------------------------------------#
# 9. Copy filtered genes but saving them in an ordered way #
#    according to fast-slow evolving                       #
#----------------------------------------------------------#
# Create objects 
name_ord_filt_genes <- rownames( ord.blengths.s2f.TN93.filt )
name_notordered     <- names( XT_aln.aln )
list_correspondance <- rep(0, length( name_ord_filt_genes ) )
# Create a list with the indexes of the genes not ordered
# corresponding to the position in ordered list
for ( i in 1:length(name_ord_filt_genes) ){
  list_correspondance[i] <- which( name_notordered == name_ord_filt_genes[i] )
}
# >> CHECK | Double check everything alright !
all.equal( name_ord_filt_genes, name_notordered[list_correspondance] ) # TRUE!
# E.g. If i = 1
k <- 1
ind.k <- which( name_notordered == name_ord_filt_genes[k] )
all.equal( name_notordered[ind.k], name_ord_filt_genes[k] )
#[1] TRUE
# E.g.2 
name_notordered[list_correspondance[1:10]]
name_ord_filt_genes[1:10]
all.equal( name_notordered[list_correspondance[1:10]], name_ord_filt_genes[1:10] )
#[1] TRUE
# >> CHECK FINISHED 

# Generate the data
# Check if directory where to copy files exist. Otherwise, create it.
if ( ! dir.exists( "filtered_genes_step2/" ) ){
  dir.create( "filtered_genes_step2/" )
}
counter = 0 
for( i in seq( 1:length( name_ord_filt_genes ) ) ){
  
  counter = counter + 1
  # Check if directory where to copy files exist. Otherwise, create it.
  if ( ! dir.exists( paste( "filtered_genes_step2/", i, "/", sep = "") ) ){
    dir.create( paste( "filtered_genes_step2/", i, "/", sep = "" ) )
  }
  # Copy files
  file.copy( paste( "baseml/", list_correspondance[i], "/",
                    "partitions12_", name_notordered[ list_correspondance[i] ],
                    ".aln", sep = "" ),
             paste( "filtered_genes_step2/", i, "/", sep = "") )
  file.copy( paste( "baseml/", list_correspondance[i], "/",
                    name_notordered[ list_correspondance[i] ], ".tree", sep = "" ),
             paste( "filtered_genes_step2/", i, "/", sep = "") )
  cat( counter, "\n" )
  write( paste( "Gene ", counter, " is: ",
                name_notordered[ list_correspondance[i] ],
                ". Copy *aln and *tree from folder ",
                list_correspondance[i], " in baseml dir! ", "\n", sep = "" ),
         "out_logs/log_06_R_copy_step2_filtered_genes.txt",
         append = TRUE, sep = "\n" )
}

#-------------------------------------------------------------------------------#
# 10. Copy filtered genes which are present in all 84 species and in an ordered #
#     way according to slow-fast evolving                                       #
#-------------------------------------------------------------------------------#
## 10.1. Prase trees from "filtered_genes_all84sp" 
# Read trees from "filtered_genes_all84sp" and save them in a vector
trees84sp.dir <- c( "filtered_genes_all84sp/" )
num.genes84sp <- length( list.files( trees84sp.dir ) )

genes84sp.ind.ordered    <- which( name_ord_filt_genes %in% list.files( trees84sp.dir ) )
name_ord_filt_genes_84sp <- name_ord_filt_genes[ genes84sp.ind.ordered ]
list_correspondance_84sp <- rep(0, length( genes84sp.ind.ordered ) )
# Create a list with the corresponding indexes of the genes ordered
# but in the not ordered list.
# Now we have 645 genes as 3 might have been trimmed due to the previous
# filterings applied !
for ( i in 1:length( genes84sp.ind.ordered ) ){
  list_correspondance_84sp[i] <- which( name_notordered == name_ord_filt_genes_84sp[i] )
}

# >> CHECK 
k <- 1
which( name_notordered == name_ord_filt_genes_84sp[k] )
#[1] 20
name_notordered[20]
#[1] "OG0002170_Telos"
genes84sp.ind.ordered[1]
#[1] 1021
name_ord_filt_genes[ genes84sp.ind.ordered[1] ]
#[1] "OG0002170_Telos"
# >> END CHECK

# Get blengths ordered object with only genes present in 84 species
ord.blengths.s2f.TN93.84sp <- ord.blengths.s2f.TN93.filt[genes84sp.ind.ordered,]
dim( ord.blengths.s2f.TN93.84sp ) # [1] 4   5
# Write table with the matrix
write.table(file = "out_RData/filtered_ordered_matrix_all84sp.csv",
            ord.blengths.s2f.TN93.84sp, quote = F, sep = ",")

## 10.2. Generate the data
# Start counter
counter <- 0 
# Check if directory where to copy files exist. Otherwise, create it.
if ( ! dir.exists( "filtered_genes_step2_all84sp/" ) ){
  dir.create( "filtered_genes_step2_all84sp/"  )
}
for( i in seq( 1:length( name_ord_filt_genes_84sp ) ) ){
  
  counter = counter + 1
  # Check if directory where to copy files exist. Otherwise, create it.
  if ( ! dir.exists( paste( "filtered_genes_step2_all84sp/", i, "/", sep = "") ) ){
    dir.create( paste( "filtered_genes_step2_all84sp/", i, "/", sep = "" ) )
  }
  # Copy files
  file.copy( paste( "baseml/", list_correspondance_84sp[i], "/",
                    "partitions12_", name_notordered[ list_correspondance_84sp[i] ],
                    ".aln", sep = "" ),
             paste( "filtered_genes_step2_all84sp/", i, "/", sep = "") )
  file.copy( paste( "baseml/", list_correspondance_84sp[i], "/",
                    name_notordered[ list_correspondance_84sp[i] ], ".tree", sep = "" ),
             paste( "filtered_genes_step2_all84sp/", i, "/", sep = "") )
  cat( counter, "\n" )
  write( paste( "Gene ", counter, " is: ",
                name_notordered[ list_correspondance_84sp[i] ],
                ". Copy *aln and *tree from folder ",
                list_correspondance_84sp[i], " in baseml dir! ", "\n", sep = "" ),
         "out_logs/log_06_R_copy_step2_filtered_all84sp_genes.txt",
         append = TRUE, sep = "\n" )
}
