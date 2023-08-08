
#This script was written by Sandra Alvarez-Carretero and adapted for the Teloschistales dataset by Theo Llewellyn

#---------------------#
#  CLEAN ENVIRONMENT  #
#---------------------#
rm( list = ls( ) )

#-------------------------#
#  SET WORKING DIRECTORY  #
#-------------------------#
setwd("/data/projects/gaya_lab/Teloschistaceae/LICHENS_MCMCTree/00_Esters2a_calibration/03_Fit_ST_to_posteriors")
# Load package to help find the path to this source file 
library(rstudioapi) 
# Set name for output tree
tname <- "LICHENS_84sp_Esters1a"

#------------#
# SCAN TREES #
#------------#
# Get trees 
tree_original     <- readLines( con = "01_trees_with_ST_calibs/Rinp/FigTree_84sp_nodelabels.tree" )
tree              <- tree_2              <- tree_original
tree_7ST         <- tree_7ST_2         <- tree_original
tree_7ST_rounded <- tree_7ST_rounded_2 <- tree_original

#------------------------#
# GET FILE WITH node.cal #
#------------------------#
# Read tree and rename cols
node.cal <- read.table( file = "01_trees_with_ST_calibs/Rinp/Node_calibrations_lichens_84sp.txt",
                        stringsAsFactors = F, sep = "|", blank.lines.skip = T, na.strings = "NAN" )
colnames( node.cal ) <- c( "name", "Uniform.cal", "Node.calib" )

ST.fit <- read.table( file = "00_fitST/Rout/ST.fitted.dists.G2.40.tsv",
                      header = TRUE, stringsAsFactors = F, sep = "\t" )
rownames( ST.fit ) <- as.numeric( gsub( pattern =  "t_n", replacement = "",
                                        x = row.names( ST.fit ) ) )

#-------------------------------#
#  REPLACE NAMES WITH node.cal  #
#-------------------------------#

# USE 7 CALIBRATIONS ONLY 
# Replace calibration names with corresponding calibration
node.labels.ST <- sort( as.numeric( rownames( ST.fit ) ), decreasing = TRUE )
for( i in node.labels.ST ){
  
  # Find calibration and replace in tree with ST-fit
  index.ST <- which( i == as.numeric( rownames( ST.fit ) ) )
  tree_7ST     <- gsub( pattern = paste0(" ",i," "),
                    x = tree_7ST, replacement = paste( " '", ST.fit[index.ST,5], "' ", sep = "" ) )
  
  # Last step -- update sp labels
  if( i == min( node.cal$Node.calib ) ){
    num_sp <- 84:1
    for( j in num_sp ){
      tree_7ST <- gsub( pattern = paste0( j,"_" ), x = tree_7ST, replacement = "" )
    }
  }
  
}


# USE 7 CALIBRATIONS ONLY - rounded
# Replace calibration names with corresponding calibration
node.labels.ST <- sort( as.numeric( rownames( ST.fit ) ), decreasing = TRUE )
for( i in node.labels.ST ){
  
  # Find calibration and replace in tree with ST-fit
  index.ST <- which( i == as.numeric( rownames( ST.fit ) ) )
  tree_7ST_rounded <- gsub( pattern = paste0(" ",i," "),
                             x = tree_7ST_rounded, replacement = paste( " '", ST.fit[index.ST,6], "' ", sep = "" ) )
  
  # Last step -- update sp labels
  if( i == min( node.cal$Node.calib ) ){
    num_sp <- 84:1
    for( j in num_sp ){
      tree_7ST_rounded <- gsub( pattern = paste0( j,"_" ), x = tree_7ST_rounded, replacement = "" )
    }
  }
  
}

#-----------------#
# WRITE NEW FILES #
#-----------------#
num_sp <- c( stringr::str_count( string = tree_original, pattern = "," ) + 1 )
phylip_header <- paste( num_sp, "  1", sep = "" )


write( x = phylip_header, file = paste( "01_trees_with_ST_calibs/Rout/",
                                        num_sp, "sp_", tname, "_7STcalib.tree", sep = "" ) )
write( x = tree_7ST, file = paste( "01_trees_with_ST_calibs/Rout/",
                                    num_sp, "sp_", tname, "_7STcalib.tree", sep = "" ),
       append = TRUE )

write( x = phylip_header, file = paste( "01_trees_with_ST_calibs/Rout/",
                                        num_sp, "sp_", tname, "_7STcalib_rounded.tree", sep = "" ) )
write( x = tree_7ST_rounded, file = paste( "01_trees_with_ST_calibs/Rout/",
                                            num_sp, "sp_", tname, "_7STcalib_rounded.tree", sep = "" ),
       append = TRUE )



## --- EXTRA 210812 | Only get tree with specific calibs --- ##

# USE 26 CALIBRATIONS ONLY 
# Replace calibration names with corresponding calibration
node.labels.ST.2 <- sort( as.numeric( rownames( ST.fit ) ), decreasing = TRUE )
check.indexes    <- sort( as.numeric( node.cal$Node.calib ), decreasing = TRUE )
for( i in node.labels.ST.2 ){
  
  # Find calibration and replace in tree with ST-fit
  if( i %in% check.indexes ){
    index.ST    <- which( i == as.numeric( rownames( ST.fit ) ) )
    tree_7ST_2 <- gsub( pattern = paste0(" ",i," "),
                         x = tree_7ST_2, replacement = paste( " '", ST.fit[index.ST,5], "' ", sep = "" ) )
  }else{ 
    tree_7ST_2 <- gsub( pattern = paste0(" ",i," "), x = tree_7ST_2, replacement = "" )
  }
  
  
  # Last step -- update sp labels
  if( i == min( node.cal$Node.calib ) ){
    num_sp <- 84:1
    for( j in num_sp ){
      tree_7ST_2 <- gsub( pattern = paste0( j,"_" ), x = tree_7ST_2, replacement = "" )
    }
  }
  
}


# USE 7 CALIBRATIONS ONLY - rounded
# Replace calibration names with corresponding calibration
node.labels.ST.2 <- sort( as.numeric( rownames( ST.fit ) ), decreasing = TRUE )
check.indexes    <- sort( as.numeric( node.cal$Node.calib ), decreasing = TRUE )
for( i in node.labels.ST.2 ){
  
  # Find calibration and replace in tree with ST-fit
  if( i %in% check.indexes ){
    index.ST <- which( i == as.numeric( rownames( ST.fit ) ) )
    tree_7ST_rounded_2 <- gsub( pattern = paste0(" ",i," "),
                               x = tree_7ST_rounded_2, replacement = paste( " '", ST.fit[index.ST,6], "' ", sep = "" ) )
  }else{
    tree_7ST_rounded_2 <- gsub( pattern = paste0(" ",i," "), x = tree_7ST_rounded_2, replacement = "" )
  }
  
  # Last step -- update sp labels
  if( i == min( node.cal$Node.calib ) ){
    num_sp <- 84:1
    for( j in num_sp ){
      tree_7ST_rounded_2 <- gsub( pattern = paste0( j,"_" ), x = tree_7ST_rounded_2, replacement = "" )
    }
  }
  
}

#-----------------#
# WRITE NEW FILES #
#-----------------#
num_sp <- c( stringr::str_count( string = tree_original, pattern = "," ) + 1 )
phylip_header <- paste( num_sp, "  1", sep = "" )


write( x = phylip_header, file = paste( "01_trees_with_ST_calibs/Rout/",
                                        num_sp, "sp_", tname, "_7STcalib_v2.tree", sep = "" ) )
write( x = tree_7ST_2, file = paste( "01_trees_with_ST_calibs/Rout/",
                                    num_sp, "sp_", tname, "_7STcalib_v2.tree", sep = "" ),
       append = TRUE )

write( x = phylip_header, file = paste( "01_trees_with_ST_calibs/Rout/",
                                        num_sp, "sp_", tname, "_7STcalib_rounded_v2.tree", sep = "" ) )
write( x = tree_7ST_rounded_2, file = paste( "01_trees_with_ST_calibs/Rout/",
                                            num_sp, "sp_", tname, "_7STcalib_rounded_v2.tree", sep = "" ),
       append = TRUE )
