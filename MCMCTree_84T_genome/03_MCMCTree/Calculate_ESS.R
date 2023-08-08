
#This script was written by Sandra Alvarez-Carretero and adapted for the Teloschistales dataset by Theo Llewellyn

#--------------#
# LOAD PACKAGE #
#--------------#
library( MCMCtreeR )
library( rstudioapi )
library( rstan )

source( "/data/projects/gaya_lab/Teloschistaceae/LICHENS_MCMCTree/00_Esters2a_calibration/Functions_plots_MCMCtreeR.R" )

#-----------------------#
# SET WORKING DIRECTORY #
#-----------------------#
#run within LICHENS_MCMCTree/00_Esters2a_calibration/01_MCMCtree_posterior
setwd("data/projects/gaya_lab/Teloschistaceae/LICHENS_MCMCTree/00_Esters2a_calibration/01_MCMCtree_posterior")

#-------------#
# DEFINE FUNS #
#-------------#
# Function to convert the matrix with the sampled values
# for each parameter in the MCMC into a 3D array readable 
# by `rstan::monitor`.
# 
# Parameters:
# x  Data frame, there are many rows as iterations and as many 
#    columns as parameters that had to be inferred during the 
#    MCMC
#
mat2arr <- function( x ){
  
  # Get rid of "Gen" column:
  ## Not needed anymore
  #x <- x[,-c(1)]
  
  # Get array format
  arr <- array( 0, dim = c( dim( x )[1], 1, dim( x )[2] ),
                dimnames = list( paste( "Iter", 1:dim( x )[1], sep = "" ),
                                 "1 chain",
                                 #paste( "Param", 1:dim( x )[2], sep = "" ) )
                                 colnames( x ) )
  )
  
  for( p in 1:dim( x )[2] ){
    arr[,1,p] <- x[,p]
  }
  
  # Return object
  return( arr )
  
}

# Function to run `rstan::monitor` and export the 
# median, min, and max values for the tail-ESS and the 
# bulk-ESS
#
# Parameters:
# x        Data frame, there are many rows as iterations and as many 
#          columns as parameters that had to be inferred during the 
#          MCMC
#
# coda_fun Boolean, TRUE if you want to also compute the ESS
#          with the coda::effectiveSize function
sum_MCMC_ESS <- function( x, coda_fun = FALSE ){
  
  # Define empty vector to report stats
  if( coda_fun == FALSE ){
    out_stats <- matrix( 0, nrow = 3, ncol = 2 )
    colnames( out_stats ) <- c( "Tail-ESS", "Bulk-ESS" )
  } else{ 
    out_stats <- matrix( 0, nrow = 3, ncol = 3 )
    colnames( out_stats ) <- c( "Tail-ESS", "Bulk-ESS", "coda-ESS" )
  }
  rownames( out_stats ) <- c( "Med.", "Min.", "Max." )
  
  # Compute stats
  stats_mcmc     <- rstan::monitor( sims = mat2arr( x ) )
  out_stats[1,1] <- median( stats_mcmc$Tail_ESS )
  out_stats[2,1] <- min( stats_mcmc$Tail_ESS )
  out_stats[3,1] <- max( stats_mcmc$Tail_ESS )
  out_stats[1,2] <- median( stats_mcmc$Bulk_ESS )
  out_stats[2,2] <- min( stats_mcmc$Bulk_ESS )
  out_stats[3,2] <- max( stats_mcmc$Bulk_ESS )
  
  if( coda_fun == TRUE ){
    ESS_coda       <- coda::effectiveSize( x = x )
    out_stats[1,3] <- median( ESS_coda )
    out_stats[2,3] <- min( ESS_coda )
    out_stats[3,3] <- max( ESS_coda )
  }
  
  # Return stats 
  return( list( tab = out_stats, stats = stats_mcmc, stats_CODA = ESS_coda ) )
  
}


# Function to load the two half of the subtrees.
# Used to get mean estimates for convergence plots
# 
# Parameters:
# mcmc    Character, path to mcmc.txt file
# subt    Character, name of the subtree
# delcol  Numeric, number of columns that are to be deleted 
#         as they do not contain divtimes. Default = 10.
load_subt <- function( mcmc1, subt, delcol = 10, perc = 0.975 ){
  
  # 1. Load files and get parameters
  cat( "Load combined mcmc.txt from ", subt, " ... ...\n" )
  run1    <- read.table( mcmc1, header = TRUE, sep = "\t" )
  
  # 2. Summarise parameters
  cat( "Generate objects with summarised estimated divergence times... ...\n")
  dim.r1   <- dim(run1)[2]
  # divtimes <- run1[,-c(1, dim(run1)[2])]
  divtimes <- run1[,-c( 1, ( dim.r1-delcol ):dim.r1 )]
  
  mean_est_divt  <- apply( X = divtimes, MARGIN = 2, FUN = mean )
  quant_est_divt <- apply( X = divtimes, MARGIN = 2, FUN = quantile, probs = c( 1-perc,perc ) )
  quant_est_divt <- t( quant_est_divt )
  test_alleq     <- all.equal( quant_est_divt[1,], quantile( divtimes[,1], probs = c( 1-perc,perc ) ) )
  if( test_alleq != TRUE ){
    stop( "There was an issue calculating quantiles!" )
  }
  
  # 3. Return object 
  cat( "\nTasks done! Return objects\n\n")
  return( list( divt = divtimes, mean_divt = mean_est_divt, quant_divt = quant_est_divt ) )
  
}


# Function to find problematic runs.
# It prints out the mean div time, qlow, and qup 
# for each node and for each run
# 
# Parameters:
# num_dirs       Numeric. Default is c(3,4,7,8,9,11,13,16). Change if needed.
# delcol         Numeric. Default is 8 so only the samples for time estimates are taken.
# name_dir_subt  Character. Name of the directory where the anlayses for this subtree ran.
# num_divt       Numeric. Number of columns of the mcmc.txt that correspond to the samples 
#                collected for the times.
# node_calib     Character. CSV file with two columns, where the first column has the names 
#                of the calibrations and in the second the corresponding node. If not available,
#                set to NULL.
# path_84sp      Character. Path to the directory that contains the 10 runs.
# tree_hyp       Character. Name of the tree hypothesis used to name the directory.
# clock          Character. GBM or ILN, default is "GBM".
# perc           Numeric. Default 0.975, change if you want different quantiles.
# out            Character. Path to where you want the output files saved. Default is 
#                current directory, change if needed.
# maintt         Boolean. TRUE if the main tree is being analysed, which is in a different 
#                directory. FALSE otherwise for the rest of tree hypotheses.
find_prob_MCMC_84sp <- function ( num_dirs = c(3,4,7,8,9,11,13,16), delcol = 8, name_dir_subt, num_divt, node_calib,
                                  path_84sp, tree_hyp, clock = "GBM", out = NULL, perc = 0.975, maintt = TRUE ){
  
  # 0. Allow for numbers not using exponentials if too low !
  options( scipen = 999 )
  
  # 1. Create global vars
  total_runs         <- length( num_dirs )
  subtree_list       <- vector( mode = "list", length = total_runs )
  subtree_meandivt   <- matrix( 0, nrow = total_runs, ncol = num_divt )
  subtree_qup        <- matrix( 0, nrow = total_runs, ncol = num_divt )
  subtree_qlow       <- matrix( 0, nrow = total_runs, ncol = num_divt )
  names( subtree_list )        <- rownames( subtree_meandivt ) <-
    rownames( subtree_qup ) <- rownames( subtree_qlow ) <- paste( "run", num_dirs, sep = "" )
  
  # 2. Get summarised data
  count     <- 0
  count_dim <- vector( "numeric", total_runs )
  for( i in num_dirs ){
      count <- count + 1
      cat( "Loading run ", i, "\n" )
      
      if ( maintt == TRUE ){
        subtree_list[[ count ]]  <- load_subt( mcmc1 = paste( path_84sp, "run", i, "/mcmctree_", clock, "/mcmc.txt", sep = "" ),
                                               subt = paste( name_dir_subt, "_run", i, sep = "" ),
                                               delcol = delcol, perc = perc )
      }else if( maintt == FALSE ){
        subtree_list[[ count ]]  <- load_subt( mcmc1 = paste( path_84sp, "mcmc", i, "/mcmctree_", clock, "/", tree_hyp, "/mcmc.txt", sep = "" ),
                                               subt = paste( name_dir_subt, "_run", i, sep = "" ),
                                               delcol = delcol, perc = perc )
      }
      
      count_dim[count]             <- dim( subtree_list[[ count ]]$divt )[1]
      subtree_meandivt[count,] <- subtree_list[[ count ]]$mean_divt
      subtree_qlow[count,]     <- subtree_list[[ count ]]$quant_divt[,1]
      subtree_qup[count,]      <- subtree_list[[ count ]]$quant_divt[,2]
  }
  colnames( subtree_meandivt ) <- colnames( subtree_qup ) <-
    colnames( subtree_qlow ) <- rownames( subtree_list[[ 1 ]]$quant_divt )
  
  # 3. Get mean data across runs. Apparently it is more accurate 
  #    than actually using meandivt per run and it matches metrics 
  #    calculated by MCMCtree :|
  mcmc_all <- matrix( 0, nrow = sum( count_dim ), ncol = num_divt )
  colnames( mcmc_all ) <- rownames( subtree_list[[ 1 ]]$quant_divt )
  start <- 1
  stop  <- 0
  cat( "Calculating mean and quantiles for all samples... ...\n\n")
  for( i in 1:length(num_dirs) ){
    # cat( "Start: ", start, "\n")
    stop <- stop + count_dim[i]
    # cat( "Stop: ", stop, "\n" )
    mcmc_all[c(start:stop),] <- matrix( unlist( subtree_list[[ i ]]$divt ), nrow = num_divt, byrow = FALSE ) 
    start <- stop + 1
  }
  
  mean_divt    <- apply( X = mcmc_all, MARGIN = 2, FUN = mean )
  mean_quants  <- apply( X = mcmc_all, MARGIN = 2, FUN = quantile, probs = c( 1-perc,perc ) )
  mean_est         <- matrix( 0, nrow = length( mean_divt ), ncol = 3 )
  mean_est_priors  <- matrix( 0, nrow = length( mean_divt ), ncol = 4 )
  mean_est[,1]     <- mean_est_priors[,1] <- round( mean_divt*100, digits = 3 )
  mean_est[,2]     <- mean_est_priors[,2] <- round( as.numeric( format( mean_quants[1,], scientific = FALSE ) )*100, digits = 2 )
  mean_est[,3]     <- mean_est_priors[,3] <- round( as.numeric( format( mean_quants[2,], scientific = FALSE ) )*100, digits = 2 )
  colnames( mean_est )        <- c( "Mean_time", "Mean_qlow", "Mean_qup" )
  colnames( mean_est_priors ) <- c( "Mean_time", "Mean_qlow", "Mean_qup", "Priors" )
  test_names <- all.equal( names( mean_divt ), colnames( mean_quants ) )
  if( test_names != TRUE ){
    stop( "Issue with names for mean divt, mean qup, and mean q low!" )
  }
  rownames( mean_est ) <- rownames( mean_est_priors ) <- names( mean_divt )
  
  # 4. Get matching nodes with calib names
  if( length( node_calib ) > 0 ) {
    match_csv <- read.table( node_calib, header = TRUE, sep = ";", stringsAsFactors = FALSE )
    ind_match <- which( as.numeric( gsub( pattern = "t_n", replacement = "", x = rownames( mean_est_priors ) ) )
                        %in% match_csv[,2] )
    for( i in ind_match ){
      tmp_ind <- which( paste( "t_n", match_csv[,2], sep = "" ) == rownames( mean_est_priors )[i] )
      rownames( mean_est_priors )[i] <- paste( "t_n", match_csv[tmp_ind,2], "_", match_csv[tmp_ind,1], sep = "" )
      mean_est_priors[i,4] <- match_csv[tmp_ind,3]
    }
  }
  
  # 5. Write separate output files for mean times, mean qup, and mean qlow for each 
  #    run
  write.table( x = subtree_meandivt, file = paste( out, "mean_divt_84sp.tsv", sep = "" ), sep = "\t", quote = FALSE )
  write.table( x = subtree_qup, file = paste( out, "mean_qup_7sp.tsv", sep = "" ), sep = "\t", quote = FALSE )
  write.table( x = subtree_qlow, file = paste( out, "mean_qlow_84sp.tsv", sep = "" ), sep = "\t", quote = FALSE )
  write.table( x = mean_est_priors, file = paste( out, "all_mean_est_84sp.tsv", sep = "" ), sep = "\t", quote = FALSE )
  
  cat( "Output files available! Check returned list with all objects generated too :) \n\n")
  return( list( tt_all = subtree_list, mean = subtree_meandivt, qup = subtree_qup, qdown = subtree_qlow,
                all_mean_est = mean_est, all_mcmc = mcmc_all ) )
  
}

# Plot convergence plot 
# Function to plot a convergence plot 
# Arguments:
#
# name_dir_subt  Character. Name of the directory where the anlayses for this subtree ran.
# mean_divt1     Character. Name of the object that contains the min divtimes for the
#                first half of runs.
# mean_divt2     Character. Name of the object that contains the min divtimes for the
#                second half of runs.
# num_runs       Integer. Number of runs.
plot_convergence <- function ( name_dir_subt, mean_divt1, mean_divt2, num_runs ){
  half <- round( num_runs/2, digits = 0 )
  tmp <- bquote( paste(  hat( italic(t) ), " | MCMC - run 1-", .(half), sep = "" ) )
  tmp2 <- bquote( paste( hat( italic(t) ), " | MCMC - run ", .(half+1), "-", .(num_runs), sep = "" ) )
  plot( x = mean_divt1, y = mean_divt2,
        main = paste( "Convergence plot - ", name_dir_subt, sep = "" ),
        xlab = tmp, ylab = tmp2 )
  abline( lm( mean_divt2~0 + mean_divt1 ),
          col="red", lty = 2 )
}


#-----------#
# LOAD DATA #
#-----------#

#-- MAIN TREE (T2) --#

#\\\\\\\\\\\\\\\\\#
# 84sp tree - GBM #
#-----------------#
if( ! dir.exists( paste("out_data/00_post_84sp_GBM", sep = "" ) ) ){
  dir.create( paste("out_data/00_post_84sp_GBM", sep = "" ) )
}
softbounds_84sp_sum  <- find_prob_MCMC_84sp( num_dirs = c(1:16), delcol = 8, name_dir_subt = "84sp",
                                             num_divt = 83, node_calib = "Calibs_nodes_84sp_formatted.csv", 
                                             tree_hyp = "Telos85T_4parts_partitions12.raxml.support_rooted_labelled_MCMCtree_calib_1a.tree", maintt = TRUE,
                                             clock = "GBM", out = "out_data/00_post_84sp_GBM/",
                                             path_84sp = paste("01_MCMCtree_posterior/",
                                                                sep = "" ),
                                             perc = 0.975 )
#> CHECK: Issues with quantiles
for ( j in 2:16 ){
  tmp.qup     <- softbounds_84sp_sum$qup[1,] - softbounds_84sp_sum$qup[j,]
  tmp.ind.qup <- which( abs( tmp.qup ) > 0.1 )
  if( length( tmp.ind.qup ) > 0 ){
    cat( "q97.5%: Check the following nodes in chain ", j, ":", names( tmp.qup[tmp.ind.qup] ), "\n" )
    cat(  "Difference: ", tmp.qup[tmp.ind.qup], "\n")
  }
  tmp.qdown    <- softbounds_84sp_sum$qdown[1,] - softbounds_84sp_sum$qdown[j,]
  tmp.ind.qdown<- which( abs( tmp.qdown ) > 0.1 )
  if( length( tmp.ind.qdown ) > 0 ){
    cat( "q2.5%: Check the following nodes in chain ", j, ":", names( tmp.qdown[tmp.ind.qdown] ), "\n" )
    cat(  "Difference: ", tmp.qdown[tmp.ind.qdown], "\n")
  }
}
#>END CHECK -- No issues!
sp84_filt_half1 <- apply( X = softbounds_84sp_sum$mean[1:8,], MARGIN = 2, FUN = mean )
sp84_filt_half2 <- apply( X = softbounds_84sp_sum$mean[9:16,], MARGIN = 2, FUN = mean )
pdf( "RPlots/Convergence_plot_84sp_GBM.pdf", paper = "a4" )
plot_convergence( name_dir_subt = "84sp", mean_divt1 = sp84_filt_half1,
                  mean_divt2 = sp84_filt_half2, num_runs = 16 )
dev.off()


#\\\\\\\\\\\\\\\\\#
# 84sp tree - ILN #
#-----------------#
if( ! dir.exists( paste("out_data/00_post_84sp_ILN", sep = "" ) ) ){
  dir.create( paste("out_data/00_post_84sp_ILN", sep = "" ) )
}
softbounds_84sp_sum_ILN  <- find_prob_MCMC_84sp( num_dirs = c(1:16), delcol = 8, name_dir_subt = "84sp",
                                             num_divt = 83, node_calib = "Calibs_nodes_84sp_formatted.csv", 
                                             tree_hyp = "Telos85T_4parts_partitions12.raxml.support_rooted_labelled_MCMCtree_calib_1a.tree", maintt = TRUE,
                                             clock = "ILN", out = "out_data/00_post_84sp_ILN/",
                                             path_84sp = paste("01_MCMCtree_posterior/",
                                                                sep = "" ),
                                             perc = 0.975 )
#> CHECK: Issues with quantiles
for ( j in 2:16 ){
  tmp.qup     <- softbounds_84sp_sum_ILN$qup[1,] - softbounds_84sp_sum_ILN$qup[j,]
  tmp.ind.qup <- which( abs( tmp.qup ) > 0.1 )
  if( length( tmp.ind.qup ) > 0 ){
    cat( "q97.5%: Check the following nodes in chain ", j, ":", names( tmp.qup[tmp.ind.qup] ), "\n" )
    cat(  "Difference: ", tmp.qup[tmp.ind.qup], "\n")
  }
  tmp.qdown    <- softbounds_84sp_sum_ILN$qdown[1,] - softbounds_84sp_sum_ILN$qdown[j,]
  tmp.ind.qdown<- which( abs( tmp.qdown ) > 0.1 )
  if( length( tmp.ind.qdown ) > 0 ){
    cat( "q2.5%: Check the following nodes in chain ", j, ":", names( tmp.qdown[tmp.ind.qdown] ), "\n" )
    cat(  "Difference: ", tmp.qdown[tmp.ind.qdown], "\n")
  }
}
#>END CHECK -- No issues!
sp84_filt_half1_ILN <- apply( X = softbounds_84sp_sum_ILN$mean[1:8,], MARGIN = 2, FUN = mean )
sp84_filt_half2_ILN <- apply( X = softbounds_84sp_sum_ILN$mean[9:16,], MARGIN = 2, FUN = mean )
pdf( "RPlots/Convergence_plot_84sp_ILN.pdf", paper = "a4" )
plot_convergence( name_dir_subt = "84sp", mean_divt1 = sp84_filt_half1_ILN,
                  mean_divt2 = sp84_filt_half2_ILN, num_runs = 16 )
dev.off()




#---------------#
# CALCULATE ESS #
#---------------#
#> ESS with RStan
# Each column is assumed to be an MCMC. Rows are iterations for parameter X
# Source explaining why it is preferable than the function in coda:
# https://nature.berkeley.edu/~pdevalpine/MCMC_comparisons/nimble_MCMC_comparisons.html

#\\\\\\\\\#
# MAIN T2 #
#---------#
ESS_s84sp        <- sum_MCMC_ESS( x = softbounds_84sp_sum$all_mcmc, coda_fun = TRUE )
dim(softbounds_84sp_sum$all_mcmc)[1]
# length = 320016
min(ESS_s84sp$stats$Rhat)
# min(Rhat) 0.9999938
ESS_s84sp$tab
#     Tail-ESS Bulk-ESS  coda-ESS
#Med.     8939     2763  7468.720
#Min.     1823     1166  1802.022
#Max.    81519    42072 85057.392


ESS_s84sp_ILN        <- sum_MCMC_ESS( x = softbounds_84sp_sum_ILN$all_mcmc, coda_fun = TRUE )
dim(softbounds_84sp_sum_ILN$all_mcmc)[1]
# length = 320016
min(ESS_s84sp_ILN$stats$Rhat)
# min(Rhat) 0.9999938
ESS_s84sp_ILN$tab
#     Tail-ESS Bulk-ESS  coda-ESS
#Med.     8939     2763  7468.720
#Min.     1823     1166  1802.022
#Max.    81519    42072 85057.392

