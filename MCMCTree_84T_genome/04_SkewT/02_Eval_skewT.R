
#This script was written by Sandra Alvarez-Carretero and adapted for the Teloschistales dataset by Theo Llewellyn

#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )

#-----------------------#
# SET WORKING DIRECTORY #
#-----------------------#
library( rstudioapi ) 
setwd("/data/projects/gaya_lab/Teloschistaceae/LICHENS_MCMCTree/00_Esters2a_calibration/03_Fit_ST_to_posteriors")

#--------------#
# LOAD LIBRARY #
#--------------#
library( sn )

#-----------#
# LOAD DATA #
#-----------#
# Load objects previously created with the fitted ST dists and 
# the posterior divtimes
load( "00_fitST/Rdata/ST.fitted.objects.RData" )
load( "00_fitST/Rdata/ST.fitted.dists.RData" )
load( "00_fitST/Rdata/post.divtimes.RData" )

# Load prior from MCMC when using ST dists as calibs
prior.ST.run1 <- read.table( "02_MCMCtree_prior_ST/MCMCtree/mcmc1/mcmc.txt",
                             header = T, sep = "\t" )
prior.ST.run2 <- read.table( "02_MCMCtree_prior_ST/MCMCtree/mcmc2/mcmc.txt",
                             header = T, sep = "\t" )
# Get parameters and divtimes
divtimes.prior.ST.run1 <- prior.ST.run1[,-c(1,(dim(prior.ST.run1)[2]) )]
divtimes.prior.ST.run2 <- prior.ST.run2[,-c(1,(dim(prior.ST.run2)[2]) )]
divtimes.prior.ST  <- rbind( divtimes.prior.ST.run1, divtimes.prior.ST.run2 )

mean_est_divt_prior.ST  <- apply( X = divtimes.prior.ST, MARGIN = 2, FUN = mean )
quant_est_divt_prior.ST <- apply( X = divtimes.prior.ST, MARGIN = 2, FUN = quantile, probs = c(0.025,0.975) )
quant_est_divt_prior.ST <- t( quant_est_divt_prior.ST )
all.equal( quant_est_divt_prior.ST[1,], quantile( divtimes.prior.ST[,1], probs = c( 0.025, 0.975 ) ) )

#-------#
# PLOTS #
#-------#

# Plot analytical ST VS prior with analytical ST VS Posterior
pdf( file = "02_MCMCtree_prior_ST/plots/00_Check_analitST_priorST_posterior_SMALL.pdf","a4")
par( mfrow = c(8,9), mai=c(0.1,0.1,0.1,0.1))
for( i in 1:length( colnames( divtimes.prior.ST ) ) ){
  
  # 1. Get a tmp_node with ST distributions
  tmp_node <- ST.fitted.objects[[ i ]]
  #tracemem( tmp_node )
  # 2. Plot fitted ST for each node previously 
  #    computed using the values sampled during the MCMC
  
  # 2.1. Find limit axis
  max_x_st            <- round( max( density( rst( n = 1000, dp = tmp_node$dp.complete ) )$x ) + 0.5 )
  max_x_chain         <- round( max( density( divtimes[,i] )$x ) + 0.5 )
  max_x_chain_priorST <- round( max( density( divtimes.prior.ST[,i] )$x ) + 0.5 )
  x_lim       <- max( max_x_chain, max_x_st, max_x_chain_priorST )
  max_y_st            <- round( max( density( rst( n = 1000, dp = tmp_node$dp.complete ) )$y ) + 0.5 )
  max_y_chain         <- round( max( density( divtimes[,i] )$y ) + 0.5 )
  max_y_chain_priorST <- round( max( density( divtimes.prior.ST[,i] )$y ) + 0.5 )
  y_lim       <- max( max_y_chain, max_y_st, max_y_chain_priorST )
  # 2.2. Plot
  plot( density( divtimes[,i], adj = 1 ),
        xlim = c( 0, x_lim ), ylim = c( 0, y_lim ), #main = colnames( divtimes )[i], cex.main = 0.6,
        main=NULL, xlab = '', ylab = '', cex.axis = 0.7, mgp = c(2.5,0,0),
        #xaxt = "n", yaxt = "n", 
        tck = 0.01, col.main = "white" )
  title( colnames( divtimes )[i], line = -2, cex.main = 0.8, adj = 0.9 )
  lines( density( divtimes.prior.ST[,i], adj = 1 ),
         col = "green" )
  curve( dst( x, xi = tmp_node$dp[1], omega = tmp_node$dp[2],
              alpha = tmp_node$dp[3], nu = tmp_node$dp[4] ),
         from = 0, to = x_lim,
         n = 1e4, add = TRUE, col = "red" )

}

plot( x = 1, y = 1, col = "white", xaxt = "n", yaxt = "n", xlab = '', ylab = '' )
info.legend <- c( "Post. ages with SB calibs",
                  "Prior ages with 26 ST calibs",
                  "Analytical ST calibrations" )
col.legend  <- c( "black", "green", "red" )
#coords.plot <- locator()
legend("center", legend = info.legend, col = col.legend,
        lty = 1, bty = 'n', cex = 0.3 )

dev.off()


# Plot analytical ST VS prior with analytical ST
pdf( file = "02_MCMCtree_prior_ST/plots/00_Check_analitST_priorST.pdf","a4")
par( mfrow = c(8,9), mai=c(0,0,0,0))
for( i in 1:length( colnames( divtimes.prior.ST ) ) ){
  
  # 1. Get a tmp_node with ST distributions
  tmp_node <- ST.fitted.objects[[ i ]]
  
  # 2. Plot fitted ST for each node previously 
  #    computed using the values sampled during the MCMC
  
  # 2.1. Find limit axis
  max_x_st            <- round( max( density( rst( n = 1000, dp = tmp_node$dp.complete ) )$x ) + 0.5 )
  max_x_chain_priorST <- round( max( density( divtimes.prior.ST[,i] )$x ) + 0.5 )
  x_lim       <- max( max_x_st, max_x_chain_priorST )
  max_y_st            <- round( max( density( rst( n = 1000, dp = tmp_node$dp.complete ) )$y ) + 0.5 )
  max_y_chain_priorST <- round( max( density( divtimes.prior.ST[,i] )$y ) + 0.5 )
  y_lim       <- max( max_y_st, max_y_chain_priorST )
  # 2.2. Plot
  plot( density( divtimes.prior.ST[,i], adj = 1 ),
        xlim = c( 0, x_lim ), ylim = c( 0, y_lim ), main = colnames( divtimes.prior.ST )[i],
        xlab = '', ylab = '', col = "green" )
  curve( dst( x, xi = tmp_node$dp[1], omega = tmp_node$dp[2],
              alpha = tmp_node$dp[3], nu = tmp_node$dp[4] ),
         from = 0, to = x_lim,
         n = 1e4, add = TRUE, col = "red" )
  
}

plot( x = 1, y = 1, col = "white" )
info.legend <- c( "Prior ages using 26 ST calibs",
                  "Analytical ST calibrations" )
col.legend  <- c( "green", "red" )
#coords.plot <- locator()
legend("topright", legend = info.legend, col = col.legend,
        lty = 1, bty = 'n', cex = 0.8 )

dev.off()

# Repeat the same but extract individual plots too
for( i in 1:length( colnames( divtimes.prior.ST ) ) ){
  
  # 1. Get a tmp_node with ST distributions
  tmp_node <- ST.fitted.objects[[ i ]]
  
  # 2. Plot fitted ST for each node previously 
  #    computed using the values sampled during the MCMC
  pdf( file = paste( "02_MCMCtree_prior_ST/plots/individual_plots/Compare_priorST.vs.analytST_", colnames( divtimes )[i], ".pdf", sep = "" ),"a4")
  # 2.1. Find limit axis
  max_x_st            <- round( max( density( rst( n = 1000, dp = tmp_node$dp.complete ) )$x ) + 0.5 )
  max_x_chain_priorST <- round( max( density( divtimes.prior.ST[,i] )$x ) + 0.5 )
  x_lim       <- max( max_x_st, max_x_chain_priorST )
  max_y_st            <- round( max( density( rst( n = 1000, dp = tmp_node$dp.complete ) )$y ) + 0.5 )
  max_y_chain_priorST <- round( max( density( divtimes.prior.ST[,i] )$y ) + 0.5 )
  y_lim       <- max( max_y_st, max_y_chain_priorST )
  # 2.2. Plot
  plot( density( divtimes.prior.ST[,i], adj = 1 ),
        xlim = c( 0, x_lim ), ylim = c( 0, y_lim ), #main = colnames( divtimes.prior.ST )[i],
        xlab = '', ylab = '', col = "green",
        main = paste( colnames( divtimes )[i], " = ",
                      "ST(", paste0( round(tmp_node$dp, 2),
                                     collapse = "," ),
                      ")", sep = "" ) )
  curve( dst( x, xi = tmp_node$dp[1], omega = tmp_node$dp[2],
              alpha = tmp_node$dp[3], nu = tmp_node$dp[4] ),
         from = 0, to = x_lim,
         n = 1e4, add = TRUE, col = "red" )
  # 2.3. Legend 
  info.legend <- c( "Est. prior times using 26 ST calibs",
                    "Analytical ST calibs" )
  col.legend  <- c( "green", "red" )
  legend( "topright", legend = info.legend, col = col.legend,
          lty = 1, bty = 'n', cex = 0.8 )
  # 2.3. Close graph
  dev.off()
  
}
