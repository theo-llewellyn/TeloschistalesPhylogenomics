
#This script was written by Sandra Alvarez-Carretero and adapted for the Teloschistales dataset by Theo Llewellyn

#--------------#
# LOAD PACKAGE #
#--------------#
library( MCMCtreeR )
library( rstudioapi )

source( "/science/projects/gaya_lab/Teloschistaceae_1/MCMCTree/LICHENS_MCMCTree/00_Esters2a_calibration/Functions_plots_MCMCtreeR.R" )

#-----------------------#
# SET WORKING DIRECTORY #
#-----------------------#
#run within LICHENS_MCMCTree/00_Esters2a_calibration/01_MCMCtree_posterior
setwd("/science/projects/gaya_lab/Teloschistaceae_1/MCMCTree/LICHENS_MCMCTree/00_Esters1a_calibration/01_MCMCtree_posterior")

#-------------#
# START TASKS #
#-------------#

# 1. Load trees and mcmc files for master tree 
## GBM - posterior
phy.post.GBM <- mcmc.post.GBM <- vector( mode = "list", 16 )
for( i in seq(1:length( phy.post.GBM ) ) ){
        
        cat( "Parsing files for run ", i, "...\n" )
        phy.post.GBM[[ i ]] <-  MCMCtreeR::readMCMCtree( paste("01_MCMCtree_posterior/","run",
                                                                i, "/mcmctree_GBM/FigTree.tre", sep = "" ),
                                                         from.file = TRUE )
        mcmc.post.GBM[[ i ]] <- read.table( paste("01_MCMCtree_posterior/","run",
                                                   i, "/mcmctree_GBM/mcmc.txt", sep = "" ),
                                            sep = "\t", header = T, stringsAsFactors = F )
}

## GBM - Prior
phy.prior.GBM <- mcmc.prior.GBM <- vector( mode = "list", 5 )
for( i in seq(1:length( phy.prior.GBM ) ) ){
        
        cat( "Parsing files for run ", i, "...\n" )
        phy.prior.GBM[[ i ]] <-  MCMCtreeR::readMCMCtree( paste("00_MCMCtree_prior/mcmc",
                                                           i, "/FigTree.tre", sep = "" ),
                                                    from.file = TRUE )
        mcmc.prior.GBM[[ i ]] <- read.table( paste("00_MCMCtree_prior/mcmc",
                                                  i, "/mcmc.txt", sep = "" ),
                                           sep = "\t", header = T, stringsAsFactors = F )
}

# 2. Create objects in APE format 
## GBM - posterior
phy.post.GBM.obj <- phy.post.GBM.edge <- mcmc.chain.GBM <- vector( mode = "list", 16 )
for( i in seq(1, length( phy.post.GBM.obj ) ) ){
        cat( "Parsing files from run ", i, "...\n" )
        phy.post.GBM.obj[[ i ]]  <- phy.post.GBM[[ i ]]$apePhy
        phy.post.GBM.edge[[ i ]] <- phy.post.GBM.obj[[ i ]]$edge
        mcmc.chain.GBM[[ i ]]    <- mcmc.post.GBM[[ i ]]
        
}

## GBM - prior 
phy.prior.GBM.obj <- phy.prior.GBM.edge <- mcmc.chain.prior.GBM <- vector( mode = "list", 5 )
for( i in seq(1, length( phy.prior.GBM.obj ) ) ){
        cat( "Parsing files from run ", i, "...\n" )
        phy.prior.GBM.obj[[ i ]]     <- phy.prior.GBM[[ i ]]$apePhy
        phy.prior.GBM.edge[[ i ]]    <- phy.prior.GBM.obj[[ i ]]$edge
        mcmc.chain.prior.GBM[[ i ]]  <- mcmc.prior.GBM[[ i ]]
        
}


# 3. Create PHY object with GBM and ILN objects 
phy.all <- vector( mode = "list", length = 21 )
names( phy.all ) <- c( paste( "GBM.post.", 1:16, sep = "" ),
                       paste( "GBM.prior.", 1:5, sep = "" ) )
j <- 0
for( i in seq(1:length( phy.all ) ) ){
        
        if( i < 17 ){
                print( i )
                cat( "Parsing GBM.post files from run ", i, "...\n" )
                phy.all[[ i ]] <- phy.post.GBM.obj[[ i ]]
        }
        else{
                j <- j + 1
                print( i )
                cat( "Parsing GBM.prior files from run ", j, "...\n" )
                phy.all[[ i ]] <- phy.prior.GBM.obj[[ j ]] 
                
        }
}


# 4. Extract ages with node age posteriors from column 2
## GBM - posterior
mcmc.node.ages.post.GBM <- all.nodes.post.GBM <- vector( mode = "list", 16 )
for( i in 1:length( mcmc.node.ages.post.GBM ) ){
        mcmc.node.ages.post.GBM[[ i ]] <- mcmc.chain.GBM[[ i ]][, 2:Ntip( phy.all[[ i ]] )] 
        all.nodes.post.GBM[[ i ]]      <- as.numeric( gsub( "t_n", "",
                                                            colnames( mcmc.node.ages.post.GBM[[ i ]] ) ) )
}

## GBM - prior
mcmc.node.ages.prior.GBM <- all.nodes.prior.GBM <- vector( mode = "list", 5 )
j <- 16
for( i in 1:5 ){
        j <- j + 1
        mcmc.node.ages.prior.GBM[[ i ]] <- mcmc.chain.prior.GBM[[ i ]][, 2:Ntip( phy.all[[ j ]] )] 
        all.nodes.prior.GBM[[ i ]]      <- as.numeric( gsub( "t_n", "",
                                                             colnames( mcmc.node.ages.prior.GBM[[ i ]] ) ) )
}

# 5. Create a vector of names for each list element as  
#    internal nodes from APE tree, using phy$edge object 
## GBM - posterior
node.ages.names.post.GBM <- vector( mode = "list", 16 )
for( i in 1:length( node.ages.names.post.GBM ) ){
        node.ages.names.post.GBM[[ i ]] <- c( Ntip( phy.all[[ i ]] ) + 1,
                                              phy.post.GBM.edge[[ i ]][ which( phy.post.GBM.edge[[ i ]][,2] > Ntip( phy.all[[ i ]] ) ), 2 ] )
        
        
}

## GBM - prior
node.ages.names.prior.GBM <- vector( mode = "list", 5 )
j <- 16
for( i in 1:length( node.ages.names.prior.GBM ) ){
        j <- j + 1
        node.ages.names.prior.GBM[[ i ]] <- c( Ntip( phy.all[[ j ]] ) + 1,
                                              phy.prior.GBM.edge[[ i ]][ which( phy.prior.GBM.edge[[ i ]][,2] > Ntip( phy.all[[ j ]] ) ), 2 ] )
        
        
}

# 6. Find where each posterior node age appears in APE 
#    edge object.
## GBM - posterior
match.nodes.post.GBM <- vector( mode = "list", 16 )
for( i in 1:length( match.nodes.post.GBM ) ){
        match.nodes.post.GBM[[ i ]] <- match( all.nodes.post.GBM[[ i ]],
                                              as.numeric( node.ages.names.post.GBM[[ i ]] ) )
        
}

## GBM - prior
match.nodes.prior.GBM <- vector( mode = "list", 5 )
for( i in 1:length( match.nodes.prior.GBM ) ){
        match.nodes.prior.GBM[[ i ]] <- match( all.nodes.prior.GBM[[ i ]],
                                              as.numeric( node.ages.names.prior.GBM[[ i ]] ) )
        
}

# 7. Create a list extracting the info from the mcmc
#    chain in APE node order 
## GBM - posterior
node.ages.post.GBM <- vector( mode = "list", 16 )
for( i in 1:length( node.ages.post.GBM ) ){
        node.ages.post.GBM[[ i ]] <- lapply( match.nodes.post.GBM[[ i ]],
                                             function( uu ) mcmc.node.ages.post.GBM[[ i ]][, uu ] )
}

## GBM - prior 
node.ages.prior.GBM <- vector( mode = "list", 5 )
for( i in 1:length( node.ages.prior.GBM ) ){
        node.ages.prior.GBM[[ i ]] <- lapply( match.nodes.prior.GBM[[ i ]],
                                             function( uu ) mcmc.node.ages.prior.GBM[[ i ]][, uu ] )
}

# 8. Name each element in list 
## GBM - posterior
for( i in 1:length( node.ages.names.post.GBM ) ){
        names( node.ages.post.GBM[[ i ]] ) <- node.ages.names.post.GBM[[ i ]]
        
}

## GBM - prior 
for( i in 1:length( node.ages.names.prior.GBM ) ){
        names( node.ages.prior.GBM[[ i ]] ) <- node.ages.names.prior.GBM[[ i ]]
        
}


# 9. Create node.ages.all object to get all GBM and ILN objects
node.ages.all <- vector( mode = "list", length = 21 )
j <- 0
for( i in 1:length( node.ages.all ) ){
        
        if( i < 17 ){
                node.ages.all[[ i ]] <- node.ages.post.GBM[[ i ]]
        }
        else{
                j <- j + 1
                node.ages.all[[ i ]] <- node.ages.prior.GBM[[ j ]]
        }
}

# 10. Get means for edges
## GBM - posterior
# Initialize matrix to get the edges of the 16 runs 
phy.matrix.post.GBM.edges <- matrix( 0, nrow = 16, ncol = length( phy.all[[1]]$edge.length ) )
for( i in 1:16 ){
        phy.matrix.post.GBM.edges[i,] <- phy.all[[ i ]]$edge.length
}
# Get mean of each edge. Also generate one mean without chains 4, 7, 8 -- according to visualization
# The latter object has been generated and added here after analysing the plots that appear later
# in the script
phy.mean.edge.post.GBM      <- apply( X = phy.matrix.post.GBM.edges, 2, mean )
# Create copy of ape object and then replace the edge.length values accordingly
phy.mean.post.GBM                   <- phy.mean.post.GBM.filt <- phy.all[[1]]
phy.mean.post.GBM$edge.length       <- phy.mean.edge.post.GBM

## GBM - prior
phy.matrix.prior.GBM.edges <- matrix( 0, nrow = 5, ncol = length( phy.all[[17]]$edge.length ) )
j <- 16
for( i in 1:5 ){
        j <- j + 1
        phy.matrix.prior.GBM.edges[i,] <- phy.all[[ j ]]$edge.length
}
# Get mean of each edge.
phy.mean.edge.prior.GBM      <- apply( X = phy.matrix.prior.GBM.edges, 2, mean )
# Create copy of ape object and then replace the edge.length values accordingly
phy.mean.prior.GBM                   <- phy.all[[17]]
phy.mean.prior.GBM$edge.length       <- phy.mean.edge.prior.GBM

# 11. Get means for ages
#------------------------------------#
# FUNCTION TO GET THE MEAN DISTANCES # 
#------------------------------------#
# Arguments:
#    len              numeric, number of internal nodes for which samples have been
#                     collected. Default: 71.
#    chains           numeric, indexes of the chains through which iterate to get the 
#                     collected samples for each age. Default: 1:10.
#    node.ages.list   list, list previously created with "chain" amount of entries 
#                     where the node ages each chain have been saved. Default: 10 chains,
#                     although it can be less, e.g., "prior.GBM".
#                     These objects "node.ages.list" could be "node.ages.post.GBM",
#                     or "node.ages.post.ILN", for instance.
#    tmp              numeric, temporary empty vector of length amount of samples 
#                     collected during the MCMCMs. Default: 20,001.
#    out              list, empty list previously created with "len" amount.
#                     E.g., "node.mean.dist.ages.post.GBM"
get_mean_dist_ages <- function( len = 71, chains = 1:10, node.ages.list, tmp ){
        
        # Create a list with 71 node entries, which we will output. We will then iterate
        # over each of the entries for each node sampled during each run. We will then
        # save the corresponding sampled vaules for this node 
        # to a the corresponding entry in the new created list.
        out          <- vector( mode = "list", length( node.ages.list[[ 1 ]] ) )
        names( out ) <- names( node.ages.list[[ 1 ]] ) 
        
        # Iterate over each node. We assume there are len = 71.
        # This will be modified if needed through the argument "len"
        for( j in 1:len ){
                
                # Reset local vars to access tmp vector
                a  <- 1
                b  <- 20001
                cte <- 20000
                
                # Iterate through every run and ten collect values for node_j
                for( i in chains ){
                        tmp[a:b] <- unlist( node.ages.list[[ i ]][j] ) # These are 20,001 values for node_j for run_i  
                        a <- b+1
                        b <- a+cte
                }
                
                # Save all the 20,001*10 values that have been sampled across the 10 runs 
                # and saved in vector "tmp"
                # Note that "out" is the list object with 71 entries, e.g., "node.mean.dist.ages.post.GBM"
                out[[ j ]] <- tmp 
                
        }
        
        # Return out
        return( out )
}

# a) Create tmp vector to allocate memory with 10*20,001 values
tmp <- rep( 0, 16*length( unlist( node.ages.post.GBM[[ 1 ]][1] ) ) )
# b) Loop from node 1 to 71 and then go through the corresponding values sampled for 
#    this node for each run. Use tmp vector to save values
node.mean.dist.ages.post.GBM <- get_mean_dist_ages( len = length( node.ages.post.GBM[[ 1 ]] ),
                                                    chains = 1:16,
                                                    node.ages.list = node.ages.post.GBM,
                                                    tmp = tmp )


# d) Now the same for GBM.prior 
tmp <- rep( 0, 5*length( unlist( node.ages.prior.GBM[[ 1 ]][1] ) ) )
node.mean.dist.ages.prior.GBM <- get_mean_dist_ages( len = length( node.ages.prior.GBM[[ 1 ]] ),
                                                     chains = 1:5,
                                                     node.ages.list = node.ages.prior.GBM,
                                                     tmp = tmp )
save( node.mean.dist.ages.prior.GBM, file = "../R_objects/node.mean.dist.ages.prior.GBM.RData" )

# TRY PLOT [[ IT WORKS! It returns the last plot ]]
# mcmc.tree.plot.RETPLOT( phy = phy.all$GBM, node.ages = node.ages.GBM, 
#                  analysis.type = "user", cex.tips = 0.2,
#                  time.correction = 100, 
#                  scale.res = c( "Eon", "Period" ), plot.type = "distributions",
#                  cex.age = 0.4, cex.labels = 0.5, relative.height = 0.08, 
#                  col.tree = "grey40", no.margin = TRUE,
#                  density.col = "blue", add.time.scale = TRUE)


#-------#
# PLOTS #
#-------#

#================================#
# 0. PLOT: Evaluate prior chains #
#================================#
## THE NEXT LINES CAN BE RUN AT ONCE
## THE OUTPUT FILE IS SAVED AS 
##    "RPlots/00_Explore_prior_runs.pdf"
# Set transparent colour
pdf("RPlots/00_Explore_prior_runs.pdf", "a4")
transp.col <- adjustcolor( col = "blue", alpha.f = 0.3 )
last.plot  <- mcmc.tree.plot.RETPLOT( phy = phy.all$GBM.prior.1, xlim.scale = c(-20,300),
                                      node.ages = node.ages.prior.GBM[[ 1 ]], 
                                      show.tip.label = TRUE,
                                      analysis.type = "user", cex.tips = 0.5,
                                      time.correction = 100, 
                                      scale.res = c( "Eon", "Period" ), plot.type = "distributions",
                                      cex.age = 0.6, cex.labels = 0.8, relative.height = 0.08, 
                                      col.tree = "grey40", no.margin = TRUE,
                                      density.col = transp.col, add.time.scale = TRUE, grey.bars = FALSE #,density.border.col = "red"
)

# Add now the other distributions with this helper function
col.mod = "grey"
transparency = 0.3 
col.plot <- c( "purple", "black", "darkgreen", "red", "darkblue", "orange", "pink",
               "lightblue", "brown", "lightgreen", "brown",
               "antiquewhite4", "grey", "yellow", "chocolate", "cyan" )
j <- 17
for( i in 2:5 ){
        j <- j + 1
        add.extra.dists( phy = phy.all[[ j ]], num.models = 1, last.plot = last.plot,
                         node.ages =  node.ages.prior.GBM[[ i ]], plot.type = "distributions",
                         time.correction = 100,
                         density.col = "white", density.border.col = col.plot[ i ],
                         distribution.height = 0.8, transparency = transparency
        )
}

# Add legend
legend( "left", legend = c( "Prior.r1", "Prior.r2", "Prior.r3",
                                "Prior.r4", "Prior.r5"),
        col = col.plot[1:5],
        lty = 1,
        bty = "n"
)
dev.off()
# All seem ok !

#====================================#
# 1. PLOT: Evaluate posterior chains #
#====================================#
## THE NEXT LINES CAN BE RUN AT ONCE
## THE OUTPUT FILE IS SAVED AS 
##    "plots/01_Explore_GBM_posterior_runs.pdf"
# Set transparent colour
pdf("RPlots/01_Explore_GBM_posterior_runs.pdf", "a4")
transp.col <- adjustcolor( col = "blue", alpha.f = 0.3 )
last.plot  <- mcmc.tree.plot.RETPLOT( phy = phy.all$GBM.post.1, xlim.scale = c(-20,300),
                                      node.ages = node.ages.post.GBM[[ 1 ]], 
                                      show.tip.label = TRUE,
                                      analysis.type = "user", cex.tips = 0.5,
                                      time.correction = 100, 
                                      scale.res = c( "Eon", "Period" ), plot.type = "distributions",
                                      cex.age = 0.6, cex.labels = 0.8, relative.height = 0.08, 
                                      col.tree = "grey40", no.margin = TRUE,
                                      density.col = transp.col, add.time.scale = TRUE, grey.bars = FALSE #,density.border.col = "red"
)

# Add now the rest of GBM distributions with this helper function
col.mod = "grey"
transparency = 0.3 
col.plot <- c( "purple", "black", "darkgreen", "red", "darkblue", "orange", "pink",
               "lightblue", "brown", "lightgreen", "brown",
               "antiquewhite4", "grey", "yellow", "chocolate", "cyan" )
for( i in 2:16 ){
        add.extra.dists( phy = phy.all[[ i ]], num.models = 1, last.plot = last.plot,
                         node.ages =  node.ages.post.GBM[[ i ]], plot.type = "distributions",
                         time.correction = 100,
                         density.col = "white", density.border.col = col.plot[ i ],
                         distribution.height = 0.8, transparency = transparency
        )
}

# Add legend
legend("left",legend = c( "GBM.post.r1", "GBM.post.r2", "GBM.post.r3", "GBM.post.r4", "GBM.post.r5",
                                "GBM.post.r6", "GBM.post.r7", "GBM.post.r8", "GBM.post.r9", "GBM.post.r10",
                                "GBM.post.r11", "GBM.post.r12", "GBM.post.r13", "GBM.post.r14",
                                "GBM.post.r15", "GBM.post.r16" ),
        col = col.plot[1:16],
        lty = 1,
        bty = "n"
)
dev.off()


#=======================================#
# 2. PLOT: Plot prior against posterior #
#=======================================#
## THE NEXT LINES CAN BE RUN AT ONCE
## THE OUTPUT FILE IS SAVED AS 
##    "plots/01_OLD-mean_filtered_postVSmean_prior.pdf"
# NOTE: We will use the mean ages of the GBM chains without 4, 7, 8 for that
# Set transparent colour
pdf("RPlots/01_mean_postVSmean_prior.GBM.pdf", "a4")
transp.col <- adjustcolor( col = "blue", alpha.f = 0.3 )
last.plot  <- mcmc.tree.plot.RETPLOT( phy = phy.mean.post.GBM, xlim.scale = c(-20,300),
                                      node.ages = node.mean.dist.ages.post.GBM, 
                                      show.tip.label = TRUE,
                                      analysis.type = "user", cex.tips = 0.5,
                                      time.correction = 100, 
                                      scale.res = c( "Eon", "Period" ), plot.type = "distributions",
                                      cex.age = 0.6, cex.labels = 0.8, relative.height = 0.08, 
                                      col.tree = "grey40", no.margin = TRUE,
                                      density.col = transp.col, add.time.scale = TRUE, grey.bars = FALSE
)

# Add now the prior distribution
col.mod = "grey"
transparency = 0.3 
col.plot <- c( "darkgreen" )
add.extra.dists( phy = phy.mean.prior.GBM, num.models = 1, last.plot = last.plot,
                 node.ages =  node.mean.dist.ages.prior.GBM, plot.type = "distributions",
                 time.correction = 100,
                 density.col = "white", density.border.col = col.plot,
                 distribution.height = 0.8, transparency = transparency
)

# Add legend
legend( "topleft", legend = c( "mean.post", "mean.prior" ),
        col = c( "blue", "darkgreen" ),
        lty = 1,
        bty = "n"
)
dev.off()
#--------------------#
# GENERATE R OBJECTS #
#--------------------#
# 1. Generate list object that can be used to plot with 
#    both posterior and prior data 

save( phy.all, file = "R_objects/phy.all.Rdata" )
save( node.ages.post.GBM, file = "R_objects/node.ages.post.GBM.Rdata" )
save( node.ages.prior.GBM, file = "R_objects/node.ages.prior.GBM.Rdata" )

# 2. Generate list object with the mean ages of the prior
save( phy.mean.prior.GBM, file = "R_objects/phy.mean.prior.GBM.RData" )

# 3. Generate list object with the mean ages of the posterior after 
#    it has been filtered (it keeps only chains without conflict)
phy.mean.post.GBM.filt.OLD <- phy.mean.post.GBM.filt
save( phy.mean.post.GBM.filt.OLD, file = "R_objects/phy.mean.post.GBM.filt.OLD.RData" )
