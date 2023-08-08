setwd("~/The Royal Botanic Gardens, Kew/Teloschistales - Documents/0_THEO_LLEWELLYN_PHD/03_RESULTS/PHYLOGENOMICS/TIME_CALIBRATION/MCMCTREE/MCMC3r/")
#set up blank dataframe to store results
BF_df <- data.frame(gene_name=character(),models=character(),logML=numeric(),SE=numeric(), bf=numeric(),logbf=numeric(),pr=numeric(),prior=numeric(),pr.ci.2.5.=numeric(),pr.ci.97.5.=numeric())

#set up counters to count number of times each model was the best model
CLK_count <- 0
ILN_count <- 0
GBM_count <- 0

#for each gene tested
for(gene in list.files(pattern="^OG*")){
  #change to that directory
  setwd(gene)
  #pull the log-marginal likelihoods for each model
  setwd("clk"); clk <- mcmc3r::stepping.stones()
  setwd("../iln"); iln <- mcmc3r::stepping.stones()
  setwd("../gbm"); gbm <- mcmc3r::stepping.stones()
  #names of model
  models <- c("CLK","ILN","GBM")
  #name of gene
  gene_name <- rep(gene,3)
  #gets log marginal likelihoods and standard errors for each
  logML <- c(clk$logml,iln$logml,gbm$logml)
  SE <- c(clk$se,iln$se,gbm$se)
  #adds logML and SE to Bayes factors, log BF, Posterior probabilities, priors and confidence intervals for Posterior probs.
  tmp <- data.frame(gene_name,models,logML,SE, mcmc3r::bayes.factors(clk, iln, gbm))
  
  #order probabilities low to high and get the name of the gene with the highest
  best_model <- tmp[order(tmp$pr)]$models[3]
  #if a particular model is the best it will add 1 to the count for that model
  if(best_model == "CLK"){
    CLK_count <- CLK_count + 1
  }
  if(best_model == "ILN"){
    ILN_count <- ILN_count + 1
  }
  if(best_model == "GBM"){
    GBM_count <- GBM_count + 1
  }
  #adds the results from that gene to dataframe
  BF_df <- rbind(BF_df,tmp)
  setwd("../..")
}

paste("CLK is best",CLK_count,"times")
paste("ILN is best",ILN_count,"times")
paste("GBM is best",GBM_count,"times")
