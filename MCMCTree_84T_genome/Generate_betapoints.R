##R-script
b = mcmc3r::make.beta(n=8, a=5, method="step-stones")
mcmc3r::make.bfctlf(b, ctlf="mcmctree.ctl", betaf="beta.txt")

dir.create("../../iln")
setwd("../../iln")
#then copy the clk .ctl file to iln and change clock model to 2
mcmc3r::make.bfctlf(b, ctlf="mcmctree.ctl", betaf="beta.txt")
dir.create("../gbm")
setwd("../gbm")
#then copy the clk .ctl file to iln and change clock model to 3
mcmc3r::make.bfctlf(b, ctlf="mcmctree.ctl", betaf="beta.txt")

setwd("../clk"); clk <- mcmc3r::stepping.stones()
setwd("../iln"); iln <- mcmc3r::stepping.stones()
setwd("../gbm"); gbm <- mcmc3r::stepping.stones()

models <- c("CLK","ILN","GBM")
logML <- c(clk$logml,iln$logml,gbm$logml)
SE <- c(clk$se,iln$se,gbm$se)
data.frame(models,logML,SE)
