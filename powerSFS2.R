# LT 13/12/2019

#
# TMB multinomial power law with random slope

# 18/12/2019
# padé approximation of power porbabilities

library(TMB)
library(dplyr)
library(tidyr)

gnomad = read.delim(file='gnomADSFS.tsv')
SFS = gnomad %>% spread(key=cat, value=obs, fill=0)

# get padé approximation coefficents
load('pade.RData')

# remove gene information: keep only counts
dat = as.matrix(SFS[,-1, drop=FALSE])

compile("powerSFS2.cpp")

dyn.load(dynlib("powerSFS2"))

f = MakeADFun(
  data=list(obs=dat, pade=coefs),
  parameters=list(b=1.8, lsigb=-3, bdevs=rep(0,nrow(dat))),
  random=c("bdevs"),
  DLL="powerSFS2", silent=FALSE)


fit = optim(f$par,f$fn,f$gr, method="BFGS", control=list(maxit=3000))

#thispars = f$env$parList()
save(f, file='fitall.RData')
