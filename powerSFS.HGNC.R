# LT 3/06/2020

# refit model using variant list from gene wich HGNC id
# use refitted pade approximation

# LT 13/12/2019

#
# TMB multinomial power law with random slope

# 18/12/2019
# padé approximation of power porbabilities

library(TMB)
library(dplyr)
library(tidyr)

gnomad = read.delim(file='gnomADSFS.3.tsv')
SFS = gnomad %>% spread(key=cat, value=obs, fill=0)

# get padé approximation coefficentsf$
load('pade2.RData')

# remove gene information: keep only counts
dat = as.matrix(SFS[,-1, drop=FALSE])

compile("powerSFS3.cpp")

dyn.load(dynlib("powerSFS3"))

f = MakeADFun(
  data=list(obs=dat, pade=coefs),
  parameters=list(b=1.8, lsigb=-2, bdevs=rep(0,nrow(dat))),
  random=c("bdevs"),
  DLL="powerSFS3", silent=FALSE)


fit = optim(f$par,f$fn,f$gr, method="BFGS", control=list(maxit=3000))
fit.nl = nlminb(f$par,f$fn,f$gr, method="BFGS", control=list(maxit=3000))
mysd = sdreport(f)

# importance sampling DID nOT WORK
g = MakeADFun(
  data=list(obs=dat, pade=coefs),
  parameters=list(b=1.8, lsigb=-2, bdevs=f$env$parList()$bdevs),
  random=c("bdevs"),
  MCcontrol = list(doMC = TRUE, seed = 9011971, n = 100),
  DLL="powerSFS3", silent=FALSE)

fit2 = optim(g$par,g$fn,g$gr, method="BFGS", control=list(maxit=3000))
#thispars = f$env$parList()
SFS$beta = f$env$parList()$bdevs
save(f, SFS, mysd, file='fitall.HGNC.RData')
