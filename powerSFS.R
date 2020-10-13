# LT 13/12/2019

#
# TMB multinomial power law with random slope

library(TMB)
library(dplyr)
library(tidyr)

gnomad = read.delim(file='gnomADSFS.tsv')
SFS = gnomad %>% spread(key=cat, value=obs, fill=0)

# remove gene information: keep only counts
dat = as.matrix(SFS[1:10,-1, drop=FALSE])

#compile("powerSFS.cpp")
#compile("powerSFS.cpp","-O1 -g",DLLFLAGS="")
compile("powerSFSpar.cpp")

#dyn.load(dynlib("powerSFS"))
dyn.load(dynlib("powerSFSpar"))

f = MakeADFun(
  data=list(obs=dat),
  parameters=list(b=1.8, lsigb=-3, bdevs=rep(0,nrow(dat))),
  random=c("bdevs"),
  DLL="powerSFS", silent=TRUE)

g = MakeADFun(
  data=list(obs=dat),
  parameters=list(b=1.8, lsigb=-3, bdevs=rep(0,nrow(dat))),
  random=c("bdevs"),
  DLL="powerSFSpar", silent=TRUE)

#fit = optim(f$par,f$fn,f$gr, method="BFGS", control=list(maxit=3000))
fit = optim(g$par,g$fn,g$gr, method="BFGS", control=list(maxit=3000))
