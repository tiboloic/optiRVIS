# LT 13/12/2019

#
# TMB multinomial power law with random slope

library(TMB)
library(dplyr)
library(tidyr)

gnomad = read.delim(file='gnomADSFS.tsv')
SFS = gnomad %>% spread(key=cat, value=obs, fill=0)

# remove gene information: keep only counts
dat = as.matrix(SFS[1:200,-1, drop=FALSE])

compile("powerSFSpar.cpp")

dyn.load(dynlib("powerSFSpar"))

p=profmem({
#system.time({
  f = MakeADFun(
  data=list(obs=dat),
  parameters=list(b=1.8, lsigb=-3, bdevs=rep(0,nrow(dat))),
#  parameters=list(b=1.77, lsigb=-32.44, bdevs=bdevs),
  random=c("bdevs"),
  DLL="powerSFSpar", silent=FALSE, atomic = TRUE)
})
sum(p$bytes, na.rm=T)

#p=profmem({
system.time({
  fit = optim(f$par,f$fn,f$gr, method="BFGS", control=list(maxit=3000))
})

dyn.unload(dynlib('powerSFSpar'))

#save.image("allfit.RData")