# LT 11/06/2020
# TMB multinomial power law with random slope without binning

library(TMB)
library(dplyr)
library(tidyr)

gnomad = read.delim(file='gnomADSFS.nobin.tsv')
SFS = gnomad %>% spread(key=cat, value=obs, fill=0)


# remove gene information: keep only counts
dat = as.matrix(SFS[,-1, drop=FALSE])

# just a plot
plot(colSums(dat), log='y')
m = lm(log(colSums(dat)) ~ log(1:1023))
lines(exp(fitted(m)), col=2, lwd=2)
# MIB1
plot(dat[which(SFS$gene == 'MIB1'),], log='xy')


# more investigation on model
sumdat = colSums(dat)
N = sum(sumdat)
mean(colSums(dat)[800:1000])

compile("powerSFS_nobin.cpp")

dyn.load(dynlib("powerSFS_nobin"))

f = MakeADFun(
  data=list(obs=dat),
  parameters=list(b=1.8, lsigb=-3, bdevs=rep(0,nrow(dat))),
  random=c("bdevs"),
  DLL="powerSFS_nobin", silent=FALSE)


fit = nlminb(f$par,f$fn,f$gr, method="BFGS", control=list(maxit=3000))

mysd = sdreport(f)

SFS$beta = f$env$parList()$bdevs

save(f, SFS, mysd, file='fit.nobin.RData')
