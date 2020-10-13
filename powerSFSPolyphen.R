# LT 08/09/2020
# use polyphen score to weight missense

# LT 19/08/2020
# use variable allele number

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

gnomad = read.delim(file='gnomADSFS.Polyphen.1.tsv')

# remove an and obs (the counts)
SFS = gnomad[,c(-3,-5)] %>% spread(key=cat, value=poly, fill=0)

# get AN: group by gene, get sum of an and total number of variants
AN = gnomad %>% group_by(gene) %>% summarize(an=sum(an)/sum(obs))

any(AN$an<16)
which(AN$an<16)
SFS = SFS[AN$an>16,]
AN=AN[AN$an>16,]

# remove gene information: keep only counts
#dat = as.matrix(SFS[,-1, drop=FALSE])
dat = as.matrix(SFS[,-1, drop=FALSE])

#compile("powerSFSan.cpp", "-O1 -g",DLLFLAGS="")
compile("powerSFSPolyphen.cpp")

dyn.load(dynlib("powerSFSPolyphen"))

f = MakeADFun(
  data=list(obs=dat, an=as.vector(AN$an, "integer")),
  parameters=list(b=1.8, lsigb=-2, bdevs=rep(0,nrow(dat))),
  random=c("bdevs"),
  DLL="powerSFSPolyphen", silent=FALSE)


#fit = optim(f$par,f$fn,f$gr, method="BFGS", control=list(maxit=3000))
fit = nlminb(f$par,f$fn,f$gr, method="BFGS", control=list(maxit=5000))
mysd = sdreport(f)

SFS$beta = f$env$parList()$bdevs
SFS$an = AN$an
save(f, SFS, mysd, file='fit.Polyphen.RData')
