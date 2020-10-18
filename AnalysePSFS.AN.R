# LT 25/08/2020

# Analyse results from varying AN

library(dplyr)

# load old SFS
load('fitall.HGNC.RData')
SFSold = SFS

load('fit.AN.RData')
SFS.an = SFS

load('fit.AN.80.RData')
SFS.an80 = SFS

# list of genes lost going from SFS.an to SFS.an80
lostgenes = !(as.character(SFS.an$gene) %in% as.character(SFS.an80$gene)) 
#SFS.missing = left_join(SFS.an80, SFS.an, by=c('gene'='gene'))
#lg = !(levels(SFS.an$gene) %in% levels(SFS.an80$gene))
summary(SFS.an[lostgenes,]$beta)
# 75% are tolerant
hist(SFS.an$beta, prob=T, col=2)
hist(SFS.an[lostgenes,]$beta, add=T, col=0,prob=T)

# look at common genes
SFS.cmon = left_join(SFS.an80, SFS.an, by=c('gene'='gene'))
plot(beta.x ~ beta.y, data=SFS.cmon, ylab='AN filtered')
abline(a=0,b=1,col=2)


# look at correlation with other metrics
cons = read.delim('data/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz')
cons$loeuf = qchisq(0.95, 2 *(cons$obs_lof+1))/2/cons$exp_lof
func = ecdf(cons$loeuf)
cons$stdloeuf = qnorm(func(cons$loeuf))
dat = inner_join(SFS.cmon, cons)

# get geVIR
gevir = read.csv('gevir.csv', skip=1)
dat = inner_join(dat, gevir, by=c(gene="gnomad_gene_name"))

dat$pSFS_percentile[order(dat$beta.x, decreasing = TRUE)] = (1:nrow(dat)) * 100 / nrow(dat)

plot(oe_mis_upper ~ beta.x, data=dat, xlab="pSFS AN filtered", ylab='LOEUF')
cor.test(dat$oe_mis_upper, dat$beta.x)
plot(oe_mis_upper ~ beta.y, data=dat, xlab="pSFS AN variable", ylab='LOEUF')
cor.test(dat$oe_mis_upper, dat$beta.y)
plot(stdloeuf ~ beta.x, data=dat, xlab="pSFS", ylab='standardized LOEUF')
cor.test(dat$oe_mis_upper, dat$beta)

# is there  is some length bias
plot(beta.x ~ log(cds_length), data=dat, xlab='cds length (log)', ylab='pSFS')
summary(lm(beta.x ~ log(cds_length), data=dat))
