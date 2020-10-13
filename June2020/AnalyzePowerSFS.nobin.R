# LT 10/06/2020
# Analyze revised powerSFS filtering on HGNC genes

library(dplyr)
library(tidyr)

load('fit.nobin.RData')


# get loeuf
cons = read.delim('../gnomad/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz')
cons$loeuf = qchisq(0.95, 2 *(cons$obs_lof+1))/2/cons$exp_lof
func = ecdf(cons$loeuf)
cons$stdloeuf = qnorm(func(cons$loeuf))
dat = inner_join(SFS, cons)

# get geVIR

gevir = read.csv('gevir.csv', skip=1)
dat = inner_join(dat, gevir, by=c(gene="gnomad_gene_name"))

dat$pSFS_percentile[order(dat$beta, decreasing = TRUE)] = (1:nrow(dat)) * 100 / nrow(dat)

# quick gnomD PLOT
plot(oe_mis_upper ~ beta, data=dat, xlab="pSFS", ylab='LOEUF')
cor.test(dat$oe_mis_upper, dat$beta)

plot(stdloeuf ~ beta, data=dat, xlab="pSFS", ylab='standardized LOEUF')
cor.test(dat$oe_mis_upper, dat$beta)

# identify some ouliers on the graph
# very intolerant in gnomad, tolerant in pSFS
dat[c(6981,12881,16226,17396,17522),]
# use threshold: beta <0, LOEUF < 0.35
set1 = subset(dat, beta < 0 & oe_lof_upper < 0.35)
# 608 genes where we do bad
# should plot the SFS

median(set1$cds_length)
#[1] 2644.5
median(dat$cds_length)
#[1] 1296


# tolerant in gnomad, very intolerant pSFS


# there  is some length bias
plot(beta ~ log(cds_length), data=dat, xlab='cds length (log)', ylab='pSFS')
summary(lm(beta ~ log(cds_length), data=dat))

m = lm(beta ~ log(cds_length), data=dat)
# quick and dirty length bias correction
dat$betacor = -residuals(m)

# pSFS percentile versus virlof percentile
plot(pSFS_percentile ~ virlof_percentile, data=dat)
lines(lowess(dat$virlof_percentile, dat$pSFS_percentile), col=2)

plot(pSFS_percentile ~ loeuf_percentile, data=dat)
lines(lowess(dat$virlof_percentile, dat$pSFS_percentile), col=2)

##############################
# plots of model fit

exactps = function(beta) {
  res = (1:1023)^-beta
  res = res/sum(res)
  res
}

SFSplot = function(i) {
  beta=f$env$parList()$b
  sig = exp(f$env$parList()$lsigb)
  plot(as.numeric(dat[i,2:1024]), log='xy', main=dat$gene[i], xlab='octave', ylab='# variants')
  lines(exactps(beta + sig * dat$beta[i]) * sum(as.numeric(dat[i,2:19])), col=2, lwd=2)
  lines(exactps(beta) * sum(dat[i,2:19]), col=1, lwd=1)
}

# plot most and least intolerant
sapply(order(dat$beta, decreasing = T)[1:20], SFSplot)

sapply(order(dat$beta, decreasing = F)[1:20], SFSplot)

# plot LOEUF intolerant, beta tolerant
sapply(match(set1$gene[order(set1$beta, decreasing = F)[1:20]], dat$gene), SFSplot)

# some LOEUF tolerant but pSFS intolerant
set2 = subset(dat, beta > 2 & oe_lof_upper > 1)
sapply(match(set2$gene[order(set2$beta, decreasing = F)[1:12]], dat$gene), SFSplot)


# plot of all mutations on arithmetic scale
beta=f$env$parList()$b
sig = exp(f$env$parList()$lsigb)
allSFS = as.matrix(dat[,2:1024])
bxs = barplot(colSums(allSFS), log='y', main='', xlab='octave', ylab='# variants')
lines(bxs, exactps(beta) * sum(allSFS), col=2, lwd=2)

bxs = barplot(colSums(allSFS), log='y', main='', xlab='octave', ylab='# variants (log)')
lines(bxs, approxps(beta, coefs) * sum(allSFS), col=2, lwd=2)

#list of most intolerant genes
cat(dat$gene[order(dat$beta, decreasing = T)[1:20]])


# load DDG2P
ddg2p = read.csv('DDG2P_11_6_2020.csv')
dat$gene[order(dat$beta, decreasing = T)[1:20]] %in% ddg2p$gene.symbol

# proportion of genes in DDG2P:
mean(dat$gene[order(dat$beta, decreasing = T)[1:20]] %in% ddg2p$gene.symbol)
# 50%

# list of most tolerant genes
cat(dat$gene[order(dat$beta, decreasing = F)[1:20]])

# plot key CHD genes
sapply(match(c('NR2F2','SMC3','BRAF','EFTUD2'), dat$gene), SFSplot)
