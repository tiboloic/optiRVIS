# LT 10/06/2020
# Analyze revised powerSFS filtering on HGNC genes

library(dplyr)
library(tidyr)

#load('fit.Polyphen.RData')
#load('fit.AN.RData')
#load('fit.AN.RData')
load('fit.AN.80.RData')

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
plot(oe_mis_upper ~ beta, data=dat, xlab="beta", ylab='LOEUF')
cor.test(dat$oe_mis_upper, dat$beta)
summary(lm(oe_mis_upper ~ beta, data=dat))
# R2 = 0.26

# logistic regression on haploinsufficient
# load gnomAD genelist
genelist = read.delim('supplementary_dataset_13_gene_lists.tsv.gz')

# keep only what is needed
HIolf = subset(genelist, gene_list %in% c('Haploinsufficient', 'Olfactory Genes'))
HIolf$Y = HIolf$gene_list == 'Haploinsufficient'
HIolf = inner_join(dat[,c('gene', 'oe_lof_upper', 'oe_mis_upper', 'gevir_percentile', 'beta', 'cds_length', 'oe_lof', 'oe_mis')], HIolf)

# fit logistic regression
m = glm(Y ~ oe_lof_upper, data = HIolf, family = binomial)
m2 = glm(Y ~ oe_lof_upper + log(cds_length), data = HIolf, family = binomial)

m3 = glm(Y ~ oe_lof_upper + log(cds_length) + gevir_percentile, data = HIolf, family = binomial)
m4 = glm(Y ~ oe_lof_upper + log(cds_length) + gevir_percentile + beta, data = HIolf, family = binomial)
m4 = glm(Y ~ log(cds_length) + qnorm(gevir_percentile/100) + beta + oe_lof_upper, data = HIolf, family = binomial)

m7 = glm(Y ~ log(cds_length) + qnorm(gevir_percentile/100) + beta + oe_lof_upper + oe_mis_upper, data = HIolf, family = binomial)
m8 = glm(Y ~ log(cds_length) + qnorm(gevir_percentile/100) + beta + oe_lof_upper + oe_mis_upper + oe_lof + oe_mis, data = HIolf, family = binomial)

library(MASS)

stepAIC(m8)
# keeps cds_length, oe_lof_upper, beta !

m.loeufbeta = glm(Y ~ oe_lof_upper + log(cds_length) + beta, data = HIolf, family = binomial)
m.beta = glm(Y ~ beta + log(cds_length), data = HIolf, family = binomial)
m.loeuf = glm(Y ~ oe_lof_upper + log(cds_length), data = HIolf, family = binomial)


library(ROCR)

ROCRperf.loeufbeta =  performance(prediction(fitted(m.loeufbeta),HIolf$Y), 'tpr','fpr')
ROCRperf.loeuf =  performance(prediction(fitted(m.loeuf),HIolf$Y), 'tpr','fpr')
ROCRperf.beta =  performance(prediction(fitted(m.beta),HIolf$Y), 'tpr','fpr')

# colorblind colorpalette
cols  = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
plot(ROCRperf.loeufbeta, lwd=2, col=cols[2])
plot(ROCRperf.loeuf, lwd=2, col=cols[4], add=T)
plot(ROCRperf.beta, lwd=2, col=cols[6], add=T)
legend('bottomright', legend = c("LOEUF + powerSFS", "LOEUF","powerSFS"),
       fill=cols[c(2,4,6)], bty='n')

# identify key genes
plot(ROCRperf.loeufbeta@y.values[[1]][150:300]~ROCRperf.loeufbeta@x.values[[1]][150:300], type='l')
lines(ROCRperf.loeuf@y.values[[1]][150:300]~ROCRperf.loeuf@x.values[[1]][150:300], col=2)
abline(v=0.027)
which.min(abs(ROCRperf.loeufbeta@x.values[[1]]-0.027))
cutoff = ROCRperf.loeufbeta@alpha.values[[1]][194]
genes.loeufbeta = HIolf[which(fitted(m.loeufbeta)>cutoff),]
genes.loeuf = HIolf[which(fitted(m.loeuf)>ROCRperf.loeuf@alpha.values[[1]][which.min(abs(ROCRperf.beta@x.values[[1]]-0.027))]),]
keygenes = subset(genes.loeufbeta, !gene %in% genes.loeuf$gene & gene_list=='Haploinsufficient')

# quick lasso
library(glmnet)
x = as.matrix(HIolf[,c('oe_lof_upper', 'oe_mis_upper', 'gevir_percentile', 'beta', 'cds_length', 'oe_lof', 'oe_mis')])
x[,'cds_length'] = log(x[,'cds_length'])
lass1 = glmnet(x,y = HIolf$Y, alpha=1, family='binomial')

plot(lass1, xvar="lambda", label=TRUE, col= cols[1:7], lwd=2)
legend("bottomright", lwd = 2, col = cols[1:7], legend = colnames(x))

cv.lass1 <- cv.glmnet(x, y = HIolf$Y, alpha=1)
plot(cv.lass1)

#plot(ROCRperf5, colorize = F, text.adj = c(-0.2,1.7))
plot(ROCRperf5, colorize = F, text.adj = c(-0.2,1.7))
plot(ROCRperf6, colorize = F, text.adj = c(-0.2,1.7), add=T, col=2)


# a quick look at recessive genes

recess = subset(genelist, gene_list %in% c('Haploinsufficient', 'Olfactory Genes', 'Autosomal Recessive'))
recess$Y = recess$gene_list == 'Autosomal Recessive'
recess = inner_join(dat[,c('gene', 'oe_lof_upper', 'oe_mis_upper', 'gevir_percentile', 'beta', 'cds_length', 'oe_lof', 'oe_mis')], recess)

x = as.matrix(recess[,c('oe_lof_upper', 'oe_mis_upper', 'gevir_percentile', 'beta', 'cds_length', 'oe_lof', 'oe_mis')])
x[,'cds_length'] = log(x[,'cds_length'])
lass1 = glmnet(x,y = recess$Y, alpha=.5, family='binomial')
plot(lass1, xvar="lambda", label=TRUE, col= cols[1:7], lwd=2)
legend("bottomright", lwd = 2, col = cols[1:7], legend = colnames(x))
cv.lass1 <- cv.glmnet(x, y = recess$Y, alpha=1)
plot(cv.lass1)

# 
# look at all genes and plot figure
SFSdeciles = quantile(SFS$beta, probs=0:10/10)
SFS$beta_decile = cut(SFS$beta, SFSdeciles, labels=1:10)
#hist(SFS$beta, breaks=SFSdeciles)

# proportion of genes by categories
