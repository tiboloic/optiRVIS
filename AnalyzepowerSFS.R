library(dplyr)
library(tidyr)

load('fit100.Rdata')
load('fit5000.Rdata')
load('fitall.RData')

# look at model
bs = f$env$parList()$b + exp(f$env$parList()$lsigb) * f$env$parList()$bdevs




# do the metric stuff
gnomad = read.delim(file='gnomADSFS.tsv')
SFS = gnomad %>% spread(key=cat, value=obs, fill=0)
#SFS = SFS[1:5000,]
SFS$beta = f$env$parList()$bdevs

# get loeuf
cons = read.delim('../gnomad/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz')

# recalculate loeuf my way
cons$loeuf = qchisq(0.95, 2 *(cons$obs_lof+1))/2/cons$exp_lof
f = ecdf(cons$loeuf)
cons$stdltloeuf = qnorm(f(cons$loeuf))
f = ecdf(cons$oe_lof_upper)
cons$stdloeuf = qnorm(f(cons$oe_lof_upper))

# some divergence above 1
plot(stdloeuf ~ stdltloeuf, data=cons)

# recalculate moeuf my way
cons$moeuf = qchisq(0.95, 2 *(cons$obs_mis+1))/2/cons$exp_mis
f = ecdf(cons$moeuf)
cons$stdltmoeuf = qnorm(f(cons$moeuf))
f = ecdf(cons$oe_mis_upper)
cons$stdmoeuf = qnorm(f(cons$oe_mis_upper))
plot(stdmoeuf ~ stdltmoeuf, data=cons)
# little divergence, only above 2

# build metric using both mis and lof
cons$toeuf = qchisq(0.95, 2 *(cons$obs_mis+cons$obs_lof+1))/2/(cons$exp_mis + cons$exp_lof)
f = ecdf(cons$toeuf)
cons$std_toeuf = qnorm(f(cons$toeuf))

plot(toeuf ~ moeuf, data=cons)
plot(std_toeuf ~ oe_mis_upper, data=cons)


#cons = cons[,c('gene', 'oe_mis', 'oe_mis_upper', 'oe_lof', 'oe_lof_upper')]

dat = inner_join(SFS, cons)

plot(oe_mis_upper ~ beta, data=dat )
cor.test(dat$oe_mis_upper, dat$beta)

plot(stdloeuf ~ beta, data=dat )
plot(std_toeuf ~ beta, data=dat )
m = lm(std_toeuf ~ beta, data=dat )
summary(lm(std_toeuf ~ beta, data=dat ))
cor.test(dat$std_toeuf, dat$beta)


# more model diagnostic
plot(beta ~ log(obs_mis + obs_lof), data=dat)


# new metric
dat$pm = dat$stdloeuf + coef(m)[2]*dat$beta

load('optiRVIS.RData')
optiRVIS = inner_join(optiRVIS, dat)
