# LT 13/10/2020

# SFS analysis of gnomad: for whole datatset instead of per gene

require(dplyr)
require(tidyr)

# load gnomad extracted variants count
gnomad = read.delim(file='gnomADSFS.worst_csq.all.tsv')
sum(gnomad$obs)
# 10 450 462 variants

# function to bin probs
getpowerps = function(beta=1, n=251496) {
   ps = (1:n)^-beta
   ps = ps/sum(ps)
   bins = floor(log(1:n,2))
   dat = data.frame(p = ps, bin=bins)
   dat = dat %>% group_by(bin) %>% summarize(p=sum(p))
   return(dat$p)
}
 
# fit powerlaw using multinomial likelihood
LL = function(pars, obs, LLonly=TRUE) {
   # calculate probabilities per bin
   ps = getpowerps(pars)
   N = sum(obs)
   pred = N * ps
   
   LL = -dmultinom(obs, prob=ps, log=TRUE)
   if (LLonly)
    return(LL)
   else
     return(ps)
 }

# one SFS per line
SFS = gnomad %>% spread(key=cat, value=obs, fill=0)

# fit model to each SFS
fit = apply(SFS, 1, function(l) {
   obs = as.numeric(l[-(1:2)])
   onefit = optimize(LL, c(1,2.5), obs=obs)
   c(beta = onefit$minimum)
})
SFS$beta = fit

save(SFS, file='fit.SFS.worst_csq.Rdata')
