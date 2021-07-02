# LT 17/10/2020

# quadratic version
# try a quadratic on log-log scale model, look at evidence for diddereing gamma across mutation class

# LT 13/10/2020

# SFS analysis of gnomad: for whole datatset instead of per gene

require(dplyr)
require(tidyr)

# load gnomad extracted variants count
gnomad = read.delim(file='gnomADSFS.worst_csq.all.tsv')
sum(gnomad$obs)
# 10 450 462 variants

# function to bin probs
getpowerps = function(beta=1, gamma = exp(-3.5), n=251496) {
   ac = 1:n
   ps = ac^-beta * (ac^gamma)^log(ac)
   ps = ps/sum(ps)
   bins = floor(log(1:n,2))
   dat = data.frame(p = ps, bin=bins)
   dat = dat %>% group_by(bin) %>% summarize(p=sum(p), .groups='drop')
   return(dat$p)
}
 
# fit powerlaw using multinomial likelihood
LL = function(pars, obs, LLonly=TRUE) {
   # calculate probabilities per bin
#   ps = getpowerps(pars[1], exp(pars[2]))
   ps = getpowerps(pars[1])
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
   cat('fitting ... \n')
   obs = as.numeric(l[-(1:2)])
   onefit = nlminb(c(beta=1.8, lgamma=-3.5), LL, obs=obs)
   onefit$par
})
SFS= cbind(SFS, t(fit))

save(SFS, file='fit.fixquad.SFS.worst_csq.Rdata')
