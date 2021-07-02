# LT 19/06/2021
# based on FIt.quad.SFS.worst_csq.R
# fit quadratic powerlaw on multiple SFSs with full model, common gamma, common beta, etc.

# LT 17/10/2020

# quadratic version
# try a quadratic on log-log scale model, look at evidence for diddereing gamma across mutation class

# LT 13/10/2020

# SFS analysis of gnomad: for whole datatset instead of per gene

require(dplyr)
require(tidyr)
require(tidyverse)

# load gnomad extracted variants count
gnomad = read.delim(file='../Oct2020/gnomADSFS.worst_csq.all.tsv')
sum(gnomad$obs)
# 10 450 462 variants

# function to bin probs
getpowerps = function(beta=1, gamma = 0, n=251496) {
   ac = 1:n
   ps = ac^-beta * (ac^gamma)^log(ac)
   ps = ps/sum(ps)
   bins = floor(log(1:n,2))
   dat = data.frame(p = ps, bin=bins)
   dat = dat %>% group_by(bin) %>% summarize(p=sum(p), .groups='drop')
   return(dat$p)
}
 
# fit second order power law using multinomial likelihood
oneLL = function(pars, obs, LLonly=TRUE) {
   # calculate probabilities per bin
   ps = getpowerps(pars[1], pars[2])
   N = sum(obs)
   pred = N * ps
   
   LL = -dmultinom(obs, prob=ps, log=TRUE)
   if (LLonly)
    return(LL)
   else
     return(ps)
 }

manyLL.full = function(pars, obs, LLonly=TRUE) {
   # pars is a matrix with two columns and as many rows as SFSs
   pars = matrix(pars, nrow(obs), 2)
   LL=0
   for (i in 1:nrow(obs))
      LL = LL + oneLL(pars[i,], obs[i,], TRUE)
   LL
}

manyLL.comon = function(pars, obs, LLonly=TRUE) {
   # pars is a vector of length 2: [beta, gamma]
   LL=0
   for (i in 1:nrow(obs))
      LL = LL + oneLL(pars, obs[i,], TRUE)
   LL
}

manyLL.gamma = function(pars, obs, LLonly=TRUE) {
   # pars is a vector of length nrow+1: [gamma, betas]
   LL=0
   for (i in 1:nrow(obs))
      LL = LL + oneLL(c(pars[i+1], pars[1]), obs[i,], TRUE)
   LL
}

# one SFS per line
SFS = gnomad %>% spread(key=cat, value=obs, fill=0)


# select 4 classes: intronic, synonymous, stop-gain, missenes
SFS2 = subset(SFS, (worst_csq %in% c('intron_variant', 'missense_variant', 'stop_gained', 'synonymous_variant')))
SFS2 = subset(SFS2, protein_coding == 'true')

# fit full model
#fit.full = nlminb(rep(c(1.8,0), each=nrow(SFS2)), manyLL.full, obs=SFS2[,-(1:2)])
fit.full = nlminb(rep(c(1.8,0), each=nrow(SFS2)), manyLL.full, obs=SFS2[,-(1:2)], control=list(trace=1))

# fit common gamma model
fit.gamma = nlminb(c(0,rep(1.8, nrow(SFS2))), manyLL.gamma, obs=SFS2[,-(1:2)], control=list(trace=1))

# fit common gamma beta model
fit.cmon = nlminb(c(0,1.8), manyLL.comon, obs=SFS2[,-(1:2)], control=list(trace=1))

# fit full model
fit.full.all = nlminb(rep(c(1.8,0), each=nrow(SFS)), manyLL.full, obs=SFS[,-(1:2)], control=list(trace=1))

# fit common gamma model
fit.gamma.all = nlminb(c(0.03,rep(1.9, nrow(SFS))), manyLL.gamma, obs=SFS[,-(1:2)], control=list(trace=1))

# fit reached max eval limit, restart from best fit param
fit.gamma.all = nlminb(fit.gamma.all$par, manyLL.gamma, obs=SFS[,-(1:2)], control=list(trace=1))
fullpar = as.vector(cbind(fit.gamma.all$par[-1], fit.gamma.all$par[1]))
fit.full.all = nlminb(fullpar, manyLL.full, obs=SFS[,-(1:2)], control=list(trace=1, eval.max=1000, iter.max=1000))
# in a couple of iterations, lost 600 units of likelihood: there is no doubt, 
# model selection is going to favor the full model

# add a fit with gamma=0
fit.gamma0.all = nlminb(fit.gamma.all$par[-1], function(par, ...) manyLL.gamma(c(0,par),...), obs=SFS[,-(1:2)], control=list(trace=1))

#save.image(file='fit.modelselection.Rdata')
f
# plot some fits, with fixed gamma and variable gamma
plotfit = function(idx) {
   n = sum(SFS[idx, -(1:2)])
   plot(0:17,SFS[idx, -(1:2)], log='y', xlab='octave', ylab='# variant (log scale)', main=paste(SFS[idx,1], '-', SFS[idx,2]))
   lines(0:17,(n*getpowerps(fit.full.all$par[idx], fit.full.all$par[idx+41])), lwd=2, col=2)
   lines(0:17,(n*getpowerps(fit.gamma.all$par[idx+1], fit.gamma.all$par[1])), lwd=2, col=3)
   lines(0:17,(n*getpowerps(fit.gamma0.all$par[idx], 0)), lwd=2, col=4)
}
# synonymous variants coding
plotfit(37)
# syn non coding
plotfit(36)
# missense coding
plotfit(17)
# missense non coding
plotfit(16)
# stop gained
plotfit(31)
#intro
plotfit(13)

sapply(1:nrow(SFS), function(i) plotfit(i))

# refit with optim, get hessian
fit.gamma.all.optim = optim(fit.gamma.all$par, manyLL.gamma, obs=SFS[,-(1:2)], control=list(trace=1), hessian=TRUE)


# hessian takes too long, try a smaller fit
SFS3 = subset(SFS, protein_coding=='true' & rowSums(SFS[,-(1:2)]) > 3000)
keep = which(SFS[,'protein_coding']=='true' & rowSums(SFS[,-(1:2)]) > 3000)

# fit common gamma
fit.gamma.SFS3.optim = optim(fit.gamma.all$par[c(1,keep+1)], manyLL.gamma, obs=SFS3[,-(1:2)], control=list(trace=1,maxit=10000), hessian=TRUE)

# try to solve hessian
gamma.SFS3.se = sqrt(diag(solve(fit.gamma.SFS3.optim$hessian)))

# order by beta
ordering = order(fit.gamma.SFS3.optim$par[-1])

res = SFS3[ordering,]
res$beta = fit.gamma.SFS3.optim$par[-1][ordering]
res$beta.se = gamma.SFS3.se[-1][ordering]
res$beta.inf = res$beta + qnorm(0.025) * res$beta.se
res$beta.sup = res$beta + qnorm(1-0.025) * res$beta.se

psf_quadratic_plot = res %>% ggplot +
   aes(x = reorder(worst_csq, beta), y = beta, ymin = beta.inf, ymax = beta.sup) + 
#   geom_hline(aes(yintercept=beta), data=res, size=1, color=variant_category_colors,show.legend=TRUE) +
   geom_pointrange(size=1) + geom_point() + xlab('vep worst consequence') + ylab('powerSFS') +
   ylim(1.7, 2.2) +
   theme_classic() + theme(legend.position="none", text=element_text(size=18)) +
   theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
   ggtitle('quadratic power law', subtitle = '')

plot(psf_quadratic_plot)