# LT 23/11/2020
# Power Analysis


## for constant.r

# function to bin probs
getpowerps = function(beta=1, n=251496) {
  ps = (1:n)^-beta
  ps = ps/sum(ps)
  bins = floor(log(1:n,2))
  dat = data.frame(p = ps, bin=bins)
  dat = dat %>% group_by(bin) %>% summarize(p=sum(p), .groups='drop')
  return(dat$p)
}

simSFS = function(p, q, beta=1, n=251496) {
  # p number of simulations
  # q number of variants per sample
  res = rmultinom(p, q, getpowerps(beta, n))
  res = t(res)
  colnames(res) = 0:(ncol(res)-1)
  rownames(res) = 1:p
  res
}

# multinomial likelihood
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

fitoneSFS = function(SFS) {
  
}

beta.syn = 1.645

deltabeta = 0.5

N = 500

ns = c(10, 20, 50, 100)

res = sapply(ns, function(n) {
  cat('n=',n,'\n')
  SFS = simSFS(N, n, beta =beta.syn+deltabeta)

  # fit model to each SFS
  fit = apply(SFS, 1, function(l) {
    obs = as.numeric(l)
    onefit = optimize(LL, c(1,5), obs=obs)
    c(beta = onefit$minimum)
  })
  fit
})

colnames(res) = ns
colMeans(res)
apply(res, 2, sd)
apply(res, 2, median)

save(res, file='PowerAnalysis.Rdata')

# produce quick figure for talk
boxplot(res, outline=FALSE, xlab="# variants per gene", ylab='intolerance')
abline(h=1.645, lwd=2, col=3)
abline(h=beta.syn+deltabeta, lwd=1)
#abline(h=1.941, lwd=2, col=2)
