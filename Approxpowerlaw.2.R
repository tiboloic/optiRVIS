# LT 1/06/2020
# redo pad'e appromixation on larger range of slope
# fit on log scale to avoild negative probabilites

# LT 181/12/2019
# approximation of pi for efficient power SFS

load('fitall.Rdata')

beta = f$env$parList()$b
sigbeta = exp(f$env$parList()$lsigb)

betamin = beta - 8 * exp(f$env$parList()$lsigb)
betamax = beta + 8 * exp(f$env$parList()$lsigb)

betas = seq(betamin, betamax, length.out = 100)

getpowerps = function(beta=1, n=251496) {
  ps = (1:n)^-beta
  ps = ps/sum(ps)
  bins = floor(log(1:n,2))
  dat = data.frame(p = ps, bin=bins)
  dat = dat %>% group_by(bin) %>% summarize(p=sum(p))
  return(dat$p)
}


dat[,1]
getpowerps(betamin)


# padé approximant
f = function(par, fp=2, fq=2, x, y) {
  p = fp
  q=fq
  if (length(par) != p+q+1) stop('wrong parameter vector\n')
  num = par[1:(p+1)]
  # by convention firt term of denominator is 1
  den = c(1,par[(p+2):length(par)])
  
  X = do.call(cbind, lapply(0:p, function(i) x^i))
  Y = do.call(cbind, lapply(0:q, function(i) x^i))


  pred = (X %*% num) / (Y %*% den)
  y-pred
}

# try a fit
# using Levenberg–Marquardt algorithm
library(minpack.lm)

x = betas
coefs= NULL
for(i in 1:18) {
  m = nls.lm(c(1,0,0,0,0,0,0), lower=NULL, upper=NULL, fn=f, fp=3, fq=3, x=x, y=dat[i,], control=list(maxiter=1000))
  plot(x, dat[i,], main=i, log='y')
  lines(x, dat[i,]-residuals(m), col=2)
  coefs = rbind(coefs, coef(m))
}  
# padé(3,3) approximant is a winner

# problem with negative probabilities, fit on the log-scale
x = betas
coefs= NULL
for(i in 1:18) {
  m = nls.lm(c(1,0,0,0,0,0,0), lower=NULL, upper=NULL, fn=f, fp=3, fq=3, x=x, y=log(dat[i,]), control=list(maxiter=1000))
  plot(x, dat[i,], main=i, log='y')
  lines(x, exp(log(dat[i,])-residuals(m)), col=2)
  coefs = rbind(coefs, coef(m))
}  
# padé(3,3) approximant is a winner

# my approximation
approxps = function(beta, coefs) {
  apply(coefs, 1, function(l) exp((l[1]+l[2]*beta+l[3]*beta^2+l[4]*beta^3)/(1+l[5]*beta+l[6]*beta^2+l[7]*beta^3)))
}

for (nbeta in rnorm(10, beta, 3*sigbeta)) {
  plot(approxps(nbeta, coefs) ~ getpowerps(nbeta), log='xy', main = nbeta)
  cat(paste('sum: ', sum(approxps(nbeta, coefs)), '\n'))
  abline(a=0, b=1)
}

#predicted on beta range
pred = sapply(betas, approxps, coefs=coefs)

save(coefs, file='pade2.RData')
write.csv(coefs, file='pade2.csv')
# ALL GOOD