# LT 181/12/2019
# approximation of pi for efficient power SFS

load('fit100.Rdata')

beta = f$env$parList()$b
sigbeta = exp(f$env$parList()$lsigb)

betamin = beta - 5 * exp(f$env$parList()$lsigb)
betamax = beta + 5 * exp(f$env$parList()$lsigb)

betas = seq(betamin, betamax, length.out = 100)

getpowerps = function(beta=1, n=251496) {
  ps = (1:n)^-beta
  ps = ps/sum(ps)
  bins = floor(log(1:n,2))
  dat = data.frame(p = ps, bin=bins)
  dat = dat %>% group_by(bin) %>% summarize(p=sum(p))
  return(dat$p)
}

dat = sapply(betas, getpowerps)

dat[,1]
getpowerps(betamin)

# fit quadratic to p1
i=2
x = 0:99
m.1 = lm(dat[i,]~ x + I(x^2))
summary(m.1)
plot(x, dat[i,])
lines(x, fitted(m.1))

# not very good for i=2, i=18 for example

# try quadratic on log
m.1 = lm(log(dat[i,])~ x + I(x^2))
summary(m.1)
plot(x, dat[i,])
lines(x, exp(fitted(m.1)), col=2)

# better for large i, but i=2 bad

# try X^3
m.1 = lm(dat[i,]~ x + I(x^2) + I(x^3))
summary(m.1)
plot(x, dat[i,])
lines(x, fitted(m.1), col=2)

# try X^3 on log
m.1 = lm(log(dat[i,])~ x + I(x^2) + I(x^3))
summary(m.1)
plot(x, dat[i,])
lines(x, exp(fitted(m.1)), col=2)

# pretty good, fit all models. worst is i=2

# try padé approximant
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
  m = nls.lm(c(1,0,0,0,0,0,0), lower=NULL, upper=NULL, fn=f, fp=3, fq=3, x=x, y=dat[i,])
  plot(x, dat[i,], main=i)
  lines(x, dat[i,]-residuals(m), col=2)
  coefs = rbind(coefs, coef(m))
}  
# padé(3,3) approximant is a winner

# my approximation
approxps = function(beta, coefs) {
  apply(coefs, 1, function(l) (l[1]+l[2]*beta+l[3]*beta^2+l[4]*beta^3)/(1+l[5]*beta+l[6]*beta^2+l[7]*beta^3))
}

for (nbeta in rnorm(10, beta, sigbeta)) {
  plot(approxps(nbeta, coefs) ~ getpowerps(nbeta), log='xy', main = nbeta)
  cat(paste('sum: ', sum(approxps(nbeta, coefs)), '\n'))
  abline(a=0, b=1)
}

save(coefs, file='pade.RData')
write.csv(coefs, file='pade.csv')
# ALL GOOD