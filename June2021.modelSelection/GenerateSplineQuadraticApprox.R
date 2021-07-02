# LT 27/06/2011

# A smooth spline approximation of the quadratic power law for fixed $\gamma$

require(mgcv)

# generate grid on range of parameters
range(fit.gamma.all$par[-1])

# numbeer of points for each dimension of the parameter space
p=200

betas = seq(1.5, 2.5, length.out = p)

n_max = 251496
ns = c(log(1:100), log(floor(exp(seq(log(101), log(n_max), length.out=(p-100))))))

getquadps = function(beta=1, gamma = 0.03135216, n=251495) {
  ac = 1:n
  ps = ac^-beta * (ac^gamma)^log(ac)
  ps = ps/sum(ps)
  bins = floor(log(1:n,2))
  dat = data.frame(p = ps, bin=bins)
  dat = dat %>% group_by(bin) %>% summarize(p=sum(p), .groups='drop')
  return(dat$p)
}

fs= sapply(betas, getquadps)

dat = data.frame(f=as.vector(fs), beta = rep(betas, each=nrow(fs)), logn=rep(log(c(2^(1:17)-1,n_max)), ncol(fs)))

m = gam(f ~ te(beta, logn), data=dat)
m = gam(f ~ te(beta, logn, k=c(20,10)), data=dat)
dat = cbind(dat, fit = fitted(m))
plot(f~logn, data=dat[1:18,], ylim=c(0,0.8))
lines(fit~logn, data=dat[1:18,], col=2, lwd=2)
points(f~logn, data=dat[100*18+1:18,])
lines(fit~logn, data=dat[100*18+1:18,], col=2, lwd=2)
points(f~logn, data=dat[199*18+1:18,])
lines(fit~logn, data=dat[199*18+1:18,], col=2, lwd=2)

# problem is predicted values go negative sometimes
# transform f using qnorm

m.trans = gam(qnorm(f) ~ te(beta, logn, k=c(150,15), bs='ts'), data=dat)
range(residuals(m.trans))
m.trans
k.check(m.trans)
range(dat$f-pnorm(fitted(m.trans)))
dat$fit = pnorm(fitted(m.trans))

m1 = gam(qnorm(f, mean=2) ~ te(beta, logn, k=c(50,10), bs='ts'), data=dat)
k.check(m1)
m2 = gam(qnorm(f, mean=2) ~ te(beta, logn, k=c(50,20), bs='ts'), data=dat)
m2 = gam(qnorm(f, mean=2) ~ te(beta, logn, k=c(50,18), bs='ts'), data=dat)
k.check(m2)
m3 = gam(qnorm(f, mean=2) ~ te(beta, logn, k=c(100,18), bs='ts'), data=dat)
k.check(m3)
range(residuals(m2))
range(dat$f-pnorm(fitted(m2), mean=2))
range(pnorm(fitted(m2), mean=2))
range(dat$f)
m4 = gam(qnorm(f, mean=2) ~ te(beta, logn, k=c(30,18), bs='ts'), data=dat)
k.check(m4)
m5 = gam(qnorm(f, mean=2) ~ te(beta, logn, k=c(70,18), bs='ts'), data=dat)
k.check(m5)
range(pnorm(fitted(m5), mean=2))
range(dat$f-pnorm(fitted(m5), mean=2))
dat$fit = pnorm(fitted(m5), mean=2)
plot(f~logn, data=dat[1:18,], ylim=c(0,0.8))
lines(fit~logn, data=dat[1:18,], col=2, lwd=2)
points(f~logn, data=dat[100*18+1:18,])
lines(fit~logn, data=dat[100*18+1:18,], col=2, lwd=2)
points(f~logn, data=dat[199*18+1:18,])
lines(fit~logn, data=dat[199*18+1:18,], col=2, lwd=2)
m6 = gam(qnorm(f, mean=2) ~ te(beta, logn, k=c(20,18), bs='cs'), data=dat)
k.check(m6)


save.image(file='SplineQuadraticApprox.RData')
