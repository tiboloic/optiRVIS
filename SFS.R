# LT 12/12/2019

# SFS analysis of gnaomad


# load gnomad extracted variants count
gnomad = read.delim(file='gnomADSFS.tsv')
sum(gnomad$obs)
# 6150902 variants (LoF + missense)

SFS = gnomad %>% group_by(cat) %>% summarize(obs = sum(obs))

# function to bin probs
getpowerps = function(beta=1, n=251496) {
   ps = (1:n)^-beta
   ps = ps/sum(ps)
   bins = floor(log(1:n,2))
   dat = data.frame(p = ps, bin=bins)
   dat = dat %>% group_by(bin) %>% summarize(p=sum(p))
   return(dat$p)
}
 
 # fit powerlaw
LL = function(pars, obs, LLonly=TRUE) {
   # calculate probabilities per bin
   ps = (1:251496)^-pars
   ps = ps/sum(ps)
   bins = floor(log(1:251496,2))
   dat = data.frame(p = ps, bin=bins)
   dat = dat %>% group_by(bin) %>% summarize(p=sum(p))
   pred = sum(obs) * dat$p
   if (LLonly)
    return(sum((obs-pred)^2))
   else
     return(pred)
 }

fit = optimize(LL, c(1,2), obs=SFS$obs)
pred = LL(fit$minimum, obs=SFS$obs, LLonly = FALSE)
 
plot(obs~cat, data=SFS)
points(0:17, pred, col=2, pch=2)

test = subset(gnomad, gene=='SMC3')
plot(obs~cat, data=test)

fitTTN = optimize(LL, c(1,2), obs=test$obs)
predTTN = LL(fitTTN$minimum, obs=test$obs, LLonly = FALSE)
meanpred = LL(fit$minimum, obs=test$obs, LLonly = FALSE)
points(0:17, meanpred, col=3, pch=3)
