
# test on folded SFS
fold = function(obs) {
  if (size(obs) %% 2 == 0) # even
  {
    mid = size(obs) / 2
    obs[1:mid] + obs[(mid+1):size(obs)]
  }
  else # odd
  {
    mid = size(obs) %/% 2
    obs[]
  }
    
}

SFSplot = function(i) {
  plot(SFSobs[i,], log='y', main=optiRVIS$gene[i])
  lines(approxps(1.7945 + 0.1186269 * optiRVIS$beta[i], coefs) * sum(SFSobs[i,]), col=2, lwd=2)
  lines(approxps(1.7945, coefs) * sum(SFSobs[i,]), col=1, lwd=1)
}

par2fit = funcion(beta) {

}

load('pade.RData')

approxps = function(beta, coefs) {
  apply(coefs, 1, function(l) (l[1]+l[2]*beta+l[3]*beta^2+l[4]*beta^3)/(1+l[5]*beta+l[6]*beta^2+l[7]*beta^3))
}

load('feb2020')
# all data
SFSobs = as.matrix(optiRVIS[,24:(24+17)])
barplot(colSums(SFSobs), log='y')
barplot(colSums(SFSobs))

plot(colSums(SFSobs), log='y')
plot(colSums(SFSobs), log='y', ylim=c(1e3,1e7))
lines(approxps(1.7945, coefs) * sum(SFSobs))

#HLA-DRB1 (outlier tolerant)
SFSplot(which(optiRVIS$gene == 'HLA-DRB1'))

# MED13 is gene with lowest LOEUF
SFSplot(which(optiRVIS$gene == 'MED13'))

# most intolerant according to powerSFS PRPF40A
SFSplot(which.max(optiRVIS$beta))

# most 10 most intolerant according to powerSFS
sapply(order(optiRVIS$beta, decreasing = T)[1:10], SFSplot)

# CLCT, DHX9 is an example where powesfs works very well

# for some genes, he fit does not make sense: MIA, UNG, SEC24D, POLR1C
# needs invetigation

# most 10 most tolerant according to powerSFS
sapply(order(optiRVIS$beta, decreasing = F)[1:10], SFSplot)

library(dplyr)
# check padé for low beta
getpowerps = function(beta=1, n=251496) {
  ps = (1:n)^-beta
  ps = ps/sum(ps)
  bins = floor(log(1:n,2))
  dat = data.frame(p = ps, bin=bins)
  dat = dat %>% group_by(bin) %>% summarize(p=sum(p))
  return(dat$p)
}

plot(getpowerps(0.9), approxps(0.9, coefs))
abline(a=0,b=1)
