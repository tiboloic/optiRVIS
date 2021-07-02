# LT 22/10/2020

#SFS subsampling formula

N = 1000
n = 10
beta=1.8
ps = 1/(1:1000)^beta
ps = ps /sum(ps)

getps = function(beta=1, n=251496) {
  ps = (1:n)^-beta
  ps/sum(ps)
}

subsampSFS = function(ps, N, n, index=1:n) {
# N size of pop, n size of sample, ps: SFS of pop  
  # to store resulting SFS
  ans = vector('integer', length=n)
  for (j in 1:n) {
    ans[j] = 0
    cat(paste0('j= '), j, '\n')
    for (i in 1:N) {
      ans[j] = ans[j] + N * ps[i] * dhyper(j, m=i, n=N-i, k=n) 
    }
  }
  ans/sum(ans)
}

# simple rescale
rescaleSFS = function(ps, N, n) {
  qs = ps[1:n]
  qs/sum(qs)
  
}

# subsample and rescale are the same only when beta=1
subsamp(ps, 1000,10)
rescaleSFS(ps, 1000, 10)

#N = 251496
N = 10000
ps = getps(1.8)
subsamp(ps, N, 7000)[1]
rescaleSFS(ps, N, 7000)[1]

# test with simulation
# a matrix with random 0/1
M = matrix(rbinom(20000, size=1, prob=0.1), ncol=100)

replicate(1000, table(rowSums(M[,sample(1:100, 10)])))
