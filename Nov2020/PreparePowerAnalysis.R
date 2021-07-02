# LT 21/11/2020

# prepare power analysis

require(tidyverse)

load("fit.AN.80.RData")

# look at how uncertainty on beta changes with #variants
tot_variants = rowSums(SFS[,2:19])
sd_betas = summary(mysd, 'random', p.val=TRUE)


# stderr beta
plot(sd_betas[,2] ~ log(tot_variants))
plot(sd_betas[,2] ~ tot_variants, log='x', xlab='# variants', ylab='stderr')

# z-score
plot(sd_betas[,3] ~ log(tot_variants))

bs = f$env$parList()$bdevs * exp(f$env$parList()$lsig) + f$env$parList()$b

# compare distribution of Lof to missense
cons = read.delim('data/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz')

down = read.delim('data/gnomad.v2.1.1.lof_metrics.downsamplings.txt.bgz')

# plot Lof, missense and synonymopus as function of sample size

p2c = down %>%
  filter(variant_type == 'LoF' & func == 'sum') %>%
  ggplot + aes(x = downsampling, y = count, linetype = obs_exp, color = variant_type) + geom_line(color=color_lof) +
  theme_classic() + scale_linetype_manual(values=linetypes, name=NULL) +
  theme(legend.justification = c(0, 1), legend.position = c(0.01, 1)) +
  scale_x_continuous(label=comma, breaks=c(seq(0, 8e4, 4e4), max(sumbounds$x))) + 
  scale_y_continuous(label=comma) + #, breaks=c(seq(0, 4e5, 1e5), max(sumbounds$y))) +
  # scale_y_continuous(label=comma, breaks=c(0, 1e5, 3e5, 4e5, sumbounds$y)) +
  xlab('Sample size') + ylab('Total number of pLoF SNVs')

satur = down %>% filter(canonical=='true' & pop == 'global') %>% group_by(downsampling) %>% 
  summarize(obs_syn=sum(obs_syn, na.rm=TRUE), obs_mis=sum(obs_mis, na.rm=TRUE), obs_lof=sum(obs_lof, na.rm=TRUE))
matplot(satur[,1], satur[,2:4])
satur[,2:4] = t(t(satur[,2:4])/apply(satur[,2:4],2,max))
matplot(satur[,1], satur[,2:4])
