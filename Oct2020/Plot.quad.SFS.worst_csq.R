# LT 14/10/2020

# plot SFS fits for synonymous, stop_gained, missense

require(tidyverse)

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

#load('Fit.quad.SFS.worst_csq.Rdata')
load('Fit.fixquad.SFS.worst_csq.Rdata')

plotSFS = function(SFS, csq_list, protein_coding = TRUE) {
  
  # filter on csq_list and protein_coding
  if (protein_coding) {
    dat = SFS %>% filter(worst_csq %in% target_variants & protein_coding=='true')
  } else {
    dat = SFS %>% filter(worst_csq %in% target_variants & protein_coding!='true')
  }
  
  # transform observed counts into obeserved proportions
  counts = dat %>% select(as.character(0:17))
  counts = counts/rowSums(counts)
  dat = dat %>% select(c('worst_csq','beta')) %>% cbind(counts) %>%
    gather(key='octave', value='obs', as.character(0:17), convert=TRUE)
  
  # get predicted proportions
  preds = t(apply(SFS[,c('beta', 'lgamma')], 1, function(l) getpowerps(l[1], exp(l[2]))))
  colnames(preds) = 0:17
  preds = as_tibble(cbind(beta=SFS$beta, preds))
  preds = preds %>% gather(key='octave', value='pred', as.character(0:17), convert = TRUE)
  
  dat = dat %>% inner_join(preds)
  #dat=subset(dat, octave > 10)
  fig = ggplot(dat) + geom_line(aes(x = octave, y = pred, color=worst_csq),lwd=1) +
    geom_point(aes(x = octave, y = obs, color=worst_csq), size=2, alpha=1) +
    theme_classic(base_size=15) + theme(legend.position = 'top') + scale_y_sqrt() +
    ggtitle(ifelse(protein_coding, 'protein-coding genes: gnomAD v2', 'non protein-coding genes')) +
    ylab('Proportion of variants')
  
  fig
}

# colors
color_syn = '#AAAAAA'
color_mis = '#FF6103'
color_lof = '#9D1309'
color_int = '#AAAAAA'
 

#variant_category_colors = c(color_int, color_syn, color_mis, color_lof)
#variant_category_colors = c(color_syn, color_mis, color_lof)

target_variants = c('stop_gained','missense_variant','synonymous_variant', 'intron_variant')
target_variants = c('stop_gained','missense_variant','synonymous_variant')

fig1a = plotSFS(SFS, target_variants, FALSE)
plot(fig1a)

fig1b = plotSFS(SFS, target_variants, protein_coding=TRUE)

fig1b = fig1b+ggtitle('Quadratic power law') + ylab('Proportion of variants') +
  xlab('Allele frequency (octave)')
plot(fig1b)
pdf(paste0('Plot.quad.SFS.worst_csq.pdf'), height=3, width=5)
print(fig1b)
dev.off()
