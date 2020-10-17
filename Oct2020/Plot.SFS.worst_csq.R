# LT 14/10/2020

# plot SFS fits for synonymous, stop_gained, missense

require(dplyr)
require(tidyr)
require(ggplot2)

# function to bin probs
getpowerps = function(beta=1, n=251496) {
  ps = (1:n)^-beta
  ps = ps/sum(ps)
  bins = floor(log(1:n,2))
  dat = data.frame(p = ps, bin=bins)
  dat = dat %>% group_by(bin) %>% summarize(p=sum(p))
  return(dat$p)
}

load('Fit.SFS.worst_csq.Rdata')

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
  preds = t(sapply(SFS$beta, getpowerps))
  colnames(preds) = 0:17
  preds = as_tibble(cbind(beta=SFS$beta, preds))
  preds = preds %>% gather(key='octave', value='pred', as.character(0:17), convert = TRUE)
  
  dat = dat %>% inner_join(preds)
  #dat=subset(dat, octave > 10)
  fig = ggplot(dat) + geom_line(aes(x = octave, y = pred, color=worst_csq),lwd=1) +
    geom_point(aes(x = octave, y = obs, color=worst_csq), size=4, alpha=0.5) +
    theme_classic() + theme(legend.position = 'top') +
    ggtitle(ifelse(protein_coding, 'protein-coding genes', 'non protein-coding genes'))
  
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
plot(fig1b)