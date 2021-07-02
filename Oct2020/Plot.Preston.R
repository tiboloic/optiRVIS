# LT 20/10/2020

# Preston plot

# linear model
getpowerps = function(beta=1, n=251496) {
  ps = (1:n)^-beta
  ps = ps/sum(ps)
  bins = floor(log(1:n,2))
  dat = data.frame(p = ps, bin=bins)
  dat = dat %>% group_by(bin) %>% summarize(p=sum(p), .groups='drop')
  return(dat$p)
}

preston_plot = function(SFS, i, model='linear') {
  

  gene = SFS[i,'gene']
  SFS_data = SFS[i,] %>% select(c('gene', as.character(0:17))) %>% cbind(getpowerps(SFS[i,'beta']))
    gather(key='octave', value='obs', as.character(0:17), convert=TRUE) %>%
    rbind(cbind(gene=gene, obs=getpowerps(SFS[i,'beta'])))
  
  
  # get predicted proportions
  pred = SFS[i, 'slope'] %>% getpowerps
  preds = t(sapply(SFS$slope, getpowerps))
  colnames(preds) = 0:17
  preds = as_tibble(cbind(beta=SFS$beta, preds))
  preds = preds %>% gather(key='octave', value='pred', as.character(0:17), convert = TRUE)

  # get average model
  dat = dat %>% inner_join(preds)

  fig = ggplot(dat) + geom_line(aes(x = octave, y = pred),lwd=1) +
    geom_point(aes(x = octave, y = obs), size=2, alpha=1) +
    theme_classic(base_size=15) + theme(legend.position = 'top') + scale_y_sqrt() +
    ggtitle(dat$gene) +
    ylab('Proportion of variants')
  
  fig
}

load('../fit.AN.80.RData')
SFSsave = SFS
b = f$env$parList()$b
lsig = f$env$parList()$lsig
SFS = SFS %>% mutate(psfs=-beta, slope = b + beta * exp(lsig))

SFS %>% filter(gene=='MIB1') %>% preston_plot %>% plot
plot(fig)

obs = SFS[i, as.character(0:17)]
pred = sum(obs) * getpowerps(SFS[i, 'beta'])
