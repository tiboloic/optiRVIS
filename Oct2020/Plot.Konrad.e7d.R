# LT 18/10/2020

# plot length boxplot per decile of intolerance metric
require(tidyverse)

gradient_colors = c('#FF2600', '#FF9300', '#FFC000', '#E1C977', '#C2D1ED')
gradient_values = c(0, 0.4, 0.6, 0.8, 1)
gradient_value_deciles = c(0, 2, 5, 7, 9)

metric_v_continuous = function(save_plot=F, metric='oe_lof_upper', y='cds_length', ylim=NULL) {

  bin = paste0(metric, '_bin')
  
  variable_label = function(variable) {
    switch(variable,
           'oe_lof_upper'='LOEUF',
           'psfs'='powerSFS',
           'cds_length'='CDS length',
           variable)
  }
  
  oe_x_label = paste(variable_label(metric), 'decile')
  
  label_function = function(x) {
    y = 10 * x + 5
    ifelse(y %% 20 == 0, paste0(y, '%'), "")
  }
  oe_x_axis = list(xlab(oe_x_label),
                   scale_x_continuous(labels=label_function, breaks=seq(-0.5, 9.5, 1), limits=c(-0.5, 9.5)))
  
  p = gene_data %>%
    filter(!is.na(.data[[bin]])) %>%
    mutate(pLI_bin = ntile(pLI, 10) - 1) %>%
    ggplot + aes_string(x = paste0(metric, '_bin'), group = paste0(metric, '_bin'), y = y) +
    geom_boxplot(aes(fill = .data[[bin]])) + oe_x_axis + scale_y_continuous() +
    scale_fill_gradientn(colors = gradient_colors, values = gradient_values, guide = FALSE) +
    theme_classic(base_size=18) + ylab(variable_label(y)) 
  
  if (!is.null(ylim)) p = p + coord_cartesian(ylim=ylim)
  if (save_plot) {
    pdf(paste0(y, metric, '.pdf'), height=3, width=5)
    print(p)
    dev.off()
  }
  return(p)
}

gene_data =   read_delim(gzfile('../data/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz'), delim = '\t',
                         col_types = do.call(cols, list())) %>%
  mutate(s = mu_lof / p)

# join with psfs
load('../fit.AN.80.RData')
#load('../fitall.HGNC.RData')
#load('../fit.Polyphen.RData')

gene_data = SFS %>% mutate(psfs = -beta) %>% select('gene', 'psfs') %>% inner_join(gene_data)
gene_data = gene_data %>% mutate(psfs_bin=as.numeric(cut(psfs, quantile(psfs, 0:10/10))) - 1)

# plot LOEUF and psfs versus CDS length
plot(metric_v_continuous(ylim=c(0,6000)))
plot(metric_v_continuous(metric='psfs', ylim=c(0,6000)))

# try on log scale
gene_data = gene_data %>% mutate(lcds = log(cds_length))
plot(metric_v_continuous(metric='oe_lof_upper',y='lcds'))
plot(metric_v_continuous(metric='psfs',y='lcds'))
          
# pSFS versus LOEUF
plot(metric_v_continuous(metric='psfs', y='oe_lof_upper'))

plot(metric_v_continuous(metric='oe_lof_upper', y='psfs_bin'))
     