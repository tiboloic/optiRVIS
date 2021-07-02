# LT 17/10/2020

# plot Konrad fig3 (gene lists per decile)

require(tidyverse)

color_benign = '#87A4DC'
color_dominant = '#F0810F'
color_recessive = '#FFBB00'
color_lof = k_lof = '#9D1309'


load_all_gene_list_data = function(list_dir = 'gene_lists/lists/') {
  all_files = list.files(list_dir, '.+tsv')
  if (length(all_files) == 0) {
    if (length(all_files) == 0) system('git clone https://github.com/macarthur-lab/gene_lists.git')
    all_files = list.files(list_dir, '.+tsv')
  }
  gene_lists = map_df(all_files, 
                      function(x) read_tsv(paste(list_dir, x, sep = '/'), 
                                           col_names = F, col_types = cols()
                      ) %>% 
                        transmute(gene = toupper(X1), gene_list=str_sub(x, end=-5)))
  return(gene_lists)
}

get_ko_gene_lists = function(list_dir = 'ko_gene_lists/') {
  gene_lists = map_df(list.files(list_dir, 'list_.+\\.tsv'), 
                      function(x) read_tsv(paste(list_dir, x, sep = '/'), 
                                           col_names = F, col_types = cols()
                      ) %>% 
                        transmute(gene = toupper(X1), gene_list=str_sub(x, 6, -5)))
  return(gene_lists)
}

gene_list_colors = c(
  'Haploinsufficient' = color_lof,
  'Essential Genes' = '#CB0000',
  'Autosomal Dominant' = color_dominant,
  'Autosomal Recessive' = color_recessive,
  'Olfactory Genes' = color_benign,  # '#598234',
  'Background' = 'lightgray',
  'Universe' = 'lightgray'
)

load_constraint_data = function(level='gene', loftee=T) {
  fname = paste0('gnomad.v2.1.1.lof_metrics.by_', level, '.txt.bgz')
  subfolder = ''
  local_name = ''
  if (!loftee) {
    subfolder = 'other_cuts/no_loftee/'
    local_name = paste0('gnomad.v2.1.1.lof_metrics.no_loftee.by_', level, '.txt.bgz')
  }
  fname = get_or_download_file(fname, subfolder, local_name)
  if (level == 'transcript') {
    col_list = list(canonical=col_logical())
  } else {
    col_list = list()
  }
  read.delim(gzfile(fname), delim = '\t',
             col_types = do.call(cols, col_list)) %>%
    mutate(s = mu_lof / p) %>%
    return
}

#constraint_metric_name = 'LOEUF'
#oe_x_label = paste(constraint_metric_name, 'decile')
#label_function = function(x) {
#  y = 10 * x + 5
#  ifelse(y %% 20 == 0, paste0(y, '%'), "")
#}
#oe_x_axis = list(xlab(oe_x_label),
#                 scale_x_continuous(labels=label_function, breaks=seq(-0.5, 9.5, 1), limits=c(-0.5, 9.5)))


######
# Konrad function for spectrum
# modified to plot different metrics
gene_list_spectrum = function(save_plot=F,
                              gene_lists_to_plot=c('Haploinsufficient', 'Autosomal Recessive', 'Olfactory Genes'),
                              gene_lists_to_skip=c(''), metric='oe_lof_upper') {
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
  
    or_genes = gene_data %>%
    filter(grepl('^OR', gene)) %>%
    transmute(gene = gene, gene_list = 'Olfactory Genes')
  
  gene_lists = load_all_gene_list_data() %>%
    bind_rows(or_genes)
  

  gene_list_spectrum_data = gene_data %>%
    left_join(gene_lists) %>%
    mutate(gene_list = if_else(grepl('haploinsufficiency', gene_list), 'Haploinsufficient', gene_list),
           presence = 1) %>%
    complete(gene_list, .data[[bin]], fill = list(presence = 0)) %>%
    count(gene_list, .data[[bin]], wt = presence) %>%
    group_by(gene_list) %>%
    mutate(prop_in_bin = n / sum(n)) %>% 
    ungroup %>%
    mutate(
      gene_list = fct_recode(gene_list, 'Autosomal Recessive' = "all_ar", 'Autosomal Dominant' = "all_ad"),
      gene_list = fct_relevel(gene_list, 'Haploinsufficient'))
  

  top_legend = max(gene_list_spectrum_data$prop_in_bin)# + ko_data$sem)
  top_legend = max(gene_list_spectrum_data$prop_in_bin)# + ko_data$sem)
  gene_list_colors_plot = gene_list_colors
  p3a = gene_list_spectrum_data %>%
    filter(gene_list %in% gene_lists_to_plot) %>%
    ggplot + aes(x = .data[[bin]], y = prop_in_bin, color = gene_list, fill = gene_list) + 
    # geom_line(lwd=1.5) + 
    geom_bar(position='dodge', stat='identity', width=0.9) + 
    theme_classic(base_size=15) + 
    scale_color_manual(values=gene_list_colors_plot, guide=F) + # name=NULL) +
    scale_fill_manual(values=gene_list_colors_plot, guide=F) + # name=NULL) +
    annotate('text', 4.5, top_legend, hjust=0.5, vjust=1, 
             label='Haploinsufficient', color=gene_list_colors['Haploinsufficient']) +
    annotate('text', 4.5, top_legend*0.88, hjust=0.5, vjust=1, 
             label='Autosomal Recessive', color=gene_list_colors['Autosomal Recessive']) +
    annotate('text', 4.5, top_legend*0.76, hjust=0.5, vjust=1, 
             label='Olfactory Genes', color=gene_list_colors['Olfactory Genes']) +
#    ylab('Percent of gene list') + oe_x_axis + scale_y_continuous(labels=percent_format(accuracy = 1)) +
    ylab('Percent of gene list') + oe_x_axis + scale_y_continuous() + coord_cartesian(ylim=c(0, 0.5))
    theme(legend.justification = c(0.5, 0.8), legend.position = c(0.5, 1))
  
  if (save_plot) {
    pdf('3a_gene_list.pdf', height = 3, width = 5)
    print(p3a)
    dev.off()
  }
  return(p3a)
}

get_proportion_in_gene_lists = function(gene_data, metric='oe_lof_upper') {
  bin = paste0(metric, '_bin')
  gene_lists = get_ko_gene_lists()
  
  gene_lists %>%
    inner_join(gene_data) %>%
    count(gene_list, .data[[bin]]) %>%
    add_count(gene_list, wt = n, name = 'nn') %>%
    mutate(prop_in_bin = n / nn,
           sem = 1.96 * sqrt(prop_in_bin * (1 - prop_in_bin) / nn)) %>%
    return
}

# mouse KO
mouse_ko_comparison = function(save_plot=F) {
  gene_data %>%
    left_join(get_ko_gene_lists() %>% filter(gene_list == 'mouse_het_lethal_genes')) %>%
    mutate(mouse_ko = !is.na(gene_list)) %>%
    glm(mouse_ko ~ oe_lof_upper + cds_length, ., family='binomial') %>%
    summary
  
  p = get_proportion_in_gene_lists(gene_data) %>%
    filter(gene_list == 'mouse_het_lethal_genes') %>%
    ggplot + aes(x = oe_lof_upper_bin, y = prop_in_bin, 
                 ymin = prop_in_bin - sem, ymax = prop_in_bin + sem) + 
    geom_bar(stat='identity', fill = color_lof, width=0.6) +
    theme_classic() + 
    ylab('Percent of mouse\nhet lethal knockout genes') + oe_x_axis + 
    scale_y_continuous(labels=percent_format(accuracy = 1))
  
  t_test_gene_list(gene_data, 'mouse_het_lethal_genes') %>% print
  
  if (save_plot) {
    pdf('3c_mouse_ko.pdf', height = 3, width = 5)
    print(p)
    dev.off()
  }
  return(p)
}

# cell KO
cell_ko_comparison = function(save_plot=F, metric='oe_lof_upper') {
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
  
  # gene_data %>%
  #   left_join(get_ko_gene_lists() %>% filter(gene_list == 'CEGv2')) %>%
  #   mutate(cell_essential = !is.na(gene_list)) %>%
  #   glm(cell_essential ~ oe_lof_upper + cds_length, ., family='binomial') %>%
  #   summary
  # gene_data %>%
  #   left_join(get_ko_gene_lists() %>% filter(gene_list == 'NEGv1')) %>%
  #   mutate(cell_non_essential = !is.na(gene_list)) %>%
  #   glm(cell_non_essential ~ oe_lof_upper + cds_length, ., family='binomial') %>%
  #   summary
  
  lists_to_plot = c('CEGv2' = color_lof,
                    'NEGv1' = color_benign)
  lists_labels = c('CEGv2' = 'Cell essential',
                   'NEGv1' = 'Cell non-essential')

  
  ko_data = get_proportion_in_gene_lists(gene_data, metric=metric) %>%
    filter(gene_list %in% names(lists_to_plot))
  
 # t_test_gene_list(gene_data, 'CEGv2') %>% print
#  t_test_gene_list(gene_data, 'NEGv1') %>% print
  
  top_legend = max(ko_data$prop_in_bin)# + ko_data$sem)
  p = ko_data %>%
    ggplot + aes(x = .data[[bin]], y = prop_in_bin, fill = gene_list,
                 ymin = prop_in_bin - sem, ymax = prop_in_bin + sem) + 
    geom_bar(stat = 'identity', position = 'dodge') + 
    #geom_bar(stat='identity', fill = color_lof, width=0.6) +
    theme_classic() + 
    ylab('Percent of essential /\n non-essential genes') + oe_x_axis + 
#    scale_y_continuous(labels=percent_format(accuracy = 1)) +
    scale_y_continuous() +
    scale_fill_manual(values=lists_to_plot, guide=F) + # labels=lists_labels, name=NULL) +
    annotate('text', 4.5, top_legend, hjust=0.5, vjust=1, 
             label='Cell essential', color=gene_list_colors['Haploinsufficient']) +
    annotate('text', 4.5, top_legend*0.88, hjust=0.5, vjust=1, 
             label='Cell non-essential', color=gene_list_colors['Olfactory Genes'])
  # theme(legend.justification = c(0.5, 0.8), legend.position = c(0.5, 1))
  
  if (save_plot) {
    pdf('3d_cell_ko.pdf', height = 3, width = 5)
    print(p)
    dev.off()
  }
  return(p)
}

# cell KO essential only 
cell_ko_essential = function(save_plot=F, metric='oe_lof_upper') {
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
  
  # gene_data %>%
  #   left_join(get_ko_gene_lists() %>% filter(gene_list == 'CEGv2')) %>%
  #   mutate(cell_essential = !is.na(gene_list)) %>%
  #   glm(cell_essential ~ oe_lof_upper + cds_length, ., family='binomial') %>%
  #   summary
  # gene_data %>%
  #   left_join(get_ko_gene_lists() %>% filter(gene_list == 'NEGv1')) %>%
  #   mutate(cell_non_essential = !is.na(gene_list)) %>%
  #   glm(cell_non_essential ~ oe_lof_upper + cds_length, ., family='binomial') %>%
  #   summary

  
  lists_to_plot = c('CEGv2' = color_lof)
  lists_labels = c('CEGv2' = 'Cell essential')
  
  ko_data = get_proportion_in_gene_lists(gene_data, metric=metric) %>%
    filter(gene_list %in% names(lists_to_plot))
  
  # t_test_gene_list(gene_data, 'CEGv2') %>% print
  #  t_test_gene_list(gene_data, 'NEGv1') %>% print
  
  top_legend = max(ko_data$prop_in_bin)# + ko_data$sem)
  p = ko_data %>%
    ggplot + aes(x = .data[[bin]], y = prop_in_bin, fill = gene_list,
                 ymin = prop_in_bin - sem, ymax = prop_in_bin + sem) + 
    #geom_bar(stat = 'identity', position = 'dodge') + 
    geom_bar(stat='identity', width=0.6) +
    theme_classic(base_size=15) + 
    ylab('Percent of essential genes') + oe_x_axis + 
    scale_y_continuous() +
    scale_fill_manual(values=lists_to_plot, guide=F) + # labels=lists_labels, name=NULL) +
    annotate('text', 4.5, top_legend, hjust=0.5, vjust=1, 
             label='Cell essential', color=gene_list_colors['Haploinsufficient']) +

  
  if (save_plot) {
    pdf('3d_cell_ko.pdf', height = 3, width = 5)
    print(p)
    dev.off()
  }
  return(p)
}

# generic gene list plotter




gene_data =   read_delim(gzfile('../data/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz'), delim = '\t',
                      col_types = do.call(cols, list())) %>%
  mutate(s = mu_lof / p)

# join with psfs
load('../fit.AN.80.RData')
#load('../fit.AN.RData')
#load('../fitall.HGNC.RData')
#load('../fit.nobin.RData')
#load('../fit.Polyphen.RData')

gene_data = SFS %>% mutate(psfs = -beta) %>% select('gene', 'psfs') %>% inner_join(gene_data)
#plot(oe_lof_upper~psfs, data=gene_data)
gene_data = gene_data %>% mutate(psfs_bin=as.numeric(cut(psfs, quantile(psfs, 0:10/10))) - 1)

plot(gene_list_spectrum())
plot(gene_list_spectrum(metric='psfs'))

#plot(gene_list_spectrum(gene_lists_to_plot = c('CEGv2_subset_universe', 'chd')))





# cell ko
plot(cell_ko_comparison())
plot(cell_ko_comparison(metric='psfs'))



# plots on small genes
gene_data = gene_data %>% mutate(cds_length_bin=as.numeric(cut(cds_length, quantile(cds_length, 0:10/10))) - 1)

gene_data = gene_data %>% filter(cds_length_bin<1)

plot(gene_list_spectrum())
plot(gene_list_spectrum(metric='psfs'))

# cell ko
plot(cell_ko_comparison())
plot(cell_ko_comparison(metric='psfs'))

# plot just cell essentials for small genes
plot(cell_ko_essential())
plot(cell_ko_essential(metric='psfs'))
