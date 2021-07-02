# LT 19/10/2020

# plot of proportion in OMIM vs metric decile
# (figure e9a in flagship paper)

load_omim_by_year_data = function() {
  omim_data = read_delim('data/forKonrad_cleaned_gene_discovery_years_2018-10-09.tsv', delim = '\t') %>%
    filter(yearDiscovered > 1990)
  
  # omim_data %>%
  #   count(yearDiscovered) %>%
  #   ggplot + aes(x = yearDiscovered, y = n) + geom_bar(stat='identity')
  
  omim_data %>%
    group_by(Ensembl.Gene.ID) %>%
    summarize(year = min(yearDiscovered), 
              discoverybyNGS = any(discoverybyNGS)
    ) %>%
    return
}

proportion_in_omim = function(save_plot=F) {
  #omim_data = load_omim_by_year_data()
  
  gene_omim_data = gene_data %>%
    mutate(in_omim = gene_id %in% omim_data$Ensembl.Gene.ID)
  
  gene_omim_data %>%
    summarize(ttest = list(t.test(oe_lof_upper ~ in_omim))) %>%
    tidy(ttest)
  
  gene_omim_data %>%
    summarize(ttest = list(glm(in_omim ~ oe_lof_upper + cds_length, family='binomial'))) %>%
    tidy(ttest)
  
  p = gene_omim_data %>%
    group_by(oe_lof_upper_bin) %>%
    summarize(num_in_omim = sum(in_omim, na.rm=T), n=n(), 
              prop_in_omim = num_in_omim / n, sem = sqrt(prop_in_omim * (1 - prop_in_omim) / n),
              prop_upper = prop_in_omim + 1.96 * sem, prop_lower = prop_in_omim - 1.96 * sem) %>%
    ggplot + aes(x = oe_lof_upper_bin, y = prop_in_omim, ymin = prop_lower, ymax = prop_upper) + 
    geom_pointrange() + scale_y_continuous(labels=percent_format(accuracy=1)) +
    theme_classic() + oe_x_axis + ylab('Proportion in OMIM')
  
  if (save_plot) {
    pdf('e9_proportion_in_omim.pdf', height=3, width=4)
    print(p)
    dev.off()
  }
  return(p)
}

metric='oe_lof_upper'
variable_label = function(variable) {
  switch(variable,
         'oe_lof_upper'='LOEUF',
         'psfs'='powerSFS',
         'cds_length'='CDS length',
         variable)
}

proportion_in_clinvar = function(metric='oe_lof_upper', save_plot=F) {
  
  oe_x_label = paste(variable_label(metric), 'decile')
  
  label_function = function(x) {
    y = 10 * x + 5
    ifelse(y %% 20 == 0, paste0(y, '%'), "")
  }
  oe_x_axis = list(xlab(oe_x_label),
                   scale_x_continuous(labels=label_function, breaks=seq(-0.5, 9.5, 1), limits=c(-0.5, 9.5)))
  
  clinvar_data = read_delim('gene_lists/lists/clinvar_path_likelypath.tsv', delim='\t', col_names=c('gene'))
  gene_clinvar_data = gene_data %>% mutate(in_clinvar = gene %in% clinvar_data$gene)
  
  p = gene_clinvar_data %>%
    group_by(.data[[paste0(metric,'_bin')]]) %>%
    summarize(num_in_clinvar = sum(in_clinvar, na.rm=T), n=n(), 
              prop_in_clinvar = num_in_clinvar / n, sem = sqrt(prop_in_clinvar * (1 - prop_in_clinvar) / n),
              prop_upper = prop_in_clinvar + 1.96 * sem, prop_lower = prop_in_clinvar - 1.96 * sem) %>%
    ggplot + aes(x = .data[[paste0(metric,'_bin')]], y = prop_in_clinvar, ymin = prop_lower, ymax = prop_upper) + 
    geom_pointrange() + scale_y_continuous() + coord_cartesian(ylim=c(0,.25)) +
    theme_classic() + oe_x_axis + ylab('Proportion in ClinVar')
  
  if (save_plot) {
    pdf('e9_proportion_in_omim.pdf', height=3, width=4)
    print(p)
    dev.off()
  }
  return(p)
}
