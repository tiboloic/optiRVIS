# LT 18/10/2020

# create gene lists in konrad's format for
# VCCRI chd gene list and DECIPHER  DDD: developmental disorder

require(tidyverse)

# DDD
ddd = read_csv('../data/DDG2P_18_10_2020.csv.gz')

ddd %>% filter(`DDD category` == 'confirmed') %>% select(`gene symbol`) %>% write_tsv(file='gene_lists/lists/ddg_conf.tsv', col_names = FALSE)

ddd %>% filter(`DDD category` == 'confirmed' & grepl('missense',`mutation consequence`)) %>% select(`gene symbol`) %>% write_tsv(file='gene_lists/lists/ddg_conf_miss.tsv', col_names = FALSE)

ddd %>% filter(grepl('missense',`mutation consequence`)) %>% select(`gene symbol`) %>% write_tsv(file='gene_lists/lists/ddg_miss.tsv', col_names = FALSE)


# CHD
chd = read_csv('../data/chdgene_table.csv')
chd %>% select(Gene) %>% write_tsv(file='gene_lists/lists/chd.tsv', col_names = FALSE)
