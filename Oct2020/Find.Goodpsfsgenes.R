
gene_data =   read_delim(gzfile('../data/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz'), delim = '\t',
                         col_types = do.call(cols, list())) %>%
  mutate(s = mu_lof / p)

# join with psfs
#load('../fit.AN.80.RData')
#load('../fitall.HGNC.RData')
load('../fit.Polyphen.RData')

gene_data = SFS %>% mutate(psfs = -beta) %>% select('gene', 'psfs') %>% inner_join(gene_data)
#plot(oe_lof_upper~psfs, data=gene_data)
#plot(log1p(oe_lof)~psfs, data=gene_data)
#lot(oe_lof~psfs, data=gene_data,log='y')
#plot(oe_mis~psfs, data=gene_data,log='y')
#plot(oe_mis_pphen~psfs, data=gene_data,log='y')

gene_data = gene_data %>% mutate(psfs_bin=as.numeric(cut(psfs, quantile(psfs, 0:10/10))) - 1)

or_genes = gene_data %>%
  filter(grepl('^OR', gene)) %>%
  transmute(gene = gene, gene_list = 'Olfactory Genes')

gene_lists = load_all_gene_list_data() %>%
  bind_rows(or_genes)


gene_list_data = gene_data %>% left_join(gene_lists) %>% 
  mutate(gene_list = if_else(grepl('haploinsufficiency', gene_list), 'Haploinsufficient', gene_list))
  

# identify ddg genes
gene_list_data %>% filter(oe_lof_upper > 1 & psfs_bin < 1 & gene_list == 'ddg_conf')
# DNMT3A   FTL     MIB1     ZMPSTE24

(gene_list_data %>% filter(oe_lof_upper > 1 & psfs_bin == 0 & gene_list == 'ddg_conf'))$gene
# ACAT1, BCKDHB, CDC6, DNMT3A, FKBP14, FTL, LZTR1, MIB1, NAGS, QDPR, SLC25A20, TWIST1, ZMPSTE24

(gene_list_data %>% filter(oe_lof_upper > 1 & psfs_bin == 0 & gene_list == 'ddg_conf_miss'))$gene
# LZTR1

(gene_list_data %>% filter(oe_lof_upper > 1 & psfs_bin == 0 & gene_list == 'ddg_miss'))$gene
#ACTA1, CFL2, LIPT1, LZTR1, TUBB3

gene_list_data %>% filter(oe_lof_upper > 1 & psfs_bin == 0 & gene_list == 'chd')
# HAND1, SMAD6



## essenial genes
gene_list_data %>% filter(oe_lof_upper > 1 & psfs_bin == 0 & gene_list == 'mgi_essential')

## essenial genes
(gene_list_data %>% filter(oe_lof_upper > 1 & psfs_bin == 0 & gene_list == 'all_ad'))$gene
(gene_list_data %>% filter(oe_lof_upper > 1 & psfs_bin == 0 & gene_list == 'all_ar'))$gene

(gene_list_data %>% filter(oe_lof_upper > 1 & psfs_bin == 0 & gene_list == 'fda_approved_drug_targets'))$gene

gene_list_data %>% filter(oe_lof_upper > 1 & psfs_bin == 0 & gene_list == 'clingen_level3_genes_2018_09_1')

gene_list_data %>% filter(oe_lof_upper > 1 & psfs_bin == 0 & gene_list == 'Haploinsufficient')

## crispr cas9 essential genes
(gene_list_data %>% filter(oe_lof_upper > 1 & oe_mis_upper >1 & psfs_bin ==0 & gene_list=='CEGv2_subset_universe'))
