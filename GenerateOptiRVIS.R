# LT 6/12/2019

# buid optiRVIS dataframe

library(dplyr)

# read Clinvar pathogenic counts from LIMBR/subRVIS analysis
load("../gnomad/intolerance/subrvis.ex.Rdata")

# object subrvis.ex contains data of CLinvar patho, missense patho and denovo
clinvar = subrvis.ex %>% group_by(gene) %>% summarize(patho=sum(patho), missense=sum(missense), denovo=sum(denovo), l=sum(l))

# load HGMD from biomart query
hgmd = read.delim(file='HGMD.biomart.tsv.gz')
hgmd = hgmd %>% group_by(Gene.name) %>% summarize(count=n())

# load gnomad extracted variants count
gnomad = read.delim(file='gnomADvariants.tsv')
gnomad %>% summarize(all = sum(FUNC) + sum(SYN))
# 8 724 064 variants as expected

# add LOEUF
cons = read.delim('../gnomad/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz')
cons = cons[,c('gene', 'oe_mis', 'oe_mis_upper', 'oe_lof', 'oe_lof_upper')]

# restrict on genes for which we have clinvar data, and add HGMD
optiRVIS = inner_join(gnomad, clinvar) %>% left_join(hgmd, by=c('gene'='Gene.name')) %>% inner_join(cons)
optiRVIS$hgmd = ifelse(is.na(optiRVIS$count), 0, optiRVIS$count)

# perform RVIS fitting
require(MASS)
m.5 = lm(AF5 ~ I(FUNC+SYN), data = optiRVIS)
#m.5 = lm(AF5 ~ l, data = optiRVIS)
optiRVIS$RVIS5 = stdres(m.5)

m.4 = lm(AF4 ~ I(FUNC+SYN), data = optiRVIS)
#m.4 = lm(AF4 ~ l, data = optiRVIS)
optiRVIS$RVIS4 = stdres(m.4)

m.3 = lm(AF3 ~ I(FUNC+SYN), data = optiRVIS)
#m.3 = lm(AF3 ~ l, data = optiRVIS)
optiRVIS$RVIS3 = stdres(m.3)

m.2 = lm(AF2 ~ I(FUNC+SYN), data = optiRVIS)
#m.2 = lm(AF2 ~ l, data = optiRVIS)
optiRVIS$RVIS2 = stdres(m.2)

m.1 = lm(FUNC ~ I(FUNC+SYN), data = optiRVIS)
#m.1 = lm(FUNC ~ l, data = optiRVIS)
optiRVIS$RVIS1 = stdres(m.1)

m.0 = lm(I(FUNC-AF5-AF4-AF3-AF2) ~ I(FUNC+SYN), data = optiRVIS)
#m.0 = lm(I(FUNC-AF5-AF4-AF3-AF2) ~ 0+I(FUNC), data = optiRVIS)
#m.0 = lm(I(FUNC-AF5-AF4-AF3-AF2) ~ l, data = optiRVIS)
optiRVIS$RVIS0 = -stdres(m.0)



save(optiRVIS, file='optiRVIS.RData')
