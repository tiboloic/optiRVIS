# LT 1/6/2020
# modify original script to keep only protein_coding genes

# LT 6/12/2019

# count Lof, missense with various AF cutoffs and all synonymous

import hail as hl

ht = hl.read_table('gs://gnomad-public/release/2.1.1/ht/exomes/gnomad.exomes.r2.1.1.sites.ht')

# first filter by most_severe_consequence, then explode ?

CSQ_LOF = [
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost"]
    
CSQ_MIS = ['missense_variant']
CSQ_SYN = ['synonymous_variant']
CSQ_ALL = CSQ_LOF + CSQ_MIS + CSQ_SYN

# functional variants
CSQ_FUNC = CSQ_LOF + CSQ_MIS

csqfunc = hl.literal(CSQ_FUNC)
selcsq = hl.literal(CSQ_ALL)

# filter PASS and AF>0 and most_severe_consequence is LoF or missense
ht = ht.filter((ht.freq[0].AF > 0) & (ht.filters.length()==0) & csqfunc.contains(ht.row.vep.most_severe_consequence))

# annotate variant with first gene matching most_severe_consequence
#ht = ht.annotate(gene = ht.vep.transcript_consequences.find(lambda x: x.consequence_terms.contains(ht.vep.most_severe_consequence)).gene_symbol)
# first protein coding gene matching worse consequence. Should I instead get all the genes and explode ?
ht = ht.annotate(gene = ht.vep.transcript_consequences.find(lambda x: x.consequence_terms.contains(ht.vep.most_severe_consequence) & (x.biotype == 'protein_coding')).gene_symbol)

ht = ht.filter(hl.is_defined(ht.gene))

# bin allele counts per octave
ht = ht.annotate(cat = hl.floor(hl.log(ht.freq[0].AC,2)))

# group by gene and summarize 
res = ht.group_by('gene','cat').aggregate(obs=hl.agg.count())

res.export('/mnt/d/optiRVIS/gnomADSFS.2.tsv')

# quality control: total number of variants shoud be equal to FUNC + SYN
ht.count()
# YES ! 8 724 064 variants
# now with the filtering: 6 027 410 variants