# LT 6/12/2019

# count Lof, missense with various AF cutoffs and all synonymous


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

# filter PASS and AF>0 and most_severe_consequence is LoF, missense or synonymous
ht = ht.filter((ht.freq[0].AF > 0) & (ht.filters.length()==0) & selcsq.contains(ht.row.vep.most_severe_consequence))

# annotate variant with first gene matching most_severe_consequence
ht = ht.annotate(gene = ht.vep.transcript_consequences.find(lambda x: x.consequence_terms.contains(ht.vep.most_severe_consequence)).gene_symbol)

# build category: 0 = synonymous, 1: functional variant AF<1e5, etc.
ht = ht.annotate(cat = hl.cond(ht.vep.most_severe_consequence == 'synonymous_variant', 0, hl.case()
.when(ht.freq[0].AF > 1e-2, 2)
.when(ht.freq[0].AF > 1e-3, 3)
.when(ht.freq[0].AF > 1e-4, 4)
.when(ht.freq[0].AF > 1e-5, 5)
.default(1)))

# group by gene and summarize 
res = ht.group_by('gene').aggregate(AF5=hl.agg.count_where(ht.cat == 5),
AF4=hl.agg.count_where(ht.cat == 4),
AF3=hl.agg.count_where(ht.cat == 3),
AF2=hl.agg.count_where(ht.cat == 2),
FUNC=hl.agg.count_where(ht.cat!=0),
SYN=hl.agg.count_where(ht.cat == 0))

res.export('/mnt/d/optiRVIS/gnomADvariants.tsv')

# quality control: total number of variants shoud be equal to FUNC + SYN
ht.count()
# YES ! 8 724 064 variants
