# LT 13/10/2020
# extract all variants to plot SFS by category of variant (Lof, missense, synonymous, intronic, etc.)


import hail as hl

ht = hl.read_table('gs://gnomad-public-requester-pays/release/2.1.1/ht/exomes/gnomad.exomes.r2.1.1.sites.ht')

# for AN cutoff
AN_ADJ_FILTER = 0.8
def get_an_filter(ht):
	samples = {'female': 57787, 'male':67961}
	return (hl.case()
	.when(ht.locus.in_autosome_or_par(), ht.freq[0].AN >= AN_ADJ_FILTER * 2 * sum(samples.values()))
	.when(ht.locus.in_x_nonpar(), ht.freq[0].AN >= AN_ADJ_FILTER * (samples['male'] + samples['female'] * 2))
	.when(ht.locus.in_y_nonpar(), ht.freq[0].AN >= AN_ADJ_FILTER * samples['male'])
	.or_missing())


# filter PASS and AF>0
ht = ht.filter((ht.freq[0].AF > 0) & (ht.filters.length()==0) & (ht.variant_type=='snv'))

# bin allele counts per octave and add worst_csq
ht = ht.annotate(cat = hl.floor(hl.log(ht.freq[0].AC,2)), protein_coding = ht.vep.transcript_consequences.any(lambda x: x.biotype == 'protein_coding'), worst_csq = ht.vep.most_severe_consequence)

# group by worst_csq, cat and summarize 
res = ht.group_by('worst_csq','protein_coding', 'cat').aggregate(obs=hl.agg.count())

#res.export('/mnt/d/optiRVIS/gnomADSFS.5.tsv')
res.export('gs://vccrieg/optiRVIS/gnomADSFS.worst_csq.all.tsv')

# quality control: all exomes variants
ht.count()
# 
