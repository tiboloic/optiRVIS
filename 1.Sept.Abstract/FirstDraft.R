# some bioconductor test to read variants

library(VariantAnnotation)
tf = TabixFile("D:/gnomad/gnomad.exomes.r2.1.1/gnomad.exomes.r2.1.1.sites.vcf.bgz")
gr <- GRanges(seqnames="1", ranges=IRanges(start=100000,
                                           end=1500000))
sp <- ScanVcfParam(which=gr, info=c("AF"), fixed="ALT")

res <- readVcf(tf, "hg19", param=sp)

library(GenomicFeatures)
library(biomaRt)

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
genes = getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position', 'ccds'),
              filters = c('chromosome_name', 'with_ccds'), values =list("1", TRUE), mart = ensembl)

genes = getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position', 'ccds'),
              filters = 'with_ccds', values =TRUE, mart = ensembl)

genes = getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position'),
              filters = 'with_hgnc', values =TRUE, mart = ensembl)

# ok let's work with CCDS genes with a HGNC symbol
genes = getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position'),
              filters = c('with_hgnc', 'with_ccds'), values =list(TRUE, TRUE), mart = ensembl)

# another way was to install package EnsDb.Hsapiens.v79
# try it
edb = EnsDb.Hsapiens.v79
genes(edb)
# is GRCh38 though
# EnsDb.Hsapiens.v75 should be HG37

# get variants per gene, with allele frequency and variant type

# create a Granges object from my gene list
genes = makeGRangesFromDataFrame(genes, keep.extra.columns = TRUE, ignore.strand = TRUE,
                                 start.field='start_position', end.field = 'end_position')

#tf = TabixFile("D:/gnomad/gnomad.exomes.r2.1.1/gnomad.exomes.r2.1.1.sites.vcf.bgz")
tf = TabixFile("D:/gnomad/gnomad.exomes.r2.1.1/filtered.vcf.bgz", index = "D:/gnomad/gnomad.exomes.r2.1.1/filtered.vcf.bgz.csi")
sp <- ScanVcfParam(which=genes[1:100], info=c("AF", "vep"), fixed="ALT")

res <- readVcf(tf, "hg19", param=sp)


sapply(info(res)$vep, function(l) strsplit(l, "|")[2])

# use sub to extract consequence
sub("([^|]+)\\|([^|]+)\\|(.+)", "\\2", info(res)$vep[[1]])
strsplit(info(res)$vep[[1]], "|", fixed=TRUE)

# try to get a list of consequences
unique(unlist(sapply(info(res)$vep, function(l) {sub("([^|]+)\\|([^|]+)\\|(.+)", "\\2", l)})))
unlist(strsplit(sub("([^|]+)\\|([^|]+)\\|(.+)", "\\2", info(res)$vep[[1000]]), "&", fixed=TRUE))

# function takes a vep annotation retuns a list of consequences
getcsqs = function (vep) {
  unlist(strsplit(sub("([^|]+)\\|([^|]+)\\|(.+)", "\\2", vep), "&", fixed=TRUE))
}


csqs = c("splice_acceptor_variant" = 3, "splice_donor_variant"=3, "stop_gained"=3, "frameshift_variant"=3,
           "stop_lost"=3, "start_lost"=3,
           "inframe_insertion"=2, "inframe_deletion"=2, "missense_variant"=2,
           "protein_altering_variant"=2, "splice_region_variant"=1, "start_retained_variant"=1,
           "stop_retained_variant" = 1, "synonymous_variant" = 1,
         "transcript_ablation" = 0, "transcript_amplification" = 0, "incomplete_terminal_codon_variant" = 0,
         "coding_sequence_variant" = 0, "mature_miRNA_variant" =0, "5_prime_UTR_variant"=0,
         "3_prime_UTR_variant" = 0, "non_coding_transcript_exon_variant" = 0, "intron_variant" = 0,
         "NMD_transcript_variant" = 0, "non_coding_transcript_variant" =0, "upstream_gene_variant" =0,
         "downstream_gene_variant" = 0, "TFBS_ablation"=0, "TFBS_amplification"=0, "TF_binding_site_variant"=0,
         "regulatory_region_ablation"=0, "regulatory_region_amplification"=0, "feature_elongation"=0,
         "regulatory_region_variant"=0, "feature_truncation"=0, "intergenic_variant"=0)
variant_class = sapply(info(res)$vep, function(onevep) max(csqs[getcsqs(onevep)], na.rm=TRUE))

# try to manage vep info properly to get the gene
strsplit(unlist(strsplit(info(res)$vep[[1]], ",", fixed=TRUE)), "|", fixed=TRUE)
matrix(unlist(strsplit(unlist(strsplit(info(res)$vep[[1]], ",", fixed=TRUE)), "|", fixed=TRUE)),nrow=67)


# master function
# extracts vep infomation, parse it, group by gene, retain worse consequence per gene
# returns dataframe, one (or zero) row per gene with columns variant_id, gene symbol, consequence
getwcsq = function(info) {
  vep = info$vep
  m = strsplit(unlist(strsplit(vep, ",", fixed=TRUE)), "|", fixed=TRUE)
  #check nb of elements should be 67 or 68 (if annotated by loftee)
  l = sapply(m, length)
  if (any(l< 67) || any(l> 68)) stop("VEP parse error")
  
  m = as.data.frame(t(sapply(m, function(vep) {
    allcsqs=unlist(strsplit(vep[2], "&", fixed=TRUE))
    c(gene=vep[4], csq = max(csqs[allcsqs]), biotype=vep[8], symbol_source = vep[25])})), stringsAsFactors=FALSE)

  # remove rows without gene, keep only protein coding HGNC genes
  m = subset(m, m$gene != "" & m$biotype == "protein_coding" & m$symbol_source == "HGNC")

  # get worst consequence per gene
  worstcsq = aggregate(m[c("gene","csq")], by=list(m$gene), max)[,-1]
  worstcsq$AF = info$AF[1]
  worstcsq
}

# row bind all the dataframes with gene and count
df = do.call("rbind", apply(info(res), 1, getwcsq))



# now count per gene
counts = table(df)
