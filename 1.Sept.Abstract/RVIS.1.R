# LT 3/09/2019
# first hopefully functional version of counting variants by class in gnomad with bioconductor


library(biomaRt)

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
# ok let's work with CCDS genes with a HGNC symbol
genes = getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position'),
              filters = c('with_hgnc', 'with_ccds'), values =list(TRUE, TRUE), mart = ensembl)
genes = makeGRangesFromDataFrame(genes, keep.extra.columns = TRUE, ignore.strand = TRUE,
                                 start.field='start_position', end.field = 'end_position')
# do it instaed by chromosome ?
# problem of genes is same variant might be counted twice:
# we want to treat each variant one once: we need non overlapping regions (which genes at not)

# gnomad VCF file
tf = TabixFile("D:/gnomad/gnomad.exomes.r2.1.1/filtered.vcf.bgz", index = "D:/gnomad/gnomad.exomes.r2.1.1/filtered.vcf.bgz.csi")

# use vcf header to get contigs
gnomadheader = scanVcfHeader(tf)
contigs = meta(gnomadheader)[[3]]
# read by batches of 100 genes
# or read gene by gene
# or chromosome by chromosome ?

# vep consequences
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

# function to parse VEP annoatation
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
  
  if (nrow(m) >0) {
    # get worst consequence per gene
   worstcsq = aggregate(m[c("gene","csq")], by=list(m$gene), max)[,-1]
   worstcsq$AF = info$AF[1]
    worstcsq
  } else c()
}

# function reads one gene from vcf and counts
countonegene = function(onegene) {
  sp <- ScanVcfParam(which=onegene, info=c("AF", "vep"), fixed="ALT")
  variants <- readVcf(tf, "hg19", param=sp)
  if (nrow(variants)> 0) {
    # row bind all the dataframes with gene and count
    df = do.call("rbind", apply(info(variants), 1, getwcsq))
    # iterate on freq cutoffs
    lapply(cutoffs, function(freq) table(subset(df, AF>freq)[,c("gene", "csq")]))
  }
}

cutoffs=c(0, 1e-6, 1e-4, 1e-2)

# do it in 100 batches
batches = floor(seq(1, length(genes) + 1, length.out=100))

res = list()
#for (i in 1:(length(batches)-1)) {
for (i in 14:(length(batches)-1)) {
cat(paste("batch ", i, "\n"))
  res[[i]] = countonegene(genes[batches[i]:(batches[i+1]-1)])
}

save(res, file="GnomadCounts.RData")
