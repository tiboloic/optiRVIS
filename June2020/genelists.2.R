# LT 11/06
# update with new score (HGNC)

# LT 7/05

# genelist analysis

# to avoid serach replace
optiRVIS=dat

# add RVIS
rvis = read.delim('GenicIntolerance_v3_12Mar16.txt')
# keep only recommended exac default: filter= 0.05%
rvis = rvis[,c(1,21)]
names(rvis) = c('gene','RVIS')
optiRVIS$RVIS = rvis[match(optiRVIS$gene, rvis$gene), 'RVIS']
## modified for gene groups

orderpercent = function(data, score, percent) {
  # data is the dataset we wnat to order
  # score is the field to order on, increasing
  # percent are the fields appearing on the axes of the plot
  # that need to be converted to percentiles
  
  # order by increasing score
  data = data[order(data[,score]),]
  data$counter = 1:nrow(data)
  
  # function to percentiize fields
  percentilize = function(data) {
    # data should be a numeric vector
    cumsum(data)/sum(data)
  }
  percent = c(percent, 'counter')
  
  for (i in 1:length(percent)) {
    field = percent[i]
    if (is.factor(data[, field]))
      data[,paste("pc", field, sep="")] = percentilize(data[, field]=='Y')
    else
      data[,paste("pc", field, sep="")] = percentilize(data[, field])
  }
  data
}

AUC = function(x,y) {
  x = c(0, x, 1)
  y = c(0, y, 1)
  sum((y[-1] - diff(y)/2) * diff(x))
}

geneGroups = function(genegroup) {
  
  cols = rainbow(6, start=0.3)
  
  # to store AUC results
  AUCs = c()
  
  pcy=paste0('pc', genegroup)
  
  counts = orderpercent(optiRVIS, "virlof_percentile", genegroup)
  plot(counts$pccounter, counts[,pcy], type="l", xlab="Percentile",
       ylab = paste0("Proportion of genes in ", genegroup), main="", lwd=2, col=cols[1])
  abline(b=1, a=0)
  AUCs["virlof"] = AUC(counts$pccounter, counts[,pcy])

  #gevir
  counts = orderpercent(optiRVIS, "gevir_percentile", genegroup)
  lines(counts$pccounter, counts[,pcy], lwd=2, col=cols[2])
  AUCs["GeVIR"] = AUC(counts$pccounter, counts[,pcy])
  
  # LOEUF
  counts = orderpercent(optiRVIS, "oe_lof_upper", genegroup)
  lines(counts$pccounter, counts[,pcy], lwd=2, col=cols[3])
  AUCs["LOEUF"] = AUC(counts$pccounter, counts[,pcy])

  # MOEUF
  counts = orderpercent(optiRVIS, "oe_mis_upper", genegroup)
  lines(counts$pccounter, counts[,pcy], lwd=2, col=cols[4])
  AUCs["MOEUF"] = AUC(counts$pccounter, counts[,pcy])  
  
  #powerSFS
  counts = orderpercent(optiRVIS, "pSFS_percentile", genegroup)
  lines(counts$pccounter, counts[,pcy], lwd=2, col=cols[5])
  AUCs["pSFS"] = AUC(counts$pccounter, counts[,pcy])
  
  #RVIS
  counts = orderpercent(optiRVIS, "RVIS", genegroup)
  lines(counts$pccounter, counts[,pcy], lwd=2, col=cols[6])
  AUCs["RVIS"] = AUC(counts$pccounter, counts[,pcy])  
  
  
  legend('bottomright', legend = c(paste0('virLOF AUC=', round(AUCs["virlof"],2)),
                                   paste0('GeVIR AUC=', round(AUCs["GeVIR"],2)),
                                   paste0('LOEUF AUC=', round(AUCs["LOEUF"],2)),
         paste0('MOEUF AUC=', round(AUCs["MOEUF"],2)),
         paste0('powerSFS AUC=', round(AUCs["pSFS"],2)),
         paste0('RVIS AUC=', round(AUCs["RVIS"],2))),
         fill=cols[1:6], bty='n')
  
  AUCs
}

geneGroups("ad_group")
geneGroups("mouse_het_lethal_group")

chdgenes = read.csv('chdgene_table.csv')

optiRVIS$chd_group = optiRVIS$gene %in% chdgenes$Gene

geneGroups("chd_group")

# add DDG@P disorders
ddg2p = read.csv('DDG2P_11_6_2020.csv')
optiRVIS$ddg = optiRVIS$gene %in% ddg2p$gene.symbol
geneGroups("ddg")

