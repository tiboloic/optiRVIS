# LT 7/05

# genelist analysis

load('feb2020')
# need to run clinvar analysis to get field pSFS

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
  
cols = rainbow(5, start=0.3)

# to store AUC results
AUCs = c()

pcy=paste0('pc', genegroup)

counts = orderpercent(optiRVIS, "RVIS3", genegroup)
plot(counts$pccounter, counts[,pcy], type="l", xlab="Percentile",
     ylab = paste0("Proportion of genes in ", genegroup), main="", lwd=2, col=cols[1])
abline(b=1, a=0)
AUCs["RVIS"] = AUC(counts$pccounter, counts[,pcy])

# MOEUF
counts = orderpercent(optiRVIS, "oe_mis_upper", genegroup)
lines(counts$pccounter, counts[,pcy], lwd=2, col=cols[2])
AUCs["MOEUF"] = AUC(counts$pccounter, counts[,pcy])

# LOEUF
counts = orderpercent(optiRVIS, "oe_lof_upper", genegroup)
lines(counts$pccounter, counts[,pcy], lwd=2, col=cols[3])
AUCs["LOEUF"] = AUC(counts$pccounter, counts[,pcy])

#powerSFS
# the scale is inverted
#optiRVIS$pSFS = -optiRVIS$beta
counts = orderpercent(optiRVIS, "pSFS", genegroup)
lines(counts$pccounter, counts[,pcy], lwd=2, col=cols[4])
AUCs["pSFS"] = AUC(counts$pccounter, counts[,pcy])

#gevir
counts = orderpercent(optiRVIS, "gevir_percentile", genegroup)
lines(counts$pccounter, counts[,pcy], lwd=2, col=cols[5])
AUCs["GeVIR"] = AUC(counts$pccounter, counts[,pcy])


legend('bottomright', legend = c(paste0('RVIS AUC=', round(AUCs["RVIS"],2)),
                                 paste0('MOEUF AUC=', round(AUCs["MOEUF"],2)),
                                 paste0('LOEUF AUC=', round(AUCs["LOEUF"],2)),
                                 paste0('powerSFS AUC=', round(AUCs["pSFS"],2)),
                                 paste0('GeVIR AUC=', round(AUCs["GeVIR"],2))),
                                 fill=cols[1:5], bty='n')

AUCs
}

geneGroups("ad_group")
geneGroups("mouse_het_lethal_group")

chdgenes = read.csv('chdgene_table.csv')

optiRVIS$chd_group = optiRVIS$gene %in% chdgenes$Gene

geneGroups("chd_group")
