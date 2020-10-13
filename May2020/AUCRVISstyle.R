# LT 6/05
# trying to finalise powerSFS results

# per base efficiency AUC RVIS style
load('feb2020')


## Final plot powerMOEUF, RVIS3, MOEUF, LOEUF
### function to order datasets for percentile plots
orderpercent = function(data, score, percent) {
  # data is the dataset we wnat to order
  # score is the field to order on, increasing
  # percent are the fields appearing on the axes of the plot
  # that need to be converted to percentiles
  
  # order by increasing score
  data = data[order(data[,score]),]
  
  # function to percentiize fields
  percentilize = function(data) {
    # data should be a numeric vector
    cumsum(data)/sum(data)
  }
  
  for (i in 1:length(percent)) {
    field = percent[i]
    data[,paste("pc", field, sep="")] = percentilize(data[, field])
  }
  data
}

AUC = function(x,y) {
  x = c(0, x, 1)
  y = c(0, y, 1)
  sum((y[-1] - diff(y)/2) * diff(x))
}


cols = rainbow(5, start=0.3)

# to store AUC results
AUCs = c()

counts = orderpercent(optiRVIS, "RVIS3", c("hgmd", "l"))
plot(pchgmd ~ pcl, data=counts, type="l", xlab="Percent Bases Covered",
     ylab = "Percent Pathogenic Variants", main="", lwd=2, col=cols[1])
abline(b=1, a=0)
AUCs["RVIS"] = AUC(counts$pcl, counts$pchgmd)

# MOEUF
counts = orderpercent(optiRVIS, "oe_mis_upper", c("hgmd", "l"))
lines(pchgmd ~ pcl, data=counts, lwd=2, col=cols[2])
AUCs["MOEUF"] = AUC(counts$pcl, counts$pchgmd)

# LOEUF
counts = orderpercent(optiRVIS, "oe_lof_upper", c("hgmd", "l"))
lines(pchgmd ~ pcl, data=counts, lwd=2, col=cols[3])
AUCs["LOEUF"] = AUC(counts$pcl, counts$pchgmd)

#powerSFS
# the scale is inverted
#optiRVIS$pSFS = -optiRVIS$beta
counts = orderpercent(optiRVIS, "pSFS", c("hgmd", "l"))
lines(pchgmd ~ pcl, data=counts, lwd=2, col=cols[4])
AUCs["pSFS"] = AUC(counts$pcl, counts$pchgmd)

#gevir
counts = orderpercent(optiRVIS, "gevir_percentile", c("hgmd", "l"))
lines(pchgmd ~ pcl, data=counts, lwd=2, col=cols[5])
AUCs["GeVIR"] = AUC(counts$pcl, counts$pchgmd)


legend('topleft', legend = c('RVIS', 'MOEUF', 'LOEUF', 'powerSFS', 'GeVIR'), fill=cols[1:6], bty='n')



                               