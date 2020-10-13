# LT 6/12/2019
#
# optiRVIS analysis
#
# original analyses prior to ABABCS conference

# LT 6/01/2020
# add powerSFS, powerLOEUF, powerMOEUF

load('optiRVIS.RData')

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
#cols = hcl.colors(5)

counts = orderpercent(optiRVIS, "RVIS1", c("patho", "l"))
plot(pcpatho ~ pcl, data=counts, type="l", xlab="Percent Bases Covered",
     ylab = "Percent Pathogenic Variants", main="", lwd=2, col=cols[1])
AUC.1 = AUC(counts$pcl, counts$pcpatho)
abline(b=1, a=0)

counts = orderpercent(optiRVIS, "RVIS5", c("patho", "l"))
lines(pcpatho ~ pcl, data=counts, lwd=2, col=cols[2])

AUC.5 = AUC(counts$pcl, counts$pcpatho)

counts = orderpercent(optiRVIS, "RVIS4", c("patho", "l"))
lines(pcpatho ~ pcl, data=counts, lwd=2, col=cols[3])
AUC.4 = AUC(counts$pcl, counts$pcpatho)

counts = orderpercent(optiRVIS, "RVIS3", c("patho", "l"))
lines(pcpatho ~ pcl, data=counts, lwd=2, col=cols[4])
AUC.3 = AUC(counts$pcl, counts$pcpatho)

counts = orderpercent(optiRVIS, "RVIS2", c("patho", "l"))
lines(pcpatho ~ pcl, data=counts, lwd=2, col=cols[5])
AUC.2 = AUC(counts$pcl, counts$pcpatho)

# LOEUF
counts = orderpercent(optiRVIS, "oe_lof_upper", c("patho", "l"))
lines(pcpatho ~ pcl, data=counts, lwd=2, col=1)

# MOEUF
counts = orderpercent(optiRVIS, "oe_mis", c("patho", "l"))
lines(pcpatho ~ pcl, data=counts, lwd=2, col=2)

# MOEUF
optiRVIS$beta = - optiRVIS$beta
counts = orderpercent(optiRVIS, "beta", c("patho", "l"))
lines(pcpatho ~ pcl, data=counts, lwd=4, col=3)

# dummy
#counts = orderpercent(optiRVIS, "patho", c("patho", "l"))
#lines(pcpatho ~ pcl, data=counts, lwd=2, col=3)

# dummy 2
optiRVIS$beta = optiRVIS$RVIS3 - 2 * optiRVIS$patho/optiRVIS$l
#counts = orderpercent(optiRVIS, "beta", c("patho", "l"))
#lines(pcpatho ~ pcl, data=counts, lwd=2, col=2)


legend('topleft', legend = c('no filtering', '0.001 %', '0.01 %', '0.1 %', '1 %'), fill=cols[1:5], bty='n')
#counts = orderpercent(optiRVIS, "RVIS0", c("patho", "l"))
#lines(pcpatho ~ pcl, data=counts, lwd=3, col=1)
#AUC.0 = AUC(counts$pcl, counts$pcpatho)

#counts = orderpercent(optiRVIS, "RVIS5", c("hgmd", "l"))
#lines(pchgmd ~ pcl, data=counts, lwd=1.5, col=6)

which.max(c(AUC.0, AUC.1, AUC.2, AUC.3, AUC.4, AUC.5))


## Final plot powerMOEUF, RVIS3, MOEUF, LOEUF

counts = orderpercent(optiRVIS, "RVIS3", c("patho", "l"))
plot(pcpatho ~ pcl, data=counts, type="l", xlab="Percent Bases Covered",
     ylab = "Percent Pathogenic Variants", main="", lwd=2, col=cols[1])
abline(b=1, a=0)

# MOEUF
counts = orderpercent(optiRVIS, "oe_mis_upper", c("patho", "l"))
lines(pcpatho ~ pcl, data=counts, lwd=2, col=cols[2])

# LOEUF
counts = orderpercent(optiRVIS, "oe_lof_upper", c("patho", "l"))
lines(pcpatho ~ pcl, data=counts, lwd=2, col=cols[3])

#powerMOEUF
counts = orderpercent(optiRVIS, "pm", c("patho", "l"))
lines(pcpatho ~ pcl, data=counts, lwd=2, col=cols[4])
#powerMOEUF
counts = orderpercent(optiRVIS, "beta", c("patho", "l"))
lines(pcpatho ~ pcl, data=counts, lwd=2, col=cols[5])
legend('topleft', legend = c('RVIS 0.1 %', 'MOEUF', 'LOEUF', 'powerMOEUF'), fill=cols[1:4], bty='n')

##
# de novo
##
counts = orderpercent(optiRVIS, "RVIS5", c("denovo", "l"))

plot(pcdenovo ~ pcl, data=counts, type="l", xlab="Percent Bases Covered",
     ylab = "Percent Pathogenic Variants", main="", lwd=1.5, lty="longdash")
abline(b=1, a=0)
AUC.5 = AUC(counts$pcl, counts$pcdenovo)

counts = orderpercent(optiRVIS, "RVIS4", c("denovo", "l"))
lines(pcdenovo ~ pcl, data=counts, lwd=1.5, col=2)
AUC.4 = AUC(counts$pcl, counts$pcdenovo)

counts = orderpercent(optiRVIS, "RVIS3", c("denovo", "l"))
lines(pcdenovo ~ pcl, data=counts, lwd=1.5, col=3)
AUC.3 = AUC(counts$pcl, counts$pcdenovo)

counts = orderpercent(optiRVIS, "RVIS2", c("denovo", "l"))
lines(pcdenovo ~ pcl, data=counts, lwd=1.5, col=4)
AUC.2 = AUC(counts$pcl, counts$pcdenovo)

counts = orderpercent(optiRVIS, "RVIS1", c("denovo", "l"))
lines(pcdenovo ~ pcl, data=counts, lwd=3, col=5)
AUC.1 = AUC(counts$pcl, counts$pcdenovo)

counts = orderpercent(optiRVIS, "RVIS0", c("denovo", "l"))
lines(pcdenovo ~ pcl, data=counts, lwd=3, col=6)
AUC.0 = AUC(counts$pcl, counts$pcdenovo)

# ideal reference: if the metrics would catch perfectly the ratio denovo/length
optiRVIS$topdenovo = -optiRVIS$denovo/optiRVIS$l
counts = orderpercent(optiRVIS, "topdenovo", c("denovo", "l"))
lines(pcdenovo ~ pcl, data=counts, lwd=2, col=3)

counts = orderpercent(optiRVIS, "pm", c("denovo", "l"))
lines(pcdenovo ~ pcl, data=counts, lwd=2, col=2)

# HGMD

counts = orderpercent(optiRVIS, "RVIS1", c("hgmd", "l"))
plot(pchgmd ~ pcl, data=counts, type="l", xlab="Percent Bases Covered",
     ylab = "Percent Pathogenic Variants", main="", lwd=2, col=cols[1])
AUC.1 = AUC(counts$pcl, counts$pcpatho)
abline(b=1, a=0)

counts = orderpercent(optiRVIS, "RVIS5", c("hgmd", "l"))
lines(pchgmd ~ pcl, data=counts, lwd=2, col=cols[2])

AUC.5 = AUC(counts$pcl, counts$pcpatho)

counts = orderpercent(optiRVIS, "RVIS4", c("hgmd", "l"))
lines(pchgmd ~ pcl, data=counts, lwd=2, col=cols[3])
AUC.4 = AUC(counts$pcl, counts$pcpatho)

counts = orderpercent(optiRVIS, "RVIS3", c("hgmd", "l"))
lines(pchgmd ~ pcl, data=counts, lwd=2, col=cols[4])
AUC.3 = AUC(counts$pcl, counts$pcpatho)

counts = orderpercent(optiRVIS, "RVIS2", c("hgmd", "l"))
lines(pchgmd ~ pcl, data=counts, lwd=2, col=cols[5])
AUC.2 = AUC(counts$pcl, counts$pcpatho)

# LOEUF
counts = orderpercent(optiRVIS, "oe_lof_upper", c("hgmd", "l"))
lines(pchgmd ~ pcl, data=counts, lwd=2, col=1)

legend('topleft', legend = c('no filtering', '0.001 %', '0.01 %', '0.1 %', '1 %'), fill=cols[1:5], bty='n')






# stupid to do benign on whole genes
# % benign
optiRVIS$benign = as.integer(optiRVIS$patho==0)
counts = orderpercent(optiRVIS, "RVIS5", c("patho", "benign"))

plot(pcpatho ~ pcbenign, data=counts, type="l", xlab="Percent benign",
     ylab = "Percent Pathogenic Variants", main="", lwd=1.5, lty="longdash")
abline(b=1, a=0)
AUC.5 = AUC(counts$pcdenovo, counts$pcl)

## SFS plot
SFS = function(beta, n=100) (1:10)^-beta / sum((1:10)^-beta) * n
barplot(SFS(1.5))
barplot(SFS(2, 50), col=2, add=T)
barplot(SFS(2, 80), col=2, add=T)

sfss = rbind(SFS(1.5), SFS(2.1, 80), SFS(2,50))
par(cex=1.5)
barplot(sfss, beside=TRUE, col=cols[1:3], space=c(0,.3), names.arg=c(1:10))
legend('topright', legend = c('neutral', 'shifted & depleted (-20%)', 'shifted & depleted (-50%)'), fill=cols[1:3], bty='n')
