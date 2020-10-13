# LT 3/09/2019

# Analyze gnomad counts
load("GnomadCounts.RData")

# fisrt get all variation without frequency cutoff
nocut = do.call(rbind, lapply(res,function(l) l[[1]]))

# counts coding variants (x axis in RVIS)
nocut = cbind(nocut, x = apply(nocut, 1, function(l) l[2]+l[3]+l[4]))

# counts missense and Lof coding variants (y axis in RVIS)
nocut = cbind(nocut, y = apply(nocut, 1, function(l) l[3]+l[4]))

# get y for other cutoffs
for (icut in 2:4) {
  dat = do.call(rbind, lapply(res,function(l) l[[icut]]))
  # counts missense and Lof coding variants (y axis in RVIS)
  target = apply(dat, 1, function(l) l[3]+l[4])
  # add column of zeros to store values
  nocut = cbind(nocut, 0)
  colnames(nocut)[ncol(nocut)] = paste("y.", icut, sep="")
  nocut[match(names(target), rownames(nocut)), ncol(nocut)] = target
}
# now for each frequency cutoff
cutoffs=c(0, 1e-6, 1e-4, 1e-2)


# buid the dataframe
counts = data.frame(gene = rownames(nocut), nocut[,5:9])

# try a fit
# on the ~ 2000 genes we have
plot(y ~ x, data = counts)
m.0 = lm(y ~ x, data = counts)
require(MASS)
score.0 = stdres(m.0)

# test performance on ClinVar
# load clinvar counts on limbr exons
# load limbr exons
load("../gnomad/intolerance/limbr.ex.Rdata")
pathopergene = aggregate(limbr.ex[,c("patho", "l")], by=list(gene=limbr.ex$gene), FUN=sum)

library(dplyr)
counts = inner_join(counts, pathopergene)

# add the scores
require(MASS)
m.1 = lm(y ~ x, data = counts)
counts$s.1 = stdres(m.1)
m.2 = lm(y.2 ~ x, data = counts)
counts$s.2 = stdres(m.2)
m.3 = lm(y.3 ~ x, data = counts)
counts$s.3 = stdres(m.3)
m.4 = lm(y.4 ~ x, data = counts)
counts$s.4 = stdres(m.4)

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


counts = orderpercent(counts, "s.1", c("patho", "l"))

plot(pcpatho ~ pcl, data=counts, type="l", xlab="Percent Bases Covered",
     ylab = "Percent Pathogenic Variants", main="", lwd=1.5, lty="longdash")
abline(b=1, a=0)
AUC.1 = AUC(counts$pcpatho, counts$pcl)

counts = orderpercent(counts, "s.2", c("patho", "l"))
lines(pcpatho ~ pcl, data=counts, lwd=1.5, col=2)
AUC.2 = AUC(counts$pcpatho, counts$pcl)

counts = orderpercent(counts, "s.3", c("patho", "l"))
lines(pcpatho ~ pcl, data=counts, lwd=1.5, col=3)
AUC.3 = AUC(counts$pcpatho, counts$pcl)

counts = orderpercent(counts, "s.4", c("patho", "l"))
lines(pcpatho ~ pcl, data=counts, lwd=1.5, col=4)
AUC.4 = AUC(counts$pcpatho, counts$pcl)
