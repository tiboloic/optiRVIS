# quick fix

sapply(seq(-1, 1,length.out = 100), function(alpha) {
optiRVIS$pm = (optiRVIS$FUNC/(optiRVIS$FUNC+optiRVIS$SYN)) + alpha * optiRVIS$beta
counts = orderpercent(optiRVIS, "pm", c("denovo", "l"))
#lines(pcdenovo ~ pcl, data=counts, lwd=2, col=cols[8])
AUC(counts$pcl, counts$pcdenovo)
})

sapply(seq(-1, 1,length.out = 100), function(alpha) {
  optiRVIS$pm = optiRVIS$stdloeuf + alpha * optiRVIS$beta
  counts = orderpercent(optiRVIS, "pm", c("hgmdpro", "l"))
  #lines(pcdenovo ~ pcl, data=counts, lwd=2, col=cols[8])
  AUC(counts$pcl, counts$pchgmdpro)
})

# 19/03/2020
# gevir
sapply(seq(2.5, 3,length.out = 100), function(alpha) {
  optiRVIS$pm = optiRVIS$gevir_percentile + alpha * optiRVIS$beta
  counts = orderpercent(optiRVIS, "pm", c("patho", "l"))
  #lines(pcpatho ~ pcl, data=counts, lwd=2, col=cols[8])
  AUC(counts$pcl, counts$pcpatho)
})

plot(stdloeuf ~ beta,data=optiRVIS)

sapply(seq(-1, 1,length.out = 100), function(alpha) {
  optiRVIS$pm = optiRVIS$gevir_percentile + alpha * optiRVIS$beta
  counts = orderpercent(optiRVIS, "pm", c("denovo", "l"))
  #lines(pcpatho ~ pcl, data=counts, lwd=2, col=cols[8])
  AUC(counts$pcl, counts$pcdenovo)
})

# RVIS + powerSFS
aucs=sapply(seq(-.1, .1,length.out = 100), function(alpha) {
  optiRVIS$pm = optiRVIS$RVIS3 + alpha * optiRVIS$beta
  counts = orderpercent(optiRVIS, "pm", c("patho", "l"))
  #lines(pcpatho ~ pcl, data=counts, lwd=2, col=cols[8])
  AUC(counts$pcl, counts$pcpatho)
})
plot(aucs)

