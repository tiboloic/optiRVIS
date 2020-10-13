# LT 19/03/2020
# lab meeting

# plot SFS versus LOEUF
plot(stdloeuf ~ beta,data=optiRVIS)

plot(stdloeuf ~ beta,data=optiRVIS, subset=ad_group=='N')
points(stdloeuf ~ beta,data=optiRVIS, subset=ad_group=='Y',col=2)
outs = subset(optiRVIS, abs(optiRVIS$beta)>3)
text(outs$beta, outs$stdloeuf, labels = outs$gene)

# plot SFS versus gevir
optiRVIS$stdgev = qnorm(optiRVIS$gevir_percentile/100)
plot(stdgev ~ beta,data=optiRVIS, subset=ad_group=='N')
points(stdgev ~ beta,data=optiRVIS, subset=ad_group=='Y',col=2)


## Final plot powerMOEUF, RVIS3, MOEUF, LOEUF

counts = orderpercent(optiRVIS, "oe_lof_upper", c("patho", "l"))
plot(pcpatho ~ pcl, data=counts, type="l", xlab="Percent Bases Covered",
     ylab = "Percent Pathogenic Variants", main="", lwd=2, col=cols[1])
abline(b=1, a=0)
AUC.RVIS = AUC(counts$pcl, counts$pcpatho)

#powerMOEUF
counts = orderpercent(optiRVIS, "pm", c("patho", "l"))
lines(pcpatho ~ pcl, data=counts, lwd=2, col=cols[4])

#gevir
counts = orderpercent(optiRVIS, "gevir_percentile", c("patho", "l"))
#lines(pcpatho ~ pcl, data=counts, lwd=2, col=cols[7])
AUC.gevir = AUC(counts$pcl, counts$pcpatho)

#gevir + power SFS
optiRVIS$gevirSFS = optiRVIS$gevir_percentile + 2.626263 * optiRVIS$beta
counts = orderpercent(optiRVIS, "gevirSFS", c("patho", "l"))
#lines(pcpatho ~ pcl, data=counts, lwd=2, col=cols[11])
AUC.pm = AUC(counts$pcl, counts$pcpatho)

#legend('topleft', legend = c('LOEUF', 'powerMOEUF', 'gevir', 'powergevir'), fill=cols[c(1,4,11,7)], bty='n')
legend('topleft', legend = c('LOEUF', 'powerMOEUF'), fill=cols[c(1,4)], bty='n')
