# LT 6/05

# length bias of various intolerance metrics

load('feb2020')
# need to run clinvar analysis to get field pSFS

## Final plot powerMOEUF, RVIS3, MOEUF, LOEUF
### function to order datasets for percentile plots
length_median_perdecile = function(data, score) {
  # data is the dataset we wnat to order
  # score is the field to order on, increasing

  
  # order by increasing score
  data = data[order(data[,score]),]
  
  # add count
  data$count = 0:(nrow(data)-1)/nrow(data)

  data$decile = floor(data$count*10)

  targetname = paste('l', score, sep='')
  res = data %>%
    group_by_at(vars('decile')) %>% 
    summarize(ltarget=median(l)) %>% ungroup
 # res[,targetname] = res$ltarget
  
#  return(res[,c(1,3)])
  return(res)  
}

scores = c('RVIS3', 'oe_mis_upper', 'oe_lof_upper', 'pSFS', 'gevir_percentile')
for (i in 1:length(scores)) {
  field = scores[i]
  if (i==1) {
    res = length_median_perdecile(optiRVIS, field)
    names(res) = c('decile', field) }
  else
    res[,field] = length_median_perdecile(optiRVIS, field)[,2]
}


cols = rainbow(5, start=0.3)
matplot(res[,-1], type="l", col=cols, lwd=2, lty=1, xlab='metric\'s decile', ylab='median gene coding length (bp)', axes=F)
axis(2)
axis(side=1, at=1:10, labels=c('0-1','1-2','2-3','3-4','4-5','5-6','6-7','7-8','8-9','9-10'))
legend('topright', legend = c('RVIS', 'MOEUF', 'LOEUF', 'powerSFS', 'GeVIR'), fill=cols[1:6], bty='n')











