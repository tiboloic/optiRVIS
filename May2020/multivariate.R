# LT 7/05

# 'multivariate' analysis
#

# function to standerdize scores
stand = function(data, score) {

  # order by increasing score
  data = data[order(data[,score]),]
  f = ecdf(data[,score])
  data[,paste("std_", score, sep="")] = qnorm(f(data[,score]))
  data
}

optiRVIS$pSFS = -optiRVIS$beta

# standardize the scores
scores = c('RVIS3', 'oe_mis_upper', 'oe_lof_upper', 'pSFS', 'gevir_percentile')

standSco = optiRVIS[, scores]
for (i in 1:length(scores)) {
  standSco = stand(standSco, scores[i])
}
std_scores = sapply(scores, function(field) paste0('std_', field))

plot(standSco[, std_scores])
                    