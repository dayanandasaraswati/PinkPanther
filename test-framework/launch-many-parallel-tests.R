## Parallel tests ##

library(multicore)
library(foreach)
library(doMC)


source('../code/test-framework/test-classification.R')


runBothSVM <- function(i){
  test.classification(10,2,c('graphsvm','plainsvm'), 'randomize', noise.gene.count=10, 
    randomize=0, 
    outfile.suffix=paste('randExp-geneStochastic-run-20,20',i,sep=''))
}

#runGraphSVM <- function(i){
#  test.classification(10,2,'graphsvm', 'randomize', noise.gene.count=10, 
#    randomize=0, outfile.suffix=paste('randExp-run',i,sep=''))
#}

registerDoMC()
foreach(i=c(1:10)) %dopar% runBothSVM(i)



