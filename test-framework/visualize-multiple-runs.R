#### 

library(ggplot2)
library(gridExtra)

plainsvmNames <- list()
graphsvmNames <- list()
for(i in c(1:10)){
  plainsvmNames[[i]] <- paste('test-output-plainsvm-randExp-geneStochastic-run-20,20',i,'.RData',sep='')
  graphsvmNames[[i]] <- paste('test-output-graphsvm-randExp-geneStochastic-run-20,20',i,'.RData',sep='')
}
 
plainsvmAUC <- list()
graphsvmAUC <- list()

for(i in c(1:10)){
  load(plainsvmNames[[i]])
  plainsvmAUC[[i]] <- apply(avg.auc,1,mean)
  rm(avg.auc)
 
  load(graphsvmNames[[i]])
  graphsvmAUC[[i]] <- apply(avg.auc, 1, mean)
  rm(avg.auc)
}

# First 10 columns of the frame are plainSVM AUCs and rest are for graphSVM
auc.list <- c(plainsvmAUC, graphsvmAUC)
auc.frame <- sapply(auc.list, function(x){return(x)})

line.types <- c(rep(1,10), rep(2,10))
color <- c(rep(1,10), rep(2,10))
leg <- c('plainSVM', rep('',9), 'graphSVM', rep('',9))

jpeg('combined-auc-10runs-randExp-geneStochastic.jpeg')
matplot(c(1:nrow(auc.frame)), auc.frame, type='l', lty=line.types,  
        col=color, xlab='Test Runs', ylab='Average AUC')
#legend(x='topleft', legend=leg, lty=line.types, col=color)
dev.off()

smooth.auc.list <- apply(auc.frame, MARGIN=2, 
                     FUN=function(x){return(loess.smooth(c(1:length(x)), x))})
jpeg('combined-auc-10runs-randExp-geneStochastic-smooth.jpeg')
plot(smooth.auc.list[[1]]$x, smooth.auc.list[[1]]$y, type='l', lty=1, ylim=c(0.5,1.0),
     xlab='Test Runs', ylab='Avg AUC', main='Points are interpolated for smoothing. Curves might be altered from original' )
for(i in c(2:10)){
  lines(smooth.auc.list[[i]]$x, smooth.auc.list[[i]]$y, lty=1)
}

for(i in c(11:20)){
  lines(smooth.auc.list[[i]]$x, smooth.auc.list[[i]]$y, lty=2, col=2)
}
legend(x='topleft', legend=c('plainsvm', 'graphsvm'), lty=c(1,2), col=c(1,2))
dev.off()


###### DIFFERENCE OF AUC GRAPH #######
diffAUC <- list()
for(i in 1:10){
  diffAUC[[i]] <- plainsvmAUC[[i]] - graphsvmAUC[[i]]
}

diffAUCFrame <- sapply(diffAUC, function(x){return(x)})

jpeg('diff-auc-10runs-randExp-geneStochastic.jpeg')
matplot(c(1:nrow(diffAUCFrame)), diffAUCFrame, type='l', 
         xlab='Test Runs', ylab='AUC of PlainSVM - AUC of GraphSVM')
#legend(x='topleft', legend=leg, lty=line.types, col=color)
dev.off()

##### Bar Graph showing the number of repeats in which PlainSVM did better
##### and the number of repeats in which GraphSVM did better, for each test run
#plainsvmWin <- apply(diffAUCFrame, 2, function(x){length(which(x>0))})
#graphsvmWin <- apply(diffAUCFrame, 2, function(x){length(which(x<=0))})

# First row is number of repeats in which PlainSVM won and second row is the
# number of repeats in which GraphSVM won for each of the ten runs
whoWon <- 
  apply(diffAUCFrame, 2, function(x){c(length(which(x>0)), length(which(x<0)))})

jpeg('whoWon-randExp-geneStochastic.jpeg')
barplot(height=whoWon, ylab='Number of wins', xlab='Test Runs', ylim=c(0,13),
        beside=TRUE, legend.text=c('plainsvm', 'graphsvm'), names.arg=c(1:10))
#legend(x='topleft', legend=c('plainsvm', 'graphsvm'))
dev.off()



