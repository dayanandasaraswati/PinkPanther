#library('igraph');
library(pathClass);
library(e1071)
library(ROCR);
library(multicore);
library(Matrix);

### COMMENTED CODE BLOCK ###
"
#set.seed(123)

# Reading the expr data and the network from files
expr = read.table('anath_random_exp.fold_10.txt', check.names=FALSE);
TEST_expr = read.table('toyTestExpr.data', check.names=FALSE);
#network = read.graph('fcNetwork.txt');
adj_matrix = read.table('fc_network_filtered.txt', check.names=FALSE);
adj_matrix = round(adj_matrix)
adj_matrix = as.matrix(adj_matrix)

# Extract portions of the data into different variables
# for convenience
newX = subset(expr, select = -disease_outcome); # without disease_outcome column
newX = as.matrix(newX);
newY = factor(expr$disease_outcome); # just the disease_outcome column
#adj_matrix = get.adjacency(network); 
#rownames(adj_matrix) = colnames(newX); # Making row,column names equal for matching
#colnames(adj_matrix) = colnames(newX);

TEST_x = subset(TEST_expr, select = -disease_outcome);
TEST_x = as.matrix(TEST_x);
TEST_y = factor(TEST_expr$disease_outcome);

# Creating the data structures required by pathClass
newMap = data.frame(probesetID=colnames(newX), graphID=colnames(adj_matrix));
newMap = as.matrix(newMap);
#newMatched = list();
#newMatched$x = newX;
#newMatched$y = newY;
#newMatched$mapping = newMap;
"
### END OF COMMENTED REGION ###


#### Randomizing Input (Reordering the samples) (Type-1)####
#### Both plain svm and network based svm shouldn't have 
#### any change in performance
# randomOrder = sample(nrow(newX))
# randomX = newX[randomOrder, ]
# randomY = newY[randomOrder]
# newX = randomX
# newY = randomY
#######################

#### Randomizing Input (No relation between X and Y) (Type-2)####
### Both plain svm and network svm must show very poor performance
# randomY = sample(newY) 
# newY = randomY
#######################

#### Randomizing Network (Degree preserving randomization) (Type-3)####
### Plain svm shouldn't have any change in performance but network
### based svm should have.
# randomNetwork = rewire(network, niter=100);
# adj_matrix = get.adjacency(randomNetwork); 
# rownames(adj_matrix) = colnames(newX); # Making row,column names equal for matching
# colnames(adj_matrix) = colnames(newX);
#######################

#### Label Permutation (Type-4) ####
#colnames(newX) = sample(colnames(newX))

######## Network-based SVM ########
nbSVM <- function(newX, newY, adj_matrix, fold, times){
ad.list <- as.adjacencyList(adj_matrix)
newCVResult.nbsvm <- crossval(newX, newY , theta.fit=fit.networkBasedSVM, 
                     folds=fold, repeats=times, parallel=TRUE, DEBUG=TRUE,
                     adjacencyList=ad.list, lambdas=10^(-3:3),
                     sd.cutoff=0.0);
#plot(newCVResult.nbsvm, fname="Result-nbsvm.pdf");                      
# newFit.nbsvm = fit.networkBasedSVM(newX, newY, DEBUG=FALSE, sd.cutoff=0.0, 
#                                 adjacencyList=ad.list, n.inner=5, lambdas=10^(-1:2));
# newPred.nbsvm = predict.networkBasedSVM(newFit.nbsvm, TEST_x, probability=TRUE, decision.values=TRUE)
}
########################

###### RRFE Algorithm ####
# newCVResult.rrfe = crossval(newX, newY, theta.fit=fit.rrfe, DEBUG=TRUE,
#                        folds=5, repeats=1, parallel=TRUE, Cs=10^(1:2));
#                        Gsub=adj_matrix, mapping=newMap, d=1/2, useAllFeatures=TRUE);
# plot(newCVResult.rrfe);                       
########################

###### Graph SVM Algorithm #####
graphSVM <- function(newX, newY, adj_matrix, betaval, fold, times){
adj_matrix = Matrix(adj_matrix)
 dk = calc.diffusionKernel(L=adj_matrix, is.adjacency=TRUE, beta=betaval);
 newCVResult.gsvm = crossval(newX, newY, theta.fit=fit.graph.svm, DEBUG=TRUE,
                            folds=fold, repeats=times, parallel=TRUE,
                            Cs=10^(-1:2), 
                            diffusionKernel=dk, useOrigMethod=TRUE);

#plot(newCVResult.gsvm, fname="result-graphsvm.pdf"); 
return(newCVResult.gsvm)
# newFit.gsvm = fit.graph.svm(newX, newY, mapping=newMap, diffusionKernel=dk);
# newPred.gsvm = predict.graphSVM(newFit.gsvm, TEST_x);
}
########################

###### Plain SVM ########
plainSVM <- function(newX, newY, kernel.type, fold, times){

## CROSS VALIDATION CODE COURTESY: Jean-Philippe Vert
# Randomly split the n samples into folds
# Returns a list of nfolds lists of indices, each corresponding to a fold
n <- length(newY)
s <- split(sample(n),rep(1:fold,length=n))

# Lists to store CV Results
ypred <- list()
ylabel <- list()

for (i in seq(fold)) {
  trainX <- newX[-s[[i]],]
  trainY <- newY[-s[[i]]]
  testX <- newX[s[[i]],]
  testY <- newY[s[[i]]]

  model <- svm(trainX, trainY, type=c('C-classification'), 
               kernel=kernel.type)

  pred <- predict(model, testX, decision.values=TRUE)
  ypred[[i]] <- attributes(pred)$decision.values
  ylabel[[i]] <- factor(testY, ordered=TRUE) # Ordered factors needed by ROCR
}

# Using ROCR to get the AUC
svm.pred <- prediction(ypred, ylabel)
svm.auc <- performance(svm.pred, measure='auc')

# Converting individual lists into y.values into single list
result <- list()
# Quite often, SVM and ROCR interpret the class labels differently. Negative
# decision values returned by SVM are compared to label 1 by ROCR which produces
# an AUC of less than 0.5. To correct for this switch of labels, the actual
# AUC is taken as 1-computedAUC.
result$auc <- sapply(svm.auc@y.values,function(x){return(ifelse(x<0.45, 1-x, x))})
cat('NOTE: Label Swap Might have been performed to get AUC values\n')
cat('AUCs from various CV Runs:\n')
cat(result$auc,'\n')
return(result)
}
########################

