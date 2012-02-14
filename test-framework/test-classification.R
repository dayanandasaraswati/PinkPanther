## Script that generates many different datasets and runs
## classification algorithms on them. 

test.classification <- function(num.tests, repeatition, 
                                classification.fn=c('graphsvm','plainsvm'),
                                exp.gen.method='synthesize',
                                mixture.type=c('snr', 'overlap', 'spread'),
                                noise.gene.count=0, randomize=0,
                                outfile.suffix=''){

  # This script generates many different synthetic datasets and runs
  # classification on them. The performance on each run is dumped into a text
  # file. It could be used for visualization. 
  #
  # Args:
  #  num.tests: Number of tests to perform
  #
  #  repeatition: Number of times to repeat each test
  #  
  #  classification.fn: Character vector. Name of the function to run. 
  #    Possible options are 'graphsvm', 'nbsvm', 'plainsvm'. The script
  #    classify.R script is loaded to locate these functions. When multiple
  #    functions are specified, each test run will generate a dataset and
  #    run both algorithms on same dataset.
  #
  #  exp.gen.method: Method to generate expr matrix. Use 'synthesize' for
  #    synthetic exp matrix or 'randomize' to randomize an actual expr 
  #    matrix. 
  #
  #  mixture.type: Varying the mixture of signal and noise components for 
  #    each of the generated datasets. Refer to generate-parameters.R file.
  #
  #  noise.gene.count: Number of noise genes to be added to dataset 
  #
  #  randomize: Integer number from 1 to 4 denoting type of randomization
  #    to perform on the expr matrix
  #    Type 1: Shuffling rows. Musn't affect classification performance
  #    Type 2: Shuffling class labels: No relation between samples and labels
  #    Type 3: Degree Preserving Randomization: Edges of the network are 
  #            shuffled in degree preserving way
  #    Type 4: Label Permutation: Gene IDs of expr matrix are permuted. 
  #            A network based classifier must get confused.
  # 
  #  outfile.suffix: Character. Suffix to be attached to the output file.
  #    All output files have prefix 'test-output' followed by the name
  #    of the classification function it represents, followed by the suffix.
  #    For each classification.fn specified, one output file will be created.
  #
  # Returns:
  #  No return value
  # 
  # Test Result Dump File Format:
  #   The output file is an RData file containing the following objects
  #   1. params: The same object as returned by generate.parameters function
  #   2. avg.auc: A matrix where rows repsent each test run and columns 
  #        represent the average AUC obtained by that algorithm in that test
  #        run. 
  #   3. classification.fn: Name of the function used for classification
  #   4. mixture.type: Mixture option used to produce the dataset
  #   5. noise.gene.count: Number of noise genes added to the dataset 

  source('~/pinpan/code/test-framework/generate-parameters.R')
  source('~/pinpan/code/classify.R')
   
  # Classification & Crossvalidation parameters
  fold = 2
  times = 2  
  betaval = 0.5
  kernel.type = 'linear'
   
  supported.fn = c('graphsvm', 'nbsvm', 'plainsvm')
  if(length(intersect(classification.fn, supported.fn))==0){
    cat('ERROR: The given classification function - ', 
         classification.fn, ' - is not supported. The following are the
         supported functions - ', supported.fn, 'Stopping.\n');
    stop()
  }

  # Generiate different parameters for different runs
  params <- list()
  exp.fold.seq <- c()
  cat('Generating parameters..\n')
  if(exp.gen.method == 'synthesize'){
    source('~/pinpan/code/create-synthetic-data/main.R')
    params <- generate.parameters(mixture.type, rate=0.1, count=num.tests) 
    if(length(params) == 0){
      cat('ERROR in generating parameters. Stopping..\n\n')
      stop()
    }
  }
  else if(exp.gen.method == 'randomize'){
    source('~/pinpan/code/create-randomized-data/main.R')
    exp.fold.seq <- seq(from=0.6, to=1.8, length.out=num.tests)
  }
  else{
    cat('Wrong value for exp.gen.method parameter. Permitted values are
         synthesize or randomize\n')
    stop()
  }

  # Creating dataframe to hold average AUC from each test run 
  # for each classification.fn
  all.avg.auc <- list()
  if('plainsvm' %in% classification.fn){
    all.avg.auc$plainsvm <- matrix(nrow=num.tests, ncol=repeatition)
  }
  if('graphsvm' %in% classification.fn){
    all.avg.auc$graphsvm <- matrix(nrow=num.tests, ncol=repeatition)
  }
  if('nbsvm' %in% classification.fn){
    all.avg.auc$nbsvm <- matrix(nrow=num.tests, ncol=repeatition)
  }

  for(i in c(1:num.tests)){
    cat('-----------------------------------------------------------\n');
    cat('TEST RUN #', i, ' \n')
    cat('--------------------\n\n')
 
    for(j in c(1:repeatition)){
      # Generate the network and expression for each parameter set
      if(exp.gen.method == 'synthesize'){
        case.mu <- params$case.mu.list[i]
        case.sd <- params$case.sd.list[i]
        control.mu <- params$control.mu.list[i]
        control.sd <- params$control.sd.list[i]

        cat('REPEAT ', j, ': Generating Synthetic Dataset..\n\n')
        dataset <- create.synthetic.data(case.mu, case.sd, control.mu, 
                                         control.sd, noise.gene.count,
                                         tofile=FALSE)
      }
      else if(exp.gen.method == 'randomize'){
        exp.fold <- exp.fold.seq[[i]]
        cat('REPEAT ', j, ': Generating Randomized Dataset with fold', 
             exp.fold,' ..\n\n')
        dataset <- create.randomized.data(exp.fold, noise.gene.count,
                                          tofile=FALSE) 
      }

      # Preparing dataset for classification
      X <- subset(dataset$exp.matrix, select = -disease_outcome)       
      Y <- factor(dataset$exp.matrix[ , 'disease_outcome'])


      # Randomizing Input if asked to
      if(randomize != 0){
        cat('RANDOMIZATION: Label Permutation (Type-4)')
        colnames(X) = sample(colnames(X))
      }

      # Creating the map between gene IDs in exp matrix and gene IDs in
      # network as required by pathClass library
      #map = data.frame(probesetID=colnames(X), 
      #        graphID=colnames(dataset$adj.matrix))
      #map = as.matrix(map)
      
      # Running the actual classification algorithm based on input selection
      if('graphsvm' %in% classification.fn){
        cat('\nREPEAT ', j, ': Performing classification using GraphSVM\n\n')
        result <- graphSVM(X, Y, dataset$adj.matrix, betaval, fold, times)
        all.avg.auc$graphsvm[i, j] <- mean(result$auc)
      }
     
      if('nbsvm' %in% classification.fn){
        cat('\nREPEAT ', j, ': Performing classification using NBSVM\n\n')
        result <- nbSVM(X, Y, dataset$adj.matrix, fold, times)
        all.avg.auc$nbsvm[i, j] <- mean(result$auc)
      }

      if('plainsvm' %in% classification.fn){
        cat('\nREPEAT ', j, ': Performing classification using PlainSVM\n\n')
        result <- plainSVM(X, Y, kernel.type, fold, times)
        all.avg.auc$plainsvm[i, j] <- mean(result$auc)  
        cat('\nAVG AUC: ', all.avg.auc$plainsvm[i, j], '\n')
      } 
    } # Repeat loop 
    cat('-----------------------------------------------------------\n');
  } # Test run loop

  # Write the test parameters and average auc values to file
  # Pardon me for the use of "tmp". I needed the variable name classification.fn
  # when saving it for compatability with other codes. FIXME!
  # In the following code classification.fn refers to a character variable.
  tmp <- classification.fn
  for(classification.fn in tmp){
    outfile <- paste('test-output-', classification.fn,'-', 
                     outfile.suffix, '.RData', sep='')
    cat('\nWriting results to ', outfile, '\n\n')
    if(randomize != 0){ 
      classification.fn <- paste(classification.fn,'-nwRandom-',randomize);
    }
    avg.auc <- all.avg.auc[[classification.fn]]
    save(exp.fold.seq, params, avg.auc, classification.fn, noise.gene.count,
         mixture.type, fold, times, file=outfile)
  }
}

## END OF FUNCTION ##
