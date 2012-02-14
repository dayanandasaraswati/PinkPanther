## Script that generates a suite of parameters that control
## the properties of datasets that are generated.

generate.parameters <- 
  function(mixture.type=c('snr', 'spread', 'overlap'), rate=0.1, count=10){

  # This function generates a suite of parameters which can be
  # fed as input to synthetize.exp function. synthesize.exp 
  # needs the mean and variance of signal and noise distributions.
  # Mean determines how distant they are from each other and variance
  # affects the spread of each curve. These two parameters control
  # the overlap of signal and noise curves and hence their separability
  # and classifiability. Quantitatively, Signal-to-Noise Ratio denotes
  # how much of signal can be retrieved from the data in the presence
  # of noise. This function provides options to generate different
  # parameters which can achieve different signal and noise characterestics.
  #
  # Args:
  #  mixture.type: character vector which is a subset of 
  #    c('snr', 'spread', 'overlap') controlling how signal and noise curves
  #    are going to mix in the synthesized data. 'snr' option changes the dataset's
  #    signal-to-noise ration which is defined as mean-of-signal/sd-of-noise.
  #    Starting from a SNR of zero, numerator and denominator are increased
  #    linearly and independently to generate various SNR values. The
  #    rate of increase is controlled by the rate parameter.
  #    'spread' option changes the signal's variance linearly whose speed
  #    is determined by rate param. Starting from 1, increase in varaince
  #    will cause the signal to be more spread out. 'overlap' causes the
  #    mean of noise to increase from zero linearly. If the signal's mean
  #    is assumed to be fixed, increase in noise's mean will move the curve
  #    away from signal, hence reducing overlap. Each of the three options
  #    may be used independently or in conjunction with the others. They
  #    will independently control their respective parameters thus producing
  #    many complicated signal and noise mixtures.
  #  TODO: Document about the maximum achievable SNR
  #   
  #  rate: Number between 0 and 1 that tells the rate of change of 
  #    parameters. Its important to visualize the signal and noise
  #    curves for the choosen mixture-model and determine the approriate
  #    rate. For example, if mixture.type is 'overlap', a high rate would
  #    quickly separate signal and noise curves thus making it trivial
  #    for classification algorithms to separate them. With 'snr' and 
  #    'overlap', separate would be even quicker than with just one of
  #    them.
  #
  #  count: Number of parameter sets to generate. 
  #
  # Returns:
  #  List containing <case.mu.list, case.sd.list, 
  #                   control.mu.list, control.sd.list>
  #  
  #  The two 'mu' lists provides mean for the cases and controls and
  #  the two 'sd' lists provides the standard deviation for cases and
  #  controls. Each list has count number of items. The i-th element
  #  in each list together forms a parameter set needed to input
  #  the synthesize.exp function.

  # Amount to increase the numerator and denominator of SNR.
  # The numerator grows twice as fast as denominator because
  # initial value of snr will be 0/1 = 0
  snr.numerator.delta <-  2
  snr.denominator.delta <- (snr.numerator.delta/10)
  # Spread and overlap increment values
  spread.delta <- 1
  overlap.delta <- 1
  
  # Signal and noise distributions are identical initially
  # case is the signal and control is the noise
  case.mu.initial <- 0
  case.sd.initial <- 1
  control.mu.initial <- 0
  control.sd.initial <- 1

  # Extracting mixture type
  if(!is.vector(mixture.type)){
    cat('mixture.type needs to be a vector with atleast one value. Stopping\n')
    return(list());
  }

  is.snr <- 'snr' %in% mixture.type
  is.spread <- 'spread' %in% mixture.type
  is.overlap <- 'overlap' %in% mixture.type

  if(!is.snr && !is.spread && !is.overlap){
    cat('Atleast one of three options for mixture.model is necessary. 
         Stopping\n') 
    return(list())
  }


  # Creating lists with some useless values
  case.mu.list <- rep(case.mu.initial, count)
  case.sd.list <- rep(case.sd.initial, count)
  control.mu.list <- rep(control.mu.initial, count)
  control.sd.list <- rep(control.sd.initial, count)
  
  # The first element in the list is always the initial
  # values. Subsequent elements represent the next parameter
  # sets. 
  i = 2
  while(i <= count){
    # Independently change the four parameters
    if(is.snr){
      case.mu.list[[i]] <- case.mu.list[[i-1]] + rate*snr.numerator.delta
      control.sd.list[[i]] <- 
        control.sd.list[[i-1]] + rate*snr.denominator.delta
    }

    if(is.spread){
      case.sd.list[[i]] <- case.sd.list[[i-1]] + rate*spread.delta
    }

    if(is.overlap){
      control.mu.list[[i]] <- control.mu.list[[i-1]] + rate*overlap.delta
    }
    i <- i + 1
  }
  
  output <- list()
  output$case.mu.list <- case.mu.list
  output$case.sd.list <- case.sd.list
  output$control.mu.list <- control.mu.list
  output$control.sd.list <- control.sd.list

  return(output)
}
## END OF FUNCTION ##
