## Script to read test results and visualize them as plots

visualize.results <- function(rdata.file=c('test-output.RData'),
                              plot.filename.suffix=''){

  # This script reads the RData file where the test results 
  # were stored and makes a bunch of plots to visualize
  # the dataset and the results. Particularly, this script
  # generates the following plots:
  # 1. Box and whisker plot for Classification AUCs written to 
  #    a file with base-name "box". (Refer documentation of 
  #    plot.filename.suffix argument for details abt base-name)
  # 2. Plots showing the signal and noise curves for first ten 
  #    params from the parameter set. The output file will have
  #    ten different plots pasted next to each other. Each plot
  #    will have both signal and noise curves together. They have
  #    a base-name "param-viz"
  #
  # Args:
  #   rdata.file: List of test output filenames
  #
  #   plot.filename.suffix: All plots are written as a png file
  #     whose filename has a base-name and a suffix separated by a
  #     hyphen. A good descriptive suffix will help user identify
  #     the plot. This has to be a valid filename string.

  library(ggplot2)
  library(gridExtra)
  
  avg.auc.list <- list()
  classification.fn.list <- list()
  # Because all rdata files contains variables with same names, the
  # variables are deleted before loading the file
  for(i in c(1:length(rdata.file))){
   load(rdata.file[[i]])
    avg.auc.list[[i]] <- avg.auc
    classification.fn.list[[i]] <- classification.fn
    rm(avg.auc)
    rm(classification.fn)
    
  }
  
  # AUC Plot
  ylabel <- paste('AUC averaged over', fold, 'repeats of CV')
  xlabel <- 'Test Runs'

  for(i in c(1:length(rdata.file))){
    png(paste('box-',classification.fn.list[[i]],'-',
              plot.filename.suffix,'.png',sep=''))
    boxplot.matrix(avg.auc.list[[i]], use.cols=FALSE, xlab=xlabel, ylab=ylabel,
      main=paste('Performance of ',classification.fn.list[[i]])) 
    dev.off()
  }
  
  fn.name.vec <- sapply(classification.fn.list,
                                     function(x){return(x)})
  # Plotting a combined AUC graph
  ylabel = paste('Average AUC')
  mean.auc.frame <- sapply(avg.auc.list, function(x){return(apply(x,1,mean))})
  jpeg(paste('combined-auc-',plot.filename.suffix,'.jpeg',sep=''))
  matplot(c(1:nrow(mean.auc.frame)), mean.auc.frame, type='l',
          xlab=xlabel, ylab=ylabel)
  legend(x="topleft", legend=fn.name.vec, lty=c(1:length(fn.name.vec)))
  dev.off()
                   
  ##### Since signal-noise plot won't change often, I am commenting it for now.
  ##### Uncomment it when needed ####
  return()

  ### UNREACHABLE CODE ###
  # Signal-Noise Plot
  plots <- list()
  for(i in c(1:min(c(dim(avg.auc)[1], 10)))){
    # Finding a good range for x-axis such that the bell shape of 
    # both curves are visible. For each curve, one standard deviation 
    # around the mean will be visible. 
    mu.1 <- params$case.mu.list[i]
    sd.1 <- params$case.sd.list[i]
    mu.2 <- params$control.mu.list[i]
    sd.2 <- params$control.sd.list[i]
    #cat(mu.1, sd.1, mu.2, sd.2, '\n') 
    left.limit <- min(c(mu.1 -2* sd.1, mu.2 - 2*sd.2))
    right.limit <- max(c(mu.1 + 2*sd.1, mu.2 + 2*sd.2))
 
    x <- seq(left.limit, right.limit, length=200)
    signal <- dnorm(x, mu.1, sd.1)
    noise <- dnorm(x, mu.2, sd.2)
 
    plot.df <- data.frame(x, signal, noise)
    xlabel = sprintf('Run #%d: Sig(%.2f,%.2f); Bck(%.2f,%.2f)',i,mu.1,sd.1,mu.2,sd.2)
    plots[[i]] <- ggplot(plot.df, aes(x = x)) +                         
                geom_line(aes(y = signal), colour = 'blue') +
                labs(x=xlabel, y='') +
                opts(axis.text.y=theme_blank()) +
                opts(axis.ticks=theme_blank()) +  
                geom_line(aes(y = noise), colour = 'red') +
                geom_area(aes(y = pmin(signal, noise)), fill = 'gray60') 
  }

  png(paste('signal-noise-',plot.filename.suffix,'.png',sep=''))
  grob = do.call(arrangeGrob, c(plots, ncol=2)) 
  grid.arrange(grob,
    main='Signal(blue), Background(red) curves for various runs')
  dev.off()
}

## END OF FUNCTON ##
