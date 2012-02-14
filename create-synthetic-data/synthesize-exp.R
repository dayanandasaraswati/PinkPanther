## Script that creates a synthetic expression matrix based on 
## gaussian distribution of values and embedding network structure
## in them.
##
## Original Script - Michael Yu (mikeyu@ucsd.edu)
## Modified by - Sanath Kumar Ramesh (skramesh@ucsd.edu)
##

########### TODO: Function description isn't good. Rewrite ########
synthesize.exp <- function(genes.in.exp, patients.count, communities,
                           case.clust, case.mu, case.sd, control.clust, 
                           control.mu, control.sd, noise.gene.count){
  # Create an expression matrix with values for cases and controls patients
  # drawn from two different normal distributions and enriched for network
  # structure. Extra genes are added to the matrix with uncorrelated 
  # expr values (noise). 
  #
  # Args:
  #   genes.in.exp: List of IDs of genes present in an actual expr matrix.
  #     The output expr matrix will contain all the genes in this list which
  #     are mapped to communities. Their expr values will be enriched for
  #     network structure in the output.
  #     
  #   patients.count: A non-zero integer denoting the number of patient 
  #     samples to be present in the output.
  #
  #   communities: A list of gene communities.  communities[[i]] is a
  #     vector of indices representing the i-th community of genes,
  #     e.g. communities[[i]] = c(1,5,6) means that the i-th community
  #     is composed of the genes having gene IDs 1, 5 and 6.
  #     The communities should be non-overlapping and have some genes in
  #     with the gene.in.exp list.
  #     
  #   case.clust: A vector of case patients. This is a vector of
  #     indices representing the cluster of case patients,
  #     e.g. case.clust = c(2,4,10) means the case cluster is composed
  #     of the patients represented by rows 2,4,10 of output expr matrix.
  #     This cluster along with control.clust should be non-overlapping
  #     and their lengths must sum up to patients.count
  #
  #   case.mu, case.sd: Mean and standard deviation for cases. Expr
  #     values for case patients are drawn from a Normal distribution
  #     with parameters case.mu and case.sd
  #   
  #   control.clust: A list of control patients. Similar to case.clust
  #
  #   control.mu, control.sd: Mean and standard deviation for controls. 
  #     Expr values for control patients are drawn from a Normal distribution
  #     with parameters control.mu and control.sd
  #    
  #   noise.gene.count: Max number of noise genes to be added to expr matrix.
  #     Among the genes that are in the communities but not in genes.in.exp
  #     list, noise.gene.count number of them will be selected if available.
  #     Otherwise, only the available genes would be selected for addition 
  #     to output exp matrix.These genes will draw their expr values from a 
  #     standard normal distribution.
  #
  #
  # Returns:
  #
  #   List containing 
  #      <exp.matrix, signal.gene.list, noise.genes.list, clust.2.communities>
  #
  #   exp.matrix: A patient by gene expression matrix where expression of 
  #     case and controls are drawn from two different distributions. Noise 
  #     genes get values from an independent distribution.
  #      
  #   signal.genes.list: Vector of character elements. Each element is the ID
  #     of signal genes in exp.matrix
  #
  #   noise.genes.list: Vector of character elements. Each element is the ID
  #     of noise genes in exp.matrix
  # 
  #   clust.2.communities: List of communities mapped to each case and control
  #     cluster. clust.2.communities[[1]] gives the list for cases and 
  #     clust.2.communities[[2]] gives the list for controls. 
 
  noise.mu <- 0
  noise.sd <- 1 
  # Making sure case.clust and control.clust are non-intersecting and
  # have full coverage of patients.
  all.patients = c(case.clust, control.clust)
  if(!identical(sort(all.patients), c(1:patients.count))){
    cat('Case and control clusters are not formed properly. Possible reasons', 
         'are \n\t 1. case.clust and control.clust is overlapping \n\t',
         '2. case.clust and control.clust do not have full coverage\n',
         'Cannot proceed any further. Stopping..\n')
    return(list())
  }

  # Making sure all gene IDs are characters. If it was an integer, 
  # named indexing of a vector using this integer would fetch
  # element at position given by that integer and not the
  # element at position with the name given by that integer.
  communities <- sapply(communities, function(x){return(as.character(x))}) 
  genes.in.exp <- as.character(genes.in.exp)

  # Selecting common genes common to genes.in.exp and communities
  # and choosing a list of noise genes from the remaining
  all.community.genes <- 
    unique(Reduce(function(x,y){return(c(x,y))}, communities))
  common.genes <- intersect(all.community.genes, genes.in.exp)
  other.genes <- setdiff(all.community.genes, common.genes)
  # Finding the number of genes available for selection as noise gene
  select.count <- 
    ifelse(noise.gene.count>length(other.genes), 
           length(other.genes), noise.gene.count)
  noise.genes.list <- sample(other.genes, size=select.count)
  
  # If there is NO overlap between the genes in expr matrix and 
  # the communities, there is no point in continuing. 
  if(length(common.genes) == 0){
    cat('ERROR: No overlap between genes in expr matrix and genes communities.
         Genes in the matrix have been randomized but none of them have been 
         scaled. Stopping.\n');
    stop()
  }

  # Pruning communities list to include only genes common to communities
  # and genes.in.exp
  communities.pruned <- sapply(communities, function(x){
      return(intersect(x, common.genes))
    });
 
  # Printing some statistics
  cat('Number of Genes in original Expr Matrix: ', 
         length(genes.in.exp), '\n')
  cat('Number of Genes in the communities: ', 
         length(all.community.genes), '\n')
  cat('Number of signal genes in final expr matrix: ',
         length(common.genes), '\n')
  cat('Number of noise genes in final expr matrix: ',
         length(noise.genes.list), '\n')

  # Randomly associate the communities to clusters
  x <- sample(1:length(communities.pruned))
  # Communitites are assigned to both cases and controls with
  # equal probabilities.
  y <- rmultinom(n=1, size=length(communities.pruned), prob=c(1,1))
  clust.2.communities.case <- x[1:y[1]]
  clust.2.communities.control <- x[(y[1]+1):length(x)]
  
  # Creating the patients by genes expression matrix
  # Columns are filled with signal genes followed by noise genes
  exp.matrix <- matrix(nrow=patients.count, 
                       ncol=(length(common.genes)+length(noise.genes.list)))
  rownames(exp.matrix) <- as.character(c(1:patients.count))
  colnames(exp.matrix) <- c(common.genes, noise.genes.list)

  # Boolean vector indicating which genes have already been assigned values.
  # This is used to keep genes found in multiple communities from
  # being assigned expr values multiple times.
  is.assigned <- mat.or.vec(length(common.genes), 1)
  names(is.assigned) <- common.genes

  # Assigning expr values for all communities in case cluster. For 
  # genes in these communities, expr values for case patients is drawn
  # from case's distribution. Rest of the patients ie. control
  # patients draw their expression from control's distributions
  for(comm.index in clust.2.communities.case){
    community <- communities.pruned[[comm.index]]
    
    # Only a subset of case patients should belong to this
    # community. They will get expr values from case's distribution.
    # Others will get expr values from control's distributions    
    case.clust.selected <- sample(case.clust, size=length(case.clust)/10)
    case.clust.unselected <- setdiff(case.clust, case.clust.selected) 
    patients.with.control.expr <- c(control.clust, case.clust.unselected) 

    for(gene in community){
      # Skipping genes which were already assigned
      if(is.assigned[gene] == 1){
        next
      }
      
      is.assigned[gene] <- 1

      exp.matrix[case.clust.selected, gene] <- 
        rnorm(length(case.clust.selected), case.mu, case.sd)
      exp.matrix[patients.with.control.expr, gene] <-
        rnorm(length(patients.with.control.expr), control.mu, control.sd)
    }
  }

  # Assign expr values for all communities in control cluster.
  # Since these communities do not contribute to the disease,
  # both case and control patients draw their expr values from
  # control's distribution. In other words, all patients draw
  # values from control's distribution
  for(comm.index in clust.2.communities.control){
    community <- communities.pruned[[comm.index]]
    
    for(gene in community){  
      # Skipping genes which were already assigned. This condition
      # might araise because of overlapping communities
      if(is.assigned[gene] == 1){
        next
      }

      is.assigned[gene] <- 1 

      exp.matrix[ , gene] <- rnorm(patients.count, control.mu, control.sd) 
    }
  }
  # Assigning expr values for noise genes
  for(gene in noise.genes.list){
    exp.matrix[ , gene] <- rnorm(patients.count, noise.mu, noise.sd);
  }
  
  # Generating the output variable
  clust.2.communities <- list()
  clust.2.communities[[1]] <- clust.2.communities.case
  clust.2.communities[[2]] <- clust.2.communities.control
  output <- list();
  output$exp.matrix <- exp.matrix
  output$signal.genes.list <- common.genes
  output$noise.genes.list <- noise.genes.list
  output$clust.2.communities <- clust.2.communities
  return(output)

" #NOTE: THIS IS A COMMENTED REGION OF CODE
  for(clust.index in 1:length(clust.2.communities)){
    comm.indices <- clust.2.communities[[clust.index]]

    # patients in this cluster
    associated.patients <- clust[[clust.index]]
    r <- length(associated.patients)
    # other patients
    other.patients <- setdiff(1:(dim(exp.matrix)[1]), associated.patients)

    for(comm.index in comm.indices){
      community <- communities[[comm.index]]

      for(gene in community){
 
        gene = as.character(gene)
        
        # Skip if this gene is not in the expression matrix or if it
        # has already been randomized. 
        if ((is.randomized[gene] == 1) || (in.matrix[gene] == 0)){
          next
        }        
        is.randomized[gene] <- 1

        gene.order <- c(gene.order, gene)
        
        z <- mean(exp.matrix[, gene])
        x <- mean(exp.matrix[associated.patients, gene])
        y <- mean(exp.matrix[other.patients, gene])
        b <- z / (y * (fold*(r/n) + (n-r)/n))
        a <- (z * fold) / (x * (fold*(r/n) + (n-r)/n))

        # Check math
        stopifnot(a > 0)
        stopifnot(b > 0)
        stopifnot(abs((x*(r/n) + y*(n-r)/n) - z) < 1e-9)
        stopifnot(abs((a*x*(r/n) + b*y*(n-r)/n) - z) < 1e-9)
        stopifnot(abs((a*x / (b*y)) - fold) < 1e-9)
          
        # Enrich levels in associated patients
        exp.matrix[associated.patients, gene] <- a * exp.matrix[associated.patients, gene]
        # Reduce levels in other patients
        exp.matrix[other.patients, gene] <- b * exp.matrix[other.patients, gene]        
      }
    }
  }
"

}
#### END OF FUNCTION ####



