## Script that randomizes expression matrices to be enriched for
## network structure.  Runs NMF on randomizations.  See the
## documentation on randomize.exp for more information.
##
## XXXXXXXXXXXXXXXXXX DOCUMENTATION ISN'T CORRECT. XXXXXXXXXXXX
## XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

randomize.exp <- function(exp.matrix.original, communities, case.clust, 
                          control.clust, fold, noise.gene.count){
  # Randomize a matrix of gene expression levels to enrich for network
  # structure.
  #
  # Args:
  #   exp.matrix: A patients-by-genes matrix of gene expression levels
  #
  #   communities: A list of gene communities.  communities[[i]] is a
  #     vector of indices representing the i-th community of genes,
  #     e.g. communities[[i]] = c(1,5,6) means that the i-th community
  #     is composed of the genes represented in columns 1,5,6 of
  #     exp.matrix.  The communities should be non-overlapping and
  #     have full coverage of the genes in exp.matrix
  #
  #   clust: A list of patient clusters.  clust[[i]] is a vector of
  #     indices representing the i-th cluster of patients,
  #     e.g. clust[[i]] = c(2,4,10) means the i-th cluster is composed
  #     of the patients represented by rows 2,4,10 of exp.matrix.  The
  #     clusters should be non-overlapping and have full coverage of
  #     the patients in exp.matrix
  #
  #   fold: The desired ratio, for any gene, the average level of that
  #     gene across patients associated with the gene's community to
  #     that across the other patients.
  #
  #
  # Returns:
  #
  #   A randomized matrix based on exp.matrix such that for every
  #   gene, the average expression is preserved while the expression
  #   is enriched in a randomly chosen cluster of patients.
  #
  #   To do the randomization, first randomly associate each community
  #   to a cluster.  Then, for every gene in every community, let
  #
  #   r : the size of the cluster associated to the gene's community
  #   n : the total number of patients
  #   z : average expression level across all patients
  #   x : average level in patients associated with the gene's community
  #   y : average level in other patients
  #
  #   Note by definition, x*(r/n) + y*(n-r)/n = z
  #
  #   Next, randomly permute (shuffle) the levels of each gene across
  #   all patients.  Then, amplify the levels in the associated
  #   patients by the factor <a>, and reduce the levels in the other
  #   patients by <b>, such that
  #
  #   1<=a and 0<b<=1
  #   ax*(r/n) + by*(n-r)/n = z
  #   (ax)/(by) = fold
  #
  #   In this way, the average expression is preserved and the desired
  #   ratio <fold> is achieved.
  #
  #   Note, <a> and <b> can be solved as
  #   b = z / (y * (fold*(r/n) + (n-r)/n))
  #   a = (z * fold) / (x * (fold*(r/n) + (n-r)/n))

  noise.mu <- 0
  noise.sd <- 1  
 
  patients.count <- nrow(exp.matrix.original)

  # Making sure case.clust and control.clust are non-intersecting and
  # have full coverage of patients.
  all.patients = c(case.clust, control.clust)
  if(!identical(sort(all.patients), c(1:patients.count))){
    cat('Case and control clusters are not formed properly. Possible reasons', 
         'are \n\t 1. case.clust and control.clust is overlapping \n\t',
         '2. case.clust and control.clust do not have full coverage\n',
         'Cannot proceed any further. Stopping..\n')
    stop()
  }
  
  # Permute each gene's levels across all patients
  for(i in 1:ncol(exp.matrix.original)){
    exp.matrix.original[,i] <- sample(exp.matrix.original[,i])
  }

  # Making sure that communities list is a list of character
  # elements. If the communities list had an integer, 
  # named indexing of a vector using this integer would fetch
  # element at position given by that integer and not the
  # element at position with the name given by that integer.
  communities = sapply(communities, function(x){return(as.character(x))}) 
  # colnames will also be used for named indexing
  colnames(exp.matrix.original) = as.character(colnames(exp.matrix.original))
  genes.in.exp <- colnames(exp.matrix.original)
                    

  # Selecting genes common to genes.in.exp and communities
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
  # and genes in expr matrix
  communities.pruned <- sapply(communities, function(x){
      return(intersect(x, common.genes))
    });

 
  # After pruning, many of the original communities might have zero size
  # or some might be really small, like one node. Removing them. 
  min.genes.in.comm = 5  
  indices <- 1:length(communities.pruned)
  to.remove <- which(sapply(indices, 
                      function(x){
                        ifelse(length(communities.pruned[[x]])<min.genes.in.comm,
                          TRUE, FALSE)
                      }))
  communities.pruned <- communities.pruned[-to.remove] 
  num.pruned.genes <- length(Reduce(
                        function(x,y){return(c(x,y))}, communities.pruned))
  # Printing some statistics
  cat('Number of Genes in original Expr Matrix: ', 
         length(genes.in.exp), '\n')
  cat('Number of Genes that will be randomized: ', 
         length(num.pruned.genes), '\n')
  cat('Number of signal genes in final expr matrix: ',
         length(common.genes), '\n')
  cat('Number of noise genes in final expr matrix: ',
         length(noise.genes.list), '\n')

  # Randomly associate the communities to clusters
  x <- sample(1:length(communities.pruned))
  # Communitites are assigned to both cases and controls with
  # equal probabilities.
  y <- rmultinom(n=1, size=length(communities.pruned), prob=c(1,1))
  #clust.2.communities.case <- x[1:y[1]]
  #clust.2.communities.control <- x[(y[1]+1):length(x)]

  # HACK:: Overriding default assignment. Assigning communities 
  # without replacement
  clust.2.communities.case <- sample(x, 2)
  x <- setdiff(x, clust.2.communities.case)
  clust.2.communities.control <- sample(x, 2)
 
  cat('$$$$ CONTROL GENES SIZE ', sum(sapply(clust.2.communities.case, 
        function(x){length(communities.pruned[[x]])})), '\n\n')
   cat('$$$$ CASE GENES SIZE ', sum(sapply(clust.2.communities.case, 
        function(x){length(communities.pruned[[x]])})), '\n\n')

  # Creating a submatrix for the noise genes and binding it to the
  # original exp matrix
  noise.matrix <- matrix(nrow=patients.count, 
                       ncol=length(noise.genes.list))
  rownames(noise.matrix) <- rownames(exp.matrix.original)
  colnames(noise.matrix) <- noise.genes.list
  exp.matrix <- cbind(exp.matrix.original, noise.matrix) 
 

  # Boolean vector indicating which genes have already been randomized.
  # This is used to keep genes found in multiple communities from
  # being randomized multiple times.
  is.randomized <- mat.or.vec(length(common.genes), 1)
  names(is.randomized) <- common.genes

  n <- nrow(exp.matrix)

  # Creating variables that are convenient to iterate on
  clust.2.communities <- list() 
  clust.2.communities[[1]] <- clust.2.communities.case
  clust.2.communities[[2]] <- clust.2.communities.control
  clust <- list()
  clust[[1]] <- case.clust
  clust[[2]] <- control.clust

  for(clust.index in 1:length(clust.2.communities)){
    comm.indices <- clust.2.communities[[clust.index]]

    # patients in this cluster
    patients.in.clust <- clust[[clust.index]]
   

    for(comm.index in comm.indices){
      community <- communities.pruned[[comm.index]]

        for(gene in community){
 
        # Skip if this gene is not in the expression matrix or if it
        # has already been randomized. 
        if(is.randomized[gene] == 1){
          next
        }
        

        # If the clust is a cluster of cases, then I'm selecting
        # only a subset of patients and associated them to this community
        if(clust.index){
          associated.patients <- sample(patients.in.clust, #size=10)
                                 size=length(patients.in.clust)/2)
        }
        else{
          associated.patients <- patients.in.clust
        }
    
        r <- length(associated.patients)
        # other patients
        other.patients <- setdiff(1:patients.count, associated.patients)
 


        is.randomized[gene] <- 1
        
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
        exp.matrix[associated.patients, gene] <- 
                   a * exp.matrix[associated.patients, gene]
        # Reduce levels in other patients
        exp.matrix[other.patients, gene] <- 
                   b * exp.matrix[other.patients, gene]        
      }
    }
  }

  # Assigning expr values for noise genes
  for(gene in noise.genes.list){
    exp.matrix[ , gene] <- rnorm(patients.count, noise.mu, noise.sd);
  }
  
  # Generating the output variable
  output <- list();
  output$exp.matrix <- exp.matrix
  output$signal.genes.list <- common.genes
  output$noise.genes.list <- noise.genes.list
  output$clust.2.communities <- clust.2.communities
  return(output)

}
#### END OF FUNCTION ####



