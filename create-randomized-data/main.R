## XXXXXXXXXXXXXXXXXX COMMENTS ARE TO BE UPDATED XXXXXXXXXXXXXXXXXXXXXXX
## XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
## Script to generate data necessary for use with classification algorithms.
## This is really a wrapper script for synthesize-exp.R script. 
##
## Written by Sanath Kumar Ramesh (skramesh@ucsd.edu) on 31-Jan-2012
##

create.randomized.data <- function(exp.fold, noise.gene.count, tofile=TRUE){

  # This script creates an expr matrix and a network and dumps them to files
  # for use with any classification algorithm.  
  #
  # The following steps are done by this script
  # 1. Loads a real network. Right now its the FunCoup's network which 
  #    Matan extracted and gave it to me in a RData file. Use a different
  #    RData file if you want other networks.
  # 
  # 2. Loads a list of genes investigated in the TCGA BRCA dataset. Again,
  #    they are loaded from a RData file, so for a new dataset, use the
  #    appropriate RData file but make sure the names are loaded in the
  #    variable exp.genes.list.
  #    
  # 3. Loads a list of communities present in the above network
  #
  # 4. Generates a synthetic expression matrix embedding the community
  #    structure. Refer synthesize-exp.R code for details about this step.
  #
  # 5. Dumps the network as an adjacent matrix and exp matrix into a file,
  #    if needed. 
  # 
  # Args:
  #   exp.fold: Non-zero numeric factor by which case's and control's 
  #     will have their expressions scaled. Refer to randomize-exp.R
  #     for detailed documentations
  #
  #   noise.gene.count: Number of noise genes to be added to exp matrix
  #
  #   tofile: Logical value. TRUE if you want the script to dump network and
  #     exp matrix to a file. FALSE if you want no file dump. This value 
  #     defaults to TRUE.
  #
  # Returns:
  #   
  #   List containing <exp.matrix, adj.matrix>
  #   
  #   exp.matrix: A patient by gene expression matrix with the disease 
  #     outcome in a column labelled 'disease_outcome'. The cases are labelled
  #     one and controls are labelled zero
  #   
  #   adj.matrix: A gene by gene matrix with boolean values that represents
  #     network as an adjacency matrix. Genes represent the vertices and a 
  #     one at position (i,j) means an edge between gene at row i and gene 
  #     at column j.

  library(Matrix)

  # Some parameters used in generation of the dataset 
  # Ratio of case and control samples. A two element vector
  # with first element for cases and second for controls.
  case.control.ratio <- c(1,3)
  
  # Load network data if it doesn't already exists
  # Since its a big file, loading it in parent environment
  # will prevent repeated loads.
  if(!exists('fc_network_pkg', envir=environment(create.randomized.data), 
              inherits=FALSE)){
    cat('Loading network..\n')
    load('~/pinpan/data/fc_network_pkg.RData', 
         envir=environment(create.randomized.data))
  }

  fc_network_pkg = get("fc_network_pkg",
                       envir=environment(create.randomized.data))

  # Extract the EntrezID list, network as adjacency matrix
  # and the community list.
  entrez <- as.character(fc_network_pkg$gvect)
  adj.matrix <- fc_network_pkg$adj.mat
  adj.matrix <- round(adj.matrix)

  # Making sure that Entrez IDs correspond to matrix indices
  stopifnot(dim(adj.matrix)[1] == length(entrez))
  rownames(adj.matrix) <- entrez
  colnames(adj.matrix) <- entrez

  # Loading the communities. The gene numbers aren't the actual Entrez ID.
  # They are indices to the Entrez vector which gives actual ID. 
  communities.as.index <- fc_network_pkg$hq.clus
  # Converting indices to actual Entrez IDs
  communities <- sapply(communities.as.index, function(x){
      return(entrez[x])
    });
 
  # Load the genes list used in creation of the exp matrix 
  if(!exists('top_exp_matrix_pos', envir=environment(create.randomized.data),
              inherits=FALSE)){
    cat('Loading Expression Matrix..\n')
    load('~/pinpan/data/read_unified_exp.RData', 
          envir=environment(create.randomized.data))
  }
  
  exp.matrix.original <- get('top_exp_matrix_pos',
                             envir=environment(create.randomized.data))
  # Making sure columns represent genes and rows represent columns
  exp.matrix.original = t(exp.matrix.original)
  # Number of patients present in the exp matrix
  patients.count <- nrow(exp.matrix.original)

  # Map gene IDs from one nomenclature to another. Use this function
  # to ensure that gene IDs used in the network and exp matrix follow
  # same nomenclature
  colnames(exp.matrix.original) <- map.geneIDs(colnames(exp.matrix.original)) 
  
  # Take the sub-matrix of the expression matrix, by filtering for
  # genes whose Hugo names are mappable to Entrez IDs
  exp.matrix.cut <- 
        exp.matrix.original[ , ! is.na(colnames(exp.matrix.original))]
  
  # Assign patients to cases and clusters randomly
  x <- sample(1:patients.count) # Permute patient indices
  # Generate case and control partitions based on the provided ratio
  y <- rmultinom(1, size=length(x), prob=case.control.ratio)
  case.clust <- x[1:y[1]]
  control.clust <- x[(y[1]+1):length(x)] 
 
  # Create the synthetic expression matrix
  cat('Creating the randomized expression matrix..\n\n')
  source('~/pinpan/code/create-randomized-data/randomize-exp.R')
  out.rand.expr <- randomize.exp(exp.matrix.original=exp.matrix.cut, 
                                 communities, case.clust, control.clust, 
                                 fold=exp.fold,noise.gene.count)

  if(length(out.rand.expr) == 0){
    cat('\n\nThere was an creating the expression matrix. Stopping\n')
    return(list())
  }

  exp.matrix <- out.rand.expr$exp.matrix
  signal.genes.list <- out.rand.expr$signal.genes.list
  noise.genes.list <- out.rand.expr$noise.genes.list
  
  # Pruning the network and matrix to retain only genes in common.
  # This includes both the signal and noise genes. Therefore, to add noise
  # genes to the network, simply use the noise.gene.count value when creating
  # the exp matrix
  all.exp.genes <- c(signal.genes.list, noise.genes.list)
  all.network.genes <- rownames(adj.matrix)
  common.genes <- intersect(all.exp.genes, all.network.genes)
  adj.matrix.pruned = adj.matrix[common.genes, common.genes]  
  exp.matrix.pruned = exp.matrix[ , common.genes] 
 
  # Adding the disease outcome to exp matrix
  disease_outcome <- rep(-1, patients.count)
  disease_outcome[case.clust] <- 1
  exp.matrix.aug <- cbind(exp.matrix.pruned, 'disease_outcome'=disease_outcome)
  
  # Writing to file if neecessary
  if(tofile){
    write.table(as.matrix(adj.matrix.pruned), file="network.txt")
    write.table(exp.matrix.aug, file="expression.txt")
  }
  
  # Packing the output variable
  output <- list()
  output$exp.matrix <- exp.matrix.aug
  output$adj.matrix <- as.matrix(adj.matrix.pruned)
 
  return(output)
}
## END OF FUNCTION ##

map.geneIDs <- function(exp.genes.list){

  # Helper function which maps geneIDs from one naming standard to another
  # Mapping HUGO names to Entrez IDs
  map <- read.table('~/pinpan/data/hugo-entrez-map.txt', na.string='null',
                    header=TRUE, sep='\t', comment.char='',
                    stringsAsFactors=FALSE)
  
  name.2.entrez.id <- map$Entrez.Gene.ID
  names(name.2.entrez.id) = map$Approved.Symbol  
 
  return(name.2.entrez.id[exp.genes.list])
 #return(exp.genes.list)
}
