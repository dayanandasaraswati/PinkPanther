## This script prepares all the data required to run the classification
## algorithm script. The following operations are performed by this script:
## 
## 1. Loads the RDump of the network from ~/pinpan/.fc_network_pkg.RData.
##    Annotates the adj.matrix with Entrez IDs and writes it to a file.
##
## 2. Dumps communities data got from the fc_network_pkg(above step),
##    to file hq_communities.R file for use by randomize-expression.R script
##
## 3. Runs the randomize-expression.R script which creates the expr-matrix
##    file.
##
## 4. Filter expr matrix and network to have only those genes present in 
##    both of them.
##
## NOTE: Run this script from ~/pinpan/data directory

library(Matrix)

#### STEP - 1
# Load data if it doesn't already exists
if(!exists("fc_network_pkg")){
  cat('Loading network data\n')
  load('~/pinpan/data/fc_network_pkg.RData')
entrez = as.character(fc_network_pkg$gvect)
adj.matrix = fc_network_pkg$adj.mat
adj.matrix = round(adj.matrix)
}

# Making sure that Entrez IDs are mapped to matrix indices properly
stopifnot(dim(adj.matrix)[1] == length(entrez))
rownames(adj.matrix) = entrez;
colnames(adj.matrix) = entrez;

#### STEP - 2
communities = fc_network_pkg$hq.clus
dump("communities", file="hq_communities.R");

#### STEP - 3
source('~/pinpan/code/create-randomized-data/randomize-expression.R')

# Expression matrix
load(file='read_unified_exp.RData')

### Use Hugo <-> Entrez ID mapping from Hugo's website (more complete)
#x <- read.table('hugo_2_entrez.txt', na.strings='null', header=TRUE,sep='\t',comment.char='', stringsAsFactors=FALSE)
x <- read.table('hugo-entrez-map.txt', na.strings='null', header=TRUE,sep='\t',comment.char='', stringsAsFactors=FALSE)
entrez.id.2.name <- x$Approved.Symbol
names(entrez.id.2.name) <- x$Entrez.Gene.ID
# Vector whose elements are Entrez IDs.  The vector's names are Hugo gene names
name.2.entrez.id <- x$Entrez.Gene.ID
names(name.2.entrez.id) <- x$Approved.Symbol

# In communities, convert Entrez IDs to Hugo names
communities <- sapply(communities, function(x){return(as.vector(entrez.id.2.name[as.character(x)]))})
# Remove genes whose Entrez IDs aren't mappable to Hugo name
communities <- sapply(communities, function(x){return(x[! is.na(x)])})

fold.list <- list()
rmatrix.list <- list()
rmatrix.ordered.list <- list()

# Generate random clustering
num.of.clusters <- 2
x <- sample(1:dim(top_exp_matrix_pos)[2])  # Permute patient indices
y <- rmultinom(1, length(x), rep(1,num.of.clusters))
random.clust <- list()
for(i in 1:length(y)){
  random.clust[[i]] <- x[1:y[i]]
  x <- x[(y[i]+1):length(x)]
}

# Map the Hugo names of expression matrix to Entrez IDs
rownames(top_exp_matrix_pos) = name.2.entrez.id[rownames(top_exp_matrix_pos)]
# exp_matrix_cut = top_exp_matrix_pos

# Take the sub-matrix of the expression matrix, by filtering for
# genes whose Hugo names are mappable to Entrez IDs
exp_matrix_cut <- top_exp_matrix_pos[! is.na(rownames(top_exp_matrix_pos)), ]

# Iterate over fold ratios
fold.vect <- as.character(10)
for(fold in fold.vect){
  print(c('fold', fold))

  # Randomize expression (Note: use transpose)
  tmp <- randomize.exp(t(exp_matrix_cut), communities, random.clust, as.numeric(fold))
  rmatrix <- tmp[[1]]
  rmatrix <- t(rmatrix)
  clust.2.communities <- tmp[[2]]
  gene.order <- tmp[[3]]
  other.genes <- setdiff(rownames(exp_matrix_cut), gene.order)

  print('done randomizing')
  
  # Order randomized matrix to look nice if printed as a heatmap
#  rmatrix.ordered <- rmatrix[c(gene.order, other.genes), Reduce(function(a,b){return(c(a,b))}, random.clust)]
    
  rmatrix.list[[fold]] <- rmatrix[c(gene.order, sample(other.genes)),] ## NOTE: Selecting only scaled genes
#  rmatrix.ordered.list[[fold]] <- rmatrix.ordered

}
print('hi')

#### STEP - 4 
# All the variables for this step come from the sourced file
for(fold in fold.vect){
  expr = t(rmatrix.list[[fold]])
  expr.matrix.genes = colnames(expr);
  network.genes = colnames(adj.matrix);
  common.genes = intersect(expr.matrix.genes, network.genes);

  expr =  expr[, common.genes]
  adj.matrix.filtered = adj.matrix[common.genes, common.genes]
  
  # Adding disease outcome to expression matrix and writing to file
  disease_outcome = rep(0, dim(expr)[1]);
  disease_outcome[random.clust[[1]]] = 1; # Assuming that cluster-1 is the case
  # NOTE: expr is appended the disease_outcome column. So, both expr and
  # expr.data.matrix contain the SAME matrix.
  expr.data.matrix = cbind(expr, 'disease_outcome' = disease_outcome);
  cat('Size of filtered expr profile is ', length(common.genes), '\n');
  write.table(expr.data.matrix, file=sprintf('sanath_random_exp.fold_%s.txt', fold))

  # Writing network to file
  # NOTE: Matrix's sparsity will be lost on writing to file. 
  # Must use sparse matrices when reading
  # XXX: ADD FOLD NUMBER TO FILENAME
  write.table(as.matrix(adj.matrix.filtered), file="fc_network_filtered.txt");

}
