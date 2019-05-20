library(igraph)

# data
head(mtcars)

# Make a correlation matrix:
mat=cor(t(mtcars[,c(1,3:6)]))

mat[mat<0.995]=0
corr <- matrix(c(1.00, 0.89, 0.67, 0.00, 0.60, 0.00, 0.89, 1.00, 0.78, 0.07, 0.50, 0.75, 0.67, 0.78, 1.00, 0.06, 0.55, 0.46, 0.00, 0.07,
                 0.06, 1.00, 0.50, 0.73, 0.60, 0.50, 0.55, 0.50, 1.00, 0.70, 0.00, 0.75, 0.46, 0.73, 0.70, 1.00), nrow=6, ncol=6, byrow=T)
row_names <- c("HGene1", "HGene2", "HGene3", "PGene1", "PGene2", "PGene3")
dimnames(corr) <- list(row_names, row_names)
# Make an Igraph object from this matrix:
network=graph_from_adjacency_matrix(corr, weighted=T, mode="undirected", diag=F)
ceb <- cluster_edge_betweenness(network) 

dendPlot(ceb, mode="hclust")