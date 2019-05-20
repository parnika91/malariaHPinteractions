# without bootstrap
# Forking: it copies the existing version of R, your entire workspace exists in each process.

# for every study individually

library(doParallel)
library(doMC)
library(foreach)
library(data.table)
library(WGCNA)
library(qgraph)
#library(profvis)

#library(Matrix)
#library(qlcMatrix)

# options(echo = TRUE)
# args <- commandArgs(TRUE)
# studyIDs <- args[1]

# study <- read.csv2(paste0("Output/ERP020067/ERP020067_filteredDisp.txt", collapse = ''), header = T, sep = '\t') 
# study <- study[,1:ncol(study)-1]
# study.na.omit <- study[,which(colSums(study) > 5)]
# study.na.omit <- study.na.omit[which(rowSums(study.na.omit) > 5),]
# study <- study.na.omit
#drop.cols <- c("study_disp", "disp_rank")
study <- read.csv2(paste0("Output/ERP020067/ERP020067_disp.txt", collapse = ''), sep = "\t", header = T)
system.time(ori_cor <- WGCNA::cor(t(study)))

  reps <- 10
  
  registerDoParallel(cores = 5)
  
  PermAsso <- function()
  {
    res <- list()
    system.time(res <- foreach(i=1:reps, .combine = '+') %dopar%
    {
      rand.col <- study[,sample(ncol(study), (ncol(study)))] %>%
        setNames(., colnames(study))
      
      (abs(WGCNA::cor(t(study), t(rand.col)) >= abs(ori_cor) +0))
    })
    system.time(red <- Reduce('+', res))
    return(res)
  }
  
  PermAssoCompiled <- compiler::cmpfun(PermAsso)
  
  outer_list <- list()
  outer_reps <- 1
  
  ptm <- proc.time() 
  for(i in 1:outer_reps)
  {
    outer_list[[i]] <- PermAssoCompiled()
  }
  proc.time() - ptm
  stopImplicitCluster()
  
  ptm1 <- proc.time()
  outer <- Reduce('+', outer_list)
  proc.time() - ptm1

pval <- outer/(reps*outer_reps)
ori_cor_copy <- ori_cor
ori_cor[which(pval>0.00)] = 0

  pdf("ERP020067_qgraph.pdf", onefile = T)
  qgraph(ori_cor, layout = "spring", minimum = 0.01, details = T)
  dev.off()

###################################################### with reduce #####################################################


library(doParallel)
#library(doMC)
library(foreach)
#library(data.table)

# options(echo = TRUE)
# args <- commandArgs(TRUE)
# studyIDs <- args[1]
# 
study <- read.csv2("/SAN/Plasmo_compare/SRAdb/Output/ERP020067/ERP020067.txt", header = T, sep = '\t')

# mat <- matrix(runif(20000, 0, 10000), byrow = T, ncol = 40)
# colnames(mat) <- c(1:ncol(mat))
# rownames(mat) <- c(1:nrow(mat))
# study <- mat
# remove the rows and cols with zero, otherwise foreach fails
study.na.omit <- study[,which(colSums(study) > 5)]
study.na.omit <- study.na.omit[which(rowSums(study.na.omit) > 5),]
study <- t(study.na.omit) # runs as columns and genes as rows

ori_cor <- cor(t(study)) # genes vs genes
#reps <- 10
reps_outer <- 10

#cl <- makeCluster(10, type = "FORK")
registerDoParallel(cores = 3)
#mcoptions <- list(preschedule=FALSE)
#.options.multicore=mcoptions

#set.seed(1234)
res <- data.frame()

ptm <- proc.time()
reduce.2 <- foreach(i=1:10) %do% # 500,000 reps
{
  print(i)
  reduce.1 <- foreach(icount(reps_outer)) %dopar% # outer reps are 100 each. Reduces list of 100 reduced matrices into 1: 10,000
  {
    # internal reps are 100 each
    #res <- foreach(icount(reps)) %do%
    #{
      rand.col <- study[,sample(colnames(study), ncol(study))]
      colnames(rand.col) <- colnames(study)
      rownames(rand.col) <- rownames(study)
      #new_cor <- cor(study, study[rand.col,]) # works, but check with tiny dataset
      
      #((abs(new_cor) >= abs(ori_cor))*1) #matrix
      #res1 <- if (res == TRUE) 1 else 0
      (abs(cor(t(study), t(rand.col))) >= abs(ori_cor) *1)
      #sh[,sort(colnames(sh))]
      #return(res)
    
    
    #Reduce('+', res) # reducing list of 100 matrices (res) into 1
  }
  Reduce('+', reduce.1) # reduces list of 10 matrices into 1
}


###################################### Feb 6 ########################################

library(doParallel)
library(foreach)

study <- read.csv2("/SAN/Plasmo_compare/SRAdb/Output/ERP020067/ERP020067.txt", header = T, sep = '\t')

# mat <- matrix(runif(20000, 0, 10000), byrow = T, ncol = 40)
# colnames(mat) <- c(1:ncol(mat))
# rownames(mat) <- c(1:nrow(mat))
# study <- mat
# remove the rows and cols with zero, otherwise foreach fails
study.na.omit <- study[,which(colSums(study) > 5)]
study.na.omit <- study.na.omit[which(rowSums(study.na.omit) > 5),]
study <- t(study.na.omit) # runs as columns and genes as rows

ori_cor <- cor(t(study)) # genes vs genes
reps_outer <- 10

registerDoParallel(cores = 3)

PermAsso <- function()
{
  res <- foreach(icount(inner_reps)) %dopar% # list of 50 large matrices
  {
    rand.col <- study[,sample(colnames(study), ncol(study))]
    colnames(rand.col) <- colnames(study)
    rownames(rand.col) <- rownames(study)
    sh <- (abs(cor(t(study), t(rand.col))) >= abs(ori_cor) *1)
    sh[,sort(colnames(sh))]
  }
  reduce.1 <- Reduce('+', res)
  return(reduce.1)
}
PermAssoCompiled <- compiler::cmpfun(PermAsso)


Outer_reduce <- matrix(rep(0,(ncol(ori_cor)**2), nrow = ncol(ori_cor)))
outer_reps = 100
inner_reps = 5

ptm <- proc.time()
for(j in 1:outer_reps)
{
  print(j)
  Outer_reduce <- Outer_reduce + PermAssoCompiled()
}
proc.time() - ptm

## partial correlation ##
# with ERP020067 #

load("Output/ERP020067/ERP020067cor.RData") # has sum, not "p-val"
pval <- Outer_reduce/100 # get "p-val"

# We can either do PCIT and then use WGCNA to make the network
# Or we can find the partial correlation and then use

# But first we need the highly ranked host-parasite pairs from "p-val" --> ori_cor

# ++++++++++++++++++++++++++++
# flattenCorrMatrix from www.sthda.com/english/wiki/correlation-matrix-an-r-function-to-do-all-you-need and 
# http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

rank.pairs <- flattenCorrMatrix(ori_cor, ERP020067)

##### change cor at high pval to zero #####
## rcorr from Hmisc gives pvalues with cor ##

ori <- ori_cor
p <- pval
dim1.gene <- rownames(ori_cor)
diag(ori_cor) <- 1
diag(pval) <- 1
ori_cor[pval > 0.05] <- 0
diag(ori_cor) <- 1

# from this improved cor matrix, either do cor2pcor, or use qgraph
# 1. using cor2pcor(cor_mat, tol = 0.05) and then maybe try a GGM here

pc <- corpcor::cor2pcor(ori_cor)
g <- graph_from_adjacency_matrix(pc)
plot.igraph(g) # not fully convinced




# 2. qgraph from http://sachaepskamp.com/files/Cookbook.html

library("qgraph")
Graph_pcor <- qgraph(ori_cor, graph = "pcor", layout = "spring") # works

Graph_pcor <- qgraph(ori_cor, graph = "pcor", layout = "spring", threshold = "bonferroni",
                     sampleSize = nrow(ori_cor), alpha = 0.05) # problem with sample_size

Graph_lasso <- qgraph(ori_cor, graph = "glasso", layout = "spring", tuning = 0.25,
                      sampleSize = 10) # does not work because not positive definite

# Using GeneNet


library(GeneNet)
true.pcor <- ggm.simulate.pcor(20, 0.1)
test.results <- ggm.list.edges(true.pcor)[1:19,]

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Rgraphviz", version = "3.8")

nlab <- LETTERS[1:20]
gr <- network.make.graph( test.results, nlab)
gr
num.nodes(gr)
edge.info(gr)
gr2 <- network.make.graph( test.results, nlab, drop.singles=TRUE)
gr2
num.nodes(gr2)
edge.info(gr2)
#plot network
#NOTE: this requires the installation of the "Rgraphviz" library
library("Rgraphviz")
plot(gr, "fdp")
plot(gr2, "fdp")

# PCIT with igraph

pc <- PCIT::pcit(ori_cor)
signif <- idx(pc)
PCIT::plotCorCoeff(ori_cor, list("PCIT Meaningful" = signif), col = c("red"))
nonsignif <- idxInvert(nrows(ori_cor), signif)
ori_cor[nonsignif] <- 0
adj <- ori_cor # or abs(cor)
library(igraph)
plot(graph_from_adjacency_matrix(adj, mode = "undirected", diag = T, weighted = T))


#### network with igraph from https://stackoverflow.com/questions/49171958/igraph-edge-width-and-color-positive-and-negative-correlation #######

set.seed(123)

cor.matrix <- matrix(runif(100, -1, 1), nrow=10)

t = which(abs(cor.matrix) > 0.6 & lower.tri(cor.matrix),arr.ind=TRUE)
t <- cbind(t, cor.matrix[which(abs(cor.matrix) > 0.6 & lower.tri(cor.matrix),arr.ind=TRUE)]) ##this adds the correlation to the graph as an edge attribute "V3"
t.graph=graph.data.frame(t,directed=F)
E(t.graph)$color <- ifelse(E(t.graph)$V3 > 0,'magenta','green') #You had this as "V3 > 0.6" which I guess works but it is more readable as 0. that way if you decide to lower the correlation threshold you do not have to change this line too.

#t.names <- colnames(cor.matrix)[as.numeric(V(t.graph)$name)]
minC <- rep(-Inf, vcount(t.graph))
maxC <- rep(Inf, vcount(t.graph))
minC[1] <- maxC[1] <- 0
l <- layout_with_fr(t.graph, minx=minC, maxx=maxC,
                    miny=minC, maxy=maxC)      
plot(t.graph, layout=l, 
     rescale=T,
     asp=0,
     edge.arrow.size=0.5, 
     vertex.label.cex=0.8, 
     vertex.label.family="Helvetica",
     vertex.label.font=2,
     #vertex.label=t.names,
     vertex.shape="circle", 
     vertex.size=3, 
     vertex.color="deepskyblue2",
     vertex.label.color="black", 
     #edge.color=E(t.graph)$color, ##do not need this since E(t.graph)$color is already defined.
     edge.width=as.integer(cut(abs(E(t.graph)$V3), breaks = 5)))


########### partial correlation in permutation ############
ERP <- read.delim("ERP020067.txt")
ERP_sh <- ERP
rownames(ERP_sh) <- sample(rownames(ERP), nrow(ERP))

lm1 <- lm(t(ERP) ~ t(ERP_sh)) # gives 0 to 0 residuals
ERP.pcor <- ppcor::pcor(t(ERP))

#ori_pcor <- cor2pcor()

##############################################################################
########### partial correlation second order using RLowPC ####################

# install package
library(devtools)
install_github('wyguo/RLowPC')
library(RLowPC)
#requires
library(corpcor)

# perm function
library(doParallel)
library(foreach)

study <- read.csv2("/SAN/Plasmo_compare/SRAdb/Output/ERP020067/ERP020067.txt", header = T, sep = '\t')

# mat <- matrix(runif(20000, 0, 10000), byrow = T, ncol = 40)
# colnames(mat) <- c(1:ncol(mat))
# rownames(mat) <- c(1:nrow(mat))
# study <- mat
# remove the rows and cols with zero, otherwise foreach fails
study.na.omit <- study[,which(colSums(study) > 5)]
study.na.omit <- study.na.omit[which(rowSums(study.na.omit) > 5),]
study <- study.na.omit # runs as columns and genes as rows
study_copy <- study

ori_cor <- cor(study) # genes vs genes
inner_reps <- 10

registerDoParallel(cores = 3)

PermAsso <- function()
{
  res <- foreach(icount(inner_reps)) %dopar% # list of 50 large matrices
  {
    rownames(study_copy) <- sample(rownames(study), nrow(study))
    (abs(cor(study, study_copy)) >= abs(ori_cor) *1)
    #sh[,sort(colnames(sh))]
  }
  reduce.1 <- Reduce('+', res)
  return(reduce.1)
}
PermAssoCompiled <- compiler::cmpfun(PermAsso)


Outer_reduce <- matrix(rep(0, ncol(ori_cor)**2), nrow = ncol(ori_cor))
outer_reps = 10

ptm <- proc.time()
for(j in 1:outer_reps)
{
  print(j)
  Outer_reduce <- Outer_reduce + PermAssoCompiled()
}
proc.time() - ptm
pval <- Outer_reduce/100 # pval is symmetric

# zero order PC

zeroPC<-function(data.exp,method='pearson'){
  message('Calulate zero order partial correlation (the general correlation)...')
  inf.cor<-cor(data.exp, method = method)
  
  ptm <- proc.time()
  for(j in 1:outer_reps)
  {
    print(j)
    Outer_reduce <- Outer_reduce + PermAssoCompiled()
  }
  proc.time() - ptm
  pval <- Outer_reduce/100
  
  inf.cor[pval > 0.05] <- 0

  diag(inf.cor)<-0
  inf.edge<-adjmatrix2edgelist(inf.cor,directed = F,order = T)
  ##calculate statistics, pval, fdr and prob
  fdr.frame<-cor2statistics(cor.vector = abs(inf.edge$weight),plot=F,verbose=F)
  inf.edge$pval = fdr.frame$pval
  inf.edge$qval = fdr.frame$qval
  inf.edge$prob = 1 - fdr.frame$lfdr
  rownames(inf.edge)<-NULL
  return(inf.edge)
}



# second order PC

secondPC<-function(data.exp,edgelist,controlist=NULL,method='pearson',progressbar=T){
  edgelist.new<-edgelist
  ####
  t1 <- Sys.time()
  ##
  inf.zeroPC<-zeroPC(data.exp=data.exp,method = method)
  
  ##calculate second order PC
  message(paste0('Calulate second order partial correlation for ',nrow(edgelist),' pairs of genes ...'))
  if(progressbar)
    pb <- txtProgressBar(min = 0, max =nrow(edgelist), style = 3)
  for(i in 1:nrow(edgelist)){
    Sys.sleep(0)
    node1<-as.vector(edgelist$from[i])
    node2<-as.vector(edgelist$to[i])
    sub.edge<-shared.neighbour(node1,node2,edgelist,verbose = F)
    if(dim(sub.edge)[1]<2){
      edgelist.new[i,3]<-inf.zeroPC[inf.zeroPC$from==node1 &inf.zeroPC$to==node2,]$weight
    } else {
      if(is.null(controlist))
        regulators<-as.vector(unique(sub.edge$from))
      if(!is.null(controlist))
        regulators<-controlist[[i]]
      if(node1 %in% regulators)
        regulators<-regulators[-which(regulators %in% node1)]
      if(node2 %in% regulators)
        regulators<-regulators[-which(regulators %in% node2)]
      targets<-c(node1,node2)
      if(length(regulators)==1){
        reg.group<-combn(x = regulators,1) } else {
          reg.group<-combn(x = regulators,2)
        }
      cor.v<-p.v<-vector()
      for(j in 1:ncol(reg.group)){
        sub.reg<-reg.group[,j]
        second.stat<-ppcor::pcor(data.exp[,c(sub.reg,targets)],method = method)
        second.stat.cor<-second.stat$estimate
        dimnames(second.stat.cor)<-list(c(sub.reg,targets),c(sub.reg,targets))
        cor.v<-c(cor.v,abs(second.stat.cor[node1,node2]))
      }
      edgelist.new[i,3]<-cor.v[abs(cor.v)==max(abs(cor.v))][1]
    }
    if(progressbar)
      setTxtProgressBar(pb, i)
  }
  if(progressbar)
    close(pb)
  edgelist.new<-edgelist.new[order(abs(edgelist.new$weight),decreasing = T),]
  fdr.frame<-cor2statistics(cor.vector = abs(edgelist.new$weight),plot=F,verbose=F)
  edgelist.new$pval = fdr.frame$pval
  edgelist.new$qval = fdr.frame$qval
  edgelist.new$prob = 1 - fdr.frame$lfdr
  rownames(edgelist.new)<-NULL
  t2 <- Sys.time()
  message(paste0('Done! Time taken:',round(as.numeric(difftime(t2,t1)),4),' ',units(difftime(t2,t1))))
  return(edgelist.new)
}


sparse.cor4 <- function(x){
  n <- nrow(x)
  cMeans <- colMeans(x)
  covmat <- (as.matrix(crossprod(x)) - n*tcrossprod(cMeans))/(n-1)
  sdvec <- sqrt(diag(covmat)) 
  cormat <- covmat/tcrossprod(sdvec)
  list(cov=covmat,cor=cormat)
}

X <- sample(0:10,1e7,replace=T,p=c(0.9,rep(0.01,10)))
x <- Matrix(X,ncol=10)
system.time(corx <- sparse.cor4(x))
system.time(corx <- cor(as.matrix(x)))

sparse.cor3 <- function(x){
  memory.limit(size=10000)
  n <- nrow(x)
  
  cMeans <- colMeans(x)
  cSums <- colSums(x)
  
  # Calculate the population covariance matrix.
  # There's no need to divide by (n-1) as the std. dev is also calculated the same way.
  # The code is optimized to minize use of memory and expensive operations
  covmat <- tcrossprod(cMeans, (-2*cSums+n*cMeans))
  crossp <- as.matrix(crossprod(x))
  covmat <- covmat+crossp
  
  sdvec <- sqrt(diag(covmat)) # standard deviations of columns
  covmat/crossprod(t(sdvec)) # correlation matrix
}



n = 400;
m = 1e6;

# Generate data
mat = as.matrix(study);
# Start timer
tic = proc.time();
# Center each variable
mat = mat - rowMeans(mat);
# Standardize each variable
mat = mat / sqrt(rowSums(mat^2));   
# Calculate correlations
cr = tcrossprod(mat);
# Stop timer
toc = proc.time();

# Show the results and the time
show(cr[1:4,1:4]);
show(toc-tic)
