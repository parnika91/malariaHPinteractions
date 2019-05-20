library(microbenchmark)

set.seed(42)
mat1 <- matrix(runif(4000, 0, 1000), ncol = 100)
set.seed(424)
mat2 <- matrix(runif(4000, 0, 1000), ncol = 100)

# corStats <- stats::cor(mat1, mat2)
#
library(HiClimR)
# corHiClimR <- HiClimR::fastCor(mat1, mat2)
# 
library(coop)
# corCoop <- coop::pcor(mat1, mat2)
# 
library(WGCNA)
# corWGCNA <- WGCNA::corFast(mat1, mat2, method = "pearson", use = "pairwise.complete.obs")
# 
library(ff)
library(propagate)
# corBigcor <- bigcor(mat1, mat2, fun = "cor", size = 2000)

bigcorr <- function(
  x, 
  y = NULL,
  fun = c("cor", "cov"), 
  size = 2000, 
  verbose = TRUE)
{
  fun <- match.arg(fun)
  if (fun == "cor") FUN <- cor else FUN <- cov
  if (fun == "cor") STR <- "Correlation" else STR <- "Covariance" 
  if (!is.null(y) & NROW(x) != NROW(y)) stop("'x' and 'y' must have compatible dimensions!")
  
  NCOL <- ncol(x)
  if (!is.null(y)) YCOL <- NCOL(y)
  
  ## calculate remainder, largest 'size'-divisible integer and block size
  REST <- NCOL %% size
  LARGE <- NCOL - REST  
  NBLOCKS <- NCOL %/% size
  
  ## preallocate square matrix of dimension
  ## ncol(x) in 'ff' single format
  if (is.null(y)) resMAT <- ff(vmode = "double", dim = c(NCOL, NCOL))  
  else resMAT <- ff(vmode = "double", dim = c(NCOL, YCOL))
  
  ## split column numbers into 'nblocks' groups + remaining block
  GROUP <- rep(1:NBLOCKS, each = size)
  if (REST > 0) GROUP <- c(GROUP, rep(NBLOCKS + 1, REST))
  SPLIT <- split(1:NCOL, GROUP)
  
  ## create all unique combinations of blocks
  COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
  COMBS <- t(apply(COMBS, 1, sort))
  COMBS <- unique(COMBS)  
  if (!is.null(y)) COMBS <- cbind(1:length(SPLIT), rep(1, length(SPLIT)))
  
  ## initiate time counter
  timeINIT <- proc.time() 
  
  ## iterate through each block combination, calculate correlation matrix
  ## between blocks and store them in the preallocated matrix on both
  ## symmetric sides of the diagonal
  for (i in 1:nrow(COMBS)) {
    COMB <- COMBS[i, ]    
    G1 <- SPLIT[[COMB[1]]]
    G2 <- SPLIT[[COMB[2]]]    
    
    ## if y = NULL
    if (is.null(y)) {
      if (verbose) message("bigcor: ", sprintf("#%d: %s of Block %s and Block %s (%s x %s) ... ", 
                                               i, STR,  COMB[1], COMB[2], length(G1),  length(G2)))
      RES <- FUN(x[, G1], x[, G2])
      resMAT[G1, G2] <- RES
      resMAT[G2, G1] <- t(RES) 
    } else ## if y = smaller matrix or vector  
    {
      if (verbose) message("bigcor: ", sprintf("#%d: %s of Block %s and 'y' (%s x %s) ... ", i, STR,  COMB[1],
                                               length(G1),  YCOL))    
      RES <- FUN(x[, G1], y)
      resMAT[G1, ] <- RES             
    }
    
    if (verbose) {
      timeNOW <- proc.time() - timeINIT
      message("bigcor: ", round(timeNOW[3], 2), " sec\n")
    }
    
    #gc()
  } 
  
  return(resMAT)
}
# for 1 matrix
mbm <- microbenchmark("corHiClimR" = {corHiClimR <- HiClimR::fastCor(mat1)}, # better than stats and coop
                      "corCoop" = {corCoop <- coop::pcor(mat1)},
                      "corWGCNA" = {corWGCNA <- WGCNA::cor(mat1)}, # corFast doesn't work
                      #"corBigcor" = {corBigcor <- bigcorr(mat1, mat1, fun = "cor", size = 100)}, # can't work without mat2 - not efficient for real data
                      "corStats" = {corStats <- stats::cor(mat1)}, # mat2
                      "corBicor" = {corBicor <- WGCNA::bicor(mat1)}) #mat2

# Unit: milliseconds
# expr       min        lq      mean   median       uq      max neval
# corHiClimR 1406.4179 1445.5099 1675.3978 1511.178 1639.434 2639.891   100
# corCoop  529.4712  540.3692  601.4164  544.977  611.233 1425.713   100
# corStats 4081.1386 4123.7648 4282.0548 4142.420 4274.770 5185.328   100
# corBicor 1120.8601 1143.7789 1298.7243 1178.990 1264.425 2306.373   100


library(ggplot2)
autoplot(mbm, title = "Benchmarking correlation computation")

library(doParallel)
library(foreach)

ori_cor <- cor(mat1, use="pairwise.complete.obs")

registerDoParallel(cores = 10)
#mcoptions <- list(preschedule=FALSE)
reps <- 10
#set.seed(1234)

#combineFun <- function() mat

PermAsso <- function(method)
{
  res <- foreach(icount(reps), .combine = "+") %dopar% # , .options.multicore=mcoptions
  {
    #rand.col <- sample(colnames(t(study)), ncol(t(study)))
    
    #new_cor <- cor(study, study[rand.col,]) # works, but check with tiny dataset
    
    #((abs(new_cor) >= abs(ori_cor))*1) #matrix
    #res1 <- if (res == TRUE) 1 else 0
    if(method == "HiClimR")
      HiClimR::fastCor(mat1)
    
    if(method == "coop")
      coop::pcor(mat1)
    
    if(method == "stats")
      stats::cor(mat1)
    
    if(method == "Bicor")
      WGCNA::bicor(mat1, use = "pairwise.complete.obs")
    
    if(method == "WGCNA")
      WGCNA::cor(mat1, mat2, use = "pairwise.complete.obs")
    
    if(method == "Bigcorr")
      bigcorr(mat1, mat2, fun = "cor", size = 100)
    
    #return(res)
  }
  return(res)
}

PermAssoCompiled <- compiler::cmpfun(PermAsso)
# ptm <- proc.time()
# res <- PermAssoCompiled("HiClimR")
# proc.time() - ptm

mbm_par <- microbenchmark("corHiClimR" = {corHiClimR <- PermAssoCompiled(method = "HiClimR")}, # better than stats and coop
                      "corCoop" = {corCoop <- PermAssoCompiled(method = "coop")},
                      #"corWGCNA" = {corWGCNA <- PermAssoCompiled(method = "WGCNA")},
                      #"corBigcor" = {corBigcor <- PermAssoCompiled(method = "Bigcorr")}, # mat2
                      "corStats" = {corStats <- PermAssoCompiled(method = "stats")}, # mat2
                      "corBicor" = {corBicor <- PermAssoCompiled(method = "Bicor")}) # mat2
