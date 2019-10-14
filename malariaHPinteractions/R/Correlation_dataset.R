library(WGCNA)
library(foreach)
library(doParallel)

studyID <- "SRP118996_SRP032775"
load(paste0("Data/", studyID, ".RData"))
study <- SRP118996_SRP032775
study <- t(study)

### correlation ####
system.time(ori_cor <- cor(study, use = 'pairwise.complete.obs'))
n <- nrow(study)

reps <- 1000

PermAsso <- function()
{
  perm <- foreach(i=1:reps, .combine = "+") %do%
    {
      #rand.col <- host.para[sample(nrow(host.para), (nrow(host.para))),]       
      (abs(cor(study, study[sample(n, n),], use = 'pairwise.complete.obs') >= abs(ori_cor) +0))
    }
  return(perm)
}

PermAssoCompiled <- compiler::cmpfun(PermAsso)

outer_reps <- 10
registerDoParallel()

ptm <- proc.time()
outer <- foreach(j=1:outer_reps, .combine = "+") %dopar%
{ 
  print(j)
  PermAssoCompiled()
}
proc.time() - ptm

save(outer, file = paste0("/short/uv67/pm4655/cor/outer_",studyID,"_",runif(1, min = 0, max = 10000),".RData", collapse = ''))
