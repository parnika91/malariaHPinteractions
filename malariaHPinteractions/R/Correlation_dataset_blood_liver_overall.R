library(WGCNA)
library(foreach)
library(doParallel)

# function to assign an RData object to a name
# https://stackoverflow.com/questions/5577221/how-can-i-load-an-object-into-a-variable-name-that-i-specify-from-an-r-data-file
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

options(echo=TRUE)
args <- commandArgs()
studyID <- args[1] # for example, for a dataset called SRP1188996.RData, write SRP118996 in the command line
load(paste0("Data/", studyID, ".RData"))
study_sampled <- study[sample(120, 120),]
study <- loadRData(studyID)
study <- t(study)

load("liver.ortho.data.RData")
load("blood.ortho.data.RData")
blood_sampled <- blood.ortho.data[,sample(214, 214)]
blood.liver <- cbind(liver.ortho.data, blood_sampled)
study <- t(blood.liver)

### correlation ####
system.time(ori_cor <- cor(study, use = 'pairwise.complete.obs'))
n <- nrow(study)

reps <- 1000

PermAsso <- function()
{
  perm <- foreach(i=1:reps, .combine = "+") %do%
    {
      blood_s <- blood.ortho.data[,sample(214, 214)]
      blood.liver <- cbind(liver.ortho.data, blood_s)
      study_s <- t(blood.liver)

      (abs(cor(study, study_s[sample(n, n),], use = 'pairwise.complete.obs') >= abs(ori_cor) +0))
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

save(outer, file = paste0("/scratch/zb88/pm4655/cor/outer_",studyID,"_",runif(1, min = 0, max = 10000),".RData", collapse = ''))
