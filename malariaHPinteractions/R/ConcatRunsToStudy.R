### Concatenate all runs in a study
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("BiocGenerics", version = "3.8")
library(rtracklayer)
library(GenomicFeatures)
library(GenomicRanges)
library(S4Vectors)

# options(echo=TRUE)
# args <- commandArgs(TRUE)
# studyIDs <- args[1]
#positive_experiments <- read.table("/SAN/Plasmo_compare/SRAdb/Input/positive_experiments.txt", header = F, sep = '\t')
#studyIDs <- positive_experiments[,1]
#studyIDs <- studyIDs[-32]
#positive_experiments <- read.delim("/SAN/Plasmo_compare/SRAdb/Input/studies.txt", header = F)
#studyIDs <- c("SRP110609", "ERP017542")
positive_experiments <- read.delim("/SAN/Plasmo_compare/SRAdb/Input/positive_experiments.txt", header = F)

ConcatRunsToStudy <- data.frame()
studyIDs <- c("Macrophage")

for(j in 1:length(studyIDs))
{
  print(studyIDs[j])
  runIDs <- read.table("/SAN/Plasmo_compare/Kai/macrophage_runs.txt", header = F, sep = ',')
  number_of_runs <- nrow(runIDs)
  if(number_of_runs == 1)
  {
    firstRun <- runIDs[1,1]
    FirstCountfile <- read.table(paste0("/SAN/Plasmo_compare/Kai/Mapping/countWithGFF3_", firstRun,".txt",collapse=''), header = T, sep = '\t')
    write.table(FirstCountfile, paste0("/SAN/Plasmo_compare/Kai/Mapping/ConcatRunsToStudy_", studyIDs[j],".txt",collapse=''), sep = '\t', row.names=F)
  }

  if(number_of_runs > 1)
    {
    firstRun <- runIDs[1,1]
    FirstCountfile <- read.table(paste0("/SAN/Plasmo_compare/Kai/Mapping/countWithGFF3_",firstRun,".txt",collapse=''), header = T, sep = '\t')
    ConcatRunsToStudy <- FirstCountfile
    colnames(ConcatRunsToStudy)[6] <- paste0(as.character(firstRun), "_",as.character(studyIDs[j]),collapse='')
    a = 2
   for(i in 2:number_of_runs)
    {
      # get runID
      runID <- runIDs[i,1]
      if(any(grepl(pattern = paste0("countWithGFF3_",runID,".txt",collapse=''), list.files())))
      {
      countfile <- read.table(paste0("/SAN/Plasmo_compare/Kai/Mapping/countWithGFF3_",runID,".txt",collapse=''), header = T, sep = '\t')
  
      ConcatRunsToStudy[1:nrow(FirstCountfile),(a+5)] <- countfile[,6]
      #ConcatRunsToStudy <- merge(ConcatRunsToStudy, countfile, by = c("seqnames", "start", "end", "width", "strand"))
      colnames(ConcatRunsToStudy)[(a+5)] <- paste0(as.character(runID), "_",as.character(studyIDs[j]),collapse='')
      a = a+1
      }
   }
   write.table(ConcatRunsToStudy, paste0("/SAN/Plasmo_compare/Kai/Mapping/ConcatRunsToStudy_", studyIDs[j],".txt", collapse=''), sep = '\t', row.names=F)
  }
}


### Concatenate all studies of the same hp pair

# #hp_pairs[1:nrow(positive_experiments)] <- unique(paste0(positive_experiments[,2], positive_experiments[,3], collapse='')) # possible hp pairs
# hp_pairs <- unique(apply(positive_experiments[,2:3], 1, function(x) paste0(x, collapse='')))
# 
# ConcatSameHPpairs <- data.frame()
# 
# for (k in 1:length(hp_pairs))
# {
#   print(hp_pairs[k])
#   relevantHPpairStudies <- positive_experiments[grep(hp_pairs[k], apply(positive_experiments[,2:3], 1, function(x) paste0(x, collapse=''))),1]
#   concatfile <- read.csv2(paste0("/SAN/Plasmo_compare/SRAdb/Analysis/ConcatRunsToStudy/ConcatRunsToStudy_",relevantHPpairStudies[1],".txt",collapse=''), sep = '\t', header=T)
#   ConcatSameHPpairs <- concatfile[,1:5] # make the first 5 columns that are common to all studies with the same hp pair
#   
#   if(length(relevantHPpairStudies) == 1)
#   {
#     concatfile <- read.csv2(paste0("/SAN/Plasmo_compare/SRAdb/Analysis/ConcatRunsToStudy/ConcatRunsToStudy_",relevantHPpairStudies,".txt",collapse=''), sep = '\t', header=T)
#     ConcatSameHPpairs <- concatfile
#   }
#   if(length(relevantHPpairStudies) > 1)
#   {
#     for(l in 1:length(relevantHPpairStudies))
#     {
#       ConcatRunsToStudy <- read.csv2(paste0("/SAN/Plasmo_compare/SRAdb/Analysis/ConcatRunsToStudy/ConcatRunsToStudy_",relevantHPpairStudies[l],".txt",collapse=''), sep = '\t', header=T)
#       
#       ConcatSameHPpairs <- cbind(ConcatSameHPpairs, ConcatRunsToStudy[,c(6:ncol(ConcatRunsToStudy))])
#     }
#   }
#   write.table(ConcatSameHPpairs, paste0("ConcatSameHPpairs_",hp_pairs[k],".txt",collapse=''), sep = '\t', row.names=F)
# }

### put gene_ids next to counts and sum up counts for same genes
### This one is for SRP029990

# merge and aggregate gene names
#positive_experiments <- read.delim("/SAN/Plasmo_compare/SRAdb/Input/studies.txt", header = F)
#studyIDs <- c("ERP107298", "SRP123516")

# host <- "mouse"
# para <- "Pyoelii"
# 
# genes <- import(paste0("/SAN/Plasmo_compare/Genomes/annotation/mousePyoelii.gtf", collapse=''), format = "gtf")
# genes <- genes[genes$type%in%"exon"]
# #genes <- genes[which(genes[,"type"] == "exon"),]
# genes.df <- as.data.frame(genes)
# genes.df.gene_name <- genes.df[,c("seqnames", "start", "end", "width", "strand", "gene_id")]

for( i in 1:length(studyIDs))
{
  print(studyIDs[i])
  study <- read.csv2(paste0("/SAN/Plasmo_compare/Kai/Mapping/ConcatRunsToStudy_", studyIDs[i],".txt", collapse=''), sep = '\t', header=T)
  
  library(rtracklayer, quietly = TRUE)
  
  #get host and parasite
  
  #host <- as.character(positive_experiments[grep(studyIDs[i],positive_experiments[,1]),2])
  #para <- as.character(positive_experiments[grep(studyIDs[i],positive_experiments[,1]),3])
  
  genes <- import(paste0("/SAN/Plasmo_compare/Genomes/annotation/",host,para,".gtf", collapse=''), format = "gtf")
  genes <- genes[genes$type%in%"exon"]
  #genes <- genes[which(genes[,"type"] == "exon"),]
  genes.df <- as.data.frame(genes)
  genes.df.gene_name <- genes.df[,c("seqnames", "start", "end", "width", "strand", "gene_id")]
  
  mergeStudy.genes.df.gene_name <- merge(study, genes.df.gene_name, by = c("seqnames", "start", "end", "width", "strand"))
  
  mergeStudy.genes.df.gene_name <- mergeStudy.genes.df.gene_name[,6:ncol(mergeStudy.genes.df.gene_name)]
  mergeStudy.genes.df.gene_name.combineGenes <- data.frame()
  mergeStudy.genes.df.gene_name.combineGenes <- aggregate(mergeStudy.genes.df.gene_name[,1] ~ gene_id, data = mergeStudy.genes.df.gene_name, sum)
  colnames(mergeStudy.genes.df.gene_name.combineGenes)[2] <- colnames(mergeStudy.genes.df.gene_name)[1]
  
  if(ncol(study) > 6)
  {
    for(k in 2:(ncol(mergeStudy.genes.df.gene_name)-1))
    {
      agg <- aggregate(mergeStudy.genes.df.gene_name[,k] ~ gene_id, data = mergeStudy.genes.df.gene_name, sum)
      mergeStudy.genes.df.gene_name.combineGenes <- merge(mergeStudy.genes.df.gene_name.combineGenes, agg, by = c("gene_id"))
      colnames(mergeStudy.genes.df.gene_name.combineGenes)[k+1] <- colnames(mergeStudy.genes.df.gene_name)[k]
    }
  }
  
  t.study <- t(mergeStudy.genes.df.gene_name.combineGenes)
  colnames(t.study) <- mergeStudy.genes.df.gene_name.combineGenes$gene_id
  t.study <- t.study[-1,]
  class(t.study) <- "numeric"
  write.table(t.study, paste0("/SAN/Plasmo_compare/Kai/Mapping/",studyIDs[i],".txt", collapse = ''), sep = '\t', row.names=T)
}


