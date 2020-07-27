### Concatenate all runs in a study
library(rtracklayer)
library(GenomicFeatures)
library(GenomicRanges)
library(S4Vectors)

options(echo=TRUE)
args <- commandArgs(TRUE)
studyIDs <- args[1]
positive_experiments <- read.delim("/SAN/Plasmo_compare/SRAdb/Input/positive_experiments.txt", sep = ",", header = F)
allHPexp <- read.delim("/SAN/Plasmo_compare/SRAdb/Output/allHPexp.txt", header = T)


################# Step 1: Bring all runs together ###################

ConcatRunsToStudy <- data.frame()

for(j in 1:length(studyIDs))
{
  print(studyIDs[j])
  runIDs <- read.table(paste0("/SAN/Plasmo_compare/SRAdb/Output/", studyIDs[j], "/runs_",studyIDs[j],".txt", collapse = ''), header = F, sep = ',')
  number_of_runs <- nrow(runIDs)
  if(number_of_runs == 1)
  {
    firstRun <- runIDs[1,1]
    FirstCountfile <- read.table(paste0("/SAN/Plasmo_compare/SRAdb/Output/", studyIDs[j], "/countWithGFF3_",firstRun,".txt",collapse=''), header = T, sep = '\t')
    write.table(FirstCountfile, paste0("/SAN/Plasmo_compare/SRAdb/Output/_", studyIDs[j],"/ConcatRunsToStudy_", studyIDs[j],".txt"), sep = '\t', row.names=F)
  }

  if(number_of_runs > 1)
    {
    firstRun <- runIDs[1,1]
    if(file.exists(paste0("/SAN/Plasmo_compare/SRAdb/Output/", studyIDs[j], "/countWithGFF3_",firstRun,".txt",collapse='')))
    {
    FirstCountfile <- read.table(paste0("/SAN/Plasmo_compare/SRAdb/Output/", studyIDs[j], "/countWithGFF3_",firstRun,".txt",collapse=''), header = T, sep = '\t')
    ConcatRunsToStudy <- FirstCountfile
    colnames(ConcatRunsToStudy)[6] <- paste0(as.character(firstRun), "_",as.character(studyIDs[j]),collapse='')
    }
    a = 2
    
   for(i in 2:number_of_runs)
    {
    runID <- runIDs[i,1]
     # get runID
      if(file.exists(paste0("/SAN/Plasmo_compare/SRAdb/Output/", studyIDs[j], "/countWithGFF3_",runID,".txt",collapse='')))
      {
      countfile <- read.table(paste0("/SAN/Plasmo_compare/SRAdb/Output/", studyIDs[j], "/countWithGFF3_",runID,".txt",collapse=''), header = T, sep = '\t')
  
      ConcatRunsToStudy[1:nrow(FirstCountfile),(a+5)] <- countfile[,6]
      #ConcatRunsToStudy <- merge(ConcatRunsToStudy, countfile, by = c("seqnames", "start", "end", "width", "strand"))
      colnames(ConcatRunsToStudy)[(a+5)] <- paste0(as.character(runID), "_",as.character(studyIDs[j]),collapse='')
      a = a+1
      }
   }
   write.table(ConcatRunsToStudy, paste0("/SAN/Plasmo_compare/SRAdb/Output/", studyIDs[j], "/ConcatRunsToStudy_", studyIDs[j],".txt", collapse=''), sep = '\t', row.names=F)
  }
}


################################# Step 2: Get gene names for all reads #################

for( i in 1:length(studyIDs))
{
  print(studyIDs[i])
  study <- read.csv2(paste0("/SAN/Plasmo_compare/SRAdb/Output/", studyIDs[i], "/ConcatRunsToStudy_", studyIDs[i],".txt", collapse=''), sep = '\t', header=T)
  
  library(rtracklayer, quietly = TRUE)
  
  #get host and parasite
  
  host <- as.character(positive_experiments[grep(studyIDs[i],positive_experiments[,1]),2])
  para <- as.character(positive_experiments[grep(studyIDs[i],positive_experiments[,1]),3])
  
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
  write.table(t.study, paste0("/SAN/Plasmo_compare/SRAdb/Output/", studyIDs[i], "/", studyIDs[i],".txt", collapse=''), sep = '\t', row.names=T)
}


#################################### Step 3: Only keep coding genes ########################################

################ host-parasite pairs ###############

require(rtracklayer)
require(dplyr)

for(i in 1:length(studyIDs))
{
  print(i)
  # get study.txt including all runs
  study <- as.data.frame(t(read.delim(paste0("Output/", studyIDs[i], "/", studyIDs[i], ".txt", collapse = ''), sep = '\t', header = T)))
  # get host parasite from allHPexp
  #hp <- as.character(unique(allHPexp[allHPexp$Study==studyIDs[i],"HostParasite"]))
  hp <- paste(as.character(positive_experiments[grep(studyIDs[i],positive_experiments[,1]),2]), as.character(positive_experiments[grep(studyIDs[i],positive_experiments[,1]),3]), sep="")
  
  if(hp == "humanPfalciparum")
  { 
    coding = as.data.frame(import("/SAN/Plasmo_compare/Genomes/annotation/humanPfalciparum.gtf")) %>%
      filter(type%in%"exon") %>%
      filter(gene_biotype%in%"protein_coding") %>%
      distinct(gene_id)
  }
  if(hp == "humanPberghei")
  { 
    coding = as.data.frame(import("/SAN/Plasmo_compare/Genomes/annotation/humanPberghei.gtf")) %>%
      filter(type%in%"exon") %>%
      filter(gene_biotype%in%"protein_coding") %>%
      distinct(gene_id) 
  }
  if(hp == "humanPvivax")
  { 
    coding = as.data.frame(import("/SAN/Plasmo_compare/Genomes/annotation/humanPvivax.gtf")) %>%
      filter(type%in%"exon") %>%
      filter(gene_biotype%in%"protein_coding") %>%
      distinct(gene_id) 
  }
  if(hp == "mousePberghei")
  { 
    coding = as.data.frame(import("/SAN/Plasmo_compare/Genomes/annotation/mousePberghei.gtf")) %>%
      filter(type%in%"exon") %>%
      filter(gene_biotype%in%"protein_coding") %>%
      distinct(gene_id) 
  }
 if(hp == "mousePyoelii")
  { 
    coding = as.data.frame(import("/SAN/Plasmo_compare/Genomes/annotation/mousePyoelii.gtf")) %>%
      filter(type%in%"exon") %>%
      filter(gene_biotype%in%"protein_coding") %>%
      distinct(gene_id) 
  }
  if(hp == "mousePchabaudi")
  { 
    coding = as.data.frame(import("/SAN/Plasmo_compare/Genomes/annotation/mousePchabaudi.gtf")) %>%
      filter(type%in%"exon") %>%
      filter(gene_biotype%in%"protein_coding") %>%
      distinct(gene_id) 
  }
  if(hp == "monkeyPcoatneyi")
  { 
    coding = as.data.frame(import("/SAN/Plasmo_compare/Genomes/annotation/monkeyPcoatneyi.gtf")) %>%
      filter(type%in%"exon") %>%
      filter(gene_biotype%in%"protein_coding") %>%
      distinct(gene_id) 
  }
  if(hp == "monkeyPcynomolgi")
  { 
    coding = as.data.frame(import("/SAN/Plasmo_compare/Genomes/annotation/monkeyPcynomolgi.gtf")) %>%
      filter(type%in%"exon") %>%
      filter(gene_biotype%in%"protein_coding") %>%
      distinct(gene_id) 
  }
  if(hp == "monkeyPknowlesi")
  { 
    coding = as.data.frame(import("/SAN/Plasmo_compare/Genomes/annotation/monkeyPknowlesi.gtf")) %>%
      filter(type%in%"exon") %>%
      filter(gene_biotype%in%"protein_coding") %>%
      distinct(gene_id) 
  }
  
  # filter study to keep only protein-coding genes
  study_coding_genes <- study %>%
    tibble::rownames_to_column('gene') %>%
    filter(rownames(study)%in%coding$gene_id) %>%
    tibble::column_to_rownames('gene')
   
  # write the table out
  write.table(study_coding_genes, paste0("/SAN/Plasmo_compare/SRAdb/Output/", studyIDs[i],"/", studyIDs[i], "_coding_genes.txt", collapse = ''), sep ='\t', row.names = T)
}

############################## Step 4: Get orthologous groups for each study ########################

parasite_orthogroups <- read.delim("/SAN/Plasmo_compare/OrthoFinder/parasite_orthogroups.txt", stringsAsFactors=FALSE) 
host_orthogroups <- read.delim("/SAN/Plasmo_compare/OrthoFinder/host_orthogroups.txt", stringsAsFactors=FALSE) 

for(i in 1:length(studyIDs))
{ print(i)
  #if the study_coding_genes exists, merge with orthogroups (join functions)
  filepath = paste0("/SAN/Plasmo_compare/SRAdb/Output/",studyIDs[i],"/",studyIDs[i],"_coding_genes.txt", collapse = "")
  if(file.exists(filepath))
  {
    # find out what host and parasite the study is
    host <- as.character(positive_experiments[grep(pattern = studyIDs[i], positive_experiments[,1]),2])
    para <- as.character(positive_experiments[grep(pattern = studyIDs[i], positive_experiments[,1]),3])

    # take the host and para orthogroups and make a df -> orthogroup | gene name 

    h.df <- data.frame(Orthogroup = host_orthogroups[,1], Org = host_orthogroups[,grep(pattern = host, colnames(host_orthogroups))])
    p.df <- data.frame(Orthogroup = parasite_orthogroups[,1], Org = parasite_orthogroups[,grep(pattern = para, colnames(parasite_orthogroups))])

    hp.df <- rbind(h.df, p.df)

    # read table
    file = read.delim(filepath, header = T) %>% tibble::rownames_to_column("Gene")

    ortho.table = merge(file, hp.df, by.x = "Gene", by.y = "Org")
    ortho.table <- data.frame(Gene = ortho.table$Gene, Orthogroup = ortho.table$Orthogroup, ortho.table[,2:(ncol(ortho.table)-1)])

    write.table(ortho.table, paste0("/SAN/Plasmo_compare/SRAdb/Output/",studyIDs[i],"/",studyIDs[i],"_orthogroups.txt", collapse = ""), sep = '\t', row.names = F)
  }

}
