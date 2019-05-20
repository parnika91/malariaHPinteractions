# Concat runs to study

ConcatRunsToStudy <- data.frame()

for(j in 1:length(studyIDs))
{
  print(studyIDs[j])
  runIDs <- read.table(paste0("/SAN/Plasmo_compare/SRAdb/Output/",studyIDs[j],"/runs_",studyIDs[j],".txt",collapse=''), header = F, sep = ',')
  number_of_runs <- nrow(runIDs)
  if(number_of_runs == 1)
  {
    firstRun <- runIDs[1,1]
    FirstCountfile <- read.table(paste0("/SAN/Plasmo_compare/SRAdb/Output/", studyIDs[j], "/countWithGFF3_", firstRun,".txt",collapse=''), header = T, sep = '\t')
    write.table(FirstCountfile, paste0("/SAN/Plasmo_compare/SRAdb/Output/", studyIDs[j], "/ConcatRunsToStudy_", studyIDs[j],".txt",collapse=''), sep = '\t', row.names=F)
  }
  
  if(number_of_runs > 1)
  {
    firstRun <- runIDs[1,1]
    FirstCountfile <- read.table(paste0("/SAN/Plasmo_compare/SRAdb/Output/",studyIDs[j],"/countWithGFF3_",firstRun,".txt",collapse=''), header = T, sep = '\t')
    ConcatRunsToStudy <- FirstCountfile
    colnames(ConcatRunsToStudy)[6] <- paste0(as.character(firstRun), "_",as.character(studyIDs[j]),collapse='')
    
    for(i in 2:number_of_runs)
    {
      # get runID
      runID <- runIDs[i,1]
      countfile <- read.table(paste0("/SAN/Plasmo_compare/SRAdb/Output/",studyIDs[j],"/countWithGFF3_",runID,".txt",collapse=''), header = T, sep = '\t')
      
      ConcatRunsToStudy[1:nrow(FirstCountfile),(i+5)] <- countfile[,6]
      #ConcatRunsToStudy <- merge(ConcatRunsToStudy, countfile, by = c("seqnames", "start", "end", "width", "strand"))
      colnames(ConcatRunsToStudy)[(i+5)] <- paste0(as.character(runID), "_",as.character(studyIDs[j]),collapse='')
    }
    write.table(ConcatRunsToStudy, paste0("/SAN/Plasmo_compare/SRAdb/Output/", studyIDs[j], "/ConcatRunsToStudy_", studyIDs[j],".txt", collapse=''), sep = '\t', row.names=F)
  }
}

# Merge and aggregate gene names


for( i in 1:length(studyIDs))
{
  study <- read.csv2(paste0("/SAN/Plasmo_compare/SRAdb/Output/", studyIDs[i], "/ConcatRunsToStudy_", studyIDs[i],".txt", collapse=''), sep = '\t', header=T)
  
  library(rtracklayer, quietly = TRUE)
  
  #get host and parasite
  
  host <- positive_experiments[grep(studyIDs[i],positive_experiments[,1]),2]
  para <- positive_experiments[grep(studyIDs[i],positive_experiments[,1]),3]
  
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
  write.table(t.study, paste0("/SAN/Plasmo_compare/SRAdb/Output/",studyIDs[i],"/",studyIDs[i],".txt", collapse = ''), sep = '\t', row.names=T)
}

# Calculate host-parasite proportions

content <- data.frame()
HostParasitePercent <- function(runs)
{
  for( i in 1:length(runs))
  {
    if(file.exists(paste0("/SAN/Plasmo_compare/SRAdb/Output/",study,"/countWithGFF3_",runs[i],".txt",collapse='')))
    {
      print(i)
      count_df <- read.csv2(file=paste0("/SAN/Plasmo_compare/SRAdb/Output/",study,"/countWithGFF3_",runs[i],".txt",collapse=''), sep='\t', header=TRUE)
      
      parasite_rows <- grep(paste0(substr(para,1,3),"_chr",collapse=''), count_df[,1])
      all_parasite <- count_df[parasite_rows,]
      all_host <- count_df[-parasite_rows,]
      parasite_count_percent <- (sum(all_parasite[,6])*100)/ (sum(all_parasite[,6]) + sum(all_host[,6]))
      host_count_percent <- (sum(all_host[,6])*100) / (sum(all_parasite[,6]) + sum(all_host[,6]))
      
      
      content[i,1] <- runs[i]
      content[i,2] <- host_count_percent
      content[i,3] <- parasite_count_percent
      
      colnames(content) <- c("Run", "Host_percent", "Parasite_percent")
    }
  }
  return(content)
}

# Get input read count for a study

study_inputReads <- data.frame()

for(j in 1:length(studyIDs))
{
  input_reads_sum <- 0
  runIDs <- read.table(paste0("/SAN/Plasmo_compare/SRAdb/Output/",studyIDs[j],"/runs_",studyIDs[j],".txt",collapse=''), header = F, sep = ',')
  number_of_runs <- nrow(runIDs)
  
  for(i in 1:number_of_runs)
  {
    # get runID
    runID <- runIDs[i,1]
    
    inputReadLine <- readLines(paste0("/SAN/Plasmo_compare/SRAdb/Output/",studyIDs[j],"/",runID,"_",studyIDs[j],".final.out",collapse=''), n=6)[6]
    inputReads <- as.numeric(strsplit(inputReadLine, split="\t")[[1]][2])
    
    input_reads_sum <- input_reads_sum + inputReads
  }
  
  study_inputReads[j,1] <- studyIDs[j]
  study_inputReads[j,2] <- input_reads_sum
  
}

# Merge host and parasite gene names and number

mergeParasiteGeneNamesAndNumber <- function(parasite_count, host_count, study)
{
  positive_experiments <- read.table("/SAN/Plasmo_compare/SRAdb/Input/positive_experiments.txt", header = F, sep = ',')
  
  require(rtracklayer)
  
  #get parasite and host
  
  para <- positive_experiments[grep(study,positive_experiments[,1]),3]
  host <- positive_experiments[grep(study,positive_experiments[,1]),2]
  
  # parasite gene merge
  p.genes <- import(paste0("/SAN/Plasmo_compare/Genomes/annotation/",para,".gtf", collapse=''), format = "gtf")
  p.genes <- p.genes[p.genes$type%in%"exon"]
  #genes <- genes[which(genes[,"type"] == "exon"),]
  p.genes.df <- as.data.frame(p.genes)
  p.genes.df.gene_name <- p.genes.df[,c("seqnames", "start", "end", "width", "strand", "gene_id")]
  
  p.mergeStudy.genes.df.gene_name <- merge(parasite_count, p.genes.df.gene_name, by = c("seqnames", "start", "end", "width", "strand"))
  
  p.mergeStudy.genes.df.gene_name <- p.mergeStudy.genes.df.gene_name[,6:ncol(p.mergeStudy.genes.df.gene_name)]
  p.mergeStudy.genes.df.gene_name.combineGenes <- data.frame()
  p.mergeStudy.genes.df.gene_name.combineGenes <- aggregate(p.mergeStudy.genes.df.gene_name[,1] ~ gene_id, data = p.mergeStudy.genes.df.gene_name, sum)
  colnames(p.mergeStudy.genes.df.gene_name.combineGenes)[2] <- "count"
  
  # host genes merge
  h.genes <- import(paste0("/SAN/Plasmo_compare/Genomes/annotation/",host,".gtf", collapse=''), format = "gtf")
  h.genes <- h.genes[h.genes$type%in%"exon"]
  #genes <- genes[which(genes[,"type"] == "exon"),]
  h.genes.df <- as.data.frame(h.genes)
  h.genes.df.gene_name <- h.genes.df[,c("seqnames", "start", "end", "width", "strand", "gene_id")]
  
  h.mergeStudy.genes.df.gene_name <- merge(host_count, h.genes.df.gene_name, by = c("seqnames", "start", "end", "width", "strand"))
  
  h.mergeStudy.genes.df.gene_name <- h.mergeStudy.genes.df.gene_name[,6:ncol(h.mergeStudy.genes.df.gene_name)]
  h.mergeStudy.genes.df.gene_name.combineGenes <- data.frame()
  h.mergeStudy.genes.df.gene_name.combineGenes <- aggregate(h.mergeStudy.genes.df.gene_name[,1] ~ gene_id, data = h.mergeStudy.genes.df.gene_name, sum)
  colnames(h.mergeStudy.genes.df.gene_name.combineGenes)[2] <- "count"
  
  
  # if(ncol(study) > 6)
  # {
  #   for(k in 2:(ncol(mergeStudy.genes.df.gene_name)-1))
  #   {
  #     agg <- aggregate(mergeStudy.genes.df.gene_name[,k] ~ gene_name, data = mergeStudy.genes.df.gene_name, sum)
  #     mergeStudy.genes.df.gene_name.combineGenes <- merge(mergeStudy.genes.df.gene_name.combineGenes, agg, by = c("gene_name"))
  #     colnames(mergeStudy.genes.df.gene_name.combineGenes)[k+1] <- "count"
  #   }
  # }
  return(list(p.mergeStudy.genes.df.gene_name.combineGenes, h.mergeStudy.genes.df.gene_name.combineGenes))
}

# MDS for a study

for(i in 4:length(studyIDs))
{
  print(studyIDs[i])
  # Read study with all reads concatenated to all gene names
  study <- read.table(paste0("/SAN/Plasmo_compare/SRAdb/Output/",studyIDs[i], "/", studyIDs[i],".txt",collapse = ''), header = T, sep = '\t')
  study_HPpercent <- read.table(paste0("/SAN/Plasmo_compare/SRAdb/Output/",studyIDs[i], "/hp_percent_", studyIDs[i],".txt",collapse = ''), header = T, sep = '\t') # get all the host and parasite percentages
  # make 1st column = row.names
  # rownames(study) <- study[,1]
  # SRP029990 <- SRP029990[,-1]
  
  # transpose
  # SRP029990.t <- t(SRP029990)
  
  # mds
  d <- 0
  if((nrow(study) != 1))
    d <- dist(study) else
      print("MDS cannot be done for study with 1 run")
  
  if(d != 0)
  {
    fit <- cmdscale(d,eig=TRUE, k=2)
    
    df <- data.frame()
    for(r in 1:nrow(study_HPpercent))
    {
      df[r,1] <- fit$points[r,1]
      df[r,2] <- fit$points[r,2]
      df[r,3] <- study_HPpercent[grep(study_HPpercent[r,1], rownames(study)),3]
      #df[r,4] <- study_HPpercent[grep(study_HPpercent[r,1], rownames(study)),1]
    }
    rownames(df) <- rownames(study)
    df[,4] <- sapply(rownames(study), function(x) strsplit(x, "_")[[1]])[1,]
    colnames(df) <- c("Dim1", "Dim2", "Parasite_percent", "run")
    
    Eigenvalues <- fit$eig
    Variance <- Eigenvalues / sum(Eigenvalues) 
    Variance1 <- 100 * signif(Variance[1], 2)
    Variance2 <- 100 * signif(Variance[2], 2)
    
    # plot: require(ggsci for colour palette)
    require(ggplot2)
    plot_MDS <- ggplot(df, aes(x = Dim1, y = Dim2, label=run)) + geom_point(aes(color = Parasite_percent), size = 5, alpha = 0.7)
    plot_MDS <- plot_MDS + geom_text(aes(label=run),hjust=0, vjust=0, size = 2) + 
      ggtitle(paste0("MDS for host and parasite genes in ",studyIDs[i],collapse = '')) + 
      theme(legend.text=element_text(size=8), legend.title = element_text(size=8)) + 
      xlab(paste0("Variance1 = ",Variance1, "%", collapse = "")) + ylab(paste0("Variance2 = ",Variance2,"%", collapse = ""))
    #+ scale_color_gradient2(breaks = c(0, 25, 50, 75, 1), labels = c("0", "25", "50", "75", "100"),low = "red", mid = "white", high = "blue", space = "Lab", na.value = "grey50", guide = "colourbar", aesthetics = "color")
    ggsave(filename = paste0("MDS_",studyIDs[i],".png", collapse= ''), path = paste0("/SAN/Plasmo_compare/SRAdb/Output/",studyIDs[i], collapse = ''), width = 60, height = 30, units = "cm", plot = print(plot_MDS),limitsize = FALSE)
  }
}
CPCOLS <- c("#1f78b4", "#33a02c", "#e31a1c")

library(ggplot2)

ggplot(iris, aes(Sepal.Length, Petal.Length)) +
      geom_point(aes(col = Species)) +
      scale_colour_manual(values = CPCOLS)

# MDS for host and parasites separately

for(i in 1:length(studyIDs))
{
  #require(scales)
  study <- studyIDs[i]
  study_host <- read.table(paste0("/SAN/Plasmo_compare/SRAdb/Output/",study,"/host_MDS_",study,".txt",collapse=''), header = T, sep = '\t')
  study_para <- read.table(paste0("/SAN/Plasmo_compare/SRAdb/Output/",study,"/parasite_MDS_",study,".txt",collapse=''), header = T, sep = '\t')
  
  study_HPpercent <- read.table(paste0("/SAN/Plasmo_compare/SRAdb/Output/",studyIDs[i], "/hp_percent_", studyIDs[i],".txt",collapse = ''), header = T, sep = '\t') # get all the host and parasite percentages
  
  d.h <- dist(t(study_host)) # euclidean distances between the rows
  fit.h <- cmdscale(d.h,eig=TRUE, k=2)
  
  df.h <- data.frame()
  for(r in 1:nrow(study_HPpercent))
  {
    df.h[r,1] <- fit.h$points[r,1]
    df.h[r,2] <- fit.h$points[r,2]
    df.h[r,3] <- study_HPpercent[grep(study_HPpercent[r,1], rownames(t(study_host))),3]
    #df[r,4] <- study_HPpercent[grep(study_HPpercent[r,1], rownames(study)),1]
  }
  rownames(df.h) <- rownames(t(study_host))
  df.h[,4] <- rownames(t(study_host))
  colnames(df.h) <- c("Dim1", "Dim2", "Parasite_percent", "run")
  
  Eigenvalues.h <- fit.h$eig
  Variance.h <- Eigenvalues.h / sum(Eigenvalues.h) 
  Variance1.h <- 100 * signif(Variance.h[1], 2)
  Variance2.h <- 100 * signif(Variance.h[2], 2)
  
  # plot: require(ggsci for colour palette)
  # require(ggsci)
  require(ggplot2)
  plot.h <- ggplot(df.h, aes(x = Dim1, y = Dim2, label=run)) + geom_point(aes(color = Parasite_percent), size = 5, alpha = 0.7)
  plot.h <- plot.h + geom_text(aes(label=run),hjust=0, vjust=0, size = 2) + ggtitle(paste0("MDS for host genes in ",study,collapse = '')) + xlab(paste0("Variance1 = ",Variance1, "%", collapse = "")) + ylab(paste0("Variance2 = ",Variance2,"%", collapse = "")) + theme(legend.text=element_text(size=8), legend.title = element_text(size=8)) + scale_color_gradient2(breaks = c(0, 25, 50, 75, 1), labels = c("0", "25", "50", "75", "100"),low = "red", mid = "white", high = "blue", space = "Lab", na.value = "grey50", guide = "colourbar", aesthetics = "color")
  ggsave(filename = paste0("MDS_host_",studyIDs[i],".pdf", collapse= ''), path = paste0("/SAN/Plasmo_compare/SRAdb/Output/",studyIDs[i], collapse = ''),plot = print(plot.h),limitsize = FALSE, width = 60, height = 30, units = "cm")
  
  d.p <- dist(t(study_para)) # euclidean distances between the rows
  fit.p <- cmdscale(d.p,eig=TRUE, k=2)
  
  df.p <- data.frame()
  for(r in 1:nrow(study_HPpercent))
  {
    df.p[r,1] <- fit.p$points[r,1]
    df.p[r,2] <- fit.p$points[r,2]
    df.p[r,3] <- study_HPpercent[grep(study_HPpercent[r,1], rownames(t(study_para))),3]
    #df[r,4] <- study_HPpercent[grep(study_HPpercent[r,1], rownames(study)),1]
  }
  rownames(df.p) <- rownames(t(study_para))
  df.p[,4] <- rownames(t(study_para))
  colnames(df.p) <- c("Dim1", "Dim2", "Parasite_percent", "run")
  
  Eigenvalues.p <- fit.p$eig
  Variance.p <- Eigenvalues.p / sum(Eigenvalues.p) 
  Variance1.p <- 100 * signif(Variance.p[1], 2)
  Variance2.p <- 100 * signif(Variance.p[2], 2)
  
  # plot: require(ggsci for colour palette)
  require(ggsci)
  plot.p <- ggplot(df.p, aes(x = Dim1, y = Dim2, label=run)) + geom_point(aes(color = Parasite_percent), size = 5, alpha = 0.7)
  plot.p <- plot.p + geom_text(aes(label=run),hjust=0, vjust=0, size = 2) + ggtitle(paste0("MDS for parasite genes in ",study,collapse = '')) + xlab(paste0("Variance1 = ",Variance1, "%", collapse = "")) + ylab(paste0("Variance2 = ",Variance2,"%", collapse = "")) + theme(legend.text=element_text(size=8), legend.title = element_text(size=8)) + scale_color_gradient2(breaks = c(0, 25, 50, 75, 1), labels = c("0", "25", "50", "75", "100"),low = "red", mid = "white", high = "blue", space = "Lab", na.value = "grey50", guide = "colourbar", aesthetics = "color")
  ggsave(filename = paste0("MDS_para_",studyIDs[i],".pdf", collapse= ''), path = paste0("/SAN/Plasmo_compare/SRAdb/Output/",studyIDs[i], collapse = ''),plot = print(plot.p),limitsize = FALSE, width = 60, height = 30, units = "cm")
  
}


# Remove genes that do not have varying expression from studies

RemoveUnvaryingGenes <- function(studyIDs)
{
  for( i in 1:length(studyIDs))
  {
    study <- read.table(paste0("/SAN/Plasmo_compare/SRAdb/Output/",studyIDs[i], "/", studyIDs[i],".txt",collapse = ''), header = T, sep = '\t')
    variance_table <- as.matrix(apply(study, 2, var))
    
    study_var <- study[,rownames(as.matrix(variance_table[variance_table > 0,]))]
    write.table(study_var, file = paste0("/SAN/Plasmo_compare/SRAdb/Output/", studyIDs[i], "/varyingGenes_", studyIDs[i], ".txt", collapse = ""), sep = '\t', row.names = T)
  }
}