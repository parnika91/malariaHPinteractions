# data("swiss")
# head(swiss)
# 
# 
# #Classical MDS
# 
# # Load required packages
library(magrittr)
library(dplyr)
library(ggpubr)
# # Cmpute MDS
# mds <- swiss %>%
#   dist() %>%          
#   cmdscale() %>%
#   as_tibble()
# colnames(mds) <- c("Dim.1", "Dim.2")
# # Plot MDS
# ggscatter(mds, x = "Dim.1", y = "Dim.2", 
#           label = rownames(swiss),
#           size = 1,
#           repel = TRUE)
# # K-means clustering
# clust <- kmeans(mds, 3)$cluster %>%
#   as.factor()
# mds <- mds %>%
#   mutate(groups = clust)
# # Plot and color by groups
# ggscatter(mds, x = "Dim.1", y = "Dim.2", 
#           label = rownames(swiss),
#           color = "groups",
#           palette = "jco",
#           size = 1, 
#           ellipse = TRUE,
#           ellipse.type = "convex",
#           repel = TRUE)
# 
# # Non-classical MDS
# library(MASS)
# mds <- swiss %>%
#   dist() %>%          
#   isoMDS() %>%
#   .$points %>%
#   as_tibble()
# colnames(mds) <- c("Dim.1", "Dim.2")
# # Plot MDS
# # Kruskal
# ggscatter(mds, x = "Dim.1", y = "Dim.2", 
#           label = rownames(swiss),
#           size = 1,
#           repel = TRUE)
# # Sammon
# # Cmpute MDS
# library(MASS)
# mds <- swiss %>%
#   dist() %>%          
#   sammon() %>%
#   .$points %>%
#   as_tibble()
# colnames(mds) <- c("Dim.1", "Dim.2")
# # Plot MDS
# ggscatter(mds, x = "Dim.1", y = "Dim.2", 
#           label = rownames(swiss),
#           size = 1,
#           repel = TRUE)

library(ggplot2)
library(reshape2)
library(gridExtra)

positive_experiments <- read.table("Input/studies.txt", header = F, sep = ',')

option(echo = TRUE)
studyIDs <- args[1]

#studyIDs <- positive_experiments[,1]

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

# separate parasite and host gene counts. MDS for only host genes and for only parasite genes. Graphs plotting Run index vs parasite genes expressed in each run in a study

# define function for aggregating gene names

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

for (i in 16:length(studyIDs))
{
  print(as.character(studyIDs[i]))
  study <- studyIDs[i]
  
  # get number of runs for studyID from runs_$study file
  runIDs <- read.table(paste0("/SAN/Plasmo_compare/SRAdb/Output/",study,"/runs_",study,".txt",collapse=''), header = F, sep = ',')
  number_of_runs <- nrow(runIDs)
  
  # for MDS
  host_MDS <- data.frame()
  para_MDS <- data.frame()
  
  # go to each run file.
  # sepatare host and parasite parts.
  # sum(gene_counts)
  host_parasite_genecount <- data.frame()
  for(j in 1:number_of_runs)
  {
    #print(j)
    # get runID
    runID <- runIDs[j,1]
    
    count_file <- read.table(paste0("/SAN/Plasmo_compare/SRAdb/Output/",study,"/countWithGFF3_",runID,".txt",collapse=''), header = T, sep = '\t')
    parasite_rows <- as.numeric(grep("P+", count_file[,1]))
    host_count <- count_file[-c(parasite_rows),]
    parasite_count <- count_file[parasite_rows,] 
    
    # get how many genes have count > 5 for parasite genes
    mergeParasiteGeneNamesAndNumberResult <- mergeParasiteGeneNamesAndNumber(parasite_count, host_count, study)
    
    para_genecount <- mergeParasiteGeneNamesAndNumberResult[[1]]
    paraGeneGreaterThanFive <- nrow(subset(para_genecount, count > 5))
    
    host_genecount <- mergeParasiteGeneNamesAndNumberResult[[2]]
    hostGeneGreaterThanFive <- nrow(subset(host_genecount, count > 5))
    
    host_sum <- sum(host_count[,6])
    parasite_sum <- sum(parasite_count[,6]) # total sum of all genes expressed
    host_parasite_genecount[j,1] <- runID
    host_parasite_genecount[j,2] <- host_sum
    host_parasite_genecount[j,3] <- parasite_sum
    host_parasite_genecount[j,4] <- paraGeneGreaterThanFive
    host_parasite_genecount[j,5] <- hostGeneGreaterThanFive
    
    # for MDS for host and parasite separately
    host_MDS[1:nrow(host_genecount),j] <- host_genecount[,2]
    para_MDS[1:nrow(para_genecount),j] <- para_genecount[,2]
    
    colnames(host_MDS)[j] <- as.character(runID)
    colnames(para_MDS)[j] <- as.character(runID)
    
  }
  
  rownames(host_MDS) <- host_genecount[,1]
  rownames(para_MDS) <- para_genecount[,1]
  
  colnames(host_parasite_genecount) <- c("Run", "Host_count", "Parasite_count", "ParaGenesGreaterThanFive", "HostGenesGreaterThanFive")
  write.table(host_parasite_genecount, paste0("/SAN/Plasmo_compare/SRAdb/Output/",study,"/host_parasite_genecount_",study,".txt",collapse=''), row.names = F, sep = '\t')
  write.table(host_MDS, paste0("/SAN/Plasmo_compare/SRAdb/Output/",study,"/host_MDS_",study,".txt",collapse=''), row.names = T, sep = '\t')
  write.table(para_MDS, paste0("/SAN/Plasmo_compare/SRAdb/Output/",study,"/parasite_MDS_",study,".txt",collapse=''), row.names = T, sep = '\t')
  
  # plot run index vs number parasite genes expressed
  
  library(reshape2)
  library(ggplot2)
  
  lessthan10e6 <- subset(host_parasite_genecount, Parasite_count < 10^6)
  require(dplyr)
  lessthan10e7 <- subset(host_parasite_genecount, between(Parasite_count, 10^6, 10^7)) # use dplyr
  morethan10e7 <- subset(host_parasite_genecount, Parasite_count > 10^7)
  
  #reshape_lessthan10e6 <- melt(lessthan10e6, id = c("Run", "Host_count", "Parasite_count"))
  #reshape_lessthan10e7 <- melt(lessthan10e7, id = c("Run", "Host_count", "Parasite_count"))
  #reshape_morethan10e7 <- melt(morethan10e7, id = c("Run", "Host_count", "Parasite_count"))
  
  
  #run_ind_vs_para_gene <- ggplot(reshaped_host_parasite_genecount, aes(x = value, fill = variable, color = variable)) + geom_density(alpha=.2)  + xlim(0,max(reshaped_host_parasite_genecount[,5])) + xlab(paste0("Number of genes for ",study, collapse = '')) + ylab("Number of runs with expression of corresponding number of genes") + theme(axis.text.y = element_blank()) + ggtitle("Number of host and parasite genes expresed") + theme_bw()
  # lessthan10e6_plot <- ggplot(lessthan10e6, aes(x = ParaGenesGreaterThanFive, y =..count.., fill = variable, color = variable)) + geom_histogram(stat="bin", alpha = 0.4) + geom_density(alpha=.2) + xlim(0,max(lessthan10e6[,4])) + xlab(paste0("Number of parasite genes with counts > 5 for ",study, collapse = '')) + ylab("Number of runs") + theme(axis.text.y = element_blank()) + ggtitle("Number of runs expressing genes with count > 5\n (total parasite gene count < 10e6)") + theme_bw()

  
  # plot for parasite gene count
  if(nrow(lessthan10e6) > 1)
    p.lessthan10e6_plot <- ggplot(lessthan10e6, aes(x = ParaGenesGreaterThanFive))  + geom_histogram(fill = "antiquewhite", color = "antiquewhite4", binwidth = 1) + scale_y_continuous(breaks = seq(0, max(lessthan10e6[,4]),1), labels = seq(0, max(lessthan10e6[,4]),1)) + xlab(paste0("Number of parasite genes with counts > 5 for ",study, collapse = '')) + ylab("Number of runs") + theme(axis.text.y = element_blank()) + ggtitle("Number of runs expressing genes with count > 5\n (total parasite gene count < 10e6)") + theme_bw()
  
  
  if(nrow(lessthan10e7) > 1)
    p.lessthan10e7_plot <- ggplot(lessthan10e7, aes(x = ParaGenesGreaterThanFive))  + geom_histogram(fill = "antiquewhite", color = "antiquewhite4", binwidth = 1) + scale_y_continuous(breaks = seq(0, max(lessthan10e7[,4]),1), labels = seq(0, max(lessthan10e7[,4]),1)) + xlab(paste0("Number of parsite genes with counts > 5 for ",study, collapse = '')) + ylab("Number of runs") + theme(axis.text.y = element_blank()) + ggtitle("Number of runs expressing genes with count > 5\n (10e6 < total parasite gene count < 10e7)") + theme_bw()
  
  
  if(nrow(morethan10e7) > 1)
    p.morethan10e7_plot <- ggplot(morethan10e7, aes(x = ParaGenesGreaterThanFive))  + geom_histogram(fill = "antiquewhite", color = "antiquewhite4", binwidth = 1) + scale_y_continuous(breaks = seq(0, max(morethan10e7[,4]),1), labels = seq(0, max(morethan10e7[,4]),1)) + xlab(paste0("Number of parasite genes with counts > 5 for ",study, collapse = '')) + ylab("Number of runs") + theme(axis.text.y = element_blank()) + ggtitle("Number of runs expressing genes with count > 5\n (total parasite gene count > 10e7)") + theme_bw()
  
  p.plots <- list(if (exists("p.lessthan10e6_plot")) p.lessthan10e6_plot else NA, if (exists("p.lessthan10e7_plot")) p.lessthan10e7_plot else NA, if (exists("p.morethan10e7_plot")) p.morethan10e7_plot else NA)
  p.plots.noNA_study <- p.plots[!is.na(p.plots)]
  
  require(gridExtra)
  #setwd("/SAN/Plasmo_comapre/SRAdb/")
  pdf(paste0("/SAN/Plasmo_compare/SRAdb/Output/", study,"/CategorisedNumberOfGenesForParasite_", study, ".pdf", collapse = ''), width = 20, height = 15)
  print(do.call(grid.arrange, p.plots.noNA_study))
  #ggsave(filename = "CategorisedNumberOfGenesAllHP.pdf", path = "/SAN/Plasmo_compare/SRAdb/Output/", plot = multiplePlotsPDF, width = 60, height = 20, units = "cm")
  dev.off()
  
  # plot for host gene count
  if(nrow(lessthan10e6) > 1)
    h.lessthan10e6_plot <- ggplot(lessthan10e6, aes(x = HostGenesGreaterThanFive))  + geom_histogram(fill = "aliceblue", color = "cadetblue", binwidth = 100) + scale_y_continuous(breaks = seq(0, max(lessthan10e6[,4]),1), labels = seq(0, max(lessthan10e6[,4]),1)) + xlab(paste0("Number of host genes with counts > 5 for ",study, collapse = '')) + ylab("Number of runs") + theme(axis.text.y = element_blank()) + ggtitle("Number of runs expressing genes with count > 5\n (total parasite gene count < 10e6)") + theme_bw()
  
  
  if(nrow(lessthan10e7) > 1)
    h.lessthan10e7_plot <- ggplot(lessthan10e7, aes(x = HostGenesGreaterThanFive))  + geom_histogram(fill = "aliceblue", color = "cadetblue", binwidth = 100) + scale_y_continuous(breaks = seq(0, max(lessthan10e7[,4]),1), labels = seq(0, max(lessthan10e7[,4]),1)) + xlab(paste0("Number of host genes with counts > 5 for ",study, collapse = '')) + ylab("Number of runs") + theme(axis.text.y = element_blank()) + ggtitle("Number of runs expressing genes with count > 5\n (10e6 < total parasite gene count < 10e7)") + theme_bw()
  
  
  if(nrow(morethan10e7) > 1)
    h.morethan10e7_plot <- ggplot(morethan10e7, aes(x = HostGenesGreaterThanFive))  + geom_histogram(fill = "aliceblue", color = "cadetblue", binwidth = 100) + scale_y_continuous(breaks = seq(0, max(morethan10e7[,4]),1), labels = seq(0, max(morethan10e7[,4]),1)) + xlab(paste0("Number of host genes with counts > 5 for ",study, collapse = '')) + ylab("Number of runs") + theme(axis.text.y = element_blank()) + ggtitle("Number of runs expressing genes with count > 5\n (total parasite gene count > 10e7)") + theme_bw()
  
  h.plots <- list(if (exists("h.lessthan10e6_plot")) h.lessthan10e6_plot else NA, if (exists("h.lessthan10e7_plot")) h.lessthan10e7_plot else NA, if (exists("h.morethan10e7_plot")) h.morethan10e7_plot else NA)
  h.plots.noNA_study <- h.plots[!is.na(h.plots)]
  
  require(gridExtra)
  #setwd("/SAN/Plasmo_comapre/SRAdb/")
  pdf(paste0("/SAN/Plasmo_compare/SRAdb/Output/", study,"/CategorisedNumberOfGenesForHost_", study, ".pdf", collapse = ''), width = 20, height = 15)
  print(do.call(grid.arrange, h.plots.noNA_study))
  #ggsave(filename = "CategorisedNumberOfGenesAllHP.pdf", path = "/SAN/Plasmo_compare/SRAdb/Output/", plot = multiplePlotsPDF, width = 60, height = 20, units = "cm")
  dev.off()
  
}


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

### get a lessthan10e6 etc plot for all studies

all_hp_genecount <- data.frame()
a = 1

for(i in 1:length(studyIDs))
{
  study <- studyIDs[i]
  study_hp_genecount <- read.csv2(paste0("/SAN/Plasmo_compare/SRAdb/Output/",study,"/host_parasite_genecount_",study,".txt",collapse=''), header = T, sep = '\t')
  n <- nrow(study_hp_genecount)
  
  #all_hp_genecount[a:(a+n-1),ncol(study_hp_genecount)] <- study_hp_genecount
  all_hp_genecount <- rbind(all_hp_genecount, study_hp_genecount)
  a = a + n
}
colnames(all_hp_genecount) <- c("Run", "Host_count", "Parasite_count", "ParaGenesGreaterThanFive", "HostGenesGreaterThanFive")
write.table(all_hp_genecount, "/SAN/Plasmo_compare/SRAdb/Output/all_hp_genecount_.txt", row.names = F, sep = '\t')

require(reshape2)
lessthan10e6 <- subset(all_hp_genecount, Parasite_count < 10^6)
reshape_lessthan10e6 <- melt(lessthan10e6, id = c("Run", "Host_count", "Parasite_count"))

#run_ind_vs_para_gene <- ggplot(reshaped_host_parasite_genecount, aes(x = value, fill = variable, color = variable)) + geom_density(alpha=.2)  + xlim(0,max(reshaped_host_parasite_genecount[,5])) + xlab(paste0("Number of genes for ",study, collapse = '')) + ylab("Number of runs with expression of corresponding number of genes") + theme(axis.text.y = element_blank()) + ggtitle("Number of host and parasite genes expresed") + theme_bw()
if(nrow(reshape_lessthan10e6) > 1)
  lessthan10e6_plot <- ggplot(reshape_lessthan10e6, aes(x = value, y =..count.., fill = variable, color = variable)) + geom_histogram(stat="bin",alpha = 0.4) + geom_density(alpha=.2) + xlim(0,max(reshape_lessthan10e6[,5])) + xlab("Number of genes with count > 5 for all studies") + ylab("Number of runs") + theme(axis.text.y = element_blank()) + ggtitle("Number of runs expressing genes with count > 5\n (total parasite gene count < 10e6)") + theme_bw()

require(dplyr)
lessthan10e7 <- subset(all_hp_genecount, between(Parasite_count, 10^6, 10^7)) # use dplyr
reshape_lessthan10e7 <- melt(lessthan10e7, id = c("Run", "Host_count", "Parasite_count"))

if(nrow(reshape_lessthan10e7) > 1)
  lessthan10e7_plot <- ggplot(reshape_lessthan10e7, aes(x = value, y =..count.., fill = variable, color = variable)) + geom_histogram(stat="bin",alpha = 0.4) + geom_density(alpha=.2) + xlim(0,max(reshape_lessthan10e7[,5])) + xlab("Number of genes with count > 5 for all studies") + ylab("Number of runs") + theme(axis.text.y = element_blank()) + ggtitle("Number of runs expressing genes with count > 5\n (10e6 < total parasite gene count < 10e7)") + theme_bw()

morethan10e7 <- subset(all_hp_genecount, Parasite_count > 10^7)
reshape_morethan10e7 <- melt(morethan10e7, id = c("Run", "Host_count", "Parasite_count"))

if(nrow(reshape_morethan10e7) > 1)
  morethan10e7_plot <- ggplot(reshape_morethan10e7, aes(x = value, y =..count.., fill = variable, color = variable))  + geom_histogram(stat="bin",alpha = 0.4) + geom_density(alpha=.2) + xlim(0,max(reshape_morethan10e7[,5])) + xlab("Number of genes with count > 5 for all studies") + ylab("Number of runs") + theme(axis.text.y = element_blank()) + ggtitle("Number of runs expressing genes with count > 5\n (total parasite gene count > 10e7)") + theme_bw()

plots <- list(if (exists("lessthan10e6_plot")) lessthan10e6_plot else NA, if (exists("lessthan10e7_plot")) lessthan10e7_plot else NA, if (exists("morethan10e7_plot")) morethan10e7_plot else NA)
plots.noNA <- plots[!is.na(plots)]

require(gridExtra)
#setwd("/SAN/Plasmo_comapre/SRAdb/")
pdf("/SAN/Plasmo_compare/SRAdb/CategorisedNumberOfGenesAllHP.pdf", width = 20, height = 15)
print(plot(do.call(grid.arrange, plots.noNA)))
#ggsave(filename = "CategorisedNumberOfGenesAllHP.pdf", path = "/SAN/Plasmo_compare/SRAdb/Output/", plot = multiplePlotsPDF, width = 60, height = 20, units = "cm")
dev.off()

## tSNE -> much slower

# colors = rainbow(length(unique(iris$Species)))
# names(colors) = unique(iris$Species)
# ecb = function(x,y){ plot(x,t='n'); text(x,labels=iris$Species, col=colors[iris$Species]) }
# tsne_iris = tsne(iris[,1:4], epoch_callback = ecb, perplexity=50)
# # compare to PCA
# dev.new()
# pca_iris = princomp(iris[,1:4])$scores[,1:2]
# plot(pca_iris, t='n')
# text(pca_iris, labels=iris$Species,col=colors[iris$Species])

# colors = rainbow(length(unique(rownames(study))))
# ecb = function(x,y){ plot(x,t='n'); text(x,labels=rownames(study), col=colors[rownames(study)]) }
# tsne_study = tsne(d, epoch_callback = ecb)

## remove genes that are not varying at all

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
