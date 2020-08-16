positive_experiments <- read.table("Input/positive_experiments.txt", header = F, sep = '\t')
studyIDs <- c("ERP109367") # leave out new monkey study and partial studies
#studyIDs <- c("ERP020067", "SRP018945", "SRP066796", "ERP106451")

# Second, we want the run IDs for each of the studies we read from the above file. So we go into the folder of each study and read the run IDs from runs_<study>.txt

getAllHPexpression <- function(study)
{ 
  study_df <- data.frame()
  
  # the following statement gives us both run ID and the study ID
  runIDs <- read.table(paste0("SRAdb/Output/",study,"/runs_",study,".Pvn.txt",collapse=''), header = F, sep = ',')
  
  # number of runs
  number_of_runs <- nrow(runIDs)
  
  # get hp_percent_<study>.txt
  hp_percent_study <- read.table(paste0("SRAdb/Output/",study,"/hp_percent_",study,".Pvn.txt",collapse=''), header = T, sep = '\t')
  
  for(i in 1:number_of_runs)
  {
    # get runID
    runID <- as.character(runIDs[i,1])
    
    # get the host and parasite expression percentages for the runID from hp_percent_study
    host_percent <- hp_percent_study[grep(runID, hp_percent_study[,1]), 2]
    parasite_percent <- hp_percent_study[grep(runID, hp_percent_study[,1]), 3]
    
    host <- "mouse"
    parasite <- "Pvinckei"
    
    count_file <- read.table(paste0("SRAdb/Output/",study,"/countWithGFF3_",runID,".txt",collapse=''), header = T, sep = '\t')
    if(file.exists(paste0("SRAdb/Output/",study,"/",runID,"_",study,".final.out",collapse='')))
    {
      inputReadLine <- readLines(paste0("SRAdb/Output/",study,"/",runID,"_",study,".final.out",collapse=''), n=6)[6]
      mapPercentLine <- readLines(paste0("SRAdb/Output/",study,"/",runID,"_",study,".final.out",collapse=''), n=10)[10]
      
      multipleReadsLine <- readLines(paste0("SRAdb/Output/",study,"/",runID,"_",study,".final.out",collapse=''), n=26)[c(24,26)]
      
      unmappedReadsLine <- readLines(paste0("SRAdb/Output/",study,"/",runID,"_",study,".final.out",collapse=''), n=31)[c(29,30,31)]
      
      inputReads <- as.numeric(strsplit(inputReadLine, split="\t")[[1]][2])
      map_percent <- as.numeric(strsplit(strsplit(mapPercentLine, split="%")[[1]][2], split="\t")[[1]][2])
      multipleReads <- sum(sapply(multipleReadsLine, function(x) as.numeric(strsplit(x, split="\t")[[1]][2])))
      unmappedReads <- sum(sapply(unmappedReadsLine, function(x) (as.numeric(strsplit(strsplit(x, split="%")[[1]][2], split="\t")[[1]][2])/100) * inputReads))
    }
    if(!file.exists(paste0("SRAdb/Output/",study,"/",runID,"_",study,".final.out",collapse='')))
    {
      inputReads <- "NA"
      map_percent <- "NA"
      multipleReads <- "NA"
      unmappedReads <- "NA"
    }
    
    parasite_rows <- as.numeric(grep("P+", count_file[,1]))
    host_count <- count_file[-c(parasite_rows),]
    parasite_count <- count_file[parasite_rows,]
    
    host_sum <- sum(host_count[,6])
    parasite_sum <- sum(parasite_count[,6])
    
    actual_host <- hp_percent_study[grep(runID, hp_percent_study[,1]),5]
    actual_para <- hp_percent_study[grep(runID, hp_percent_study[,1]),6]
    
    #mergeParasiteGeneNamesAndNumberResult <- mergeParasiteGeneNamesAndNumber(parasite_count, host_count, study)
    #para_genecount <- mergeParasiteGeneNamesAndNumberResult[[1]] # dataframe of host part
    #host_genecount <- mergeParasiteGeneNamesAndNumberResult[[2]] # dataframe of parasite part
    
    #parasite_uniquegenes <- nrow(subset(para_genecount, count > 0))
    #host_uniquegenes <- nrow(subset(host_genecount, count > 0))
    
    
    study_df[i,1] <- runID
    study_df[i,2] <- study
    study_df[i,3] <- host
    study_df[i,4] <- host_percent
    study_df[i,5] <- parasite
    study_df[i,6] <- parasite_percent
    study_df[i,7] <- sum(count_file[,6]) # sum(counts)
    study_df[i,8] <- inputReads # number of input reads
    study_df[i,9] <- map_percent # map_percent
    study_df[i,10] <- paste0(host,parasite, collapse='')
    study_df[i,11] <- host_sum
    study_df[i,12] <- parasite_sum
    study_df[i,13] <- multipleReads
    study_df[i,14] <- unmappedReads
    study_df[i,15] <- actual_para # actual proportion of parasites, taking into account unique map%
    study_df[i,16] <- actual_host
    # study_df[i,17] <- parasite_uniquegenes
    # study_df[i,18] <- host_uniquegenes
    
    
  }
  return (study_df)
}
studyIDs <- c("ERP109367")
allHPexpmPchv1 <- data.frame()
for(j in 1:length(studyIDs))
{
  print(studyIDs[j])
  allHPexpmPchv1 <- rbind(allHPexpmPchv1, getAllHPexpression(studyIDs[j]))
}

colnames(allHPexpmPchv1) <- c("RunID", "Study", "Host", "Host_percent", "Parasite", "Parasite_percent", "CountSum", "InputReads",  "MapPercent" ,"HostParasite", "HostSum", "ParasiteSum", "MultipleMapReads", "UnmappedReads", "ActualParasite", "ActualHost")
write.table(allHPexpmPchv, "allHPexpmPchv.txt", sep = '\t', row.names = F)
save(allHPexpmPchv, file = "allHPexpmPchv.RData")
# Plot

require(ggplot2)

# col=c("mousePyoelii"="red", "monkeyPcynomolgi"="blue", "humanPfalciparum"="green")
# shape=c("mousePyoelii"=0, "monkeyPcynomolgi"=1)

allHPexpression_ggplot <- ggplot(allHPexp, aes(as.numeric(ActualParasite), ActualHost, col=HostParasite, shape=HostParasite))

allHPexpression_gg <- allHPexpression_ggplot + geom_point(alpha = 0.65) + geom_hline(yintercept=50, color="brown") + geom_vline(xintercept=50, color="brown") + theme_bw() +   xlim(0,100) + ylim(0,100) + ggtitle("Host expression vs parasite expression (uniquely mapped)") +   xlab("Parasite gene expression percentage") + ylab("Host gene expression percentage") +theme(axis.text=element_text(size=12), axis.title=element_text(size=14), plot.title = element_text(size=18)) + labs(color='Host-parasite\npairs') + labs(shape='Host-parasite\npairs') + scale_shape_manual(values = 0:7) + scale_color_manual(values = col)
allHPexpression_gg
ggsave("UniquelyMappedReads.png")

# For this we require two columns of allHPexpression - Study and Parasite_percent

# parasiteProportion_ggplot <- ggplot(allHPexp, aes(x = as.factor(Study), y = Parasite_percent, col=HostParasite, fill=HostParasite, shape=HostParasite))
# 
# parasiteProportion_gg <- parasiteProportion_ggplot + geom_violin(scale="width", alpha=0.3) +  ylim(0,100)+ geom_jitter(height = 0, width = 0.1, size = 2.5) + ggtitle("Parasite distribution in studies") + xlab("SRA study ID") + ylab("Parasite gene expression percentage") + theme(axis.text=element_text(size=14), axis.title=element_text(size=16), plot.title = element_text(size=20), axis.text.x = element_text(angle = 90, vjust = 1, size = 7)) + theme_classic() + labs(color='Host-parasite\npairs') + labs(shape='Host-parasite\npairs') +labs(fill='Host-parasite\npairs') +theme(strip.background = element_rect(fill=c("aliceblue"), colour = "lightblue"), strip.text.x = element_text(size = 12, colour = "black")) +facet_wrap(~Host, scales = "free") + coord_flip() + scale_shape_manual(values = 0:7)
# 
# parasiteProportion_gg
# ggsave("ParasiteProportion.png")

cols <- data.frame(HostParasite = c("humanPfalciparum", "humanPvivax", "humanPberghei", "mousePberghei", "mousePchabaudi", "mousePyoelii", "monkeyPcoatneyi", "monkeyPcynomolgi"),
                   colours = c("chocolate1", "hotpink2", "brown3", "chartreuse3", "darkgreen", "lightgreen", "dodgerblue3", "navy"))

allHPexp <- merge(allHPexp, cols, by = "HostParasite")

require(ggplot2)
# col <- as.character(allHPexp$colours)
# names(col) <- as.character(allHPexp$HostParasite)
# 
# parasiteProportion_ggplot <- ggplot(allHPexp, aes(x = as.factor(Study),y = Parasite_percent, colour = HostParasite, fill = HostParasite, shape=HostParasite))
# 
# parasiteProportion_gg <- parasiteProportion_ggplot + geom_violin(scale="width", alpha=0.3) +  ylim(0,100)+ geom_jitter(height = 0, width = 0.1, size = 2.5) + scale_color_manual(values = col) + scale_fill_manual(values = col) + ggtitle("Parasite distribution in studies") + xlab("SRA study ID") + ylab("Parasite gene expression percentage") + theme(axis.text=element_text(size=14), axis.title=element_text(size=16), plot.title = element_text(size=20), axis.text.x = element_text(angle = 90, vjust = 1, size = 7)) + theme_classic() + labs(color='Host-parasite\npairs') + labs(shape='Host-parasite\npairs') +labs(fill='Host-parasite\npairs') +theme(strip.background = element_rect(fill=c("aliceblue"), colour = "lightblue"), strip.text.x = element_text(size = 12, colour = "black")) +facet_wrap(~Host, scales = "free") + coord_flip() + scale_shape_manual(values = 0:7)
# 
# parasiteProportion_gg
# ggsave("ParasiteProportion.png")

## 
col <- as.character(allHPexp$colours)
names(col) <- as.character(allHPexp$HostParasite)
require(ggplot2)
parasiteProportion_ggplot <- ggplot(allHPexp, aes(x = as.factor(Study),y = as.numeric(as.character(Parasite_percent)), colour = HostParasite, fill = HostParasite))

parasiteProportion_gg <- parasiteProportion_ggplot  + ylim(0,100) + geom_jitter(height = 0, width = 0.1, size = 1.5, alpha = 0.6) + scale_color_manual(values = col) + scale_fill_manual(values = col) + ggtitle("Parasite distribution in studies") + xlab("SRA study ID") + ylab("Parasite gene expression percentage") + theme(axis.text=element_text(size=14), axis.title=element_text(size=16), plot.title = element_text(size=20), axis.text.x = element_text(angle = 90, vjust = 1, size = 7)) + theme_classic() + labs(color='Host-parasite\npairs') +labs(fill='Host-parasite\npairs') +theme(strip.background = element_rect(fill=c("aliceblue"), colour = "lightblue"), strip.text.x = element_text(size = 12, colour = "black")) +facet_wrap(~Host, scales = "free") + coord_flip()

parasiteProportion_gg
ggsave("ParasiteProportion.png")
##

parasiteProportionCluster_ggplot <- ggplot(allHPexp, aes(x=Parasite_percent, y=log10(CountSum), col=HostParasite, shape=HostParasite))

parasiteProportionCluster_gg <- parasiteProportionCluster_ggplot + geom_point() + xlim(0,100) +
  labs(color='Host-parasite\npairs') + labs(shape='Host-parasite\npairs') + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), plot.title = element_text(size=18)) + theme_bw() + scale_shape_manual(values = 0:7)

parasiteProportionCluster_gg

require(ggplot2)
hp_parasiteProportion_ggplot <- ggplot(allHPexp, aes(log10(InputReads), ActualParasite, shape=HostParasite, color=Study))

hp_parasiteProportion_gg <- hp_parasiteProportion_ggplot + geom_point(alpha=0.7) + ylim(0,100) +
  ylab("Parasite proportion uniquely mapped") + xlab("Number of input reads in a sample (run)") + 
  ggtitle("Parasite proportion relative to number of input reads for a host-parasite pair") + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), plot.title = element_text(size=18)) + theme_light() + scale_shape_manual(values = 0:7)

hp_parasiteProportion_gg

hp_parasiteProportion_parasite_ggplot <- ggplot(allHPexp, aes(log10(InputReads), ActualParasite, shape=Parasite, color=Parasite))

hp_parasiteProportion_parasite_gg <- hp_parasiteProportion_parasite_ggplot + 
  geom_point(alpha=0.7) + ylim(0,100) + ylab("Parasite proportion uniquely mapped") + xlab("Number of input reads in a sample (run)") + 
  ggtitle("Parasite proportion relative to number of input reads for each parasite") + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), plot.title = element_text(size=18)) + theme_minimal() + scale_shape_manual(values = 0:7)

hp_parasiteProportion_parasite_gg

Summary_table <- data.frame()

for(i in 1:length(studyIDs))
{
  total_runs <- nrow(read.table(paste0("/SAN/Plasmo_compare/SRAdb/Output/",studyIDs[i],"/runs_",studyIDs[i],".txt",collapse=''), header = F, sep = ','))
  table_study <- t(read.table(paste0("/SAN/Plasmo_compare/SRAdb/Output/",studyIDs[i],"/table_",studyIDs[i],".txt",collapse=''), header = T, sep = '\t'))
  
  Summary_table[i,1] <- studyIDs[i]
  Summary_table[i,2] <- total_runs
  Summary_table[i,3] <- (as.numeric(table_study[3])/total_runs)*100
  Summary_table[i,4] <- (as.numeric(table_study[4])/total_runs)*100
  Summary_table[i,5] <- (as.numeric(table_study[5])/total_runs)*100
  Summary_table[i,6] <- (as.numeric(table_study[6])/total_runs)*100
}

colnames(Summary_table) <- c("Study", "TotalRunsInStudy", "PercentageRunsIn70_30", "PercentageRunsIn90_10", "PercentageRunsIn99_1", "WorseThan99_1")

# for Latex:
# require(stargazer) # Hlavac, Marek (2018). stargazer: Well-Formatted Regression and Summary Statistics Tables.
# table_gg <- stargazer(Summary_table, rownames = FALSE, summary = FALSE, align = TRUE)
require(ggpubr)
ggtexttable(Summary_table, rows = NULL, theme = ttheme("mBlue"))
# require(pander)
# table_gg <- knitr::kable(Summary_table)
# table_gg
write.table(Summary_table, "Summary_table.txt", sep='\t', row.names = F)


### suitability ###
suitability <- list()

number_of_df <- levels(as.factor(allExpression$HostParasite))

for(hp in number_of_df)
{
  # first get all runs for hp and remove rows with NA (for Pyoelii)
  hp.fulldf <- allExpression[allExpression$HostParasite==hp,]
  hp.fulldf <- hp.fulldf[,-c(13:14)]
  hp.fulldf <- na.omit(hp.fulldf)
  
  # select runs which have mapped above 65%
  hp.df <- hp.fulldf[as.numeric(hp.fulldf$MapPercent)>65,]
  
  # select runs with parasite reads greater than 10e6
  hp.df <- hp.df[as.numeric(hp.df$ParasiteSum)>1e6,]
  
  # select runs with host reads greater than 10e6
  hp.df <- hp.df[as.numeric(hp.df$HostSum)>1e6,]
  
  suitability[[hp]] <- hp.df
}

save(suitability, file = "suitability.RData")

Pfal<-allExpression[allExpression$Parasite=="Pfalciparum",]
Pber<-allExpression[allExpression$Parasite=="Pberghei",]
Pviv<-allExpression[allExpression$Parasite=="Pvivax",]
Pcha<-allExpression[allExpression$Parasite=="Pchabaudi",]
Pcyn<-allExpression[allExpression$Parasite=="Pcynomolgi",]
Pcoa<-allExpression[allExpression$Parasite=="Pcoatneyi",]
Pyoe<-allExpression[allExpression$Parasite=="Pyoelii",]
hPfal<-allExpression[allExpression$HostParasite=="humanPfalciparum",]
hPber<-allExpression[allExpression$HostParasite=="humanPberghei",]
hPviv<-allExpression[allExpression$HostParasite=="humanPvivax",]
mPcha<-allExpression[allExpression$HostParasite=="mousePchabaudi",]
mPyoe<-allExpression[allExpression$HostParasite=="mousePyoelii",]
mPber<-allExpression[allExpression$HostParasite=="mousePberghei",]
moPcyn<-allExpression[allExpression$HostParasite=="monkeyPcynomolgi",]
moPcoa<-allExpression[allExpression$HostParasite=="monkeyPcoatneyi",]

parasiteReads <- data.frame(Parasite = c("Pfalciparum", "Pberghei", "Pvivax", "Pchabaudi", "Pyoelii", "Pcynomolgi", "Pcoatneyi"), MappedReads10E9 = c(sum(Pfal$CountSum/10^9),sum(Pber$CountSum/10^9),sum(Pviv$CountSum/10^9),sum(Pcha$CountSum/10^9),sum(Pyoe$CountSum/10^9),sum(Pcyn$CountSum/10^9),sum(Pcoa$CountSum/10^9)))
hostparasiteReads <- data.frame(HostParasite = c("humanPfalciparum", "humanPberghei", "humanPvivax", "mousePchabaudi", "mousePberghei", "mousePyoelii", "monkeyPcynomolgi", "monkeyPcoatneyi"), MappedReads10E9 = c(sum(hPfal$CountSum/10^9),sum(hPber$CountSum/10^9),sum(hPviv$CountSum/10^9),sum(mPcha$CountSum/10^9),sum(mPber$CountSum/10^9),sum(mPyoe$CountSum/10^9),sum(moPcyn$CountSum/10^9),sum(moPcoa$CountSum/10^9)))
write.table(parasiteReads, "NumberOfReadsMappedOntoEachParasite.txt", sep = '\t', row.names = F)
write.table(hostparasiteReads, "NumberOfReadsMappedOntoEachHostParasite.txt", sep = '\t', row.names = F)

positive_experiments <- read.table("/SAN/Plasmo_compare/SRAdb/Input/positive_experiments.txt", header = F, sep = '\t')
require(rtracklayer)

for(s in 1:length(studyIDs))
{
  study <- studyIDs[s]
  print(study)
  
  para <- positive_experiments[grep(study,positive_experiments[,1]),3]
  host <- positive_experiments[grep(study,positive_experiments[,1]),2]
  # the following statement gives us both run ID and the study ID
  runIDs <- read.table(paste0("SRAdb/Output/",study,"/runs_",study,".txt",collapse=''), header = F, sep = ',')
  
  # number of runs
  number_of_runs <- nrow(runIDs)
  
  p.genes <- import(paste0("/SAN/Plasmo_compare/Genomes/annotation/",para,".gtf", collapse=''), format = "gtf")
  p.genes <- p.genes[p.genes$type%in%"exon"]
  #genes <- genes[which(genes[,"type"] == "exon"),]
  p.genes.df <- as.data.frame(p.genes)
  p.genes.df.gene_name <- p.genes.df[,c("seqnames", "start", "end", "width", "strand", "gene_id")]
  
  h.genes <- import(paste0("/SAN/Plasmo_compare/Genomes/annotation/",host,".gtf", collapse=''), format = "gtf")
  h.genes <- h.genes[h.genes$type%in%"exon"]
  #genes <- genes[which(genes[,"type"] == "exon"),]
  h.genes.df <- as.data.frame(h.genes)
  h.genes.df.gene_name <- h.genes.df[,c("seqnames", "start", "end", "width", "strand", "gene_id")]
  
  
  for(i in 1:number_of_runs)
  {
    print(i)
    # get runID
    runID <- as.character(runIDs[i,1])
    
    count_file <- read.table(paste0("/SAN/Plasmo_compare/SRAdb/Output/",study,"/countWithGFF3_",runID,".txt",collapse=''), header = T, sep = '\t')
    parasite_rows <- as.numeric(grep("P+", count_file[,1]))
    host_count <- count_file[-c(parasite_rows),]
    parasite_count <- count_file[parasite_rows,] 
    
    # parasite gene merge
    
    p.mergeStudy.genes.df.gene_name <- merge(parasite_count, p.genes.df.gene_name, by = c("seqnames", "start", "end", "width", "strand"))
    
    p.mergeStudy.genes.df.gene_name <- p.mergeStudy.genes.df.gene_name[,6:ncol(p.mergeStudy.genes.df.gene_name)]
    p.mergeStudy.genes.df.gene_name.combineGenes <- data.frame()
    p.mergeStudy.genes.df.gene_name.combineGenes <- aggregate(p.mergeStudy.genes.df.gene_name[,1] ~ gene_id, data = p.mergeStudy.genes.df.gene_name, sum)
    colnames(p.mergeStudy.genes.df.gene_name.combineGenes)[2] <- "count"
    
    # host genes merge
    
    
    h.mergeStudy.genes.df.gene_name <- merge(host_count, h.genes.df.gene_name, by = c("seqnames", "start", "end", "width", "strand"))
    
    h.mergeStudy.genes.df.gene_name <- h.mergeStudy.genes.df.gene_name[,6:ncol(h.mergeStudy.genes.df.gene_name)]
    h.mergeStudy.genes.df.gene_name.combineGenes <- data.frame()
    h.mergeStudy.genes.df.gene_name.combineGenes <- aggregate(h.mergeStudy.genes.df.gene_name[,1] ~ gene_id, data = h.mergeStudy.genes.df.gene_name, sum)
    colnames(h.mergeStudy.genes.df.gene_name.combineGenes)[2] <- "count"
    
    parasite_uniquegenes <- nrow(subset(p.mergeStudy.genes.df.gene_name.combineGenes, count > 0))
    host_uniquegenes <- nrow(subset(h.mergeStudy.genes.df.gene_name.combineGenes, count > 0))
    
    allHPexp[allHPexp$RunID==runID,17] <- parasite_uniquegenes
    allHPexp[allHPexp$RunID==runID,18] <- host_uniquegenes
  }
}

colnames(allHPexp)[17] <- "NumberOfParaGenes"
colnames(allHPexp)[18] <- "NumberOfHostGenes"

write.table(allHPexp, "allHPexp.txt", sep = '\t', row.names = F)
save(allHPexp, file = "allHPexp.RData")

#### parasite genes vs parasite reads ####
col <- as.character(allHPexp$colours)
names(col) <- as.character(allHPexp$HostParasite)

#### blood ##
blood <- allHPexp[allHPexp$Tissue=="blood",]
blood$ProtCodGenesHostPerMaxValue_b <- blood$NumberProtCodGenesHost / max(blood$NumberProtCodGenesHost)
blood$ProtCodGenesParaPerMaxValue_b <- blood$NumberProtCodGenesPara / max(blood$NumberProtCodGenesPara)

paragenes_ggplot <- ggplot(blood, aes(x = log10(ProteinCodPara), 
                                         y = ProtCodGenesParaPerMaxValue_b,
                                         colour = HostParasite, 
                                         fill = HostParasite))
paragenes_gg <- paragenes_ggplot + geom_point(alpha = 0.7) +
  geom_jitter() + 
  ggtitle("Number of reads mapping onto parasite genes") + 
  xlab("(Log10) Number of reads mapping onto parasite transcriptome") + 
  ylab("Number of genes in blood samples") + 
  theme_bw() + scale_color_manual(values = col) + scale_fill_manual(values = col) + 
  theme(axis.text=element_text(size=22), axis.title=element_text(size=22), 
        plot.title = element_text(size=24))
paragenes_gg
ggsave("paragenes_blood_30July2020.png", width = 35, height = 30, units = "cm")

#### host genes vs host reads ####
hostgenes_ggplot <- ggplot(blood, aes(x = log10(ProteinCodHost), 
                                         y = ProtCodGenesHostPerMaxValue_b,
                                         colour = HostParasite, 
                                         fill = HostParasite))
hostgenes_gg <- hostgenes_ggplot + geom_point(alpha = 0.7) + 
  ggtitle("Number of reads mapping onto host genes") +
  xlab("(Log10) Number of reads mapping onto host transcriptome") + 
  ylab("Number of genes represented in blood samples") + 
  theme_bw()  + scale_color_manual(values = col) + scale_fill_manual(values = col) +
  theme(axis.text=element_text(size=18), axis.title=element_text(size=22), 
        plot.title = element_text(size=24))
hostgenes_gg
ggsave("hostgenes_blood_30July2020.png", width = 30, height = 30, units = "cm")

#liver 
liver <- allHPexp[allHPexp$Tissue=="liver",]
liver$ProtCodGenesHostPerMaxValue_l <- liver$NumberProtCodGenesHost / max(liver$NumberProtCodGenesHost)
liver$ProtCodGenesParaPerMaxValue_l <- liver$NumberProtCodGenesPara / max(liver$NumberProtCodGenesPara)

paragenes_ggplot <- ggplot(liver, aes(x = log10(ProteinCodPara), 
                                         y =  ProtCodGenesParaPerMaxValue_l,
                                         colour = HostParasite, 
                                         fill = HostParasite))
paragenes_gg <- paragenes_ggplot + geom_point(alpha = 0.7) +
  geom_jitter() + 
  ggtitle("Number of reads mapping onto parasite genes") + 
  xlab("(Log10) Number of reads mapping onto parasite transcriptome") + 
  ylab("Number of genes in liver samples") + 
  theme_bw() + scale_color_manual(values = col) + scale_fill_manual(values = col) + 
  theme(axis.text=element_text(size=22), axis.title=element_text(size=22), 
        plot.title = element_text(size=24))
paragenes_gg
ggsave(plot = paragenes_gg, "paragenes_liver_30July2020.png", width = 35, height = 30, units = "cm")

#### host genes vs host reads ####
hostgenes_ggplot <- ggplot(liver, aes(x = log10(ProteinCodHost), 
                                         y = ProtCodGenesHostPerMaxValue_l,
                                         colour = HostParasite, 
                                         fill = HostParasite))
hostgenes_gg <- hostgenes_ggplot + geom_point(alpha = 0.7) + 
  ggtitle("Number of reads mapping onto host genes") +
  xlab("(Log10) Number of reads mapping onto host transcriptome") + 
  ylab("Number of genes in liver samples") + 
  theme_bw()  + scale_color_manual(values = col) + scale_fill_manual(values = col) +
  theme(axis.text=element_text(size=18), axis.title=element_text(size=22), 
        plot.title = element_text(size=24))
hostgenes_gg
ggsave(plot = hostgenes_gg, "hostgenes_liver_30July2020.png", width = 30, height = 30, units = "cm")


#### overall host and para gene vs read #####

paragenes_ggplot <- ggplot(allHPexp, aes(x = log10(ProteinCodPara), 
                                      y =  allHPexp$NumberProtCodGenesPara / max(allHPexp$NumberProtCodGenesPara),
                                      colour = HostParasite, 
                                      fill = HostParasite))
paragenes_gg <- paragenes_ggplot + geom_point(alpha = 0.7) +
  geom_jitter() + 
  ggtitle("Number of reads mapping onto parasite genes") + 
  xlab("(Log10) Number of reads mapping onto parasite transcriptome") + 
  ylab("Number of genes in a sample") + 
  theme_bw() + scale_color_manual(values = col) + scale_fill_manual(values = col) + 
  theme(axis.text=element_text(size=22), axis.title=element_text(size=22), 
        plot.title = element_text(size=24))
paragenes_gg
ggsave(plot = paragenes_gg, "paragenes_all_30July2020.png", width = 35, height = 30, units = "cm")

#### host genes vs host reads ####

hostgenes_ggplot <- ggplot(allHPexp, aes(x = log10(ProteinCodHost), 
                                      y = allHPexp$NumberProtCodGenesHost / max(allHPexp$NumberProtCodGenesHost),
                                      colour = HostParasite, 
                                      fill = HostParasite))
hostgenes_gg <- hostgenes_ggplot + geom_point(alpha = 0.7) + 
  ggtitle("Number of reads mapping onto host genes") +
  xlab("(Log10) Number of reads mapping onto host transcriptome") + 
  ylab("Number of genes in a sample") + 
  theme_bw()  + scale_color_manual(values = col) + scale_fill_manual(values = col) +
  theme(axis.text=element_text(size=18), axis.title=element_text(size=22), 
        plot.title = element_text(size=24))
hostgenes_gg
ggsave(plot = hostgenes_gg, "hostgenes_all_30July2020.png", width = 30, height = 30, units = "cm")


ggplot(allHPexp, aes(x = log10(NumberOfParaGenes), y = log10(NumberOfHostGenes))) + geom_point(alpha = 0.7) # Most runs have host genes > 10000 and para genes > 5000(!)

### suitability analysis ###
hp <- allHPexp[allHPexp$HostParasite=="monkeyPcynomolgi",]

firstthreshold <- hp[hp$NumberOfHostGenes>=3000 & hp$NumberOfParaGenes>= 100 & hp$HostSum >=100000 & hp$ParasiteSum >= 10000 & as.numeric(as.character(hp$MapPercent)) >= 65,]
nrow(firstthreshold)
(studiesfirstthreshold <- unique(as.character(firstthreshold$Study)))

secondthreshold <- hp[hp$NumberOfHostGenes>=3000 & hp$NumberOfParaGenes>= 1000 & hp$HostSum >=100000 & hp$ParasiteSum >= 10000 & as.numeric(as.character(hp$MapPercent)) >= 65,]
nrow(secondthreshold)
(studiessecondthreshold <- unique(as.character(secondthreshold$Study)))

thirdthreshold <- hp[hp$NumberOfHostGenes>=10000 & hp$NumberOfParaGenes>= 3000 & hp$HostSum >=100000 & hp$ParasiteSum >= 10000 & as.numeric(as.character(hp$MapPercent)) >= 65,]
nrow(thirdthreshold)
(studiesthirdthreshold <- unique(as.character(thirdthreshold$Study)))

# find chosen runs per study based on prc and hrc
positive_experiments <- read.table("Input/positive_experiments.txt", header = F, sep = '\t')
studyIDs <- positive_experiments[,1] # leave out new monkey study and partial studies

for(i in 1:length(studyIDs))
{
  print(i)
  study <- allHPexp[allHPexp$Study==studyIDs[i],]
  
  threshold <- study[study$NumberOfHostGenes>=3000 & study$NumberOfParaGenes>= 500 & study$HostSum >=1000000 & study$ParasiteSum >= 100000,]
  
  if(nrow(threshold) > 1)
  {
    runsthreshold <- unique(as.character(threshold$RunID))
    
    useful_runs <- write.table(runsthreshold, paste0("Output/", studyIDs[i], "/", studyIDs[i], "_runs_above_threshold.txt", collapse = ''))
    df <- read.table(paste0("Output/", studyIDs[i], "/", studyIDs[i], "_filteredDisp.txt", collapse = ''), header = T, sep = '\t')
    
    runs_threshold <-  sapply(runsthreshold, function(x) paste(x, studyIDs[i], sep = '_'))
    df_useful_runs <- df[,runs_threshold]
  }
  else
    df_useful_runs <- "none"
  
  write.table(df_useful_runs, paste0("Output/", studyIDs[i], "/", studyIDs[i], "_threshold.txt", collapse = ''), sep = '\t')
}

## countSum vs host gene count and countSum vs parasite gene count #
allHPexp_mapthreshold <- allHPexp[allHPexp$MapPercent>=65,]
col <- as.character(allHPexp_mapthreshold$colours)
names(col) <- as.character(allHPexp_mapthreshold$HostParasite)

ipparagenes_ggplot <- ggplot(allHPexp, aes(x = log10(CountSum), y = log10(NumberOfParaGenes), colour = HostParasite, fill = HostParasite))
ipparagenes_gg <- ipparagenes_ggplot + geom_point(alpha = 0.7) + ggtitle("Number of parasite genes vs total input reads in each run") + xlab("(Log10) Number of input reads in the sequenced sample") + ylab("(Log10) Number of parasite genes in the run") + theme_bw() + scale_color_manual(values = col) + scale_fill_manual(values = col) + theme(axis.text=element_text(size=14), axis.title=element_text(size=16), plot.title = element_text(size=20))
ipparagenes_gg
ggsave("inputreads_paragenes.png")


iphostgenes_ggplot <- ggplot(allHPexp_mapthreshold, aes(x = log10(CountSum), y = log10(NumberOfHostGenes), colour = HostParasite, fill = HostParasite))
iphostgenes_gg <- iphostgenes_ggplot + geom_point(alpha = 0.7) + ggtitle("Number of host genes vs total input reads in each run") + xlab("(Log10) Number of input reads in the sequenced sample") + ylab("(Log10) Number of host genes in the run") + theme_bw()  + scale_color_manual(values = col) + scale_fill_manual(values = col) + theme(axis.text=element_text(size=14), axis.title=element_text(size=16), plot.title = element_text(size=20))
iphostgenes_gg
ggsave("inputreads_hostgenes.png")


