load("~/Documents/Data/overall_bipartite_with_6_datasets.RData")
Barseq20200228 <- read.csv("~/Downloads/Barseq20200228.csv", stringsAsFactors=FALSE)
parasite_orthogroups <- read.delim("~/Downloads/parasite_orthogroups.txt", stringsAsFactors=FALSE)

para_genes <- sapply(overall[,1], function(x) strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\1 \\2', x), ' ')[[1]][2])
para_df <- as.data.frame(para_genes)
para_ <- table(para_df)
para <- parasite_orthogroups[,c(1,4)]
h <- merge(para, Barseq20200228, by.x = "Pb_g", by.y = "gene")
i <- merge(h, para_, by.x = "Orthogroup", by.y = "para_df")

phenotype <-c("Essential","Dispensable", "Fast", "Insufficient data", "Slow")
color.codes<-as.character(c("#3399FF", "#FF00F0", "#FF0F00", 
                            "#FFF000", "#5efa0b00"))
color.names<-c("blue", "purple", "red", "yellow", "green")
df2=data.frame(phenotype, color.codes, color.names); df2

df <-left_join(i,df2, by = "phenotype"); df 

pdf("PlasmoGEM.pdf")
ggplot(data=df, aes(Freq, Confidence, colour = phenotype))+ 
  geom_point(alpha = 0.4) +
  scale_colour_manual(values=setNames(color.codes, phenotype))+
  ggtitle("Parasite genes degree centrality vs rel growth rate (PlasmoGEM)")
dev.off()


# for gene pairs with abs(cor) >= 0.9
overall0.9 <- overall[which(abs(overall$cor) >= 0.9),]
para_genes <- sapply(overall0.9[,1], function(x) strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\1 \\2', x), ' ')[[1]][2])

para_df <- as.data.frame(para_genes)
para_ <- table(para_df)
para <- parasite_orthogroups[,c(1,4)]
h <- merge(para, Barseq20200228, by.x = "Pb_g", by.y = "gene")
i <- merge(h, para_, by.x = "Orthogroup", by.y = "para_df")

phenotype <-c("Essential","Dispensable", "Fast", "Insufficient data", "Slow")
color.codes<-as.character(c("#3399FF", "#FF00F0", "#FF0F00", 
                            "#FFF000", "#5efa0b00"))
color.names<-c("blue", "purple", "red", "yellow", "green")
df2=data.frame(phenotype, color.codes, color.names); df2

df <-left_join(i,df2, by = "phenotype"); df 

pdf("PlasmoGEM_0.9.pdf")
ggplot(data=df, aes(Freq, Confidence, colour = phenotype))+ 
  geom_point(alpha = 0.4) +
  scale_colour_manual(values=setNames(color.codes, phenotype)) +
  ggtitle("Parasite genes degree centrality vs confidence (PlasmoGEM)")
dev.off()