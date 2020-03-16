library(dplyr)

load("~/Documents/Data/overall_bipartite_with_6_datasets.RData")
Barseq20200228 <- read.csv("~/Downloads/Barseq20200305.csv", stringsAsFactors=FALSE)
parasite_orthogroups <- read.delim("~/Documents/Data/parasite_orthogroups.txt", stringsAsFactors=FALSE)

overall_0.9 <- overall[abs(overall$cor) >= 0.95,]
para_genes <- sapply(overall_0.9[,1], function(x) strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\1 \\2', x), ' ')[[1]][2])
para_df <- as.data.frame(para_genes)
colnames(para_df) <- "Orthogroup"
#para_ <- table(para_df)
para <- parasite_orthogroups[,c(1,4)]
colnames(para)[2] <- "gene_id"
#h <- merge(para, Barseq20200228, by.x = "Pb_g", by.y = "gene")
i.p <- inner_join(para_df, para, by.x = "para_genes", by.y = "Orthogroup")

# write.table(i.p, "parasite0.9_details_plasmogem.txt", sep = '\t', row.names = F)
host_orthogroups <- read.delim("~/Documents/Data/host_orthogroups.txt", stringsAsFactors=FALSE)
host_genes <- sapply(overall_0.9[,1], function(x) strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\1 \\2', x), ' ')[[1]][1])
host_df <- as.data.frame(host_genes)
colnames(host_df) <- "Orthogroup"
host <- host_orthogroups[,c(1,3)]
colnames(host)[2] <- "gene_id"
i.h <- inner_join(host_df, host)

#BiocManager::install("rtracklayer")
gtf <- import("/home/parnika/Documents/Data/Kai/mousePberghei.gtf", format = "gtf")
gtf.df <- as.data.frame(gtf)
gtf.df <- gtf.df[which(gtf.df$type=="exon"),]
gtf.df.genename <- gtf.df[,c("gene_id", "gene_name")]
gtf.host <- gtf.df.genename[grep(pattern = "ENS", gtf.df.genename$gene_id),]
n <- gtf.host[!duplicated(gtf.host), ]

gtf.para <- gtf.df.genename[grep(pattern = "P", gtf.df.genename$gene_id),]
gtf.para$gene_name <- gtf.para$gene_id

i.h.gn <- left_join(i.h, n)
i.p.gn <- i.p
i.p.gn$gene_name <- i.p.gn$gene_id

j <- cbind(i.h.gn, i.p.gn)
write.table(j, "host_para_overall0.95cor_interactions.txt", sep = '\t', row.names = F)

# add gene product for plasmodium
para_gp <- Barseq20200228[,c(1,4)]
colnames(para_gp)[1] <- "gene_name"
i.p.gn.gp <- left_join(i.p.gn, para_gp)

library(readxl)
mouse_overall0_9_genes_uniprot <- read_excel("Downloads/mouse_overall0.9_genes_uniprot.xlsx")[,c(1,6)]
host_gp <- mouse_overall0_9_genes_uniprot[!duplicated(mouse_overall0_9_genes_uniprot$gene_id),]

i.h.gn.gp <- left_join(i.h.gn, host_gp)
j <- cbind(i.h.gn.gp, i.p.gn.gp)
write.table(j, "host_para_overall0.95cor_interactions.txt", sep = '\t', row.names = F)


# phenotype tests
# 
i.p.agg <- i.p %>% count(gene_id)
h <- merge(i.p.agg, Barseq20200228, by.x = "gene_id", by.y = "gene")

phenotype <-c("Essential","Dispensible", "Fast", "Insufficient Data", "Slow")
color.codes<-as.character(c("#3399FF", "#FF00F0", "#Fd1F65",
                            "#FFF000", "#5efa0b00"))
color.names<-c("blue", "purple", "pink", "yellow", "green")
df2=data.frame(phenotype, color.codes, color.names); df2

df <- left_join(h,df2); df
v <- rep("#FF00F0", 15)
v[14] <- "Fd1F65"

df$color.codes <- v

x <- rep("purple", 15)
x[14] <- "pink"
df$color.names <- x

df <- df[,c(1,2,7,11,12,13)]

ggplot(data=df, aes(x=n, y=Confidence, colour = phenotype))+
  geom_point(alpha = 0.4) +
  scale_colour_manual(values=setNames(color.codes, phenotype))

png("Pb_cor0.9_PlasmoGEM.png")
ggplot(data=df, aes(x=n, y=Confidence, colour = phenotype, label=gene_id))+
  geom_point(alpha = 0.4) +
  geom_text(aes(label=df$gene_id),hjust=0, vjust=0, size = 2) +
  ggtitle("Plasmodium berghei gene interactors in correlation >= 0.9 (PlasmoGEM)") +
  xlab("node degree")
dev.off()
  