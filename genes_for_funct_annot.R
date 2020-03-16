library(rtracklayer)
setwd("~/")

gtf <- import("~/Documents/Data/Kai/mousePberghei.gtf", format = "gtf")
gtf.df <- as.data.frame(gtf)
gtf.df <- gtf.df[gtf.df$type%in%"exon",]
rm(gtf)
gtf.gene.names <- gtf.df[,c("gene_id", "gene_name")]
gtf.gene.names <- gtf.gene.names[!duplicated(gtf.gene.names$gene_name),]

load("~/Documents/Data/overall_bipartite_with_6_datasets.RData")
parasite_orthogroups <- read.delim("~/Documents/Data/parasite_orthogroups.txt", stringsAsFactors=FALSE)
host_orthogroups <- read.delim("~/Documents/Data/host_orthogroups.txt", stringsAsFactors=FALSE)

overall <- overall[which(abs(overall$cor) >= 0.95),]

para_genes <- sapply(overall[,1], function(x) strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\1 \\2', x), ' ')[[1]][2])
para_df <- as.data.frame(para_genes)
para <- parasite_orthogroups[,c(1,4)]
colnames(para_df) <- "Orthogroup"
ip <- inner_join(para_df, para)

host_genes <- sapply(overall[,1], function(x) strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\1 \\2', x), ' ')[[1]][1])
host_df <- as.data.frame(host_genes)
host <- host_orthogroups[,c(1,3)]
colnames(host_df) <- "Orthogroup"
ih <- inner_join(host_df, host)

colnames(ih)[2] <- "gene_id"
host.gene.names <- inner_join(ih, gtf.gene.names)

hp <- cbind(host.gene.names, ip)
hp <- hp[!duplicated(hp),]

write.table(hp, "hp_0.85interactions.txt", sep = "\t", row.names=F)
