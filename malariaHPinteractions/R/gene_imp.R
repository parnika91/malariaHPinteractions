library(dplyr)
library(igraph)
library(ggplot2)
library(lrtest)
library(fitdistrplus)

load("~/Documents/Data/overall_bipartite_with_6_datasets.RData")
Barseq20200228 <- read.csv("~/Downloads/Barseq20200228.csv", stringsAsFactors=FALSE)
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

write.table(i.p, "parasite0.9_details_plasmogem.txt", sep = '\t', row.names = F)
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

load("df_concat_allhosts/cor/overall_6_datasets_para_edges.RData")
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
i.p.agg <- i.p %>% count(Orthogroup, gene_id)
h <- merge(i.p.agg, Barseq20200228, by.x = "gene_id", by.y = "gene")

phenotype.y <-c("Essential","Dispensable", "Fast", "Insufficient data", "Slow")
color.codes<-as.character(c("#3399FF", "#FF00F0", "#Fd1F65",
                            "#FFF000", "#5efa0b00"))
color.names<-c("blue", "purple", "pink", "yellow", "green")
df2=data.frame(phenotype.y, color.codes, color.names); df2

df <- left_join(rgr_mis_allgenes_clean,df2); df
# v <- rep("#FF00F0", 15)
# v[14] <- "Fd1F65"
# 
# df$color.codes <- v
# 
# x <- rep("purple", 15)
# x[14] <- "pink"
# df$color.names <- x
# 
# df <- df[,c(1,2,7,11,12,13)]

ggplot(data=df, aes(x=pb_pp_ec, y=RGR, colour = phenotype.y))+
  geom_point(alpha = 0.8) +
  scale_colour_manual(values=unique(as.character(df$color.codes))) +
  theme_bw() + 
  ggtitle("Relative growth rate vs eigen centrality") +
  xlab("Eigenvector centrality") + 
  ylab("Relative Growth Rate")
ggsave("RGR_vs_ec_pb_pp_ec.png")

n_ess <- df[which(df$phenotype == "Essential"), "n"]
n_dis <- df[which(df$phenotype == "Dispensable"), "n"]

png("Pb_cor0.9_PlasmoGEM.png")
ggplot(data=df, aes(x=n, y=Confidence, colour = phenotype, label=gene_id))+
  geom_point(alpha = 0.4) +
  # geom_text(aes(label=df$gene_id),hjust=0, vjust=0, size = 2) +
  ggtitle("Plasmodium berghei gene interactors in correlation >= 0.9 (PlasmoGEM)") +
  xlab("node degree")
dev.off()


# betweenness and closness
host_genes <- sapply(overall_0.9[,1], function(x) strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\1 \\2', x), ' ')[[1]][1])
hp <- data.frame(h = host_genes, p = para_genes, cor = overall[,2])

# degree 
load("overall_6_datasets_para_edges.RData")
para_genes <- c(as.character(para[,1]), as.character(para[,2]))
para_genes <- as.data.frame(para_genes)
colnames(para_genes) <- "para"
para_deg <- para_genes %>% count(para)
save(para_deg, file = "SameSize_para_para_degree.RData")

ig <-graph_from_data_frame(hp, directed = F)
bw <- betweenness(ig, directed=F, weights=NA)
save(bw, file = "bw_ERP004598_all_bipartite.RData")
bw.df <- as.data.frame(bw) %>% tibble::rownames_to_column("Orthogroup")
f <- inner_join(df, bw.df)

cl <- closeness(ig)
save(cl, file = "cl_ERP004598_all_bipartite.RData")
cl.df <- as.data.frame(cl) %>% tibble::rownames_to_column("Orthogroup")
df <- inner_join(f, cl.df)

save(df,file = "degree_bw_cl_ERP004598_all_bipartite_edges.RData")
write.table(df, "degree_bw_cl_ERP004598_all_bipartite_edges.txt", row.names = F, sep = '\t')

Disp <- df[which(df$phenotype == "Dispensable"),]

lm_eqn <- function(df, y, x){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

################# plots #######################
## disp_ov_bip: RGR vs log10(n) ##
png("disp_ov_bip_RGR_vs_log10degree.png", width = 40, height = 30, res = 300, units = "cm")
ggplot(disp_ov_bip, aes(log10(n), Relative.Growth.Rate)) + 
  geom_point(alpha = 0.7) +
  geom_smooth(method='lm', formula= y~x) +
  ggtitle("RGR vs log10(node degree) in dispensable genes\n Bipartite edges (ISIGEM 0)") +
  xlab("log10 node degree") +
  ylab("Relative growth rate") +
  theme_bw() + theme(text = element_text(size = 20)) +
  geom_text(size = 10, x = 2.5, y = 0.65, label = lm_eqn(disp_ov_bip, disp_ov_bip$Relative.Growth.Rate, log10(disp_ov_bip$n)), parse = T)
dev.off()
##
## disp_ov_bip: RGR vs log10(bw+1) ##
png("disp_ov_bip_RGR_vs_log10bw.png", width = 40, height = 30, res = 300, units = "cm")
ggplot(disp_ov_bip, aes(log10(bw+1), Relative.Growth.Rate)) + 
  geom_point(alpha = 0.7) +
  geom_smooth(method='lm', formula= y~x) +
  ggtitle("RGR vs log10(betweenness) in dispensable genes\n Bipartite edges (ISIGEM 0)") +
  xlab("log10 betweenness") +
  ylab("Relative growth rate") +
  theme_bw() + theme(text = element_text(size = 20)) +
  geom_text(size = 10, x = 2.5, y = 0.65, label = lm_eqn(disp_ov_bip, disp_ov_bip$Relative.Growth.Rate, log10(disp_ov_bip$bw+1)), parse = T)
dev.off()
##

## disp_ov_bip: RGR vs -log(cl) ##
png("disp_ov_bip_RGR_vs_log10cl.png", width = 40, height = 30, res = 300, units = "cm")
ggplot(disp_ov_bip, aes(-log10(cl), Relative.Growth.Rate)) + 
  geom_point(alpha = 0.7) +
  geom_smooth(method='lm', formula= y~x) +
  ggtitle("RGR vs -log10(closeness) in dispensable genes\n Bipartite edges (ISIGEM 0)") +
  xlab("-log10 closeness") +
  ylab("Relative growth rate") +
  theme_bw() + theme(text = element_text(size = 20)) +
  geom_text(size = 10, x = 4.6, y = 0.65, label = lm_eqn(disp_ov_bip, disp_ov_bip$Relative.Growth.Rate, -log10(disp_ov_bip$cl)), parse = T)
dev.off()
##

## disp_ov_para: RGR vs n ##
png("disp_ov_para_RGR_vs_log10degree.png", width = 40, height = 30, res = 300, units = "cm")
ggplot(disp_ov_para, aes(log10(n), Relative.Growth.Rate)) + 
  geom_point(alpha = 0.7) +
  geom_smooth(method='lm', formula= y~x) +
  ggtitle("RGR vs log10(node degree) in dispensable genes\n para-para edges (ISIGEM 0)") +
  xlab("log10 node degree") +
  ylab("Relative growth rate") +
  theme_bw() + theme(text = element_text(size = 20)) +
  geom_text(size = 10, x = 2.5, y = 0.65, label = lm_eqn(disp_ov_para, disp_ov_para$Relative.Growth.Rate, log10(disp_ov_para$n)), parse = T)
dev.off()
##

## disp_ov_para: RGR vs log10(bw+1) ##
png("disp_ov_para_RGR_vs_log10bw.png", width = 40, height = 30, res = 300, units = "cm")
ggplot(disp_ov_para, aes(log10(bw+1), Relative.Growth.Rate)) + 
  geom_point(alpha = 0.7) +
  geom_smooth(method='lm', formula= y~x) +
  ggtitle("RGR vs log10(betweenness) in dispensable genes\n para-para edges (ISIGEM 0)") +
  xlab("log10 betweenness") +
  ylab("Relative growth rate") +
  theme_bw() + theme(text = element_text(size = 20)) +
  geom_text(size = 10, x = 2.5, y = 0.65, label = lm_eqn(disp_ov_para, disp_ov_para$Relative.Growth.Rate, log10(disp_ov_para$bw+1)), parse = T)
dev.off()
##

## disp_ov_para: RGR vs -log(cl) ##
png("disp_ov_para_RGR_vs_log10cl.png", width = 40, height = 30, res = 300, units = "cm")
ggplot(disp_ov_para, aes(-log10(cl), Relative.Growth.Rate)) + 
  geom_point(alpha = 0.7) +
  geom_smooth(method='lm', formula= y~x) +
  ggtitle("RGR vs -log10(closeness) in dispensable genes\n para-para edges (ISIGEM 0)") +
  xlab("-log10 closeness") +
  ylab("Relative growth rate") +
  theme_bw() + theme(text = element_text(size = 20)) +
  geom_text(size = 10, x = 3.75, y = 0.65, label = lm_eqn(disp_ov_para, disp_ov_para$Relative.Growth.Rate, -log10(disp_ov_para$cl)), parse = T)
dev.off()
##

## disp_ss_bip: RGR vs log10(n) ##
png("disp_ss_bip_RGR_vs_log10degree.png", width = 40, height = 30, res = 300, units = "cm")
ggplot(disp_ss_bip, aes(log10(n), Relative.Growth.Rate)) + 
  geom_point(alpha = 0.7) +
  geom_smooth(method='lm', formula= y~x) +
  ggtitle("RGR vs log10(node degree) in dispensable genes\nBipartite edges (ISIGEM 0) smaller dataset") +
  xlab("log10 node degree") +
  ylab("Relative growth rate") +
  theme_bw() + theme(text = element_text(size = 20)) +
  geom_text(size = 10, x = 2.5, y = 0.65, label = lm_eqn(disp_ss_bip, disp_ss_bip$Relative.Growth.Rate, log10(disp_ss_bip$n)), parse = T)
dev.off()
##
## disp_ss_bip: RGR vs log10(bw) ##
png("disp_ss_bip_RGR_vs_log10bw.png", width = 40, height = 30, res = 300, units = "cm")
ggplot(disp_ss_bip, aes(log10(bw+1), Relative.Growth.Rate)) + 
  geom_point(alpha = 0.7) +
  geom_smooth(method='lm', formula= y~x) +
  ggtitle("RGR vs log10(betweenness) in dispensable genes\nBipartite edges (ISIGEM 0) smaller dataset") +
  xlab("log10 betweenness") +
  ylab("Relative growth rate") +
  theme_bw() + theme(text = element_text(size = 20)) +
  geom_text(size = 10, x = 2.5, y = 0.65, label = lm_eqn(disp_ss_bip, disp_ss_bip$Relative.Growth.Rate, log10(disp_ss_bip$bw+1)), parse = T)
dev.off()
##
## disp_ss_bip: RGR vs log10(cl) ##
png("disp_ss_bip_RGR_vs_log10cl.png", width = 40, height = 30, res = 300, units = "cm")
ggplot(disp_ss_bip, aes(-log10(cl), Relative.Growth.Rate)) + 
  geom_point(alpha = 0.7) +
  geom_smooth(method='lm', formula= y~x) +
  ggtitle("RGR vs log10(closeness) in dispensable genes\nBipartite edges (ISIGEM 0) smaller dataset") +
  xlab("-log10 closeness") +
  ylab("Relative growth rate") +
  theme_bw() + theme(text = element_text(size = 20)) +
  geom_text(size = 10, x = 6.5, y = 0.65, label = lm_eqn(disp_ss_bip, disp_ss_bip$Relative.Growth.Rate, -log10(disp_ss_bip$cl)), parse = T)
dev.off()
##
## disp_ss_para: RGR vs log10(n) ##
png("disp_ss_para_RGR_vs_log10degree.png", width = 40, height = 30, res = 300, units = "cm")
ggplot(disp_ss_para, aes(log10(n), Relative.Growth.Rate)) + 
  geom_point(alpha = 0.7) +
  geom_smooth(method='lm', formula= y~x) +
  ggtitle("RGR vs log10(node degree) in dispensable genes\npara-para edges (ISIGEM 0) smaller dataset") +
  xlab("log10 node degree") +
  ylab("Relative growth rate") +
  theme_bw() + theme(text = element_text(size = 20)) +
  geom_text(size = 10, x = 2.5, y = 0.65, label = lm_eqn(disp_ss_para, disp_ss_para$Relative.Growth.Rate, log10(disp_ss_para$n)), parse = T)
dev.off()
##

## disp_ss_para: RGR vs log10(bw) ##
png("disp_ss_para_RGR_vs_log10bw.png", width = 40, height = 30, res = 300, units = "cm")
ggplot(disp_ss_para, aes(log10(bw+1), Relative.Growth.Rate)) + 
  geom_point(alpha = 0.7) +
  geom_smooth(method='lm', formula= y~x) +
  ggtitle("RGR vs log10(betweenness) in dispensable genes\npara-para edges (ISIGEM 0) smaller dataset") +
  xlab("log10 betweenness") +
  ylab("Relative growth rate") +
  theme_bw() + theme(text = element_text(size = 20)) +
  geom_text(size = 10, x = 2.5, y = 0.65, label = lm_eqn(disp_ss_para, disp_ss_para$Relative.Growth.Rate, log10(disp_ss_para$bw+1)), parse = T)
dev.off()
##

## disp_ss_para: RGR vs log10(cl) ##
png("disp_ss_para_RGR_vs_log10cl.png", width = 40, height = 30, res = 300, units = "cm")
ggplot(disp_ss_para, aes(-log10(cl), Relative.Growth.Rate)) + 
  geom_point(alpha = 0.7) +
  geom_smooth(method='lm', formula= y~x) +
  ggtitle("RGR vs log10(closeness) in dispensable genes\npara-para edges (ISIGEM 0) smaller dataset") +
  xlab("-log10 closeness") +
  ylab("Relative growth rate") +
  theme_bw() + theme(text = element_text(size = 20)) +
  geom_text(size = 10, x = 3.8, y = 0.65, label = lm_eqn(disp_ss_para, disp_ss_para$Relative.Growth.Rate, -log10(disp_ss_para$cl)), parse = T)
dev.off()
##
######################################
load("df_concat_allhosts/cor/overall_6_datasets_para_edges.RData")
Barseq20200228 <- read.csv("Barseq20200228.csv", stringsAsFactors=FALSE)
parasite_orthogroups <- read.delim("/SAN/Plasmo_compare/OrthoFinder/parasite_orthogroups.txt", stringsAsFactors=FALSE)
para_genes <- c()
para_genes <- c(c(as.character(para[,1])), c(as.character(para[,2])))
para_df <- as.data.frame(para_genes)
colnames(para_df) <- "Orthogroup"
#para_ <- table(para_df)
para_ <- parasite_orthogroups[,c(1,4)]
colnames(para_)[2] <- "gene_id"
#h <- merge(para, Barseq20200228, by.x = "Pb_g", by.y = "gene")
i.p <- inner_join(para_df, para_, by.x = "para_genes", by.y = "Orthogroup")

i.p.agg <- i.p %>% count(Orthogroup, gene_id)
h <- merge(i.p.agg, Barseq20200228, by.x = "gene_id", by.y = "gene")

phenotype <-c("Essential","Dispensable", "Fast", "Insufficient data", "Slow")
color.codes<-as.character(c("#3399FF", "#FF00F0", "#Fd1F65",
                            "#FFF000", "#5efa0b00"))
color.names<-c("blue", "purple", "pink", "yellow", "green")
df2=data.frame(phenotype, color.codes, color.names); df2

df <- left_join(h,df2); df

ig <-graph_from_data_frame(para[,1:3], directed = F)
bw <- betweenness(ig, directed=F, weights=NA)
save(bw, file = "bw_ERP004598_all_para_edges.RData")
bw.df <- as.data.frame(bw) %>% tibble::rownames_to_column("Orthogroup")
f <- inner_join(df, bw.df)

cl <- closeness(ig, weight=NA)
save(cl, file = "cl_ERP004598_all_para_edges.RData")
cl.df <- as.data.frame(cl) %>% tibble::rownames_to_column("Orthogroup")
df <- inner_join(f, cl.df)
save(df,file = "degree_bw_cl_ERP004598_all_para_edges.RData")
write.table(df, "degree_bw_cl_ERP004598_all_para_edges.txt", row.names = F, sep = '\t')

Disp <- df[which(df$phenotype == "Dispensable"),]

lm_eqn <- function(df, y, x){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

degree_g <- ggplot(Disp, aes(n, Relative.Growth.Rate)) + 
  geom_point(alpha = 0.7) +
  geom_smooth(method='lm', formula= y~x) +
  geom_text(x = 2000, y = 0.65, label = lm_eqn(Disp, Disp$Relative.Growth.Rate, Disp$n), parse = T)

bw_g <- ggplot(Disp, aes(log10(bw), Relative.Growth.Rate)) + 
  geom_point(alpha = 0.7) +
  geom_smooth(method='lm', formula= y~x) +
  geom_text(x = 2, y = 0.65, label = lm_eqn(Disp, Disp$Relative.Growth.Rate, Disp$bw), parse = T)

cl_g <- ggplot(Disp, aes(cl, Relative.Growth.Rate)) + 
  geom_point(alpha = 0.7) +
  geom_smooth(method='lm', formula= y~x) +
  geom_text(x = 4.575, y = 0.65, label = lm_eqn(Disp, Disp$Relative.Growth.Rate, Disp$cl), parse = T)

############ get pcor, pval, scor, pval ##################
# load df
cor.test(df$n, df$Relative.Growth.Rate)
cor.test(df$n, df$Relative.Growth.Rate, method = "spearman")
cor.test(df$bw, df$Relative.Growth.Rate)
cor.test(df$bw, df$Relative.Growth.Rate,method = "spearman")
cor.test(df$cl, df$Relative.Growth.Rate)
cor.test(df$cl, df$Relative.Growth.Rate, method = "spearman")
disp <- df[df$phenotype%in%"Dispensable",]
cor.test(disp$n, disp$Relative.Growth.Rate)
cor.test(disp$n, disp$Relative.Growth.Rate, method = "spearman")
cor.test(disp$bw, disp$Relative.Growth.Rate)
cor.test(disp$bw, disp$Relative.Growth.Rate, method = "spearman")
cor.test(disp$cl, disp$Relative.Growth.Rate)
cor.test(disp$cl, disp$Relative.Growth.Rate, method = "spearman")
ess <- df[df$phenotype%in%"Essential",]
cor.test(ess$n, ess$Relative.Growth.Rate)
cor.test(ess$n, ess$Relative.Growth.Rate,  method = "spearman")
cor.test(ess$bw, ess$Relative.Growth.Rate)
cor.test(ess$bw, ess$Relative.Growth.Rate , method = "spearman")
cor.test(ess$cl, ess$Relative.Growth.Rate)
cor.test(ess$cl, ess$Relative.Growth.Rate , method = "spearman")
essl <- df[(df$phenotype%in%"Essential" | df$phenotype%in%"Slow"),]
cor.test(essl$n, essl$Relative.Growth.Rate)
cor.test(essl$n, essl$Relative.Growth.Rate, method = "spearman")
cor.test(essl$bw, essl$Relative.Growth.Rate)
cor.test(essl$bw, essl$Relative.Growth.Rate, method = "spearman")
cor.test(essl$cl, essl$Relative.Growth.Rate)
cor.test(essl$cl, essl$Relative.Growth.Rate , method = "spearman")

##################################################

# lm_eqn <- function(df, y, x, z){
#   m <- lm(y ~ x + z, df);
#   eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
#                    list(a = format(unname(coef(m)[1]), digits = 2),
#                         b = format(unname(coef(m)[2]), digits = 2),
#                         r2 = format(summary(m)$r.squared, digits = 3)))
#   as.character(as.expression(eq));
# }
# 
# degree_g <- ggplot(Disp, aes(n, Relative.Growth.Rate)) + 
#   geom_point(alpha = 0.7) +
#   geom_smooth(method='lm', formula= y~x) +
#   geom_text(x = 2000, y = 0.65, label = lm_eqn(Disp, Disp$Relative.Growth.Rate, Disp$n), parse = T)
# 
install.packages("devtools")
devtools::install_github("cardiomoon/ggiraphExtra")

require(ggplot2)
install.packages("stargazer") #Produces easy to read regression results (similar to what you get in SPSS)
install.packages("effects") #We will use this to create our interactions 
#install.packages("ggplot2")
install.packages("psych")
m.add <- lm(formula = Relative.Growth.Rate ~ cl+n, data = df)
m.int <- lm(formula = Relative.Growth.Rate ~ cl*n, data = df)

library(stargazer)
stargazer(m.add, m.int,type="text", 
          column.labels = c("Main Effects", "Interaction"), 
          intercept.bottom = FALSE, 
          single.row=FALSE,     
          notes.append = FALSE, 
          header=FALSE) 

ggplot(df,aes(y=Relative.Growth.Rate,x=-log10(cl),color=n))+geom_point()+stat_smooth(method="lm",se=FALSE)
summary(m.int)
vcov(m.int)
describe(df)
install.packages("sjPlot")
plot_model(m.int, type = "int", mdrt.values = "meansd")

#check correlation of original RGR and fitted RGR

oriRGR <- na.omit(df$Relative.Growth.Rate)
cor(oriRGR, fitted(m.int))


# likelihood ratio test
df <- df[df$Relative.Growth.Rate <= 1,]
require(lmtest)
m.add <- glm(formula = Relative.Growth.Rate ~ cl+n+bw, data = na.omit(df))
m.int <- glm(formula = Relative.Growth.Rate ~ cl*n*bw, data = na.omit(df))
lrtest(m.int, m.add)

require(fitdistrplus)
fitdistrplus::descdist(as.numeric(na.omit(df$Relative.Growth.Rate)))


require(betareg)
df <- df[which(df$Relative.Growth.Rate<=1),]
mod1 <- betareg(Relative.Growth.Rate ~  n*bw, data = na.omit(df))
