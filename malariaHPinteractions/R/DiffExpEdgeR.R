###################################################
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("edgeR")
# BiocManager::install("org.Mm.eg.db")
# install.packages("ggplot2")
# install.packages("gridExtra")
# install.packages("plotly")
# suppressPackageStartupMessages(library("plotly"))

# x <- read.delim("Macrophage_coding_genes.txt")

# Conditions <- read.csv("~/Documents/Data/Kai/Conditions_edited_to_E.csv", stringsAsFactors=FALSE)
# Conditions$Trt_time <- paste(Conditions$treatment, Conditions$time)

# CountSamples <- as.matrix(as.character(rownames(y$samples)))
# CountSamples[] <- apply(CountSamples, 1, as.character)
# CountSamples <- as.matrix(sapply(CountSamples, function(x) strsplit(x, split = "_")[[1]][1]))
# #saved this and cleaned the file. Added 0 to S3E03 etc and loaded it back again.
# CountSamples <- as.matrix(sapply(CountSamples, function(x) substring(x, 2)))
# 
# colnames(CountSamples) <- "sample"
# CountSampleCond <- merge(CountSamples, Conditions)
# 
# CountSampleCond <- CountSampleCond[,c(1,7)]
# for(i in 1:nrow(CountSampleCond))
# {
#   CountSampleCond[i,3] <- colnames(y$counts)[grep(CountSampleCond[i,1], colnames(y$counts))]
# }

# SampleCond <- data.frame(CountfileSample = colnames(x))
# for(i in 1:nrow(Conditions))
# {
#   a <- grep(pattern = Conditions[i,1], SampleCond[,1])
#   SampleCond[a,2] <- Conditions[i,1]
#   SampleCond[a,3] <- Conditions[i,7]
# }
# write.table(SampleCond, "SampleCond.txt", row.names = F)
# edited where necessary (1E11 is not in Conditions) and loaded back
# SampleCond <- read.csv("~/Documents/Data/Kai/SampleCond.txt", sep="", stringsAsFactors=FALSE)

# SampleCond <- na.omit(SampleCond)
# x <- x[,SampleCond$CountfileSample]
# colnames(x) <- SampleCond[,2]
# 
# write.table(x, "Macrophage_RNAseqReadCount.txt", row.names = T, sep = '\t')

######## Analysis starts here #######

library(edgeR)
library(org.Mm.eg.db)
library(pheatmap)
library(rtracklayer)
library(dplyr)
library(ggplot2)
library(ggrepel)
setwd("~/Documents/Data/Kai/")

SampleCond <- read.csv("~/Documents/Data/Kai/SampleCond.txt", sep="", stringsAsFactors=FALSE)
countfile <- read.delim("~/Documents/Data/Kai/Macrophage_RNAseqReadCount.txt", stringsAsFactors=FALSE)
countfile <- countfile[-grep(rownames(countfile), pattern = "PB"),]

group <- SampleCond[,3]

# remove space so I can use the names later to makeContrasts
group <- gsub(" ", "_", group)

# make the count file a DGEList file
y <- DGEList(counts=countfile, group=group)

# remove low expression genes so that FDR depends on fewer genes
keep <- filterByExpr(y) 
y <- y[keep, , keep.lib.sizes=FALSE]

# before filtering, nrow(y) = 26961
# after filtering, nrow(y) = 13019

# normalise gene counts for lib size and RNA composition effect
y <- calcNormFactors(y)

# make design matrix. We will do many different models
# If the desired condition is not automatically made the baseline
# when using ~group, do it manually using relevel
design <- model.matrix(~0+group)

y$samples$group <- relevel(y$samples$group, ref="untr_4h")
design1 <- model.matrix(~group)

# get estimation of dispersion
y <- estimateDisp(y, design = design)

# fit negative binomial distribution 
fit <- glmQLFit(y, design = design)

# test hypothesis: gene different in any of the groups
qlf <- glmQLFTest(fit, coef = 2:16)
# or, to give a threshold of fold change,
# qlf <- glmTreat(fit, lfc = 2)
#res <- topTags(qlf, n=nrow(qlf$table))
# de.genes <- rownames(topTags(qlf, n=100)$table)
diff <- topTags(qlf, n = nrow(qlf$table))$table
write.table(diff, "diff_genes_all_conditions.txt", sep = '\t', row.names = T)

# heatmap
logcpm <- cpm(y, log=TRUE)
colnames(logcpm) <- y$samples$group

############## One-time plots and analysis ##########
# of samples
# sampleDists <- as.matrix(dist(t(logcpm)))
# colnames(sampleDists) <- rownames(sampleDists) <- y$samples$group

# pdf("macrophage_RNAseq_sample_heatmap.pdf", onefile = T)#,  width = 50, height  = 50, unit = "cm", res = 450
# p <- pheatmap(sampleDists, fontsize = 3)
# dev.off()

# png("macrophage_RNAseq_sample_heatmap.png", width = 50, height  = 50, unit = "cm", res = 450)# onefile = T,  
# pheatmap(sampleDists, fontsize = 5)
# dev.off()

# use normalised y
# pdf("macrophage_RNAseq_genes_heatmap.pdf", onefile = T)
# p <- pheatmap(logcpm, fontsize = 1)
# dev.off()

# selected samples
# pdf("macrophage_RNAseq_genes_heatmap_selected.pdf", onefile = T)
# p <- pheatmap(logcpm[,c(12:15,21,99)], fontsize = 10)
# dev.off()

# top genes
# diff <- topTags(qlf, n = 20)$table
# res_genes <- rownames(diff)
# res_logcpm <- data.frame()
# for(i in 1:length(res_genes))
# {
#   if(grep(rownames(logcpm), pattern = res_genes[i]))
#   {
#     a <- grep(rownames(logcpm), pattern = res_genes[i])
#     res_logcpm[i,1:ncol(logcpm)] <- logcpm[a,]
#     rownames(res_logcpm)[i] <- res_genes[i]
#   }
# }
# colnames(res_logcpm) <- colnames(logcpm)
# 
# gene_name <- c()
# 
# for(i in 1:length(res_genes))
# {
#   if(grep(rownames(logcpm), pattern = res_genes[i]))
#   {
#     a <- grep(rownames(logcpm), pattern = res_genes[i])
#     res_logcpm[i,1:ncol(logcpm)] <- logcpm[a,]
#     rownames(res_logcpm)[i] <- res_genes[i]
#     
#     if(!is.na(gtf.df[grep(pattern = res_genes[i], gtf.df$gene_id)[1], "gene_name"])) {
#       gene_name[i] <- gtf.df[grep(pattern = res_genes[i], gtf.df$gene_id)[1], "gene_name"]
#     } else {
#       gene_name[i] <- rownames(res_logcpm)[i]
#     }
#   }
# }
# colnames(res_logcpm) <- colnames(logcpm)
# res_logcpm$gene_name <- gene_name
# 
# pdf("macrophage_RNAseq_top20genes_heatmap.pdf", onefile = T)
# p <- pheatmap(res_logcpm[,c(1:140)], fontsize_row = 5, fontsize_col = 2.5, labels_row = gene_name)
# dev.off()

# MDS
# pdf("macrophage_RNAseq_sample_MDS.pdf", onefile = T)#,  width = 50, height  = 50, unit = "cm", res = 450
# plotMDS(sampleDists, cex = 0.3)
# dev.off()

# MA plot
# pdf("macrophage_RNAseq_sample_MA.pdf", onefile = T)
# plotSmear(y, de.tags = de.genes)
# dev.off()

# get gene names from ensembl gene IDs
gtf <- import("mousePberghei.gtf", format = "gtf")
gtf <- gtf[gtf$type%in%"exon"]
gtf <- gtf[gtf$gene_biotype%in%"protein_coding"]
gtf.df <- as.data.frame(gtf)
# gtf.h <- gtf[grep(pattern = "ENS", gtf$gene_id),]

############ many contrasts #############

# plot heatmap for each condition and save table
Diff_gene_condition <- function(ref, trt, time)
{
  if(time == 24)
    qlf <- glmQLFTest(fit, contrast=my.contrasts.24[,paste0(ref, "_vs_", trt, collapse = '')])
  if(time == 4)
    qlf <- glmQLFTest(fit, contrast=my.contrasts.4[,paste0(ref, "_vs_", trt, collapse = '')])
  if(time == "inter")
    qlf <- glmQLFTest(fit, contrast=my.contrasts.inter[,paste0(ref, "_vs_", trt, collapse = '')])
  if(time == "overall")
    qlf <- glmQLFTest(fit, contrast=my.contrasts.single[,paste0(ref, "_vs_", trt, collapse = '')])
  
  diff <- topTags(qlf, n = 50)$table
  res_genes <- rownames(diff)
  
  res_logcpm <- data.frame()
  gene_name <- c()
  
  if(time == "overall")
    logcpm = logcpm2
  
  for(i in 1:length(res_genes))
  {
    if(grep(rownames(logcpm), pattern = res_genes[i]))
    {
      a <- grep(rownames(logcpm), pattern = res_genes[i])
      res_logcpm[i,1:ncol(logcpm)] <- logcpm[a,]
      rownames(res_logcpm)[i] <- res_genes[i]
      
      if(!is.na(gtf.df[grep(pattern = res_genes[i], gtf.df$gene_id)[1], "gene_name"])) {
        gene_name[i] <- gtf.df[grep(pattern = res_genes[i], gtf.df$gene_id)[1], "gene_name"]
      } else {
        gene_name[i] <- rownames(res_logcpm)[i]
      }
    }
  }
  colnames(res_logcpm) <- colnames(logcpm)
  res_logcpm$gene_name <- gene_name
  
  pdf(paste0("macrophage_RNAseq_top50_", ref, "_vs_", trt, "_heatmap.pdf", collapse = ''), onefile = T)
  p <- pheatmap(res_logcpm[,c(which(colnames(res_logcpm)==ref),which(colnames(res_logcpm)==trt))], 
                fontsize_row = 5, fontsize_col = 5, border_color = NA,
                main = paste0("Top 50 genes ",ref, " vs ", trt, collapse = ''), 
                labels_row = res_logcpm$gene_name)
  dev.off()
  
  # save the entire table
  diff <- topTags(qlf, n = nrow(qlf$table))$table
  
  host_gene_name_table <- na.omit(data.frame(Ensembl_geneID = gtf$gene_id[grep(pattern = "ENS", gtf$gene_id)], Gene_name = gtf$gene_name[grep(pattern = "ENS", gtf$gene_id)]))
  para_gene_name_table <- na.omit(data.frame(Ensembl_geneID = gtf$gene_id[grep(pattern = "PB", gtf$gene_id)], Gene_name = gtf$gene_id[grep(pattern = "PB", gtf$gene_id)]))
  gene_name_table <- rbind(host_gene_name_table, para_gene_name_table)
  gene_name_table <- unique(gene_name_table)
  
  diff$Ensembl_geneID <- rownames(diff)
  ens.df <- as.data.frame(diff$Ensembl_geneID)
  colnames(ens.df) <- "Ensembl_geneID"
  
  diff_names <- inner_join(ens.df, gene_name_table)
  diff$Gene_name <- diff_names$Gene_name
  diff <- diff[,c(6,7,1,2,3,4,5)]
  write.table(diff, paste0(ref,"_vs_", trt, ".txt", collapse = ''), row.names = F)
}

# 4 hours: untr@4 vs others@4
y$samples$group <- relevel(y$samples$group, ref="untr_4h")
design1 <- model.matrix(~0+group)
colnames(design1) <- substring(colnames(design1), 6)
fit <- glmQLFit(y, design1)

my.contrasts.4 <- makeContrasts(
  untr_4h_vs_SPZhi_4h = SPZhi_4h - untr_4h,
  untr_4h_vs_iRBChi_4h = iRBChi_4h - untr_4h,
  untr_4h_vs_iRBClo_4h = iRBClo_4h - untr_4h,
  untr_4h_vs_RBC_4h = RBC_4h - untr_4h,
  untr_4h_vs_SPZlo_4h = SPZlo_4h - untr_4h,
  untr_4h_vs_LPS_4h = LPS_4h - untr_4h,
  untr_4h_vs_PBS_4h = PBS_4h - untr_4h,
  
  LPS_4h_vs_SPZhi_4h =  SPZhi_4h - LPS_4h,
  LPS_4h_vs_SPZlo_4h = SPZlo_4h - LPS_4h,
  LPS_4h_vs_iRBChi_4h = iRBChi_4h - LPS_4h,
  LPS_4h_vs_iRBClo_4h = iRBClo_4h - LPS_4h,
  LPS_4h_vs_PBS_4h = PBS_4h - LPS_4h,
  
  PBS_4h_vs_SPZhi_4h = SPZhi_4h - PBS_4h,
  PBS_4h_vs_SPZlo_4h = SPZlo_4h - PBS_4h,
  PBS_4h_vs_iRBChi_4h = iRBChi_4h - PBS_4h,
  PBS_4h_vs_iRBClo_4h = iRBClo_4h - PBS_4h,
  PBS_4h_vs_LPS_4h = LPS_4h - PBS_4h,
  levels=design1
  )

Diff_gene_condition(ref = "PBS_4h", trt = "LPS_4h", time = 4)

# 24 hours: untr@24 vs others@24
y$samples$group <- relevel(y$samples$group, ref="untr_24h")
design1 <- model.matrix(~0+group)
colnames(design1) <- substring(colnames(design1), 6)
fit <- glmQLFit(y, design1)

my.contrasts.24 <- makeContrasts(
  untr_24h_vs_SPZhi_24h = SPZhi_24h - untr_24h,
  untr_24h_vs_iRBChi_24h = iRBChi_24h - untr_24h,
  untr_24h_vs_iRBClo_24h = iRBClo_24h - untr_24h,
  untr_24h_vs_RBC_24h = RBC_24h - untr_24h,
  untr_24h_vs_SPZlo_24h = SPZlo_24h - untr_24h,
  untr_24h_vs_LPS_24h = LPS_24h - untr_24h,
  untr_24h_vs_PBS_24h = PBS_24h - untr_24h,
  
  LPS_24h_vs_SPZhi_24h = SPZhi_24h - LPS_24h,
  LPS_24h_vs_SPZlo_24h = SPZlo_24h - LPS_24h,
  LPS_24h_vs_iRBChi_24h = iRBChi_24h - LPS_24h,
  LPS_24h_vs_iRBClo_24h = iRBClo_24h - LPS_24h,
  LPS_24h_vs_PBS_24h = PBS_24h - LPS_24h,
  
  PBS_24h_vs_SPZhi_24h = SPZhi_24h - PBS_24h,
    PBS_24h_vs_SPZlo_24h = SPZlo_24h - PBS_24h,
    PBS_24h_vs_iRBChi_24h = iRBChi_24h - PBS_24h,
    PBS_24h_vs_iRBClo_24h = iRBClo_24h - PBS_24h,
    PBS_24h_vs_LPS_24h = LPS_24h - PBS_24h,
    levels=design1
  )

Diff_gene_condition(ref = "LPS_24h", trt = "SPZlo_24h", time = 24)

# 4vs24 hours: 4 vs 24

y$samples$group <- relevel(y$samples$group, ref="untr_4h")
design1 <- model.matrix(~0+group)
colnames(design1) <- substring(colnames(design1), 6)
fit <- glmQLFit(y, design1)

my.contrasts.inter <- makeContrasts(
  untr_4h_vs_untr_24h = untr_4h - untr_24h,
  SPZhi_4h_vs_SPZhi_24h = SPZhi_4h - SPZhi_24h,
  iRBChi_4h_vs_iRBChi_24h = iRBChi_4h - iRBChi_24h,
  iRBClo_4h_vs_iRBClo_24h = iRBClo_4h - iRBClo_24h,
  RBC_4h_vs_RBC_24h = RBC_4h - RBC_24h,
  SPZlo_4h_vs_SPZlo_24h = SPZlo_4h - SPZlo_24h,
  LPS_4h_vs_LPS_24h = LPS_4h - LPS_24h,
  PBS_4h_vs_PBS_24h = PBS_4h - PBS_24h,
  levels=design1
)

Diff_gene_condition(ref = "PBS_4h", trt = "PBS_24h", time = "inter")

# overall contrasts

y2$samples$group <- relevel(y2$samples$group, ref="untr")
design2 <- model.matrix(~0+group2)
colnames(design2) <- substring(colnames(design2), 7)
fit <- glmQLFit(y2, design2)

my.contrasts.overall <- makeContrasts(
  untr_vs_SPZhi = SPZhi - untr,
  untr_vs_iRBChi = iRBChi - untr,
  untr_vs_iRBClo = iRBClo - untr,
  untr_vs_RBC = RBC - untr,
  untr_vs_SPZlo = SPZlo - untr,
  untr_vs_LPS = LPS - untr,
  untr_vs_PBS = PBS - untr,
  
  LPS_vs_SPZhi = SPZhi - LPS,
  LPS_vs_SPZlo = SPZlo - LPS,
  LPS_vs_iRBChi = iRBChi - LPS,
  LPS_vs_iRBClo = iRBClo - LPS,
  LPS_vs_PBS = PBS - LPS,
  
  PBS_vs_SPZhi = SPZhi - PBS,
  PBS_vs_SPZlo = SPZlo - PBS,
  PBS_vs_iRBChi = iRBChi - PBS,
  PBS_vs_iRBClo = iRBClo - PBS,
  PBS_vs_LPS = LPS - PBS,
  levels=design2
)

my.contrasts.single <- makeContrasts(
  LPS_vs_SPZhi = (SPZhi_24h - SPZhi_4h) - (LPS_24h - LPS_4h),
  levels=design1
)
Diff_gene_condition(ref = "LPS", trt = "SPZhi", time = "single")

############ Convert IDs using biomaRt #############

library(biomaRt)
  mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
genes <- getBM(
  filters="ensembl_gene_id",
  attributes=c("ensembl_gene_id", "external_gene_name"),
  values=res_genes,
  mart=mart)

###################################

qlf <- glmQLFTest(fit, coef=2)
go <- goana(qlf, species="Mm")
topGO(go, sort="up")
keg <- kegga(qlf, species="Mm")
topKEGG(keg, sort="up")

############ diff_df between conditions 1 and 2 #############

untr4_vs_conditions4 <- c("untr_4h_vs_SPZhi_4h",
                  "untr_4h_vs_iRBChi_4h",
                  "untr_4h_vs_iRBClo_4h",
                  "untr_4h_vs_RBC_4h",
                  "untr_4h_vs_SPZlo_4h",
                  "untr_4h_vs_LPS_4h",
                  "untr_4h_vs_PBS_4h",
                  "LPS_4h_vs_SPZhi_4h",
                  "LPS_4h_vs_SPZlo_4h",
                  "LPS_4h_vs_iRBChi_4h",
                  "LPS_4h_vs_iRBClo_4h",
                  "LPS_4h_vs_PBS_4h",
                  "PBS_4h_vs_SPZhi_4h",
                  "PBS_4h_vs_SPZlo_4h",
                  "PBS_4h_vs_iRBChi_4h",
                  "PBS_4h_vs_iRBClo_4h",
                  "PBS_4h_vs_LPS_4h")

untr24_vs_conditions24 <- c("untr_24h_vs_SPZhi_24h",
                   "untr_24h_vs_iRBChi_24h",
                   "untr_24h_vs_iRBClo_24h",
                   "untr_24h_vs_RBC_24h",
                   "untr_24h_vs_SPZlo_24h",
                   "untr_24h_vs_LPS_24h",
                   "untr_24h_vs_PBS_24h",
                   "LPS_24h_vs_SPZhi_24h",
                   "LPS_24h_vs_SPZlo_24h",
                   "LPS_24h_vs_iRBChi_24h",
                   "LPS_24h_vs_iRBClo_24h",
                   "LPS_24h_vs_PBS_24h",
                   "PBS_24h_vs_SPZhi_24h",
                   "PBS_24h_vs_SPZlo_24h",
                   "PBS_24h_vs_iRBChi_24h",
                   "PBS_24h_vs_iRBClo_24h",
                   "PBS_24h_vs_LPS_24h")

condition4_vs_condition24 <- c("untr_4h_vs_untr_24h",
                     "SPZhi_4h_vs_SPZhi_24h",
                     "iRBChi_4h_vs_iRBChi_24h",
                     "iRBClo_4h_vs_iRBClo_24h",
                     "RBC_4h_vs_RBC_24h",
                     "SPZlo_4h_vs_SPZlo_24h",
                     "LPS_4h_vs_LPS_24h",
                     "PBS_4h_vs_PBS_24h") 

volcano_plots <- function(conditions)
{
  if(conditions == "untr4_vs_conditions4")
    conds = untr4_vs_conditions4
  if(conditions == "untr24_vs_conditions24")
    conds = untr24_vs_conditions24
  if(conditions == "condition4_vs_condition24")
    conds = condition4_vs_condition24
  
  for(cond in conds)
  {
    DGE <- read.csv(paste0(conditions,"/", cond, ".txt"), sep="", stringsAsFactors=FALSE)
    diff_df <- DGE[c("Gene_name", "logFC", "FDR")]
    # diff_df <- diff_df[-grep(diff_df$Gene_name, pattern = "PB"),]
    
    # add a grouping column; default value is "not significant"
    diff_df["group"] <- "NotSignificant"
    
    # change the grouping for the entries with significance but not a large enough Fold change
    diff_df[which(diff_df['FDR'] < 0.05 & abs(diff_df['logFC']) < 1.5 ),"group"] <- "Significant"
    
    # change the grouping for the entries a large enough Fold change but not a low enough p value
    diff_df[which(diff_df['FDR'] > 0.05 & abs(diff_df['logFC']) > 1.5 ),"group"] <- "logFC"
    
    # change the grouping for the entries with both significance and large enough logFC change
    diff_df[which(diff_df['FDR'] < 0.05 & abs(diff_df['logFC']) > 1.5 ),"group"] <- "Significant&logFC"
    
    diff_df <- as.data.frame(diff_df)
    diff_df <-  diff_df[order(-abs(diff_df$logFC), diff_df$FDR),]
    
    # top_peaks <- diff_df[with(diff_df, order(logFC, FDR)),][1:5,]
    # top_peaks <- rbind(top_peaks, diff_df[with(diff_df, order(-logFC, FDR)),][1:5,])

    p <- ggplot(diff_df, aes(x = logFC, y = -log10(FDR), colour = group)) +
    geom_point() +
    ggtitle(cond) +
    geom_text_repel(
      data = head(diff_df, 100),
      aes(label = Gene_name),
      size = 2,
      colour = "black",
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines"))
      
      # a <- list()
      # for (i in seq_len(nrow(top_peaks))) {
      #   m <- top_peaks[i, ]
      #   a[[i]] <- list(   
      #     x = m[["logFC"]],
      #     y = -log10(m[["FDR"]]),
      #     text = m[["Ensembl_geneID"]],
      #     xref = "x",
      #     yref = "y",
      #     showarrow = TRUE,
      #     arrowhead = 0.5,
      #     ax = 20,
      #     ay = -40
      #   )
      # }
      # 
      # 
      # # make the Plot.ly plot
      # p <- plot_ly(data = diff_df, x = diff_df[,2], y = -log10(diff_df[,3]), text = diff_df[,1], mode = "markers", color = group) %>% 
      #   layout(title ="Volcano Plot") %>%
      #   layout(annotations = a)
      # p
    ggsave(paste0(conditions,"/", cond,".png"), plot = p)
  }
}

volcano_plots(conditions = "condition4_vs_condition24")