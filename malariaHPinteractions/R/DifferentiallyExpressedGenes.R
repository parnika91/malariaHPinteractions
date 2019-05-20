# Find differentially expressed genes in my studies
# Parnika Mukherjee
# Dec 6, 2018

rm(list = ls())

# load packages
# EdgeR
# library("BiocManager", lib.loc="/localstorage/parnika/R/x86_64-pc-linux-gnu-library/3.5")
# BiocManager::install("edgeR")
library(edgeR)
library(tidyverse)

#######################  Demo  ###########################

# load dataset 
# RNAseq_counts <- read.delim(file = "http://supraHex.r-forge.r-project.org/RNAseq_counts.txt", header = T, row.names = 1)
# # df of 62757 obs and 6 variables. Has counts from 3 different cell lines with controls and treatment and their gene names
# 
# RNAseq_geneinfo <- read.delim(file = "http://supraHex.r-forge.r-project.org/RNAseq_geneinfo.txt", header = T, row.names = 1)
# # has gene type, gene name, description, entrez ID and gene length
# 
# # Create DGEList object (edgeR's container for RNA-seq count data)
# d_obj <- DGEList(counts = RNAseq_counts, genes = RNAseq_geneinfo)
# 
# str(d_obj)
# # Formal class 'DGEList' [package "edgeR"] with 1 slot
# # ..@ .Data:List of 3, named $counts, $samples and $genes
# # .. ..$ : int [1:62757, 1:6] 2104 0 2493 506 957 3 13134 4776 1563 2182 ...
# # .. .. ..- attr(*, "dimnames")=List of 2
# # .. .. .. ..$ : chr [1:62757] "ENSG00000000003" "ENSG00000000005" "ENSG00000000419" "ENSG00000000457" ...
# # .. .. .. ..$ : chr [1:6] "Cellline1_CON" "Cellline1_DEX" "Cellline2_CON" "Cellline2_DEX" ...
# # .. ..$ :'data.frame':	6 obs. of  3 variables:
# #   .. .. ..$ group       : Factor w/ 1 level "1": 1 1 1 1 1 1
# # .. .. ..$ lib.size    : num [1:6] 60786378 50620475 57843782 58620216 57637493 ...
# # .. .. ..$ norm.factors: num [1:6] 1 1 1 1 1 1
# # .. ..$ :'data.frame':	62757 obs. of  5 variables:
# #   .. .. ..$ GeneType    : Factor w/ 35 levels "3prime_overlapping_ncrna",..: 18 18 18 18 18 18 18 18 18 18 ...
# # .. .. ..$ GeneName    : Factor w/ 57341 levels "5_8S_rRNA","5S_rRNA",..: 54647 53904 14067 50048 8220 15887 9676 16299 16593 26291 ...
# # .. .. ..$ Description : Factor w/ 35130 levels "","1-acylglycerol-3-phosphate O-acyltransferase 1",..: 31210 31059 6193 27613 3665 7859 4753 8312 8886 18311 ...
# # .. .. ..$ EntrezGeneID: int [1:62757] 7105 64102 8813 57147 55732 2268 3075 2519 2729 4800 ...
# # .. .. ..$ GeneLength  : int [1:62757] 12882 15083 23688 44636 192073 23213 95626 16879 119629 27031 ...
# 
# # In edgeR it is recommended to remove features without at least 1 read/count per million (cpm) in n of the samples, where n is the size of the smallest
# # group of replicates. In this case, n = 3, but we set it to be 6 as the strictest filtering. <------------------------------- did not understand why 6
# 
# cpms <- edgeR::cpm(d_obj$counts)
# keep <- rowSums(cpms >= 1) >= 6
# d <- d_obj[keep,]
# summary(d)
# # Length Class      Mode   
# # counts  81522  -none-     numeric
# # samples     3  data.frame list   
# # genes       5  data.frame list   
# 
# # Reset the lib sizes
# d$samples$lib.size <- colSums(d$counts)
# 
# # Estimate normalisation factors
# d <- calcNormFactors(d)
# 
# # Define design matrix
# cnames <- colnames(d)
# cellline <- gsub("_.*", "", cnames, perl = T)
# treatment <- gsub(".*_", "", cnames, perl = T)
# targets <- data.frame(sample = cnames, cellline = cellline, treatment = treatment)
# design <- model.matrix(~ cellline + treatment, targets)
# design
# 
# # Inspect relationship between samples
# plotMDS(d, labels = colnames(d), col = c("red", "darkgreen", "blue")[factor(targets$cellline)], xlim = c(-2,2), ylim = c(-2,2))
# 
# # This plot clearly shows the variances between celllines, calling for paired design test
# # See dispersion and GLM fitting
# # Estimate overall dispersion to get an idea of the overall level of biological variability
# d <- estimateGLMCommonDisp(d, design, verbose = T)
# summary(d)
# #                   Length Class      Mode   
# # counts            81522  -none-     numeric
# # samples               3  data.frame list   
# # genes                 5  data.frame list   
# # common.dispersion     1  -none-     numeric
# # AveLogCPM         13587  -none-     numeric
# 
# # Estimate dispersion values relative to design matrix
# d2 <- estimateGLMTrendedDisp(d, design)
# d2 <- estimateGLMTagwiseDisp(d2, design)
# 
# # Given design matrix and dispersion estimates, fit a GLM to each feature
# f <- glmFit(d2, design)
# 
# # PErform likelihood test
# # Specify difference of interest: DEX vs CON
# contrasts <- rbind(c(0,0,0,1))
# 
# # Prepare results
# logFC <- matrix(nrow=nrow(d2), ncol = 2)
# PValue <- matrix(nrow=nrow(d2), ncol = 2)
# FDR <- matrix(nrow = nrow(d2), ncol = 2)
# tmp <- c("DEX", "CON")
# colnames(logFC) <- paste0(tmp, "_logFC")
# colnames(PValue) <- paste0(tmp, "_PValue")
# colnames(FDR) <-  paste0(tmp, "_FDR")
# rownames(logFC) <-  rownames(PValue) <-  rownames(FDR) <- rownames(d2)
# 
# # Perform the test, calculating p-values and FDR
# for(i in 1:nrow(contrasts))
# {
#   lrt <- glmLRT(f, contrast = contrasts[i,])
#   tt <- topTags(lrt, n = nrow(d2), adjust.method = "BH", sort.by = "none")
#   logFC[,i] <- tt$table$logFC
#   PValue[,i] <- tt$table$PValue
#   FDR <- tt$table$FDR
# }
# 
# # MA plots for RNA-seq data
# lrt <- glmLRT(f, contrast = contrasts[1,])
# # lrt <-  glmLRT(f, coef = 4)
# summary(de <- decideTestsDGE(lrt, adjust.method = "BH", p.value = 0.05))
# 
# detags <- rownames(d2)[as.logical(de)]
# plotSmear(lrt, de.tags = detags)
# 
# # Output edgeR results
# # log counts per million
# cpms <- edgeR::cpm(d2, log = T, prior.count = 2)
# colnames(cpms) <- paste0(colnames(d2), '_CPM')
# 
# # log ratio between treatment and control
# # odd_indexes <- seq(1, ncol(cpms), 2)
# # even_indexes <- seq(2, ncol(cpms), 2)
# # logFC_cpm <- cpms[,even_indexes] - cpms[,odd_indexes]
# # colnames(logFC_cpm) <- gsub("_CPM", "_CPM_logFC", colnames(logFC), perl = T)
# 
# # Write out the file
# out <- data.frame(EnsemblGeneID = rownames(d2$genes), d2$genes, cpms, logFC, PValue, FDR)
# write.table(out, file = "RNAseq_edgeR.txt", col.names = T, row.names = F, sep = '\t', quote = F)


#####################  My study  #########################

# load dataset
study <- t(read.delim("SRAdb/Output/ERP106451/ERP106451.txt"))
genes <- as_tibble(rownames(study))

# Create DGEList object (edgeR's container for RNA-seq count data)
d_obj <- DGEList(counts = study, genes = genes)
summary(d_obj)
#         Length Class      Mode   
# counts  116221 -none-     numeric
# samples      3 data.frame list   
# genes        1 data.frame list 

str(d_obj)
# Formal class 'DGEList' [package "edgeR"] with 1 slot
# ..@ .Data:List of 3
# .. ..$ : int [1:16603, 1:7] 135 752 54188 0 0 68 0 4226 109 0 ...
# .. .. ..- attr(*, "dimnames")=List of 2
# .. .. .. ..$ : chr [1:16603] "A1BG" "A1CF" "A2M" "A2ML1" ...
# .. .. .. ..$ : chr [1:7] "ERR1744069_ERP020067" "ERR1744068_ERP020067" "ERR1744072_ERP020067" "ERR1744067_ERP020067" ...
# .. ..$ :'data.frame':	7 obs. of  3 variables:
#   .. .. ..$ group       : Factor w/ 1 level "1": 1 1 1 1 1 1 1
# .. .. ..$ lib.size    : num [1:7] 1.68e+07 1.29e+07 6.88e+06 1.21e+08 4.03e+06 ...
# .. .. ..$ norm.factors: num [1:7] 1 1 1 1 1 1 1
# .. ..$ :'data.frame':	16603 obs. of  1 variable:
#   .. .. ..$ value: chr [1:16603] "A1BG" "A1CF" "A2M" "A2ML1" ...

# In edgeR it is recommended to remove features without at least 1 read/count per million (cpm) in n of the samples, where n is the size of the smallest
# group of replicates.
cpms <- edgeR::cpm(d_obj)
summary(cpms)
# ERR1744069_ERP020067 ERR1744068_ERP020067 ERR1744072_ERP020067 ERR1744067_ERP020067 ERR1744071_ERP020067 ERR1744066_ERP020067 ERR1744070_ERP020067
# Min.   :    0.000    Min.   :    0.00     Min.   :     0.0     Min.   :    0.00     Min.   :     0.0     Min.   :    0.00     Min.   :    0.000   
# 1st Qu.:    0.000    1st Qu.:    0.00     1st Qu.:     0.0     1st Qu.:    0.00     1st Qu.:     0.0     1st Qu.:    0.00     1st Qu.:    0.000   
# Median :    1.663    Median :    0.00     Median :     0.0     Median :    0.63     Median :     0.0     Median :    0.95     Median :    0.746   
# Mean   :   60.230    Mean   :   60.23     Mean   :    60.2     Mean   :   60.23     Mean   :    60.2     Mean   :   60.23     Mean   :   60.230   
# 3rd Qu.:   20.963    3rd Qu.:    0.15     3rd Qu.:     0.0     3rd Qu.:   12.37     3rd Qu.:     0.0     3rd Qu.:   15.64     3rd Qu.:   18.413   
# Max.   :28776.684    Max.   :68706.15     Max.   :483529.4     Max.   :64778.08     Max.   :570980.8     Max.   :61865.89     Max.   :26849.363 

countCheck <- cpms > 1
keep <- which(rowSums(countCheck) >= 2)
d_obj <- d_obj[keep,]
str(d_obj)
# Formal class 'DGEList' [package "edgeR"] with 1 slot
# ..@ .Data:List of 3
# .. ..$ : int [1:527, 1:7] 1438 1710 8988 126364 107636 4577 20612 513 286 226 ...
# .. .. ..- attr(*, "dimnames")=List of 2
# .. .. .. ..$ : chr [1:527] "ABCB6" "ABHD4" "ACADVL" "ACTB" ...
# .. .. .. ..$ : chr [1:7] "ERR1744069_ERP020067" "ERR1744068_ERP020067" "ERR1744072_ERP020067" "ERR1744067_ERP020067" ...
# .. ..$ :'data.frame':	7 obs. of  3 variables:
#   .. .. ..$ group       : Factor w/ 1 level "1": 1 1 1 1 1 1 1
# .. .. ..$ lib.size    : num [1:7] 1.68e+07 1.29e+07 6.88e+06 1.21e+08 4.03e+06 ...
# .. .. ..$ norm.factors: num [1:7] 1 1 1 1 1 1 1
# .. ..$ :'data.frame':	527 obs. of  1 variable:
#   .. .. ..$ value: chr [1:527] "ABCB6" "ABHD4" "ACADVL" "ACTB" ...

summary(cpm(d_obj))
# ERR1744069_ERP020067 ERR1744068_ERP020067 ERR1744072_ERP020067 ERR1744067_ERP020067 ERR1744071_ERP020067 ERR1744066_ERP020067 ERR1744070_ERP020067
# Min.   :    3.029    Min.   :    1.004    Min.   :     1.16    Min.   :    1.066    Min.   :     1.24    Min.   :    1.821    Min.   :    1.493   
# 1st Qu.:   59.236    1st Qu.:    6.874    1st Qu.:     5.89    1st Qu.:   50.392    1st Qu.:     6.20    1st Qu.:   61.034    1st Qu.:   64.632   
# Median :  148.164    Median :  118.093    Median :    21.51    Median :  135.098    Median :    18.60    Median :  156.260    Median :  181.020   
# Mean   :  483.709    Mean   :  558.665    Mean   :   625.92    Mean   :  513.091    Mean   :   467.93    Mean   :  549.480    Mean   :  540.718   
# 3rd Qu.:  411.414    3rd Qu.:  437.192    3rd Qu.:    86.49    3rd Qu.:  462.831    3rd Qu.:    68.45    3rd Qu.:  419.205    3rd Qu.:  440.045   
# Max.   :18072.982    Max.   :22352.234    Max.   :202505.08    Max.   :22404.190    Max.   :118470.71    Max.   :18595.477    Max.   :18639.288

# Normalisation
d_obj <- calcNormFactors(d_obj, method = "TMM")

# Exploration
plotMDS(d_obj) # Clear distinction between samples! Similar to before, but not exactly the same. Same data points are together.

# Set up model
replicates <- colnames(study)
groups <- c("1","4","2","3","2","3","1")
design <- model.matrix(~ 0 + groups)

# estimate dispersion
d_obj1 <- estimateGLMCommonDisp(d_obj, design = design)
d_obj1 <- estimateGLMTrendedDisp(d_obj1, design = design)
d_obj1 <- estimateGLMTagwiseDisp(d_obj1, design = design)

plotBCV(d_obj1)

# Differential expression
fit <- glmFit(d_obj1, design)
lrt <- glmLRT(fit, coef = 4)

edgeR_res <- topTags(lrt)

deGenes <- decideTestsDGE(lrt, p = 0.001)
deGenes <- rownames(lrt)[as.logical(deGenes)]
plotSmear(lrt, de.tags = deGenes)
abline(h = c(-1,1), col = 2)
