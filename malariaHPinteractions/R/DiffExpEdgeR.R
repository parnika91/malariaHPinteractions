if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")
library(edgeR)
BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
install.packages("ggplot2")
install.packages("gridExtra")
install.packages("plotly")
suppressPackageStartupMessages(library("plotly"))

setwd("Documents/Data/Kai/")
x <- read.delim("Macrophage_coding_genes.txt")

y <- DGEList(counts=x)
Conditions <- read.csv("~/Documents/Data/Kai/Conditions.csv", stringsAsFactors=FALSE)
Conditions$Trt_time <- paste(Conditions$treatment, Conditions$time)


group <- c(rep(1, 94), rep(2, 48), rep(3, 5))
y <- DGEList(counts=x, group=group)

keep <- filterByExpr(y) 
y_keep <- y[keep, , keep.lib.sizes=FALSE]

# nrow(y) = 26961
# nrow(y_keep) = 7527

y_keep <- calcNormFactors(y_keep)
y_keep <- estimateDisp(y_keep)
et <- exactTest(y_keep)
topTags(et)

fit <- glmQLFit(y_keep)
tr <- glmLRT(fit, coef=2)
topTags(tr)

############ Convert IDs #############

library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- getBM(
  filters="ensembl_gene_id",
  attributes=c("ensembl_gene_id", "entrezgene"),
  values=ensembl.genes,
  mart=mart)

###################################

qlf <- glmQLFTest(fit, coef=2)
go <- goana(qlf, species="Mm")
topGO(go, sort="up")
keg <- kegga(qlf, species="Mm")
topKEGG(keg, sort="up")

logcpm <- cpm(y_keep, log=TRUE)

## diff_df between conditions 1 and 2 ##

diff_df <- DGE_macro[c("gene", "logFC", "FDR")]

# preview the dataset; data required for the plot
head(diff_df)
# add a grouping column; default value is "not significant"
diff_df["group"] <- "NotSignificant"

# for our plot, we want to highlight 
# FDR < 0.05 (significance level)
# Fold Change > 1.5

# change the grouping for the entries with significance but not a large enough Fold change
diff_df[which(diff_df['FDR'] < 0.05 & abs(diff_df['logFC']) < 1.5 ),"group"] <- "Significant"

# change the grouping for the entries a large enough Fold change but not a low enough p value
diff_df[which(diff_df['FDR'] > 0.05 & abs(diff_df['logFC']) > 1.5 ),"group"] <- "logFC"

# change the grouping for the entries with both significance and large enough logFC change
diff_df[which(diff_df['FDR'] < 0.05 & abs(diff_df['logFC']) > 1.5 ),"group"] <- "Significant&logFC"


# # Find and label the top peaks..
# top_peaks <- diff_df[with(diff_df, order(logFC, FDR)),][1:5,]
# top_peaks <- rbind(top_peaks, diff_df[with(diff_df, order(-logFC, FDR)),][1:5,])


# Add gene labels to the plot
# Single Gene Annotation example
# m <- diff_df[with(diff_df, order(logFC, FDR)),][1,]
# a <- list(
#   x = m[["logFC"]],
#   y = -log10(m[["FDR"]]),
#   text = m[["external_gene_name"]],
#   xref = "x",
#   yref = "y",
#   showarrow = TRUE,
#   arrowhead = 7,
#   ax = 20,
#   ay = -40
# )

# Add gene labels for all of the top genes we found
# here we are creating an empty list, and filling it with entries for each row in the dataframe
# each list entry is another list with named items that will be used by Plot.ly
# a <- list()
# for (i in seq_len(nrow(top_peaks))) {
#   m <- top_peaks[i, ]
#   a[[i]] <- list(
#     x = m[["logFC"]],
#     y = -log10(m[["FDR"]]),
#     text = m[["gene"]],
#     xref = "x",
#     yref = "y",
#     showarrow = TRUE,
#     arrowhead = 0.5,
#     ax = 20,
#     ay = -40
#   )
# }


# make the Plot.ly plot
# p <- plot_ly(data = diff_df, x = logFC, y = -log10(FDR), text = external_gene_name, mode = "markers", color = group) %>% 
#   layout(title ="Condition 1 vs 2") %>%
#   layout(annotations = a)
# p
diff_df <- as.data.frame(diff_df)
diff_df <-  diff_df[order(-abs(diff_df$logFC), diff_df$FDR),]

p <- ggplot(diff_df, aes(x = logFC, y = -log10(FDR), colour = group)) + 
  geom_point() + 
  ggtitle("Condition 1 vs 2") + 
  geom_text_repel(
    data = head(diff_df, 10),
    aes(label = gene),
    size = 2,
    colour = "black",
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )
ggsave("macrophage1vs2.png", plot = p)
