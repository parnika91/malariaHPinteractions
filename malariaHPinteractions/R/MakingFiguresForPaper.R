library(pheatmap)

# Pheatmap for log10 bipartite edges with clusters in breaks
load("/home/parnika/Documents/Data/b_edges_overlap_raw_numbers_matrix_all_datasets.RData")
load("/home/parnika/Documents/Data/anno.RData")
# remove "b_edges" from each rowname
rownames(mat1) <- sapply(rownames(mat1), function(x) paste(strsplit(x, split = "_")[[1]][1], 
                        strsplit(x, split = "_")[[1]][2], sep = '_'))
colnames(mat1) <- sapply(colnames(mat1), function(x) paste(strsplit(x, split = "_")[[1]][1], 
                                                           strsplit(x, split = "_")[[1]][2], sep = '_'))
pdf("pheatmap_b_edges_overlap_log10_clustered_divided.pdf")
pheatmap::pheatmap(log10(mat1+1), annotation_row = anno, fontsize = 8, main = "Log10 (Overlapping bipartite edges)",
                   cluster_cols = T, cluster_rows = T, cutree_cols = 3, cutree_rows = 3)
dev.off()

# relationship between #samples and #bipartite edges ---> DOES NOT REALLY WORK
o <- read.delim("~/Documents/Data/b_edges_overlap_signif_all_datasets.txt", stringsAsFactors=FALSE)
samplesize_bipsize <- data.frame(Dataset = o$Set2[1:27], OwnBipEdges = o$Set2_edges.[1:27])
samplesize_bipsize$Samples <- c(42, 48, 53, 23, 48, 80, 17, 143, 319, 
                                58, 118, 317, 8, 13, 53, 9, 13, 122, 
                                14, 14, 60, 14, 16, 26, 4, 6, 23)
plot(x = samplesize_bipsize$Samples, y = samplesize_bipsize$OwnBipEdges)
library(ggplot2)
gg <- ggplot(samplesize_bipsize, aes(x = log10(Samples), y = log10(OwnBipEdges))) + geom_point()
ggsave("gg.png")

# topGO stuff