##script to do GO plots with GOplot

#install.packages("GOplot")
library(GOplot)
library(tidyverse)
library(ggplot2)
library(circlize)

p_ov_BP <- read.csv("~/p_OG_topGO_BP_liver_overall_para_result_KS.txt")

#data(EC)

# EC$genelist
# EC$genes
# EC$process
# EC$david
# circ <- circle_dat(EC$david, EC$genelist)
# chord_EC <- chord_dat(circ, EC$genes, EC$process)
# GOChord(chord_EC, gene.order = 'none')

expanded_p_ov_BP <- p_ov_BP %>% 
  mutate(category = rep("BP", nrow(.))) %>% 
  mutate(Genes = strsplit(as.character(GenesForGOterm), ",")) %>% 
  unnest(Genes) %>% 
  select(category, GO.ID, Term, Significant, KS, Genes) %>% 
  rename(ID = GO.ID, term = Term, count = Significant, genes = Genes, adj_pval = KS) %>% 
  filter(adj_pval <= 0.05)

genes.of.interest <- tail(names(sort(table(expanded_p_ov_BP$genes))), n = 15)
chord <- chord_dat(expanded_p_ov_BP,
                   genes = genes.of.interest)
# #chord <- chord_dat(data = circ, genes = genes, process = process)
# which(rowMeans(chord) == 0)
# which(colMeans(chord) == 0)
# 
# #Eliminate genes related to no processes
# if(length(which(rowMeans(chord) == 0) != 0)){
#   paste("one or more genes are not related any of the processes")
#   chord <- chord[-which(rowMeans(chord) == 0),]
# }else {
#   paste("all genes are related to at least one process")
# }
# 
# #Eliminate processes genes that contain no genes
# if(length(which(colMeans(chord) == 0) != 0)){
#   paste("one or more processes are not related any of the selected genes")
#   chord <- chord[,-which(colMeans(chord) == 0)]
# }else {
#   paste("all processes have at least one related gene")
# }
# 
# png("GOchord_liver_ov_p_BP.png", height = 50, width = 50, units = "cm", res = 300)
# GOChord(chord, 
#         gene.order = 'none', 
#         space = 0.02, 
#         gene.size = 5, 
#         gene.space = 0.25, 
#         border.size = 0,
#         process.label = 10#,
#         #ribbon.col = "blue"
#         )
# dev.off()
# 
# install.packages("circlize")
# library(circlize)
# 
# numbers <- sample(c(1:1000), 100, replace = T)
# data <- matrix( numbers, ncol=5)
# rownames(data) <- paste0("orig-", seq(1,20))
# colnames(data) <- paste0("dest-", seq(1,5))

# Load the circlize library
library(circlize)

# Make the circular plot
gene_col = rep("grey", length(genes.of.interest))
term_col = rand_color(length(unique(expanded_p_ov_BP$term)))
grid.col = setNames(c(term_col, gene_col), c(unique(expanded_p_ov_BP$term), genes.of.interest))
png("GOchord_liver_ov_p_BP.png", height = 70, width = 70, units = "cm", res = 300)
circos.clear()
chordDiagram(t(chord), 
             transparency = 0.5, 
             grid.col = grid.col, 
             annotationTrack = c("grid"))
circos.track(track.index = 1, panel.fun = function(x, y)
{
  circos.text(CELL_META$xcenter, 
              CELL_META$ylim[1], 
              CELL_META$sector.index, cex = 1,
              facing = "clockwise", 
              niceFacing = TRUE, 
              adj = c(-0.1, 0.5))
  },
bg.border = 0.2)
dev.off()




### dot plot
GO_processed <- p_ov_BP %>% 
  filter(KS <= 0.05) %>% 
  select(GO.ID, Term, Significant, KS, GenesForGOterm)
sorted_GO <- GO_processed[order(GO_processed$KS),]

p <- ggplot(sorted_GO) +
  geom_point(pch = 21,
             aes(x = Significant, y = reorder(Term, Significant), fill = KS, size = Significant)) + #aes(x = Significant, y = reorder(Term, Significant), fill = KS, size = Significant)
  labs(size = "Count", colour = "adj. p-value") +
  scale_fill_gradient(low = "blue", 
                        high = "red") +
  xlab("Number of genes in GO term") +
  ylab("GO term") +
  ggtitle("Number of genes in significantly enriched GO terms") +
  theme_bw()
ggsave(plot = p, "Dotplot_GO_count_pval.png")
