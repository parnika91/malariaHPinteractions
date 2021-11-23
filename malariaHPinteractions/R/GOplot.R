##script to do GO plots with GOplot

#install.packages("GOplot")
library(GOplot)
library(tidyverse)
library(ggplot2)
library(circlize)

GOfile <- read.csv("~/p_OG_topGO_BP_liver_overall_para_result_KS.txt")
GOfile <- read.delim("~/p_OG_topGO_BP_liver_core_para_result_KS.txt")
GOfile <- read.delim("~/Mmus_topGO_BP_liver_overall_host_result_KS.txt")
GOfile <- read.delim("~/Mmus_topGO_BP_liver_core_host_result_KS.txt")

expanded <- GOfile %>% 
  mutate(category = rep("BP", nrow(.))) %>% 
  mutate(Genes = strsplit(as.character(GenesForGOterm), ",")) %>% 
  unnest(Genes) %>% 
  select(category, GO.ID, Term, Significant, KS, Genes) %>% 
  rename(ID = GO.ID, term = Term, count = Significant, genes = Genes, adj_pval = KS) %>% 
  filter(adj_pval <= 0.05)

genes.of.interest <- tail(names(sort(table(expanded$genes))), n = 15)
chord <- chord_dat(expanded,
                   genes = genes.of.interest)

gene_col = rep("grey", length(genes.of.interest))
term_col = rand_color(length(unique(expanded$term)))
grid.col = setNames(c(term_col, gene_col), c(unique(expanded$term), genes.of.interest))
png("GOchord_liver_core_m_BP.png", height = 70, width = 70, units = "cm", res = 300)
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

GOfile_MF <- read.delim("~/p_OG_topGO_MF_liver_overall_para_result_KS.txt")
GOfile_MF <- read.delim("~/p_OG_topGO_MF_liver_core_para_result_KS.txt")
GOfile_MF <- read.delim("~/Mmus_topGO_MF_liver_overall_host_result_KS.txt")
GOfile_MF <- read.delim("~/Mmus_topGO_MF_liver_core_host_result_KS.txt")


GO_processed <- GOfile %>% 
  filter(KS <= 0.05) %>% 
  select(GO.ID, Term, Significant, KS, GenesForGOterm) %>% 
  filter(GenesForGOterm != "")
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
ggsave(plot = p, "Dotplot_GO_BP_liver_core_m.png")


GO_processed_MF <- GOfile_MF %>% 
  filter(KS <= 0.05) %>% 
  select(GO.ID, Term, Significant, KS, GenesForGOterm) %>% 
  filter(GenesForGOterm != "")
sorted_GO_MF <- GO_processed_MF[order(GO_processed_MF$KS),]

p <- ggplot(sorted_GO_MF) +
  geom_point(pch = 21,
             aes(x = Significant, y = reorder(Term, Significant), fill = KS, size = Significant)) + #aes(x = Significant, y = reorder(Term, Significant), fill = KS, size = Significant)
  labs(size = "Count", colour = "adj. p-value") +
  scale_fill_gradient(low = "yellow", 
                      high = "orange") +
  xlab("Number of genes in GO term") +
  ylab("GO term") +
  ggtitle("Number of genes in significantly enriched GO terms") +
  theme_bw()
ggsave(plot = p, "Dotplot_GO_MF_liver_core_m.png")
