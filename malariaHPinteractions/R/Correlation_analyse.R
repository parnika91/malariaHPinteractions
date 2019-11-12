# This script is to accummulate all the outer files
# and analyse them: first sum them up all,
# then find the p-values
# third, get gene pairs with p-values == 0
# fourth, get this from the other study and finc common pairs

library(WGCNA)
library(reshape2)
library(ggplot2)
library(dplyr)

sum <- matrix(rep(0, 17991*17991), nrow = 17991)

files <- grep(pattern = "outer_", list.files())
file.names <- list.files()[files]

for(i in 1:length(file.names))
{
  print(i)
  #setwd("/short/uv67/pm4655/cor/")
  if(grepl(pattern = ".RData", file.names[i]))
  {
    load(file.names[i])
    sum <- sum + outer 
  }
}

all_genes2 <- rownames(sum)

load("/SAN/Plasmo_compare/SRAdb/Output/ERP106451/cor/ERP106451.ortho.data.allruns.RData")
study <- t(ERP106451.ortho.data.allruns)
ori_cor <- cor(study, use = "pairwise.complete.obs")
ori_cor[lower.tri(ori_cor, diag = T)] <- NA
cor_melt <- melt(ori_cor)
colnames(cor_melt) <- c("gene1", "gene2", "cor")

pval <- sum
pval[lower.tri(pval, diag = T)] <- NA
pval_melt <- melt(pval)
colnames(pval_melt) <- c("gene1", "gene2", "permute_score")

pval_cor <- cbind(cor_melt, pval_melt$permute_score)
colnames(pval_cor)[4] <- "permute_score"

#Do this step in next analysis:
pval_cor_na.omit <- na.omit(pval_cor)
pval0 <- pval_cor_na.omit[pval_cor_na.omit$permute_score==0,]

hg_col1 <- as.character(pval0[grep(pattern = "h_OG", pval0$gene1),1])
unique_hg_col1 <- unique(hg_col1)          
hg_col2 <- as.character(pval0[grep(pattern = "h_OG", pval0$gene2),2])
unique_hg_col2 <- unique(hg_col2)
unique_hg <- unique(c(unique_hg_col1), c(unique_hg_col2))

pg_col1 <- as.character(pval0[grep(pattern = "p_OG", pval0$gene1),1])
unique_pg_col1 <- unique(pg_col1)  
pg_col2 <- as.character(pval0[grep(pattern = "p_OG", pval0$gene2),2])
unique_pg_col2 <- unique(pg_col2)
unique_pg <- unique(c(unique_pg_col1), c(unique_pg_col2))

bipartite <- pval0[(grepl(pattern = "h_OG", pval0$gene1) & grepl(pattern = "p_OG", pval0$gene2)) | (grepl(pattern = "p_OG", pval0$gene1) & grepl(pattern = "h_OG", pval0$gene2)),]

b_hg_col1 <- as.character(bipartite[grep(pattern = "h_OG", bipartite$gene1),1])
b_unique_hg_col1 <- unique(b_hg_col1)
b_hg_col2 <- as.character(bipartite[grep(pattern = "h_OG", bipartite$gene2),2])
b_unique_hg_col2 <- unique(b_hg_col2)
b_unique_hg <- unique(c(b_unique_hg_col1), c(b_unique_hg_col2))

b_pg_col1 <- as.character(bipartite[grep(pattern = "p_OG", bipartite$gene1),1])
b_unique_pg_col1 <- unique(b_pg_col1)
b_pg_col2 <- as.character(bipartite[grep(pattern = "p_OG", bipartite$gene2),2])
b_unique_pg_col2 <- unique(b_pg_col2)
b_unique_pg <- unique(c(b_unique_pg_col1), c(b_unique_pg_col2))

# bip <- rbind(ERP106451_bipartite[,1:2], SRP118996_bipartite[,1:2])
# bip[] <- lapply(bip, as.character)
# bip_dupl <- unique(bip[(duplicated(bip) | duplicated(bip, fromLast = TRUE)),])

load("ERP106451_bipartite.RData")
load("ERP106451_allruns_bipartite.RData")
load("ERP106451_stringent15_bipartite.RData")

ERP106451_allruns_bip_names_concat <- apply(ERP106451_allruns_bipartite[,c("gene1","gene2") ] , 1 , paste , collapse = "_" )
ERP106451_stringent15_bip_names_concat <- apply(ERP106451_stringent15_bipartite[,c("gene1","gene2") ] , 1 , paste , collapse = "_" )
ERP106451_bip_names_concat <- apply(ERP106451_bipartite[,c("gene1","gene2") ] , 1 , paste , collapse = "_" )

ERP106451_intersection_all <- list(allruns = ERP106451_allruns_bip_names_concat, stringent = ERP106451_stringent15_bip_names_concat, normal = ERP106451_bip_names_concat)

library(VennDiagram)
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")
venn.diagram(
        x = ERP106451_intersection_all,
        category.names = c("allruns" , "stringent" , "normal"),
        filename = 'ERP106451_venn_diagram.png',
        output=TRUE,
        
        # Output features
        imagetype="png" ,
        height = 480 , 
        width = 480 , 
        resolution = 300,
        compression = "lzw",
        
        # Circles
        lwd = 2,
        lty = 'blank',
        fill = myCol,
        
        # Numbers
        cex = .3,
        fontfamily = "sans",
        
        # Set names
        cat.cex = 0.6,
        cat.default.pos = "outer",
        cat.pos = c(-27, 27, 135),
        cat.dist = c(0.055, 0.055, 0.085),
        cat.fontfamily = "sans",
        rotation = 1
)


#library(grid)
#png(file="SRP118996_intersect.png", width = 15, height = 15, units = "cm", res = 450) # or other device; , onefile = F for pdf()
#upset(fromList(SRP118996_upset), sets = names(SRP118996_upset), order.by = "freq",  mainbar.y.label = "Genes pairs in intersections", sets.x.label = "Genes pairs per dataset", empty.intersections = "on", main.bar.color = "darkblue", sets.bar.color=c("maroon"), matrix.color="darkgreen", text.scale = c(1.2, 1, 1, 1, 0.8, 0.8))
#grid.text("SRP118996 datasets",x = 0.65, y=0.95, gp=gpar(fontsize=10))
#dev.off()

hg_common <- unique(bip_dupl$gene1)
pg_common <- unique(bip_dupl$gene2)

# join two studies correlation values #

bip_cor1 <- inner_join(bip_dupl, SRP032775_bipartite, by = c("gene1", "gene2"))
bip_cor12 <- inner_join(bip_cor1, SRP118996_bipartite, by = c("gene1", "gene2"))
bip_cor12 <- bip_cor12[,c(1,2,3,5)]

common_all <- inner_join(SRP032775_pval0, SRP118996_pval0, by = c("gene1", "gene2"))
common_hh <- common_all[which(grepl(pattern= "h_OG", common_all$gene1) & grepl(pattern = "h_OG", common_all$gene2)),]
common_pp <- common_all[which(grepl(pattern= "p_OG", common_all$gene1) & grepl(pattern = "p_OG", common_all$gene2)),]

png("SRP118996_SRP032775_permsc_vs_cor_paraprop.png")
plot(x = pval_cor_na.omit$cor, y = pval_cor_na.omit$permute_score, 
     main = "SRP032775 and SRP118996: Distribution of permutation score and corr coef", 
     xlab = "Correlation coefficient [-1, 1]", ylab = "Permutation score of gene pair",
     colour = "#E7B800", cex = 0.2)
#plot(density(pval_cor_na.omit$cor))
dev.off()

pv = pval_cor_na.omit[sample(nrow(pval_cor_na.omit), size = nrow(pval_cor_na.omit)/1000),]

gg_pvsc <- ggplot(pval_cor_na.omit, aes(x = cor, y = permute_score)) + geom_point(alpha = 0.7) + theme_bw() +
             ggtitle("ERP106451_allruns: Distribution of permutation score and correlation coef") +
                 xlab("Correlation coefficient [-1, 1]") + ylab("Permutation score of gene pair")

ggsave(plot = gg_pvsc, filename = "ERP106451_allruns.png")

### to plot perm_sc == 0 vs number of perms

sum2 <- matrix(rep(0, 17991*17991), nrow = 17991)
sum2[lower.tri(sum2, diag = T)] <- NA
df <- data.frame(Perms = seq(10000, 100000, 10000))
a = 1
for(i in 1:length(file.names))
{
  print(i)
  if(grepl(pattern = ".RData", file.names[i]))
  {
  load(file.names[i])
  #outer[lower.tri(outer, diag = T)] <- NA
  sum2 <- outer + sum2
  #sum2[lower.tri(sum2, diag = T)] <- NA
  sum2_melt <- melt(sum2)
  
  colnames(sum2_melt) <- c("gene1", "gene2", "permute_score")

  pval_cor <- cbind(cor_melt, sum2_melt$permute_score)
  colnames(pval_cor)[4] <- "permute_score"

  pval_cor_na.omit <- na.omit(pval_cor)
  pval0 <- pval_cor_na.omit[pval_cor_na.omit$permute_score==0,]

  cor0.9 <- nrow(pval0[abs(pval0$cor) >= 0.9,])
  cor0.8 <- nrow(pval0[(abs(pval0$cor) >= 0.8 & abs(pval0$cor) < 0.9),])
  cor0.7 <- nrow(pval0[(abs(pval0$cor) >= 0.7 & abs(pval0$cor) < 0.8),])
  cor0.6 <- nrow(pval0[(abs(pval0$cor) >= 0.6 & abs(pval0$cor) < 0.7),])
  cor0.5 <- nrow(pval0[(abs(pval0$cor) >= 0.5 & abs(pval0$cor) < 0.6),])

  df[a,2] <- cor0.9
  df[a,3] <- cor0.8
  df[a,4] <- cor0.7
  df[a,5] <- cor0.6
  df[a,6] <- cor0.5
  
  a = a + 1}
}

save(df, file = "SRP118996_hp_12Sep2019df.RData")

colnames(df) <- c("Permutations", "Cor0.9", "Cor0.8", "Cor0.7", "Cor0.6", "Cor0.5")
df_melt <- melt(df, id.vars = "Permutations")

png("SRP032775_hp_100k_12Sep2019.png")
ggplot(df_melt, aes(x = Permutations, y = log10(value), colour = variable)) + geom_line() +
       scale_x_continuous(breaks = seq(10000, 100000, 20000), limits = c(10000, 100000)) + 
       ggtitle("SRP032775: #edges with perm_score zero; para read prop >5%") + 
       ylab("log10 of number of edges with perm_score == 0")
dev.off()


# For every study, read the outer files, melt them, read the different kinds of edges and nodes
# count these and put then in the table
# Things to count:
# 1. Number of hh edges, pp edges, all edges and bipartite edges
# 2. Number of h nodes, p nodes, h nodes in bipartite edges and p nodes in bipartite edges
# The idea is to pool the numbers together and plot as different facets
# So, for a single study, a matrix will have 
# 10 row to represent 100000 permutations
# cols for 1 and 2., 8 in total
# In ggplot, I would reshape then by the edge type and node type -> 2 groups of facets
# 0.9 10k each is a separate matrix. then rbind them.
# 0.8 10k
# 0.7 10k
# 0.6 10k
# 0.5 10k


number_of_edges_and_nodes <- function(pcc)
{
  sum2 <- matrix(rep(0, 17991*17991), nrow = 17991)
  sum2[lower.tri(sum2, diag = T)] <- NA
  df <- data.frame(Perms = seq(10000, 100000, 10000))
  a = 1
  
  for(i in 1:length(file.names))
  {
    print(i)
    if(grepl(pattern = ".RData", file.names[i]))
    {
      load(file.names[i])
      sum2 <- outer + sum2
      sum2_melt <- melt(sum2)
      
      colnames(sum2_melt) <- c("gene1", "gene2", "permute_score")
      
      pval_cor <- cbind(cor_melt, sum2_melt$permute_score)
      colnames(pval_cor)[4] <- "permute_score"
      
      pval_cor_na.omit <- na.omit(pval_cor)
      pval0 <- pval_cor_na.omit[pval_cor_na.omit$permute_score==0,]
      
      if(pcc >= 0.9){ r1 = 0.9; r2 = 1.0 }
      if(pcc >= 0.8 & pcc < 0.9){ r1 = 0.8; r2 = 0.9 }
      if(pcc >= 0.7 & pcc < 0.8){ r1 = 0.7; r2 = 0.8 }
      if(pcc >= 0.6 & pcc < 0.7){ r1 = 0.6; r2 = 0.7 }
      if(pcc >= 0.5 & pcc < 0.6){ r1 = 0.5; r2 = 0.6 }
      
      all.edges <- pval0[(abs(pval0$cor) >= r1 & abs(pval0$cor) < r2),]
      all.edges.number <- nrow(all.edges)
      hh.edges <- all.edges[(grepl(pattern = "h_OG", as.character(all.edges$gene1)) & grepl(pattern = "h_OG", as.character(all.edges$gene2))),]
      hh.edges.number <- nrow(hh.edges)
      pp.edges <- all.edges[(grepl(pattern = "p_OG", as.character(all.edges$gene1)) & grepl(pattern = "p_OG", as.character(all.edges$gene2))),]
      pp.edges.number <- nrow(pp.edges)
      hp.edges <- all.edges[((grepl(pattern = "h_OG", as.character(all.edges$gene1)) & grepl(pattern = "p_OG", as.character(all.edges$gene2))) | (grepl(pattern = "p_OG", as.character(all.edges$gene1)) & grepl(pattern = "h_OG", as.character(all.edges$gene2)))),]
      hp.edges.number <- nrow(hp.edges)
      
      hh.h.nodes <- unique(c(unique(as.character(hh.edges$gene1))), c(unique(as.character(hh.edges$gene2))))
      pp.p.nodes <- unique(c(unique(as.character(pp.edges$gene1))), c(unique(as.character(pp.edges$gene2))))
      
      hp.h.nodes.col1 <- unique(as.character(hp.edges$gene1))
      # if(length(hp.h.nodes.col1) > 0 & length(hp.h.nodes.col2) > 0)
      #   hp.h.nodes <- unique(hp.h.nodes.col1, hp.h.nodes.col2)
      # if(length(hp.h.nodes.col1) > 0 & length(hp.h.nodes.col2) == 0)
      #   hp.h.nodes <- unique(hp.h.nodes.col1)
      # if(length(hp.h.nodes.col2) > 0 & length(hp.h.nodes.col1) == 0)
      hp.h.nodes <- length(hp.h.nodes.col1)
      
      hp.p.nodes.col1 <- unique(as.character(hp.edges$gene2))
      # hp.p.nodes.col2 <- unique(as.character(grep(pattern = "p_OG",hp.edges$gene2)))
      # if(length(hp.p.nodes.col1) > 0 & length(hp.p.nodes.col2) > 0)
      #   hp.p.nodes <- unique(hp.p.nodes.col1, hp.p.nodes.col2)
      # if(length(hp.p.nodes.col1) > 0 & length(hp.p.nodes.col2) == 0)
      #   hp.p.nodes <- unique(hp.p.nodes.col1)
      # if(length(hp.p.nodes.col2) > 0 & length(hp.p.nodes.col1) == 0)
      hp.p.nodes <- length(hp.p.nodes.col1)
      
      df[a,2] <- pcc
      df[a,3] <- all.edges.number
      df[a,4] <- hh.edges.number
      df[a,5] <- pp.edges.number
      df[a,6] <- hp.edges.number
      df[a,7] <- length(hh.h.nodes)
      df[a,8] <- length(pp.p.nodes)
      df[a,9] <- hp.h.nodes
      df[a,10] <- hp.p.nodes
      
      a = a+1
    }
  }
  return(df)
}
nodes_and_edges <- rbind(number_of_edges_and_nodes(0.9), 
                             number_of_edges_and_nodes(0.8), 
                             number_of_edges_and_nodes(0.7),
                             number_of_edges_and_nodes(0.6),
                             number_of_edges_and_nodes(0.5))

save(nodes_and_edges, file = "ERP106451_allruns_nodes_and_edges.RData")
colnames(nodes_and_edges) <- c("Perms", "cor", "all.edges", "hh.edges", "pp.edges", "hp.edges", "hh.h.nodes", "pp.p.nodes", "hp.h.nodes", "hp.p.nodes")
df_individual_melt <- melt(nodes_and_edges, id = c("Perms", "cor"))
colnames(df_individual_melt) <- c("Perms", "cor", "edge.node.type", "count")

edge.node.plot <- ggplot(df_individual_melt, aes(x = Perms, y = count, colour = factor(cor))) +
                         geom_line() + scale_x_continuous(breaks = seq(10000, 100000, 25000), limits = c(10000, 100000))+
			 theme_bw() +
                         facet_wrap(. ~ edge.node.type, scales="free_y") + ggtitle("ERP106451_allruns: Number of edges/nodes with permutation score 0")

ggsave(plot = edge.node.plot, filename = "ERP106451_allruns_no_para_prop_edge.node.plot.png")

edges.df <- df_individual_melt[grep(pattern = "edge", df_individual_melt$edge.node.type),]
nodes.df <- df_individual_melt[grep(pattern = "node", df_individual_melt$edge.node.type),]

edge.plot <- ggplot(edges.df, aes(x = Perms, y = log10(count), colour = factor(cor))) +
                         geom_line() + scale_x_continuous(breaks = seq(10000, 100000, 25000), limits = c(10000, 100000))+
                         theme_bw() +
                         facet_wrap(. ~ edge.node.type, scales="free_y") + ggtitle("ERP106451_allruns: Number of edges with permutation score 0") +
                         xlab("Permutations")

ggsave(plot = edge.plot, filename = "ERP106451_allruns_no_para_prop_edge.plot.png")

node.plot <- ggplot(nodes.df, aes(x = Perms, y = log10(count), colour = factor(cor))) +               
                         geom_line() + scale_x_continuous(breaks = seq(10000, 100000, 25000), limits = c(10000, 100000))+
                         theme_bw() +
                         facet_wrap(. ~ edge.node.type, scales="free_y") + ggtitle("ERP106451_allruns: Number of nodes with permutation score 0") +
                         xlab("Permutations")
 

ggsave(plot = node.plot, filename = "ERP106451_allruns_no_para_prop_node.plot.png")












# find common edges and nodes between the two studies at every 10000 step

# https://stackoverflow.com/questions/5577221/how-can-i-load-an-object-into-a-variable-name-that-i-specify-from-an-r-data-file
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

common_edges_and_nodes <- function(pcc)
{
  # monkey init
  mon_sum <- matrix(rep(0, 17991*17991), nrow = 17991)
  mon_sum[lower.tri(mon_sum, diag = T)] <- NA
  # human init
  hum_sum <- matrix(rep(0, 17991*17991), nrow = 17991)
  hum_sum[lower.tri(hum_sum, diag = T)] <- NA
  
  df <- data.frame(Permutations = seq(10000, 100000, 10000))
  a = 1

  for(i in 1:length(mon_file.names))
  {
    print(i)
    
    mon_outer = loadRData(paste0(mon_prefix, mon_file.names[i], collapse = ''))
    mon_sum <- mon_outer + mon_sum
    mon_sum_melt <- melt(mon_sum)
    colnames(mon_sum_melt) <- c("gene1", "gene2", "permute_score")
    mon_pval_cor <- cbind(mon_cor_melt, mon_sum_melt$permute_score)
    colnames(mon_pval_cor)[4] <- "permute_score"
    mon_pval_cor_na.omit <- na.omit(mon_pval_cor)
    mon_pval0 <- mon_pval_cor_na.omit[mon_pval_cor_na.omit$permute_score==0,]
    ####
    
    hum_outer = loadRData(paste0(hum_prefix, hum_file.names[i], collapse = ''))
    hum_sum <- hum_outer + hum_sum
    hum_sum_melt <- melt(hum_sum)
    colnames(hum_sum_melt) <- c("gene1", "gene2", "permute_score")
    hum_pval_cor <- cbind(hum_cor_melt, hum_sum_melt$permute_score)
    colnames(hum_pval_cor)[4] <- "permute_score"
    hum_pval_cor_na.omit <- na.omit(hum_pval_cor)
    hum_pval0 <- hum_pval_cor_na.omit[hum_pval_cor_na.omit$permute_score==0,]
    ####
    
    if(pcc >= 0.9){ r1 = 0.9; r2 = 1.0 }
    if(pcc >= 0.8 & pcc < 0.9){ r1 = 0.8; r2 = 0.9 }
    if(pcc >= 0.7 & pcc < 0.8){ r1 = 0.7; r2 = 0.8 }
    if(pcc >= 0.6 & pcc < 0.7){ r1 = 0.6; r2 = 0.7 }
    if(pcc >= 0.5 & pcc < 0.6){ r1 = 0.5; r2 = 0.6 }
    
    ####
    mon.all.edges <- mon_pval0[(abs(mon_pval0$cor) >= r1 & abs(mon_pval0$cor) < r2),]
    mon.hh.edges <- mon.all.edges[(grepl(pattern = "h_OG", as.character(mon.all.edges$gene1)) & grepl(pattern = "h_OG", as.character(mon.all.edges$gene2))),]
    mon.pp.edges <- mon.all.edges[(grepl(pattern = "p_OG", as.character(mon.all.edges$gene1)) & grepl(pattern = "p_OG", as.character(mon.all.edges$gene2))),]
    mon.hp.edges <- mon.all.edges[((grepl(pattern = "h_OG", as.character(mon.all.edges$gene1)) & grepl(pattern = "p_OG", as.character(mon.all.edges$gene2))) | (grepl(pattern = "p_OG", as.character(mon.all.edges$gene1)) & grepl(pattern = "h_OG", as.character(mon.all.edges$gene2)))),]
    
    mon.hh.h.nodes <- unique(c(unique(as.character(mon.hh.edges$gene1))), c(unique(as.character(mon.hh.edges$gene2))))
    mon.pp.p.nodes <- unique(c(unique(as.character(mon.pp.edges$gene1))), c(unique(as.character(mon.pp.edges$gene2))))
    
    mon.hp.h.nodes.col1 <- unique(as.character(mon.hp.edges$gene1))
    mon.hp.p.nodes.col1 <- unique(as.character(mon.hp.edges$gene2))
    
    ####
    hum.all.edges <- hum_pval0[(abs(hum_pval0$cor) >= r1 & abs(hum_pval0$cor) < r2),]
    hum.hh.edges <- hum.all.edges[(grepl(pattern = "h_OG", as.character(hum.all.edges$gene1)) & grepl(pattern = "h_OG", as.character(hum.all.edges$gene2))),]
    hum.pp.edges <- hum.all.edges[(grepl(pattern = "p_OG", as.character(hum.all.edges$gene1)) & grepl(pattern = "p_OG", as.character(hum.all.edges$gene2))),]
    hum.hp.edges <- hum.all.edges[((grepl(pattern = "h_OG", as.character(hum.all.edges$gene1)) & grepl(pattern = "p_OG", as.character(hum.all.edges$gene2))) | (grepl(pattern = "p_OG", as.character(hum.all.edges$gene1)) & grepl(pattern = "h_OG", as.character(hum.all.edges$gene2)))),]
    
    hum.hh.h.nodes <- unique(c(unique(as.character(hum.hh.edges$gene1))), c(unique(as.character(hum.hh.edges$gene2))))
    hum.pp.p.nodes <- unique(c(unique(as.character(hum.pp.edges$gene1))), c(unique(as.character(hum.pp.edges$gene2))))
    
    hum.hp.h.nodes.col1 <- unique(as.character(hum.hp.edges$gene1))
    hum.hp.p.nodes.col1 <- unique(as.character(hum.hp.edges$gene2))
    
    ####
    common.all.edges <- inner_join(mon.all.edges, hum.all.edges, by = c("gene1", "gene2"))
    common.hh.edges <- inner_join(mon.hh.edges, hum.hh.edges, by = c("gene1", "gene2"))
    common.pp.edges <- inner_join(mon.pp.edges, hum.pp.edges, by = c("gene1", "gene2"))
    common.hp.edges <- inner_join(mon.hp.edges, hum.hp.edges, by = c("gene1", "gene2"))
    common.hh.h.nodes <- unique(c(mon.hh.h.nodes, hum.hh.h.nodes))
    common.pp.p.nodes <- unique(c(mon.pp.p.nodes, hum.pp.p.nodes))
    common.hp.h.nodes <- unique(c(mon.hp.h.nodes.col1, hum.hp.h.nodes.col1))
    common.hp.p.nodes <- unique(c(mon.hp.p.nodes.col1, hum.hp.p.nodes.col1))
    
    df[a,2] <- pcc
    df[a,3] <- nrow(common.all.edges)
    df[a,4] <- nrow(common.hh.edges)
    df[a,5] <- nrow(common.pp.edges)
    df[a,6] <- nrow(common.hp.edges)
    df[a,7] <- length(common.hh.h.nodes)
    df[a,8] <- length(common.pp.p.nodes)
    df[a,9] <- length(common.hp.h.nodes)
    df[a,10] <- length(common.hp.p.nodes)

    a = a+1
  }
  return(df)
}

df_common_study <- rbind(common_edges_and_nodes(0.9), 
                         common_edges_and_nodes(0.8), 
                         common_edges_and_nodes(0.7),
                         common_edges_and_nodes(0.6),
                         common_edges_and_nodes(0.5))

save(df_common_study, file = "df_common_study.RData")
colnames(df_common_study) <- c("Perms", "cor", "all.edges", "hh.edges", "pp.edges", "hp.edges", "hh.h.nodes", "pp.p.nodes", "hp.h.nodes", "hp.p.nodes")
df_common_melt <- melt(df_common_study, id = c("Perms", "cor"))
colnames(df_common_melt) <- c("Perms", "cor", "edge.node.type", "count")

#edge.node.plot <- ggplot(df_individual_melt, aes(x = Perms, y = count, colour = factor(cor))) +
#                         geom_line() + scale_x_continuous(breaks = seq(10000, 100000, 25000), limits = c(10000, 100000))+
#                         theme_bw() +
#                         facet_wrap(. ~ edge.node.type, scales="free_y") + ggtitle("SRP118996: Number of edges/nodes with permute score == 0")

#ggsave(plot = edge.node.plot, filename = "edge.node.plot.png")

com.edges.df <- df_common_melt[grep(pattern = "edge", df_common_melt$edge.node.type),]
com.nodes.df <- df_common_melt[grep(pattern = "node", df_common_melt$edge.node.type),]

com.edge.plot <- ggplot(com.edges.df, aes(x = Perms, y = log10(count), colour = factor(cor))) +
                         geom_line() + scale_x_continuous(breaks = seq(10000, 100000, 25000), limits = c(10000, 100000))+
                         theme_bw() +
                         facet_wrap(. ~ edge.node.type, scales="free_y") + ggtitle("ERP106451 and SRP118996: Common edges with permutation score 0") +
                         xlab("Permutations")

ggsave(plot = com.edge.plot, filename = "com.edge.plot.png")

com.node.plot <- ggplot(com.nodes.df, aes(x = Perms, y = log10(count), colour = factor(cor))) +
                         geom_line() + scale_x_continuous(breaks = seq(10000, 100000, 25000), limits = c(10000, 100000))+
                         theme_bw() +
                         facet_wrap(. ~ edge.node.type, scales="free_y") + ggtitle("ERP106451 and SRP118996: Common nodes with permutation score 0") +
                         xlab("Permutations")
 

ggsave(plot = com.node.plot, filename = "com.node.plot.png")

