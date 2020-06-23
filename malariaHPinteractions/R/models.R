# all modelling stuff
library(dplyr)
library(igraph)
library(betareg)
library(ggeffects)
library(ggplot2)
library(reshape2)

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
############################## data pre-processing

# PARASITE - PARASITE EDGES

pb_pp <- loadRData("df_concat_allhosts/cor/pberall_6_datasets_para_edges2.RData")
pb_pp <- loadRData("ERP004598/cor/ERP004598_all_para.RData")

ig_pb_pp <- graph_from_data_frame(pb_pp[,1:3], directed = F)
ig_pb_pp <- set_edge_attr(ig_pb_pp, "weight", index = E(ig_pb_pp), pb_pp[,3])
ig_pb_pp <- graph_from_data_frame(pb_pp[,1:3], directed = F)
ig_pb_pp <- set_edge_attr(ig_pb_pp, "weight", index = E(ig_pb_pp), pb_pp[,3])

# pbERALL NW
# get node properties with edge weights
pb_pp_bw_w <- betweenness(ig_pb_pp, v = V(ig_pb_pp), directed = FALSE, weights = abs(E(ig_pb_pp)$weight))
pb_pp_cl_w <- closeness(ig_pb_pp, vids = V(ig_pb_pp), weights = abs(E(ig_pb_pp)$weight))
pb_pp_dg_w <- pb_pp[,1:3] %>%
  setNames(., c("para1", "para2", "cor")) %>%
  melt(., "cor") %>%
  select(value, cor) %>%
  group_by(value) %>%
  summarise(degree = sum(abs(cor))) %>%
  tibble::deframe()
save.image("pberall_weighted_properties_para_corrected.RData")
# node properties without edge weights
pb_pp_bw <- betweenness(ig_pb_pp, v = V(ig_pb_pp), directed = FALSE)
pb_pp_cl <- closeness(ig_pb_pp, vids = V(ig_pb_pp))
pb_pp_dg <- degree(ig_pb_pp, v = V(ig_pb_pp))
save.image("pberall_unweighted_properties_para_corrected.RData")

# PBERGHEI NW
# get node properties with edge weights
pb_pp_bw_w <- betweenness(ig_pb_pp, v = V(ig_pb_pp), directed = FALSE, weights = abs(E(ig_pb_pp)$weight))
pb_pp_cl_w <- closeness(ig_pb_pp, vids = V(ig_pb_pp), weights = abs(E(ig_pb_pp)$weight))
pb_pp_dg_w <- pb_pp[,1:3] %>%
  setNames(., c("para1", "para2", "cor")) %>%
  melt(., "cor") %>%
  select(value, cor) %>%
  group_by(value) %>%
  summarise(degree = sum(abs(cor))) %>%
  tibble::deframe()
  save.image("Pfal_weighted_properties_para.RData")
  # node properties without edge weights
ig_pf_pp <- graph_from_data_frame(pf_pp[,1:3], directed = F)
pb_pp_bw <- betweenness(ig_pb_pp, v = V(ig_pb_pp), directed = FALSE)
pb_pp_cl <- closeness(ig_pb_pp, vids = V(ig_pb_pp))
pb_pp_dg <- degree(ig_pb_pp, v = V(ig_pb_pp))
save.image("Pfal_unweighted_properties_para.RData")
# PFal


# HOST - PARASITE EDGES (ONLY PARASITE GENES)

pb_bp <- loadRData("df_concat_allhosts/cor/df_concat_allhosts2_all_bipartite.RData")
pb_bp <- loadRData("ERP004598/cor/ERP004598_all_bipartite.RData")

ig_pb_bp <- graph_from_data_frame(pb_bp[,1:3], directed = F)
ig_pb_bp <- set_edge_attr(ig_pb_bp, "weight", index = E(ig_pb_bp), pb_bp[,3])
ig_pb_bp <- graph_from_data_frame(pb_bp[,1:3], directed = F)
ig_pb_bp <- set_edge_attr(ig_pb_bp, "weight", index = E(ig_pb_bp), pb_bp[,3])

# pbERALL NW
# get node properties with edge weights


pb_bp_bw_w <- betweenness(ig_pb_bp, v = V(ig_pb_bp), directed = FALSE, weights = abs(E(ig_pb_bp)$weight))
pb_bp_cl_w <- closeness(ig_pb_bp, vids = V(ig_pb_bp), weights = abs(E(ig_pb_bp)$weight))
pb_bp_dg_w <- pb_bp[,1:3] %>%
  setNames(., c("host", "para", "cor")) %>%
  select(para, cor) %>%
  group_by(para) %>%
  summarise(degree = sum(abs(cor))) %>%
  tibble::deframe()
  save.image("pberall_weighted_properties_bipartite_corrected.RData")

# node properties without edge weights
pb_bp_bw <- betweenness(ig_pb_bp, v = V(ig_pb_bp), directed = FALSE)
pb_bp_cl <- closeness(ig_pb_bp, vids = V(ig_pb_bp))
pb_bp_dg <- degree(ig_pb_bp, v = V(ig_pb_bp))
save.image("pberall_unweighted_properties_bipartite.RData")


pb_bp_bw <- betweenness(ig_pb_bp, v = V(ig_pb_bp), directed = FALSE)
pb_bp_cl <- closeness(ig_pb_bp, vids = V(ig_pb_bp))
pb_bp_dg <- degree(ig_pb_bp, v = V(ig_pb_bp))

save.image("Pberghei_unweighted_properties_bipartite.RData")
# PBERGHEI NW
# get node properties with edge weights
pb_bp_bw_w <- betweenness(ig_pb_bp, v = V(ig_pb_bp), directed = FALSE, weights = abs(E(ig_pb_bp)$weight))
pb_bp_cl_w <- closeness(ig_pb_bp, vids = V(ig_pb_bp), weights = abs(E(ig_pb_bp)$weight))
pb_bp_dg_w <- pb_bp[,1:3] %>%
  setNames(., c("host", "para", "cor")) %>%
  select(para, cor) %>%
  group_by(para) %>%
  summarise(degree = sum(abs(cor))) %>%
  tibble::deframe()library(dplyr)
library(igraph)
library(betareg)
library(ggeffects)
library(ggplot2)
library(reshape2)

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
#
 save.image("Pberghei_weighted_properties_bipartite_corrected.RData")
# node properties without edge weights

#Pfal
 
library(dplyr)
library(igraph)
library(betareg)
library(ggeffects)
library(ggplot2)
library(reshape2)

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
#
pf_pp <- loadRData("DRP000987/cor/DRP000987_str_para.RData")
ig_pf_pp <- graph_from_data_frame(pf_pp[,1:3], directed = F)

# node properties without edge weights
pf_pp_bw <- betweenness(ig_pf_pp, v = V(ig_pf_pp), directed = FALSE)
pf_pp_cl <- closeness(ig_pf_pp, vids = V(ig_pf_pp))
pf_pp_dg <- degree(ig_pf_pp, v = V(ig_pf_pp))
save.image("Pfal_unweighted_properties_para.RData")

ig_pf_pp <- set_edge_attr(ig_pf_pp, "weight", index = E(ig_pf_pp), pf_pp[,3])

pf_pp_bw_w <- betweenness(ig_pf_pp, v = V(ig_pf_pp), directed = FALSE, weights = abs(E(ig_pf_pp)$weight))
pf_pp_cl_w <- closeness(ig_pf_pp, vids = V(ig_pf_pp), weights = abs(E(ig_pf_pp)$weight))
pf_pp_dg_w <- pf_pp[,1:3] %>%
  setNames(., c("para1", "para2", "cor")) %>%
  melt(., "cor") %>%
  select(value, cor) %>%
  group_by(value) %>%
  summarise(degree = sum(abs(cor))) %>%
  tibble::deframe()
save.image("Pfal_weighted_properties_para.RData")
rm(list = list(pf_pp, ig_pf_pp))

pf_bp <- loadRData("DRP000987/cor/DRP000987_str_bipartite.RData")
ig_pf_bp <- graph_from_data_frame(pf_bp[,1:3], directed = F)

pf_bp_bw <- betweenness(ig_pf_bp, v = V(ig_pf_bp), directed = FALSE)
pf_bp_cl <- closeness(ig_pf_bp, vids = V(ig_pf_bp))
pf_bp_dg <- degree(ig_pf_bp, v = V(ig_pf_bp))
save.image("Pfal_unweighted_properties_bipartite.RData")

# get node properties with edge weights
ig_pf_bp <- set_edge_attr(ig_pf_bp, "weight", index = E(ig_pf_bp), pf_bp[,3])

pf_bp_bw_w <- betweenness(ig_pf_bp, v = V(ig_pf_bp), directed = FALSE, weights = abs(E(ig_pf_bp)$weight))
pf_bp_cl_w <- closeness(ig_pf_bp, vids = V(ig_pf_bp), weights = abs(E(ig_pf_bp)$weight))
pf_bp_dg_w <- pf_bp[,1:3] %>%
  setNames(., c("host", "para", "cor")) %>%
  select(para, cor) %>%
  group_by(para) %>%
  summarise(degree = sum(abs(cor))) %>%
  tibble::deframe()
save.image("Pfal_weighted_properties_bipartite.RData")


## PlasmoGEM dataset
Barseq20200228 <- read.csv("Barseq20200228.csv", stringsAsFactors=FALSE)
parasite_orthogroups <- read.delim("/SAN/Plasmo_compare/OrthoFinder/parasite_orthogroups.txt", stringsAsFactors=FALSE)
para <- parasite_orthogroups[,c(1,4)]
colnames(para)[2] <- "gene"
RGR_OG <- inner_join(Barseq20200228, para, by.x="gene", by.y= "Pb_g")
# WEIGHTED PLASMOGEM
# inner_join with the parasite genes from orthogroups with all 8 parasite vectors
# pick only parasite genes from bw objects
pb_bp_bw_w_p <- pb_bp_bw_w[grep(pattern = "p_OG", names(pb_bp_bw_w))]
pb_bp_bw_w_p <- pb_bp_bw_w[grep(pattern = "p_OG", names(pb_bp_bw_w))]
pf_bp_bw_w_p <- pf_bp_bw_w[grep(pattern = "p_OG", names(pf_bp_bw_w))]


pb_pp_bw_w_p <- pb_pp_bw_w[grep(pattern = "p_OG", names(pb_pp_bw_w))]
pb_pp_bw_w_p <- pb_pp_bw_w[grep(pattern = "p_OG", names(pb_pp_bw_w))]
pf_pp_bw_w_p <- pf_pp_bw_w[grep(pattern = "p_OG", names(pf_pp_bw_w))]

pb_bp_bw_w_df <- as.data.frame(pb_bp_bw_w_p) %>% tibble::rownames_to_column("Orthogroup")
pb_bp_bw_w_df <- as.data.frame(pb_bp_bw_w_p) %>% tibble::rownames_to_column("Orthogroup")
pb_pp_bw_w_df <- as.data.frame(pb_pp_bw_w_p) %>% tibble::rownames_to_column("Orthogroup")
pb_pp_bw_w_df <- as.data.frame(pb_pp_bw_w_p) %>% tibble::rownames_to_column("Orthogroup")

pf_bp_bw_w_df <- as.data.frame(pf_bp_bw_w_p) %>% tibble::rownames_to_column("Orthogroup")
pf_pp_bw_w_df <- as.data.frame(pf_pp_bw_w_p) %>% tibble::rownames_to_column("Orthogroup")

pb_bp_dg_w_df <- as.data.frame(pb_bp_dg_w) %>% tibble::rownames_to_column("Orthogroup")
pb_bp_dg_w_df <- as.data.frame(pb_bp_dg_w) %>% tibble::rownames_to_column("Orthogroup")
pb_pp_dg_w_df <- as.data.frame(pb_pp_dg_w) %>% tibble::rownames_to_column("Orthogroup")
pb_pp_dg_w_df <- as.data.frame(pb_pp_dg_w) %>% tibble::rownames_to_column("Orthogroup")

pf_bp_dg_w_df <- as.data.frame(pf_bp_dg_w) %>% tibble::rownames_to_column("Orthogroup")
pf_pp_dg_w_df <- as.data.frame(pf_pp_dg_w) %>% tibble::rownames_to_column("Orthogroup")

#join.all <- plyr::join_all(list(pb_bp_bw_w_df, pb_bp_bw_w_df, pb_pp_bw_w_df, pb_pp_bw_w_df, pb_bp_dg_w_df, pb_bp_dg_w_df, pb_pp_dg_w_df, pb_pp_dg_w_df), by='Orthogroup', type='inner')

# don't just join but also keep all the genes
join.all.full_w <- plyr::join_all(list(pb_bp_bw_w_df, pb_bp_bw_w_df, 
                                       pb_pp_bw_w_df, pb_pp_bw_w_df, 
                                       pb_bp_dg_w_df, pb_bp_dg_w_df, 
                                       pb_pp_dg_w_df, pb_pp_dg_w_df,
                                       pf_pp_dg_w_df, pf_bp_df_w_df,
                                       pf_pp_bw_w_df, pf_bp_bw_w_df), 
                                  by='Orthogroup', type='full')


join.all.full_w <- plyr::join_all(list(join.all.full_w, pf_pp_dg_w_df, pf_bp_dg_w_df,
                                  pf_pp_bw_w_df, pf_bp_bw_w_df), by='Orthogroup', type='full')

join.all.full_w[is.na(join.all.full_w)] <- 0
RGR_w <- inner_join(join.all.full_w, RGR_OG)

RGR_w <- na.omit(RGR_w)
RGR <- RGR_w$Relative.Growth.Rate
RGR[which(RGR >= 1)] <- 0.99999999
#RGR[which(RGR == 0)] <- 0.00000001
range(RGR)
#[1] 0.03307655 1.00000000
RGR_w$Relative.Growth.Rate <- RGR
# UNWEIGHTED PLASMOGEM
# inner_join with the parasite genes from orthogroups with all 8 parasite vectors

Barseq20200228 <- read.csv("Barseq20200228.csv", stringsAsFactors=FALSE)
parasite_orthogroups <- read.delim("/SAN/Plasmo_compare/OrthoFinder/parasite_orthogroups.txt", stringsAsFactors=FALSE)
para <- parasite_orthogroups[,c(1,4)]
colnames(para)[2] <- "gene"
RGR_OG <- inner_join(Barseq20200228, para, by.x="gene", by.y= "Pb_g")

pb_bp_bw_p <- pb_bp_bw[grep(pattern = "p_OG", names(pb_bp_bw))]
pb_bp_bw_p <- pb_bp_bw[grep(pattern = "p_OG", names(pb_bp_bw))]

pb_pp_bw_p <- pb_pp_bw[grep(pattern = "p_OG", names(pb_pp_bw))]
pb_pp_bw_p <- pb_pp_bw[grep(pattern = "p_OG", names(pb_pp_bw))]

pb_bp_bw_df <- as.data.frame(pb_bp_bw_p) %>% tibble::rownames_to_column("Orthogroup")
pb_bp_bw_df <- as.data.frame(pb_bp_bw_p) %>% tibble::rownames_to_column("Orthogroup")
pb_pp_bw_df <- as.data.frame(pb_pp_bw_p) %>% tibble::rownames_to_column("Orthogroup")
pb_pp_bw_df <- as.data.frame(pb_pp_bw_p) %>% tibble::rownames_to_column("Orthogroup")

pb_bp_dg_p <- pb_bp_dg[grep(pattern = "p_OG", names(pb_bp_dg))]
pb_bp_dg_df <- as.data.frame(pb_bp_dg_p) %>% tibble::rownames_to_column("Orthogroup")
pb_bp_dg_p <- pb_bp_dg[grep(pattern = "p_OG", names(pb_bp_dg))]
pb_bp_dg_df <- as.data.frame(pb_bp_dg_p) %>% tibble::rownames_to_column("Orthogroup")
pb_pp_dg_df <- as.data.frame(pb_pp_dg) %>% tibble::rownames_to_column("Orthogroup")
pb_pp_dg_df <- as.data.frame(pb_pp_dg) %>% tibble::rownames_to_column("Orthogroup")

pf_bp_bw_p <- pf_bp_bw[grep(pattern = "p_OG", names(pf_bp_bw))]
pf_pp_bw_p <- pf_pp_bw[grep(pattern = "p_OG", names(pf_pp_bw))]
pf_bp_bw_df <- as.data.frame(pf_bp_bw_p) %>% tibble::rownames_to_column("Orthogroup")
pf_pp_bw_df <- as.data.frame(pf_pp_bw_p) %>% tibble::rownames_to_column("Orthogroup")
pf_bp_dg_p <- pf_bp_dg[grep(pattern = "p_OG", names(pf_bp_dg))]
pf_bp_dg_df <- as.data.frame(pf_bp_dg_p) %>% tibble::rownames_to_column("Orthogroup")
pf_pp_dg_df <- as.data.frame(pf_pp_dg) %>% tibble::rownames_to_column("Orthogroup")


#join.all.unweighted <- plyr::join_all(list(pb_bp_bw_df, pb_bp_bw_df, pb_pp_bw_df, pb_pp_bw_df, pb_bp_dg_df, pb_bp_dg_df, pb_pp_dg_df, pb_pp_dg_df), by='Orthogroup', type='inner')

join.all.full_un<- plyr::join_all(list(pb_bp_bw_df, pb_bp_bw_df, pb_pp_bw_df, pb_pp_bw_df, pb_bp_dg_df, pb_bp_dg_df, pb_pp_dg_df, pb_pp_dg_df), by='Orthogroup', type='full')

join.all.full_un <- plyr::join_all(list(join.all.full_un, pf_pp_dg_df, pf_bp_dg_df,
                                       pf_pp_bw_df, pf_bp_bw_df), by='Orthogroup', type='full')

join.all.full_un[is.na(join.all.full_un)] <- 0
RGR_un <- inner_join(join.all.full_un, RGR_OG)

RGR_un <- na.omit(RGR_un)
RGR <- RGR_un$Relative.Growth.Rate
RGR[which(RGR >=1)] <- 0.99999999
range(RGR)
#[1] 0.03307655 1.00000000
RGR_un$Relative.Growth.Rate <- RGR


#################### models

pb_bp_dg_m <- betareg(Relative.Growth.Rate ~ pb_bp_dg, data = RGR_un) #works
pb_bp_dg_bw_m <- betareg(Relative.Growth.Rate ~ pb_bp_dg + pb_bp_bw_p, data = RGR_un) #does not work
pb_bp_pp_dg_m <- betareg(Relative.Growth.Rate ~ pb_bp_dg + pb_pp_dg, data = RGR_un) # works
pb_pp_dg_m <- betareg(Relative.Growth.Rate ~ pb_pp_dg, data = RGR_un)
pb_pp_dg_bw_m <- betareg(Relative.Growth.Rate ~ pb_pp_dg + pb_pp_bw_p, data = RGR_un)


pb_bp_dg_m <- betareg(Relative.Growth.Rate ~ pb_bp_dg, data = RGR_un) # works - positive effect
pb_pb_bp_dg_m <- betareg(Relative.Growth.Rate ~ pb_bp_dg + pb_bp_dg, data = RGR_un) # works
pb_pp_dg_m <- betareg(Relative.Growth.Rate ~ pb_pp_dg, data = RGR_un) # works very well
pb_pb_pp_dg_m <- betareg(Relative.Growth.Rate ~ pb_pp_dg + pb_pp_dg, data = RGR_un) # works very well no warning

pb_bp_dg_bw_m <- betareg(Relative.Growth.Rate ~ pb_bp_dg + pb_bp_bw_p, data = RGR_un) # opposite effects
pb_pb_bp_dg_bw_m <- betareg(Relative.Growth.Rate ~ pb_bp_dg + pb_bp_dg + 
                              pb_bp_bw_p + pb_bp_bw_p, data = RGR_un)# pberghei effects are opposite
pb_pp_dg_bw_m <- betareg(Relative.Growth.Rate ~ pb_pp_dg + pb_pp_bw_p, data = RGR_un) # works
pb_pb_pp_dg_bw_m <- betareg(Relative.Growth.Rate ~ pb_pp_dg + pb_pp_dg + 
                              pb_pp_bw_p + pb_pp_bw_p, data = RGR_un)# all good


pb_bp_pp_dg_m <- betareg(Relative.Growth.Rate ~ pb_bp_dg + pb_pp_dg, data = RGR_un) # works all good

# all degrees
pb_pb_bp_pp_dg_m <- betareg(Relative.Growth.Rate ~ pb_bp_dg + pb_pp_dg + 
                              pb_bp_dg + pb_pp_dg, data = RGR_un)
#pp_dg, data = RGR_un)# pb degree has opposite effect
# all betweenness and degres
pb_pb_bp_pp_dg_bw_m <- betareg(Relative.Growth.Rate ~ pb_bp_dg + pb_pp_dg +
                                 pb_bp_dg + pb_pp_dg + 
                                 pb_bp_bw_p + pb_pp_bw_p + 
                                 pb_bp_bw_p + pb_pp_bw_p, data = RGR_un) # pberall bip hs opposite effects
  


## weighted

pb_bp_dg_w_m <- betareg(Relative.Growth.Rate ~ pb_bp_dg_w, data = RGR_w) #has positive effect
pb_bp_dg_bw_w_m <- betareg(Relative.Growth.Rate ~ pb_bp_dg_w + pb_bp_bw_w_p, data = RGR_w) #both have positive effect
pb_bp_pp_dg_w_m <- betareg(Relative.Growth.Rate ~ pb_bp_dg_w + pb_pp_dg_w, data = RGR_w) # works:both have negative effect
pb_bp_pp_dg_bw_w_m <- betareg(Relative.Growth.Rate ~ pb_bp_dg_w + pb_pp_dg_w + 
                              pb_bp_bw_w_p + pb_pp_bw_w_p, data = RGR_w)# pp bw has negative effect


pb_bp_dg_w_m <- betareg(Relative.Growth.Rate ~ pb_bp_dg_w, data = RGR_w) # works - positive effect
pb_pb_bp_dg_w_m <- betareg(Relative.Growth.Rate ~ pb_bp_dg_w + pb_bp_dg_w, data = RGR_w) # both positive effects
pb_pp_dg_w_m <- betareg(Relative.Growth.Rate ~ pb_pp_dg_w, data = RGR_w) # works very well
pb_pb_pp_dg_w_m <- betareg(Relative.Growth.Rate ~ pb_pp_dg_w + pb_pp_dg_w, data = RGR_w) # works very well no warning


pb_bp_pp_dg_bw_w_m <- betareg(Relative.Growth.Rate ~ pb_bp_dg_w + pb_bp_bw_w_p + pb_pp_dg_w + pb_pp_bw_w_p, data = RGR_w)
pb_bp_dg_bw_w_m <- betareg(Relative.Growth.Rate ~ pb_bp_dg_w + pb_bp_bw_w_p, data = RGR_w) # opposite effects
pb_pb_bp_dg_bw_w_m <- betareg(Relative.Growth.Rate ~ pb_bp_dg_w + pb_bp_dg_w + 
                              pb_bp_bw_w_p + pb_bp_bw_w_p, data = RGR_w)# pberghei effects are opposite
pb_pp_dg_bw_w_m <- betareg(Relative.Growth.Rate ~ pb_pp_dg_w + pb_pp_bw_w_p, data = RGR_w) # works
pb_pb_pp_dg_bw_w_m <- betareg(Relative.Growth.Rate ~ pb_pp_dg_w + pb_pp_dg_w + 
                              pb_pp_bw_w_p + pb_pp_bw_w_p, data = RGR_w)# all good
pb_pb_bp_pp_dg_w_m <- betareg(Relative.Growth.Rate ~ pb_bp_dg_w + pb_pp_dg_w + pb_bp_dg_w + pb_pp_dg_w, data = RGR_w)
pb_pb_bp_pp_dg_bw_w_m <- betareg(Relative.Growth.Rate ~ pb_bp_dg_w + pb_pp_dg_w + pb_bp_dg_w + pb_pp_dg_w + pb_bp_bw_w_p + pb_pp_bw_w + pb_bp_bw_w_p + pb_pp_bw_w_p, data = RGR_w)


pb_bp_pp_dg_bw_m <- betareg(Relative.Growth.Rate ~ pb_bp_dg_p + pb_pp_dg + 
                              pb_bp_bw_p + pb_pp_bw_p, data = RGR_un)# works

pb_bp_pp_dg_bw_m <- betareg(Relative.Growth.Rate ~ pb_bp_dg_p + pb_pp_dg + 
                              pb_bp_bw_p + pb_pp_bw_p, data = RGR_un)# pberghei bipartite had opposite effects

pf_bp_pp_dg_bw_m <- betareg(Relative.Growth.Rate ~ pf_bp_dg_p + pf_pp_dg + 
                              pf_bp_bw_p + pf_pp_bw_p, data = RGR_un)# pberghei bipartite had opposite effects


# lower the AIC/BIC value, the better
pb_ef <- ggeffect(pb_bp_pp_dg_bw_m, terms = c("pb_pp_dg"))
  
  ggplot(pb_ef, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)
ggsave("pb_ef_pp_dg_bw.png")

pb_ef <- ggeffect(pb_bp_pp_dg_bw_m, terms = c("pb_pp_dg"))

ggplot(pb_ef, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)
ggsave("pb_ef_pp.png")

pf_ef <- ggeffect(pf_bp_pp_dg_bw_m, terms = c("pf_pp_dg"))

ggplot(pf_ef, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)
ggsave("pf_ef_pp.png")

##### compute unweighted stats for full networks of pberall, Pb and Pf

# load full networks

pb <- loadRData("/SAN/Plasmo_compare/SRAdb/Output/df_concat_allhosts/cor/df_concat_allhosts_all_na.omit.RData")
pb <-  pb[pb$permute_score==0,]

pb <- loadRData("/SAN/Plasmo_compare/SRAdb/Output/ERP004598/cor/ERP004598_all_na.omit.RData")
pb <-  pb[pb$permute_score==0,]

pf <- loadRData("/SAN/Plasmo_compare/SRAdb/Output/DRP000987/cor/DRP000987_str_na.omit.RData")
pf <-  pf[pf$permute_score==0,]

# make unweighted igraphs for all
pb_ig <- graph_from_data_frame(pb[,1:3], directed = F)

pb_all_bw <- betweenness(pb_ig, v = V(pb_ig), directed = FALSE)
pb_all_cl <- closeness(pb_ig, vids = V(pb_ig))
pb_all_dg <- degree(pb_ig, v = V(pb_ig))

pb_ig <- graph_from_data_frame(pb[,1:3], directed = F)

pb_all_bw <- betweenness(pb_ig, v = V(pb_ig), directed = FALSE)
pb_all_cl <- closeness(pb_ig, vids = V(pb_ig))
pb_all_dg <- degree(pb_ig, v = V(pb_ig))

pf_ig <- graph_from_data_frame(pf[,1:3], directed = F)

pf_all_bw <- betweenness(pf_ig, v = V(pf_ig), directed = FALSE)
pf_all_cl <- closeness(pf_ig, vids = V(pf_ig))
pf_all_dg <- degree(pf_ig, v = V(pf_ig))

save.image("Unweighted_full_network_stats_pb_pb_pf.RData")

load("Unweighted_deg_bw_pberall_Pberghei_Pfal_joined_RGR.RData")

ov_all_bw_p <- ov_all_bw[grep(pattern = "p_OG", names(ov_all_bw))]
ov_all_cl_p <- ov_all_cl[grep(pattern = "p_OG", names(ov_all_cl))]
ov_all_dg_p <- ov_all_dg[grep(pattern = "p_OG", names(ov_all_dg))]

pb_all_bw_p <- pb_all_bw[grep(pattern = "p_OG", names(pb_all_bw))]
pb_all_cl_p <- pb_all_cl[grep(pattern = "p_OG", names(pb_all_cl))]
pb_all_dg_p <- pb_all_dg[grep(pattern = "p_OG", names(pb_all_dg))]

pf_all_bw_p <- pf_all_bw[grep(pattern = "p_OG", names(pf_all_bw))]
pf_all_cl_p <- pf_all_cl[grep(pattern = "p_OG", names(pf_all_cl))]
pf_all_dg_p <- pf_all_dg[grep(pattern = "p_OG", names(pf_all_dg))]

ov_all_bw_df <- as.data.frame(ov_all_bw_p) %>% 
  tibble::rownames_to_column("Orthogroup")
ov_all_cl_df <- as.data.frame(ov_all_cl_p) %>% 
  tibble::rownames_to_column("Orthogroup")
ov_all_dg_df <- as.data.frame(ov_all_dg_p) %>% 
  tibble::rownames_to_column("Orthogroup")

pb_all_bw_df <- as.data.frame(pb_all_bw_p) %>% 
  tibble::rownames_to_column("Orthogroup")
pb_all_cl_df <- as.data.frame(pb_all_cl_p) %>% 
  tibble::rownames_to_column("Orthogroup")
pb_all_dg_df <- as.data.frame(pb_all_dg_p) %>% 
  tibble::rownames_to_column("Orthogroup")

pf_all_bw_df <- as.data.frame(pf_all_bw_p) %>% 
  tibble::rownames_to_column("Orthogroup")
pf_all_cl_df <- as.data.frame(pf_all_cl_p) %>% 
  tibble::rownames_to_column("Orthogroup")
pf_all_dg_df <- as.data.frame(pf_all_dg_p) %>% 
  tibble::rownames_to_column("Orthogroup")

ov_pb_pf.full_un <- plyr::join_all(list(ov_all_bw_df, ov_all_cl_df, ov_all_dg_df, 
                                        pb_all_bw_df, pb_all_cl_df, pb_all_dg_df,
                                        pf_all_bw_df, pf_all_cl_df, pf_all_dg_df), by='Orthogroup', type='full')

ov_pb_pf.full_un[is.na(ov_pb_pf.full_un)] <- 0

Barseq20200228 <- read.csv("Barseq20200228.csv", stringsAsFactors=FALSE)
Pfal_ess <- read.csv("Pfal_ess.csv", stringsAsFactors=FALSE)
Pfal_ess$Gene_ID <- substring(Pfal_ess$Gene_ID, 1, 13)
parasite_orthogroups <- read.delim("/SAN/Plasmo_compare/OrthoFinder/parasite_orthogroups.txt", stringsAsFactors=FALSE)
para <- parasite_orthogroups[,c(1,4)]
para_pf <- parasite_orthogroups[,c(1,5)]
colnames(para)[2] <- "gene"
colnames(para_pf)[2] <- "Gene_ID"

RGR_OG <- inner_join(Barseq20200228, para, by.x="gene", by.y= "Pb_g")
Pfal_ess_ortho <- inner_join(Pfal_ess, para_pf)
RGR_un <- inner_join(pb_pb_pf.full_un, RGR_OG)
RGR_Ov_Pb_Pf <- inner_join(RGR_j, Pfal_ess_ortho, b = "Orthogroup")

RGR_un <- na.omit(RGR_un)
RGR <- RGR_un$Relative.Growth.Rate
RGR[which(RGR >=1)] <- 0.99999999
range(RGR)
#[1] 0.03307655 1.00000000
RGR_un$Relative.Growth.Rate <- RGR

MIS <- RGR_Ov_Pb_Pf$MIS
MIS[which(MIS >=1)] <- 0.99999999
RGR_Ov_Pb_Pf$MIS <- MIS


ov_bp_cl_p <- ov_bp_cl[grep(pattern = "p_OG", names(ov_bp_cl))]
ov_pp_cl_p <- ov_pp_cl[grep(pattern = "p_OG", names(ov_pp_cl))]

pb_bp_cl_p <- pb_bp_cl[grep(pattern = "p_OG", names(pb_bp_cl))]
pb_pp_cl_p <- pb_pp_cl[grep(pattern = "p_OG", names(pb_pp_cl))]

pf_bp_cl_p <- pf_bp_cl[grep(pattern = "p_OG", names(pf_bp_cl))]
pf_pp_cl_p <- pf_pp_cl[grep(pattern = "p_OG", names(pf_pp_cl))]

ov_bp_cl_df <- as.data.frame(ov_bp_cl_p) %>% 
  tibble::rownames_to_column("Orthogroup")
ov_pp_cl_df <- as.data.frame(ov_pp_cl_p) %>% 
  tibble::rownames_to_column("Orthogroup")

pb_bp_cl_df <- as.data.frame(pb_bp_cl_p) %>% 
  tibble::rownames_to_column("Orthogroup")
pb_pp_cl_df <- as.data.frame(pb_pp_cl_p) %>% 
  tibble::rownames_to_column("Orthogroup")

pf_bp_cl_df <- as.data.frame(pf_bp_cl_p) %>% 
  tibble::rownames_to_column("Orthogroup")
pf_pp_cl_df <- as.data.frame(pf_pp_cl_p) %>% 
  tibble::rownames_to_column("Orthogroup")

join.cl<-plyr::join_all(list(ov_bp_cl_df, ov_pp_cl_df,
                      pb_bp_cl_df, pb_pp_cl_df, pf_bp_cl_df, pf_pp_cl_df), by = "Orthogroup", type = "full")
join.cl[is.na(join.cl)] <- 0

RGR_k <- inner_join(RGR_un, RGR_unw, by = "Orthogroup")

RGR_j <- inner_join(join.cl, RGR_k, by = "Orthogroup")

# models
RGR_j <- na.omit(RGR_j)
RGR <- RGR_j$Relative.Growth.Rate.x
RGR[which(RGR >=1)] <- 0.99999999
range(RGR)
#[1] 0.03307655 1.00000000
RGR_j$Relative.Growth.Rate <- RGR

ov_pp_m <- betareg(data = RGR_j, Relative.Growth.Rate ~ ov_pp_dg + ov_pp_bw_p + ov_pp_cl_p) # does not converge
ov_pp_m <- betareg(data = RGR_j, Relative.Growth.Rate ~ ov_pp_dg + ov_pp_bw_p) # awesome

pb_pp_m <- betareg(data = RGR_j, Relative.Growth.Rate ~ pb_pp_dg + pb_pp_bw_p + pb_pp_cl_p) # bad model
pb_pp_m <- betareg(data = RGR_j, Relative.Growth.Rate ~ pb_pp_dg + pb_pp_bw_p) # bw notsig

pf_pp_m <- betareg(data = RGR_j, Relative.Growth.Rate ~ pf_pp_dg + pf_pp_bw_p + pf_pp_cl_p) # cl has -ve effect
pf_pp_m <- betareg(data = RGR_j, Relative.Growth.Rate ~ pf_pp_dg + pf_pp_bw_p) # bw not sig

ov_all_m <- betareg(data = RGR_j, Relative.Growth.Rate ~ ov_all_dg + ov_all_bw_p + ov_all_cl_p) # bw not sig
pb_all_m <- betareg(data = RGR_j, Relative.Growth.Rate ~ pb_all_dg_p.x + pb_all_bw_p.x + pb_all_cl_p.x) # bw not sig
pf_all_m <- betareg(data = RGR_j, Relative.Growth.Rate ~ pf_all_dg_p.x + pf_all_bw_p.x + pf_all_cl_p.x) # horrible models


pb_all_ef <- ggeffect(pb_all_m, terms = c("pb_all_dg_p.x"))

ggplot(pb_all_ef, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)
ggsave("pb_all_ef.png")

pb_pp_ef <- ggeffect(pb_pp_m, terms = c("pb_pp_dg"))

ggplot(pb_pp_ef, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)
ggsave("pb_pp_ef.png")

ov_all_ef <- ggeffect(ov_all_m, terms = c("ov_all_dg_p"))

ggplot(ov_all_ef, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)
ggsave("ov_all_ef.png")

ov_pp_ef <- ggeffect(ov_pp_m, terms = c("ov_pp_dg"))

ggplot(ov_pp_ef, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)
ggsave("ov_pp_ef.png")

# new stats for bipartite network: clustering coefficient, density, size
bip <- loadRData("df_concat_allhosts/cor/df_concat_allhosts_all_bipartite.RData")

bip_ig <- graph_from_data_frame(ov, directed = F)

# size is 8 million
# clustering coef == transitivity
bip_ig_tr <- transitivity(bip_ig, type = "global", weights = NULL)

# density
bip_ig_dens <- edge_density(bip_ig, loops = FALSE)

# new stats for full nw ov, pb, pf : eigen centrality, k-shells, expected force
# 1. Ov

ov <- loadRData("df_concat_allhosts/cor/df_concat_allhosts_all_na.omit.RData") %>%
  filter(permute_score == 0)

ov_ig <- graph_from_data_frame(ov, directed = F)

# a. eigen centrality

ov_ig_ec <- eigen_centrality(ov_ig, directed = FALSE)

# b. k-shells or k-core
ov_ig_kc <- coreness(ov_ig, mode = c("all"))

ov_ig_ec_p <- ov_ig_ec$vector[grep(pattern = "p_OG", names(ov_ig_ec$vector))]
ov_ig_ec_df <- as.data.frame(ov_ig_ec_p) %>% 
  tibble::rownames_to_column("Orthogroup")

ov_ig_kc_p <- ov_ig_kc[grep(pattern = "p_OG", names(ov_ig_kc))]
ov_ig_kc_df <- as.data.frame(ov_ig_kc_p) %>% 
  tibble::rownames_to_column("Orthogroup")

# c. expected force


# 2. Pb
pb <- loadRData("ERP004598/cor/ERP004598_all_na.omit.RData") %>%
  filter(permute_score == 0)

pb_ig <- graph_from_data_frame(pb, directed = F)

# a. eigen centrality

pb_ig_ec <- eigen_centrality(pb_ig, directed = FALSE)

# b. k-shells or k-core
pb_ig_kc <- coreness(pb_ig, mode = c("all"))

pb_ig_ec_p <- pb_ig_ec$vector[grep(pattern = "p_OG", names(pb_ig_ec$vector))]
pb_ig_ec_df <- as.data.frame(pb_ig_ec_p) %>% 
  tibble::rownames_to_column("Orthogroup")

pb_ig_kc_p <- pb_ig_kc[grep(pattern = "p_OG", names(pb_ig_kc))]
pb_ig_kc_df <- as.data.frame(pb_ig_kc_p) %>% 
  tibble::rownames_to_column("Orthogroup")


# c. expected force

# 3. pf
pf <- loadRData("DRP000987/cor/DRP000987_str_na.omit.RData") %>%
  filter(permute_score == 0)

pf_ig <- graph_from_data_frame(pf, directed = F)

# a. eigen centrality

pf_ig_ec <- eigen_centrality(pf_ig, directed = FALSE)

# b. k-shells or k-core
pf_ig_kc <- coreness(pf_ig, mode = c("all"))

# c. expected force

pf_ig_ec_p <- pf_ig_ec$vector[grep(pattern = "p_OG", names(pf_ig_ec$vector))]
pf_ig_ec_df <- as.data.frame(pf_ig_ec_p) %>% 
  tibble::rownames_to_column("Orthogroup")

pf_ig_kc_p <- pf_ig_kc[grep(pattern = "p_OG", names(pf_ig_kc))]
pf_ig_kc_df <- as.data.frame(pf_ig_kc_p) %>% 
  tibble::rownames_to_column("Orthogroup")

join <-plyr::join_all(list(ov_ig_ec_df, ov_ig_kc_df, 
                             pb_ig_ec_df, pb_ig_kc_df,
                             pf_ig_ec_df, pf_ig_kc_df), by = "Orthogroup", type = "full")
join[is.na(join)] <- 0

RGR_MIS <- inner_join(join, RGR_j, by = "Orthogroup")


# new stats for para nw ov, pb, pf : eigen centrality, k-shells, expected force
#1. Ov

ov_p <- loadRData("df_concat_allhosts/df_concat_allhosts_all_para.RData")
ov_ig_p <- graph_from_data_frame(ov_p, directed = F)

# a. eig_en centrality

ov_ig_p_ec <- eigen_centrality(ov_ig_p, directed = FALSE)
ov_ig_p_ec <- ov_ig_p_ec$vector

# b. k-shells or k-core
ov_ig_p_kc <- coreness(ov_ig_p, mode = c("all"))

ov_ig_p_ec_p <- ov_ig_p_ec[grep(pattern = "p_OG", names(ov_ig_p_ec))]
ov_ig_p_ec_df <- as.data.frame(ov_ig_p_ec_p) %>% 
  tibble::rownames_to_column("Orthogroup")

ov_ig_p_kc_p <- ov_ig_p_kc[grep(pattern = "p_OG", names(ov_ig_p_kc))]
ov_ig_p_kc_df <- as.data.frame(ov_ig_p_kc_p) %>% 
  tibble::rownames_to_column("Orthogroup")

# 2. Pb
pb_p <- loadRData("ERP004598/cor/ERP004598_all_para_edges.RData")

pb_ig_p <- graph_from_data_frame(pb_p, directed = F)

# a. eigen centrality

pb_ig_p_ec <- eigen_centrality(pb_ig_p, directed = FALSE)
pb_ig_p_ec <- pb_ig_p_ec$vector

# b. k-shells or k-core
pb_ig_p_kc <- coreness(pb_ig_p, mode = c("all"))

pb_ig_p_ec_p <- pb_ig_p_ec[grep(pattern = "p_OG", names(pb_ig_p_ec))]
pb_ig_p_ec_df <- as.data.frame(pb_ig_p_ec_p) %>% 
  tibble::rownames_to_column("Orthogroup")

pb_ig_p_kc_p <- pb_ig_p_kc[grep(pattern = "p_OG", names(pb_ig_p_kc))]
pb_ig_p_kc_df <- as.data.frame(pb_ig_p_kc_p) %>% 
  tibble::rownames_to_column("Orthogroup")

# 3. pf
pf_p <- loadRData("DRP000987/cor/DRP000987_str_para.RData")

pf_ig_p <- graph_from_data_frame(pf_p, directed = F)

# a. eigen centrality

pf_ig_p_ec <- eigen_centrality(pf_ig_p, directed = FALSE)
pf_ig_p_ec <- pf_ig_p_ec$vector

# b. k-shells or k-core
pf_ig_p_kc <- coreness(pf_ig_p, mode = c("all"))

pf_ig_p_ec_p <- pf_ig_p_ec[grep(pattern = "p_OG", names(pf_ig_p_ec))]
pf_ig_p_ec_df <- as.data.frame(pf_ig_p_ec_p) %>% 
  tibble::rownames_to_column("Orthogroup")

pf_ig_p_kc_p <- pf_ig_p_kc[grep(pattern = "p_OG", names(pf_ig_p_kc))]
pf_ig_p_kc_df <- as.data.frame(pf_ig_p_kc_p) %>% 
  tibble::rownames_to_column("Orthogroup")


join_p <-plyr::join_all(list(ov_ig_p_ec_df, ov_ig_p_kc_df, 
                             pb_ig_p_ec_df, pb_ig_p_kc_df,
                             pf_ig_p_ec_df, pf_ig_p_kc_df), by = "Orthogroup", type = "full")
join_p[is.na(join_p)] <- 0

RGR_MIS <- inner_join(join_p, RGR_Ov_Pb_Pf, by = "Orthogroup")

# cleaning up RGR_MIS
n <- grep(pattern = ".x", colnames(RGR_MIS))
RGR_MIS <- RGR_MIS[,c(-n)]
n <- c(54,55,56,57,58,59,60,61,62,63,64,39,44,45,48)
RGR_MIS <- RGR_MIS[,c(-n)]

colnames(RGR_MIS)[2:7] <- c("ov_all_ec", "ov_all_kc", "pb_all_ec", "pb_all_c", "pf_all_ec", "pf_all_kc")
colnames(RGR_MIS)[8:13] <- c("ov_pp_ec", "ov_pp_kc", "pb_pp_ec", "pb_pp_kc", "pf_pp_ec", "pf_pp_kc")
colnames(RGR_MIS)[14:15] <- c("ov_bp_cl","ov_pp_cl")
colnames(RGR_MIS)[16:19] <- c("pb_bp_cl", "pb_pp_cl", "pf_bp_cl", "pf_pp_cl")
colnames(RGR_MIS)[20:22] <- c("ov_bp_bw", "pb_bp_bw", "pb_pp_bw")
colnames(RGR_MIS)[23:24] <- c("ov_bp_dg", "pb_bp_dg")
colnames(RGR_MIS)[27:30] <- c("pf_bp_dg", "pf_pp_bw", "pf_bp_bw", "ov_pp_bw")
colnames(RGR_MIS)[32:37] <- c("pb_all_bw", "pb_all_cl", "pb_all_dg", "pf_all_bw", "pf_all_cl", "pf_all_dg")

RGR_MIS_cleaned <- RGR_MIS
save(RGR_MIS_cleaned, file = "Unweighted_RGR_MIS_cleaned.RData")

join_ov_full <-plyr::join_all(list(ov_all_bw_df, ov_all_cl_df, 
                             ov_all_dg_df), by = "Orthogroup", type = "full")
join_ov_full[is.na(join_ov_full)] <- 0

RGR_MIS_cleaned <- inner_join(RGR_MIS_cleaned, join_ov_full)
#models

ov_pp_m <- betareg(data = RGR_j, Relative.Growth.Rate ~ ov_pp_dg + ov_pp_bw_p) # awesome


# ERP106451 - full nw

pfE <- loadRData("ERP106451/cor/ERP106451_int_na.omit.RData") %>%
  filter(permute_score == 0)

pfE_ig <- graph_from_data_frame(pfE, directed = F)

# a. eigen centrality

pfE_ig_ec <- eigen_centrality(pfE_ig, directed = FALSE)

# b. k-shells or k-core
pfE_ig_kc <- coreness(pfE_ig, mode = c("all"))
# c. degree
pfE_ig_dg <- degree(pfE_ig, v = V(pfE_ig))
# d. betweenness
pfE_ig_bw <- betweenness(pfE_ig, directed = FALSE, v = V(pfE_ig))
# e. closeness
pfE_ig_cl <- closeness(pfE_ig, vids = V(pfE_ig))

pfE_ig_dg_p <- pfE_ig_dg[grep(pattern = "p_OG", names(pfE_ig_dg))]
pfE_ig_dg_df <- as.data.frame(pfE_ig_dg_p) %>% 
  tibble::rownames_to_column("Orthogroup")

pfE_ig_bw_p <- pfE_ig_bw[grep(pattern = "p_OG", names(pfE_ig_bw))]
pfE_ig_bw_df <- as.data.frame(pfE_ig_bw_p) %>% 
  tibble::rownames_to_column("Orthogroup")

pfE_ig_cl_p <- pfE_ig_cl[grep(pattern = "p_OG", names(pfE_ig_cl))]
pfE_ig_cl_df <- as.data.frame(pfE_ig_cl_p) %>% 
  tibble::rownames_to_column("Orthogroup")

pfE_ig_ec_p <- pfE_ig_ec$vector[grep(pattern = "p_OG", names(pfE_ig_ec$vector))]
pfE_ig_ec_df <- as.data.frame(pfE_ig_ec_p) %>% 
  tibble::rownames_to_column("Orthogroup")

pfE_ig_kc_p <- pfE_ig_kc[grep(pattern = "p_OG", names(pfE_ig_kc))]
pfE_ig_kc_df <- as.data.frame(pfE_ig_kc_p) %>% 
  tibble::rownames_to_column("Orthogroup")

# ERP106451 - para nw

pfE_p <- loadRData("ERP106451/cor/ERP106451_int_para.RData")

pfE_p_ig <- graph_from_data_frame(pfE_p, directed = F)

# a. ep_igen centrality

pfE_p_ig_ec <- eigen_centrality(pfE_p_ig, directed = FALSE)

# b. k-shells or k-core
pfE_p_ig_kc <- coreness(pfE_p_ig, mode = c("all"))
# c. degree
pfE_p_ig_dg <- degree(pfE_p_ig, v = V(pfE_p_ig))
# d. betweenness
pfE_p_ig_bw <- betweenness(pfE_p_ig, directed = FALSE, v = V(pfE_p_ig))
# e. closeness
pfE_p_ig_cl <- closeness(pfE_p_ig, vids = V(pfE_p_ig))

pfE_p_ig_dg_p <- pfE_p_ig_dg[grep(pattern = "p_OG", names(pfE_p_ig_dg))]
pfE_p_ig_dg_df <- as.data.frame(pfE_p_ig_dg_p) %>% 
  tibble::rownames_to_column("Orthogroup")

pfE_p_ig_bw_p <- pfE_p_ig_bw[grep(pattern = "p_OG", names(pfE_p_ig_bw))]
pfE_p_ig_bw_df <- as.data.frame(pfE_p_ig_bw_p) %>% 
  tibble::rownames_to_column("Orthogroup")

pfE_p_ig_cl_p <- pfE_p_ig_cl[grep(pattern = "p_OG", names(pfE_p_ig_cl))]
pfE_p_ig_cl_df <- as.data.frame(pfE_p_ig_cl_p) %>% 
  tibble::rownames_to_column("Orthogroup")

pfE_p_ig_ec_p <- pfE_p_ig_ec$vector[grep(pattern = "p_OG", names(pfE_p_ig_ec$vector))]
pfE_p_ig_ec_df <- as.data.frame(pfE_p_ig_ec_p) %>% 
  tibble::rownames_to_column("Orthogroup")

pfE_p_ig_kc_p <- pfE_p_ig_kc[grep(pattern = "p_OG", names(pfE_p_ig_kc))]
pfE_p_ig_kc_df <- as.data.frame(pfE_p_ig_kc_p) %>% 
  tibble::rownames_to_column("Orthogroup")