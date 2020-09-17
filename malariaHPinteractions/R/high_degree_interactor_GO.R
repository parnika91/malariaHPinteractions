# script to see what GO terms are enriched among the interactors of hogh degree genes
# 1. take the highest degree host genes, for each gene, find all parasite interactors, do GO term analysis for these interactors. same for top degree parasite genes
# 2. take all parasite genes in a GO term, get all host interctors, do GO terms analysis on these host genes - eg, what host GO terms are parasite pathogenesis associated with

library(dplyr)
library(igraph)
library(betareg)
library(ggeffects)
library(ggplot2)
library(reshape2)
library(stringr)

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
pOG <- read.delim("~/Documents/Data/parasite_orthogroups.txt", 
                                   stringsAsFactors=FALSE)
hOG <- read.delim("~/Documents/Data/host_orthogroups.txt", 
                                   stringsAsFactors=FALSE)
pb_bp <- loadRData("ERP004598/cor/ERP004598_all_bipartite.RData")#; colnames(pb_bp)[1] <- "host"; colnames(pb_bp)[2] <- "para"

d <- data.frame(h = as.character(pb_bp[,1]), p = as.character(pb_bp[,2]))
ig <- graph_from_data_frame(d, directed = F)

dg <- degree(ig, v = V(ig), loops = F, normalized = F)

dg_p <- dg[grep(pattern = "p_OG", names(dg))]
dg_p <- as.data.frame(dg_p) %>%
  tibble::rownames_to_column("Orthogroup")

dg_h <- dg[grep(pattern = "h_OG", names(dg))]
dg_h <- as.data.frame(dg_h) %>%
  tibble::rownames_to_column("Orthogroup")

para_high_dg <- dg_p[which(dg_p[,2]>=100),]

high_degree_parasite_gene <- list()

for(i in 1:length(para_high_dg))
{
	p_ortho <- para_high_dg[i]
	pg <- pOG[grep(pattern = p_ortho, pOG$Orthogroup), "Pb_g"]
	df <- data.frame()
	a = 1

	# get all host interactors
	interactors <- pb_bp[grep(pattern = pg, pb_bp[,2]),1]

	# GO function for host
	allres <- hostGOenr(interactors)%>%
				filter(KS <= 0.05)
	#allres$Study <- rep("ERP004598", nrow(allres))
	high_degree_parasite_gene[[i]] <- allres
	names(high_degree_parasite_gene)[i] <- paste0(pg, "_", p_ortho)
}

host_high_dg <- dg_h[which(dg_h[,2]>=100),]

high_degree_host_gene <- list()

for(j in 1:length(host_high_dg))
{
	h_ortho <- host_high_dg[j]
	hg <- hOG[grep(pattern = h_ortho, hOG$Orthogroup), "m_g"]
	df <- data.frame()
	a = 1

	# get all host interactors
	interactors <- pb_bp[grep(pattern = pg, pb_bp[,1]),2]

	# GO function for host
	allres <- paraGOenr(interactors)%>%
				filter(KS <= 0.05)
	#allres$Study <- rep("ERP004598", nrow(allres))
	high_degree_host_gene[[j]] <- allres
	names(high_degree_host_gene)[j] <- paste0(hg, "_", h_ortho)
}

save(high_degree_parasite_gene, file = "high_degree_parasite_gene.RData")
save(high_degree_host_gene, file = "high_degree_host_gene.RData")