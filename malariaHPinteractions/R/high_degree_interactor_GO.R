# script to see what GO terms are enriched among the interactors of hogh degree genes
# 1. take the highest degree host genes, for each gene, find all parasite interactors, do GO term analysis for these interactors. same for top degree parasite genes
# 2. take all parasite genes in a GO term, get all host interctors, do GO terms analysis on these host genes - eg, what host GO terms are parasite pathogenesis associated with
library(topGO)
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

hostGOenr <- function(host_genes)
{
	host_orthogroups <- read.delim("/SAN/Plasmo_compare/OrthoFinder/host_orthogroups.txt", stringsAsFactors=FALSE)
	host_genes <- data.frame(Orthogroup = host_genes)
    #colnames(host_genes)[1] <- "Orthogroup"
    host_in <- inner_join(host_genes, host_orthogroups)
    host_in <- host_in[,c(1,3)]
    host_in <- unique(as.character(host_in[,2]))
    
    host_bg <- host_orthogroups[,3]
    h_geneList <- as.integer(host_bg %in% host_in)
    names(h_geneList) <- host_bg
    
    topDiffGenes <- function(allScore) 
    {
      return(allScore == 1)
    }
    x <- topDiffGenes(h_geneList)
    
    hGOdata <- new("topGOdata",
                   ontology = "BP",
                   allGenes = h_geneList,
                   nodeSize = 10,
                   annotationFun = annFUN.org,
                   geneSelectionFun = topDiffGenes,
                   mapping = "org.Mm.eg",
                   ID = "ensembl")
    
    resultKS=runTest(hGOdata, algorithm='weight01', statistic='KS')
    allGO=usedGO(hGOdata)
    all_res=GenTable(hGOdata, KS=resultKS, orderBy="KS", topNodes=length(allGO))
    
    GenesForGOterm <- c()
    myterms = all_res$GO.ID
    mygenes <- genesInTerm(hGOdata, myterms)
    for (i in 1:length(myterms))
    {
      myterm <- mygenes[myterms[i]][[1]]
      mygenesforterm <- myterm[which(myterm %in% host_in == TRUE)]
      mygenesforterm <- paste(mygenesforterm, collapse=',')
      GenesForGOterm[i] <- mygenesforterm
    }
    
    all_res$GenesForGOterm <- GenesForGOterm
    return(all_res)
}

paraGOenr <- function(para_genes)
{
	para_genes <- data.frame(Orthogroup = para_genes)
	pOG <- read.delim("/SAN/Plasmo_compare/OrthoFinder/parasite_orthogroups.txt", 
                                   stringsAsFactors=FALSE)
	geneID2GO <- readMappings(file = "p_OG_GOterms.txt") # 3659 genes
    para_in <- unique(as.character(para_genes[,1]))
    para_bg <- names(geneID2GO)
    geneList = factor(as.integer(para_bg %in% para_in))
    names(geneList)= para_bg

    GOdata <- new("topGOdata",
                  ontology = "BP",
                  allGenes = geneList,
                  annot = annFUN.gene2GO,
                  gene2GO = geneID2GO,
                  nodeSize = 1)
    resultKS=runTest(GOdata, algorithm='weight01', statistic='KS') 
    allGO=usedGO(GOdata)
    all_res=GenTable(GOdata, KS=resultKS, orderBy="KS", topNodes=length(allGO))
    
    GenesForGOterm <- c()
    myterms = all_res$GO.ID
    mygenes <- genesInTerm(GOdata, myterms)
    for (i in 1:length(myterms))
    {
      myterm <- mygenes[myterms[i]][[1]]
      mygenesforterm <- myterm[which(myterm %in% para_in == TRUE)]
      mygenesforterm <- sapply(mygenesforterm, 
                               function(x) pOG[grep(pattern = x, pOG$Orthogroup),"Pb_g"])
      mygenesforterm <- paste(mygenesforterm, collapse=',')
      GenesForGOterm[i] <- mygenesforterm
    }

    all_res$GenesForGOterm <- GenesForGOterm
    return(all_res)
}

pOG <- read.delim("/SAN/Plasmo_compare/OrthoFinder/parasite_orthogroups.txt", 
                                   stringsAsFactors=FALSE)
hOG <- read.delim("/SAN/Plasmo_compare/OrthoFinder/host_orthogroups.txt", 
                                  stringsAsFactors=FALSE)

### blood

bp <- loadRData("blood_core_edges.RData")#; colnames(bp)[1] <- "host"; colnames(bp)[2] <- "para"

d <- data.frame(h = as.character(bp[,1]), p = as.character(bp[,2]))
ig <- graph_from_data_frame(d, directed = F)

dg <- degree(ig, v = V(ig), loops = F, normalized = F)

dg_p <- dg[grep(pattern = "p_OG", names(dg))]
dg_p <- as.data.frame(dg_p) %>%
  tibble::rownames_to_column("Orthogroup")

dg_h <- dg[grep(pattern = "h_OG", names(dg))]
dg_h <- as.data.frame(dg_h) %>%
  tibble::rownames_to_column("Orthogroup")

para_high_dg <- dg_p[which(dg_p[,2]>=70),]

high_degree_parasite_gene <- list()

for(i in 1:nrow(para_high_dg))
{
	p_ortho <- para_high_dg[i,1]
	pg <- pOG[grep(pattern = p_ortho, pOG$Orthogroup), "Pberghei"]
	# get all host interactors
	interactors <- bp[grep(pattern = p_ortho, bp[,2]),1]

	# GO function for host
	allres <- hostGOenr(interactors)%>%
				filter(KS <= 0.05)
	#allres$Study <- rep("ERP004598", nrow(allres))
	high_degree_parasite_gene[[i]] <- allres
	names(high_degree_parasite_gene)[i] <- paste0(pg, "_", p_ortho)
}

host_high_dg <- dg_h[which(dg_h[,2]>=100),]

high_degree_host_gene <- list()

for(j in 1:nrow(host_high_dg))
{
	h_ortho <- host_high_dg[j,1]
	hg <- hOG[grep(pattern = h_ortho, hOG$Orthogroup), "mouse"]
	# get all host interactors
	interactors <- bp[grep(pattern = h_ortho, bp[,1]),2]

	# GO function for host
	allres <- paraGOenr(interactors)%>%
				filter(KS <= 0.05)
	#allres$Study <- rep("ERP004598", nrow(allres))
	high_degree_host_gene[[j]] <- allres
	names(high_degree_host_gene)[j] <- paste0(hg, "_", h_ortho)
}

save(high_degree_parasite_gene, file = "high_degree_parasite_gene.RData")
save(high_degree_host_gene, file = "high_degree_host_gene.RData")

high_degree_parasite_gene_edges <- data.frame()
a = 1
for(i in 1:length(high_degree_parasite_gene))
{
	rows <- nrow(high_degree_parasite_gene[[i]])
	high_degree_parasite_gene_edges[a:(a+rows-1),1] <- rep(str_sub(names(high_degree_parasite_gene)[i], 1, 13), rows)
	high_degree_parasite_gene_edges[a:(a+rows-1),2] <- high_degree_parasite_gene[[i]]$Term

	a = a+rows
}
colnames(high_degree_parasite_gene_edges) <- c("Parasite_gene", "Host_GO_term")

high_degree_host_gene_edges <- data.frame()
a = 1
for(i in 1:length(high_degree_host_gene))
{
	rows <- nrow(high_degree_host_gene[[i]])
	high_degree_host_gene_edges[a:(a+rows-1),1] <- rep(str_sub(names(high_degree_host_gene)[i], 1, -13), rows)
	high_degree_host_gene_edges[a:(a+rows-1),2] <- high_degree_host_gene[[i]]$Term

	a = a+rows
}
colnames(high_degree_host_gene_edges) <- c("Host_gene", "Parasite_GO_term")

save(high_degree_parasite_gene_edges, file = "high_degree_parasite_gene_edges.RData")
save(high_degree_host_gene_edges, file = "high_degree_host_gene_edges.RData")


### liver


bp <- loadRData("liver_core_edges.RData")#; colnames(bp)[1] <- "host"; colnames(bp)[2] <- "para"

d <- data.frame(h = as.character(bp[,1]), p = as.character(bp[,2]))
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
	pg <- pOG[grep(pattern = p_ortho, pOG$Orthogroup), "Pberghei"]
	df <- data.frame()
	a = 1

	# get all host interactors
	interactors <- bp[grep(pattern = pg, bp[,2]),1]

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
	hg <- hOG[grep(pattern = h_ortho, hOG$Orthogroup), "mouse"]
	df <- data.frame()
	a = 1

	# get all host interactors
	interactors <- bp[grep(pattern = pg, bp[,1]),2]

	# GO function for host
	allres <- paraGOenr(interactors)%>%
				filter(KS <= 0.05)
	#allres$Study <- rep("ERP004598", nrow(allres))
	high_degree_host_gene[[j]] <- allres
	names(high_degree_host_gene)[j] <- paste0(hg, "_", h_ortho)
}

save(high_degree_parasite_gene, file = "high_degree_parasite_gene.RData")
save(high_degree_host_gene, file = "high_degree_host_gene.RData")