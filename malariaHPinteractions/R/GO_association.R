# script to find associations between genes
# From the parasite side:
# get all the enriched GO terms, eg. pathogenesis
# get all parasite genes in there, find out all the host interactors, do GO term analysis of host GO terms
# gives us what host GO terms pathogenesis is associated with

# From ov

library(topGO)
library(dplyr)
library(Rgraphviz)
library(biomaRt)
library(org.Mm.eg.db)
library(ggplot2)


library(dplyr)
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

hostGOenr <- function(host_genes)
{
	host_orthogroups <- read.delim("~/Downloads/host_orthogroups.txt", stringsAsFactors=FALSE)
    colnames(host_genes)[1] <- "Orthogroup"
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
                   ontology = GeneOnt,
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
	geneID2GO <- readMappings(file = "topGO/p_OG_GOterms.txt") # 3659 genes
    para_in <- unique(as.character(para_genes[,1]))
    para_bg <- names(geneID2GO)
    geneList = factor(as.integer(para_bg %in% para_in))
    names(geneList)= para_bg

    GOdata <- new("topGOdata",
                  ontology = GeneOnt,
                  allGenes = geneList,
                  annot = annFUN.gene2GO,
                  gene2GO = geneID2GO,
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

pOG <- read.delim("~/Documents/Data/parasite_orthogroups.txt", 
                                   stringsAsFactors=FALSE)
hOG <- read.delim("~/Documents/Data/host_orthogroups.txt", 


p_GO <- read.table(paste0("p_OG_topGO_BP_overall_addblood_para_result.txt", collapse = ""), header = T, sep = '\t') %>%
				filter(KS <= 0.05)
h_GO <- read.table(paste0("Mmus_topGO_BP_overall_addblood_host_result.txt", collapse = ''), header = T, sep = '\t') %>%
					filter(KS <= 0.05)
ov_ad_bp <- loadRData("overall_addblood/cor/overall_addblood_all_bipartite.RData"); colnames(pb_bp)[1] <- "host"; colnames(pb_bp)[2] <- "para"


GO_asso_from_parasite <- list()
for(i in 1:nrow(p_GO))
{
	genes_in_GO_term <- p_GO[i,"GenesForGOterm"]

	# Find their interactors
	hg_vector <- c()
	for(j in 1:length(genes_in_GO_term))
	{
		pg <- genes_in_GO_term[j]
		p_ortho <- pOG[grep(pattern = pg, pOG$Pb_g), "Orthogroup"]

		# find its interactors
		hg_vector <- c(hg_vector, ov_ad_bp[grep(pattern = p_ortho, ov_ad_bp$para), "host"])
	}

	# Go analysis of host genes
	allres <- hostGOenr(hg_vector)%>%
		filter(KS <= 0.05)

	GO_asso_from_parasite[[j]] <- allres
	names(GO_asso_from_parasite)[j] <- paste0(p_GO[i,"GO.ID"], "_", p_GO[i,"Term"])
}

GO_asso_from_host <- list()
for(i in 1:nrow(h_GO))
{
	hgenes_in_GO_term <- h_GO[i,"GenesForGOterm"]

	# Find their interactors
	pg_vector <- c()
	for(j in 1:length(hgenes_in_GO_term))
	{
		hg <- hgenes_in_GO_term[j]
		h_ortho <- hOG[grep(pattern = hg, hOG$m_g), "Orthogroup"]

		# find its interactors
		pg_vector <- c(pg_vector, ov_ad_bp[grep(pattern = h_ortho, ov_ad_bp$host), "para"])

	}

	# Go analysis of host genes
	allres <- paraGOenr(pg_vector)%>%
		filter(KS <= 0.05)
		
	GO_asso_from_host[[j]] <- allres
	names(GO_asso_from_host)[j] <- paste0(h_GO[i,"GO.ID"], "_", h_GO[i,"Term"])
}

save(GO_asso_from_host, file = "GO_asso_from_host.RData")
save(GO_asso_from_parasite, file = "GO_asso_from_parasite.RData")