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

hostGOenr <- function(host_genes, GO)
{
	host_orthogroups <- read.delim("host_orthogroups.txt", stringsAsFactors=FALSE)
	host_genes <- data.frame(Orthogroup = host_genes)
    colnames(host_genes)[1] <- "Orthogroup"
    host_in <- inner_join(host_genes, host_orthogroups)
    host_in <- host_in[,c(1,3)]
    host_in <- unique(as.character(host_in[,2]))
    host_bg <- host_orthogroups[,3]
    h_geneList <- as.integer(host_bg %in% host_in)
    names(h_geneList) <- host_bg

    # universe for core blood 
    # host_uni <- loadRData("host_universe_for_blood_core.RData")
    # host_bg <- host_orthogroups$Orthogroup
    # h_bg <- intersect(host_uni, host_bg)
    # h_bg <- host_orthogroups[host_orthogroups$Orthogroup%in%h_bg,"mouse"]
    # host_in <- host_genes[,1]
    # host_in <- host_orthogroups[host_orthogroups$Orthogroup%in%host_in,"mouse"]
    # h_geneList <- as.factor(as.integer(h_bg %in% host_in))
    # h_genenames <- sapply(h_bg, function(x) host_orthogroups[grep(pattern = x, host_orthogroups$Orthogroup), 3])
    # names(h_geneList) <- h_genenames
    # 
    
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
    
    resultKS=runTest(hGOdata, algorithm='weight01', statistic='Fisher')
    allGO=usedGO(hGOdata)
    all_res=GenTable(hGOdata, Fisher=resultKS, orderBy="Fisher", topNodes=length(allGO))
    #par(cex = 0.4)
    # Plotting results
    #showSigOfNodes(hGOdata, score(resultKS), firstSigNodes = 5, useInfo ='all')
    # red: significant, yellow = not sig
    #printGraph(hGOdata, resultKS, firstSigNodes = length(allGO), fn.prefix = paste0("tGO_", GO,"_host_BP", collapse = ''), useInfo = "all", pdfSW = TRUE)
    
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

paraGOenr <- function(para_genes, GO)
{
	para_genes <- data.frame(Orthogroup = para_genes)
	pOG <- read.delim("parasite_orthogroups.txt", 
                                   stringsAsFactors=FALSE)
	geneID2GO <- readMappings(file = "~/topGO/p_OG_GOterms.txt") # 3659 genes
  all_res <- data.frame(GO.ID = "GO:000000", Term = "xyz", Annotated = 0, Significant = 0, Expected = 0, Fisher = 1, GenesForGOterm = "")
	#normal:
    para_in <- unique(as.character(para_genes[,1]))
    para_bg <- names(geneID2GO)
    geneList = factor(as.integer(para_bg %in% para_in))

    if(nlevels(factor(as.integer(para_bg %in% para_in))) == 2)
    {
    names(geneList)= para_bg

    # universe for core network: 
    # para_uni <- loadRData("para_universe_for_blood_core.RData")
    # para_in <- unique(as.character(para_genes[,1]))
    # bg <- intersect(para_uni, names(geneID2GO))
    # geneList = factor(as.integer(bg %in% para_in))
    # names(geneList)= bg
    # 

    GOdata <- new("topGOdata",
                  ontology = "BP",
                  allGenes = geneList,
                  annot = annFUN.gene2GO,
                  gene2GO = geneID2GO,
                  nodeSize = 1)
    resultKS=runTest(GOdata, algorithm='weight01', statistic='Fisher') 
    allGO=usedGO(GOdata)
    all_res=GenTable(GOdata, Fisher=resultKS, orderBy="Fisher", topNodes=length(allGO))
    #par(cex = 0.4)
    # Plotting results
    #showSigOfNodes(GOdata, score(resultKS), firstSigNodes = length(allGO), useInfo ='all')
    # red: significant, yellow = not sig
    #printGraph(GOdata, resultKS, firstSigNodes = 5, fn.prefix = paste0("tGO_", GO,"_parasite_BP", collapse = ''), useInfo = "all", pdfSW = TRUE)
    
    
    GenesForGOterm <- c()
    myterms = all_res$GO.ID
    mygenes <- genesInTerm(GOdata, myterms)
    for (i in 1:length(myterms))
    {
      myterm <- mygenes[myterms[i]][[1]]
      mygenesforterm <- myterm[which(myterm %in% para_in == TRUE)]
      mygenesforterm <- sapply(mygenesforterm, 
                               function(x) pOG[grep(pattern = x, pOG$Orthogroup),"Pberghei"])
      mygenesforterm <- paste(mygenesforterm, collapse=',')
      GenesForGOterm[i] <- mygenesforterm
    }

    all_res$GenesForGOterm <- GenesForGOterm
  }
    return(all_res)
  
}

pOG <- read.delim("parasite_orthogroups.txt", 
                                   stringsAsFactors=FALSE)
hOG <- read.delim("host_orthogroups.txt", stringsAsFactors=FALSE)


p_GO <- read.table(paste0("~/p_OG_topGO_BP_blood_overall_para_result.txt", collapse = ""), header = T, sep = '\t') %>%
				dplyr::filter(KS <= 0.05)
h_GO <- read.delim(paste0("~/Mmus_topGO_BP_blood_overall_host_result.txt", collapse = ''), header = T, sep = '\t', stringsAsFactors=FALSE) %>%
					dplyr::filter(KS <= 0.05)
#ov_ad_bp <- loadRData("blood_core_network_nomo.RData"); colnames(ov_ad_bp)[1] <- "host"; colnames(ov_ad_bp)[2] <- "para"
ov_ad_bp <- loadRData("blood_all_bipartite.RData"); colnames(ov_ad_bp)[1] <- "host"; colnames(ov_ad_bp)[2] <- "para"


GO_asso_from_parasite <- list()
for(i in 1:nrow(p_GO))
{
	genes_in_GO_term <- as.character(p_GO[i,"GenesForGOterm"])
	if(genes_in_GO_term != "")
	{
	  genes_in_GO_term <- strsplit(genes_in_GO_term, split = ",")[[1]]

	# Find their interactors
  	hg_vector <- c()
  	for(j in 1:length(genes_in_GO_term))
  	{
  		pg <- genes_in_GO_term[j]
  		p_ortho <- pOG[grep(pattern = pg, pOG$Pberghei), "Orthogroup"]
  
  		# find its interactors
  		hg_vector <- c(hg_vector, as.character(ov_ad_bp[grep(pattern = p_ortho, ov_ad_bp$para), "host"]))
  	}
	

	# Go analysis of host genes
	allres <- hostGOenr(hg_vector, GO = p_GO[i,1])%>%
		dplyr::filter(Fisher <= 0.05)

	GO_asso_from_parasite[[i]] <- allres
	names(GO_asso_from_parasite)[i] <- paste0(p_GO[i,"GO.ID"], "_", p_GO[i,"Term"])
	}
}
save(GO_asso_from_parasite, file = "GO_asso_from_parasite.RData")


GO_asso_from_host <- list()
for(i in 1:nrow(h_GO))
{
  print(i)
	hgenes_in_GO_term <- as.character(h_GO[i,"GenesForGOterm"])
	if(hgenes_in_GO_term != "")
	{ 
	hgenes_in_GO_term <- strsplit(hgenes_in_GO_term, split = ",")[[1]]


	# Find their interactors
	pg_vector <- c()
	for(j in 1:length(hgenes_in_GO_term))
	{
		hg <- as.character(hgenes_in_GO_term[j])
	 
  		h_ortho <- hOG[grep(pattern = hg, hOG$mouse), "Orthogroup"]
  
  		# find its interactors
  		pg_vector <- c(pg_vector, as.character(ov_ad_bp[grep(pattern = h_ortho, ov_ad_bp$host), "para"]))

  	}

  	# Go analysis of host genes
  	allres <- paraGOenr(pg_vector, GO = h_GO[i,1])%>%
  		filter(Fisher <= 0.05)
  		
  	GO_asso_from_host[[i]] <- allres
  	names(GO_asso_from_host)[i] <- paste0(h_GO[i,"GO.ID"], "_", h_GO[i,"Term"])
	}
}

save(GO_asso_from_host, file = "GO_asso_from_host.RData")

# make GO term edges to vis as networks

GO_asso_from_parasite_edges <- data.frame()
a = 1
for(i in 1:length(GO_asso_from_parasite))
{
	rows <- nrow(GO_asso_from_parasite[[i]])
	if(length(rows) != 0)
	{GO_asso_from_parasite_edges[a:(a+rows-1),1] <- rep(strsplit(names(GO_asso_from_parasite)[i], split = "_")[[1]][2], rows)
	GO_asso_from_parasite_edges[a:(a+rows-1),2] <- GO_asso_from_parasite[[i]]$Term

	a = a+rows}
}

g = 1
GO_asso_list = list()
for(k in 1:length(GO_asso_from_host))
{
  if(!is.null(GO_asso_from_host[[k]]))
    {
      GO_asso_list[[g]] <- GO_asso_from_host[[k]]
      names(GO_asso_list)[[g]] <- names(GO_asso_from_host)[[k]]
      g = g+1
    }
}

GO_asso_from_host_edges <- data.frame()
a = 1
for(i in 1:length(GO_asso_list))
{
  if(nrow(GO_asso_list[[i]]) > 0)
  {
   rows <- nrow(GO_asso_list[[i]])
	if(!is.null(rows))
	{
	GO_asso_from_host_edges[a:(a+rows-1),1] <- rep(strsplit(names(GO_asso_list)[i], split = "_")[[1]][2], rows)
	GO_asso_from_host_edges[a:(a+rows-1),2] <- GO_asso_list[[i]]$Term

	a = a+rows
	}
}
}

adhesion_edges <- data.frame()
a = 1
for(i in 1:length(adhesion))
{
  if(nrow(adhesion[[i]]) > 0 | !is.null(adhesion[[i]]))
  {
    rows <- nrow(adhesion[[i]])
    if(!is.null(rows))
    {
      adhesion_edges[a:(a+rows-1),1] <- rep(strsplit(names(adhesion)[i], split = "_")[[1]][2], rows)
      adhesion_edges[a:(a+rows-1),2] <- adhesion[[i]]$Term
      
      a = a+rows
    }
  }
}

GO_asso_from_parasite_edges[,1] <- paste("p_", GO_asso_from_parasite_edges[,1], sep ="")
GO_asso_from_parasite_edges[,2] <- paste("h_", GO_asso_from_parasite_edges[,2], sep ="")

GO_asso_from_host_edges[,1] <- paste("h_", GO_asso_from_host_edges[,1], sep ="")
GO_asso_from_host_edges[,2] <- paste("p_", GO_asso_from_host_edges[,2], sep ="")

save(GO_asso_from_parasite_edges, file = "GO_asso_from_parasite_edges.RData")
save(GO_asso_from_host_edges, file = "GO_asso_from_host_edges.RData")
save(adhesion_edges, file = "adhesion_edges.RData")

GO_asso_overall <- full_join(GO_asso_overall_from_host_edges, GO_asso_overall_from_parasite_edges)

############################################################################################################################################


#### liver

liver_bp <- loadRData("liver.int.overall/cor/liver.int.overall_bipartite.RData"); colnames(liver_bp)[1] <- "host"; colnames(liver_bp)[2] <- "para"
l_p_GO <- read.table(paste0("p_OG_topGO_BP_liver.int.overall_para_result.txt", collapse = ""), header = T, sep = '\t') %>%
				filter(KS <= 0.05)
l_h_GO <- read.table(paste0("Mmus_topGO_BP_liver.int.overall_host_result.txt", collapse = ''), header = T, sep = '\t') %>%
					filter(KS <= 0.05)

liver_GO_asso_from_parasite <- list()
for(i in 1:nrow(l_p_GO))
{
	genes_in_GO_term <- l_p_GO[i,"GenesForGOterm"]
	genes_in_GO_term <- strsplit(genes_in_GO_term, split = ",")[[1]]

	# Find their interactors
	if(length(genes_in_GO_term) > 0)
	{
		hg_vector <- c()
	for(j in 1:length(genes_in_GO_term))
	{
		pg <- genes_in_GO_term[j]
		p_ortho <- pOG[grep(pattern = pg, pOG$Pberghei), "Orthogroup"]

		# find its interactors
		hg_vector <- c(hg_vector, as.character(liver_bp[grep(pattern = p_ortho, liver_bp$para), "host"]))
		
	}


	# Go analysis of host genes
	allres <- hostGOenr(hg_vector)%>%
		filter(KS <= 0.05)

	liver_GO_asso_from_parasite[[i]] <- allres
	names(liver_GO_asso_from_parasite)[i] <- paste0(l_p_GO[i,"GO.ID"], "_", l_p_GO[i,"Term"])
}
}

liver_GO_asso_from_host <- list()
for(i in 1:nrow(l_h_GO))
{
	hgenes_in_GO_term <- l_h_GO[i,"GenesForGOterm"]
	hgenes_in_GO_term <- strsplit(hgenes_in_GO_term, split = ",")[[1]]


	# Find their interactors
	if(length(genes_in_GO_term) > 0)
	{
	pg_vector <- c()
	for(j in 1:length(hgenes_in_GO_term))
	{
		hg <- hgenes_in_GO_term[j]
		h_ortho <- hOG[grep(pattern = hg, hOG$mouse), "Orthogroup"]

		# find its interactors
		pg_vector <- c(pg_vector, as.character(liver_bp[grep(pattern = h_ortho, liver_bp$host), "para"]))

	}

	# Go analysis of host genes
	allres <- paraGOenr(pg_vector)%>%
		filter(KS <= 0.05)
		
	liver_GO_asso_from_host[[i]] <- allres
	names(liver_GO_asso_from_host)[i] <- paste0(l_h_GO[i,"GO.ID"], "_", l_h_GO[i,"Term"])
}
}

save(liver_GO_asso_from_host, file = "liver_GO_asso_from_host.RData")
save(liver_GO_asso_from_parasite, file = "liver_GO_asso_from_parasite.RData")

# make GO term edges to vis as networks

liver_GO_asso_from_parasite_edges <- data.frame()
a = 1
for(i in 1:length(liver_GO_asso_from_parasite))
{
	if(!is.null(nrow(liver_GO_asso_from_parasite[[i]]))) 
	{
		rows <- nrow(liver_GO_asso_from_parasite[[i]])
		liver_GO_asso_from_parasite_edges[a:(a+rows-1),1] <- rep(strsplit(names(liver_GO_asso_from_parasite)[i], split = "_")[[1]][2], rows)
		liver_GO_asso_from_parasite_edges[a:(a+rows-1),2] <- liver_GO_asso_from_parasite[[i]]$Term

		a = a+rows
	}
}

liver_GO_asso_from_host_edges <- data.frame()
a = 1
for(i in 1:length(liver_GO_asso_from_host))
{
	if(!is.null(nrow(liver_GO_asso_from_parasite[[i]]))) 
	{
	rows <- nrow(liver_GO_asso_from_host[i])
	liver_GO_asso_from_host_edges[a:(a+rows-1),1] <- rep(strsplit(names(liver_GO_asso_from_host)[i], split = "_")[[1]][2], rows)
	liver_GO_asso_from_host_edges[a:(a+rows-1),2] <- liver_GO_asso_from_host[[i]]$Term

	a = a+rows
}
}

save(liver_GO_asso_from_parasite_edges, file = "liver_GO_asso_from_parasite_edges.RData")
save(liver_GO_asso_from_host_edges, file = "liver_GO_asso_from_host_edges.RData")


######## for paper1, I want to show the parasite GO terms connected to host GO term cell-cell adhesion with cadherin
# 1. cell-cell adhesion by cadherin
cell_adh <- GO_asso_from_host[[93]]
cell_adh <- cell_adh[which(cell_adh$GenesForGOterm != ""),]

# # get all the genes, see how many unique genes are here
# cadh_genes <- paste(cell_adh$GenesForGOterm, collapse = ",")
# cadh_genes <- strsplit(cadh_genes, split = ',')[[1]]
# cadh_genes <- unique(cadh_genes)
# 
# ### get the interactor in overall_addblood
# host_cell_adh <- h_GO[which(h_GO$GO.ID == "GO:0044331"),]
# host_cadh_genes <- paste(host_cell_adh$GenesForGOterm, collapse = ",")
# host_cadh_genes <- strsplit(host_cadh_genes, split = ',')[[1]]
# 
# host_cadh_ortho <- sapply(host_cadh_genes, function(x) hOG[grep(pattern = x, hOG$mouse),1])
# 
# host_cadh_interactors <- data.frame()
# a = 1
# 
# for(i in 1:length(host_cadh_ortho))
# {
#   para_ortho <- ov_ad_bp[grep(pattern = host_cadh_ortho[i], ov_ad_bp[,1]),2]
#   host_cadh_interactors[a:(a+(length(para_ortho)-1)),1] <- rep(host_cadh_ortho[i], length(para_ortho))
#   host_cadh_interactors[a:(a+(length(para_ortho)-1)),2] <- para_ortho
#   a = a+length(para_ortho)
# }
# colnames(host_cadh_interactors) <- c("host", "para")
# hg <- sapply(host_cadh_interactors[,1], function(x) hOG[grep(pattern = x, hOG$Orthogroup), "mouse"])
# pg <- sapply(host_cadh_interactors[,2], function(x) pOG[grep(pattern = x, pOG$Orthogroup), "Pberghei"])
# 
# host_cadh_interactors <- data.frame(host = hg, para = pg)

cell_adh_interactions <- data.frame(host = rep("GO:0044331 Cell-cell adhesion with cadherin", nrow(cell_adh)),
                                    para = apply(cell_adh[,c(1,2)], 1, paste, collapse = " "))
write.csv2(cell_adh_interactions, "cell_adh_interactions.csv", row.names = F, quote = F)

patho <- GO_asso_from_parasite[[9]]
patho <- patho[which(patho$GenesForGOterm != ""),]
patho_interactions <- data.frame(host = rep("GO:0009405_pathogenesis", nrow(patho)),
                                    para = apply(patho[,c(1,2)], 1, paste, collapse = " "))
write.csv2(patho_interactions, "patho_interactions.csv", row.names = F, quote = F)
