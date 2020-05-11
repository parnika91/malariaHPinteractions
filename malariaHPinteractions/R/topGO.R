#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
  
#BiocManager::install("topGO")
library(topGO)
library(dplyr)
library(Rgraphviz)
library(biomaRt)
library(org.Mm.eg.db)
library(ggplot2)
#library(org.Pf.plasmo.db)
#In the first step a convenientRobject of class topGOdata is created 
#containing all the information requiredfor  the  remaining  two  steps.


##Prepping gene sets
# load("ERP106451_SRP118996/ERP106451_SRP118996_bipartite.RData")
# h_OG <- data.frame(host_genes = unique(as.character(common_individual_bipartite_common[,1])))
# p_OG <- data.frame(para_genes = unique(as.character(common_individual_bipartite_common[,2])))
# write.table(h_OG, "ERP106451_SRP118996_concat_indi_bipartite_host_genes.txt",
#             sep = "\t", quote = F, row.names = F)
# write.table(p_OG, "ERP106451_SRP118996_concat_indi_bipartite_para_genes.txt",
#             sep = "\t", quote = F, row.names = F)
# 
## Parasite
#paraGO <- readGAF("/home/parnika/Downloads/Pberghei.gaf")
# para <- read.delim("~/topGO/Pb_annot_biomaRt.txt", header=T,
#                    stringsAsFactors=FALSE)
# parasite_orthogroups <- read.delim("~/Documents/Data/parasite_orthogroups.txt", stringsAsFactors=FALSE)
# 
# para_annot <- para[,c(1,5)]
# para_annot[,1] <- sapply(para_annot[,1], function(x)
#   strsplit(x, split = "\\.")[[1]][[1]])
# colnames(para_annot) <- c("Gene", "GO")
# para_annot <- aggregate(para_annot$GO,para_annot['Gene'],paste,collapse=', ')
# para_annot <- para_annot[which(as.character(para_annot$GO) != ""),]
# colnames(para_annot) <- NULL
# write.table(para_annot, "Pb_annot_biomaRt.txt", row.names = F, sep = "\t", quote = F)
# paste(Pb_annot[,1], "0", sep = '')
# 1. Backgound genes and their annotations
# Found via GeneDB redirecting to

# combining ensembland geneDB Pb GO annotations
# para_annot[,1] <- paste(para_annot[,1], "0", sep = '')
# colnames(para_annot) <- c("Gene", "GO")
# colnames(Pb_annot) <- c("Gene", "GO")
# 
# both <- full_join(Pb_annot, para_annot, by ="Gene")
# 
# save(both, file = "combinedGOensemblGeneDBberghei.rds")
# #vec <- strsplit(paste(both$GO.x, both$GO.y, sep = ', ')[2], split = ", ")[[1]]
# GOterms <- paste(both[,2], both[,3], sep = ', ')
# GOtermsvec <- sapply(GOterms, function(x) strsplit(x, split = ', ')[[1]])
# GOtermsvec <- lapply(GOtermsvec, function(x) x[which(x!="NA")]) #d <- d[!is.na(d)]
# GOtermsvec <- sapply(GOtermsvec, function(x) unique(x))
# GOtermsvec <- lapply(GOtermsvec, function(x) paste(x, collapse = ', '))
# 
# both$GO <- GOtermsvec
# both <- both[,c(1,4)]
# both$GO <- vapply(both$GO, paste, collapse = ", ", character(1L))
# colnames(both) <- NULL
# write.table(as.data.frame(both), "Pb_annot_ens_genedb.txt", sep = '\t', row.names = F, quote = F)

# Index of ftp://ftp.sanger.ac.uk/pub/genedb/releases/latest/Pberghei/
# bip_studies <- c("SameSize_bipartite", "df_concat_allhosts_bipartite",
#                  "ERP004598_all_bipartite", "SRP118996_bipartite",
#                  "ERP106451_bipartite", "ERP106451_SRP118996_concat_bipartite",
#                  "ERP106451_SRP118996_concat_indi_bipartite")
bip_studies <- "df_concat_allhosts_bipartite"
geneont <- c("BP", "CC", "MF")

for(m in 1:length(bip_studies))
{
  for(n in 1:length(geneont))
  {
    study <- bip_studies[m]
    # study is followed by "_para_genes.txt" or "_host_genes.txt"
    # study can be
    # 1. SameSize_bipartite or SameSize
    # 2. df_concat_allhosts_bipartite or df_concat_allhosts
    # 3. ERP004598_all_bipartite or ERP004598_all
    # 4. SRP118996_bipartite
    # 5. ERP106451_bipartite
    # 6. ERP106451_SRP118996_concat_bipartite
    # 7. ERP106451_SRP118996_concat_indi_bipartite
    GeneOnt <- geneont[n]
    
    parasite_orthogroups <- read.delim("~/Documents/Data/parasite_orthogroups.txt",
                                       stringsAsFactors=FALSE) %>%
      dplyr::as_tibble() %>%
      dplyr::select(Orthogroup, Pb_g) %>%
      mutate(Pb_g = paste(Pb_g, "0", sep = ""))
    
    geneID2GO <- readMappings(file = "topGO/Pb_annot_ens_genedb.txt") # 3659 genes
    
    # 2. Make list of interesting genes. They don't need p.values
    # Pvalues would be required only to categorise interesting genes from the universe
    # Here we achieve this by 0/1 using %in%
    para_genes <- read.csv(paste0("~/Documents/Data/", study ,"_para_genes.txt",
                                  collapse = ""), sep="", stringsAsFactors=FALSE)
    colnames(para_genes)[1] <- "Orthogroup"
    para_in <- inner_join(para_genes, parasite_orthogroups)
    #para_in <- para_in[,c(1,4)]
    para_in <- unique(as.character(para_in[,2]))
    # para_in <- paste(para_in, "0", sep = "")
    # 
    # # 3. Interesting genes
    # # para_in <- inner_join(ERP004598_all_bipartite_paragenes, parasite_orthogroups)
    # # para_in <- para_in[,c(1,4)]
    # # para_in <- unique(as.character(para_in[,2]))
    # # para_in <- as.data.frame(para_in)
    # # para_in$p.value <- rep(0.000001, times = nrow(para_in))
    # # para_int <- as.character(para_in[,1])
    # # para_in <- para_in$p.value
    # # para_names <- para_int
    # # names(para_in) <- para_in$para_names
    # # Make topGOdata
    # 
    #3. To know which genes are interesting in the universe, we do %in% with background genes
    para_bg <- as.data.frame(names(geneID2GO))
    colnames(para_bg) <- "Pb_g"
    para_annot <- parasite_orthogroups[,2]
    para_bg <- unlist(c(inner_join(para_bg, para_annot)))
    
    # universe containing annotated genes out of 4010 orthogroups
    # para_bg <- gsub("\"", "", para_bg, fixed = T)
    geneList = factor(as.integer(para_bg %in% para_in))
    names(geneList)= para_bg

    GOdata <- new("topGOdata",
                  ontology = GeneOnt,
                  allGenes = geneList,
                  annot = annFUN.gene2GO,
                  #geneSelectionFun = function(x)x==1,
                  gene2GO = geneID2GO,
                  nodeSize = 1)
    #Expected: Under random chance, number of genes that would be expected 
    # to be significantly DE and annotated with that term
    # The column Expected represents the expected number of interesting genes mapped to the 
    # GO term if the interesting genes were randomly distributed over all GO terms.
    
    # sg <- sigGenes(GOdata)
    resultKS=runTest(GOdata, algorithm='weight01', statistic='KS') 
    #resultFisher=runTest(GOdata, algorithm='weight01', statistic='Fisher') 
    allGO=usedGO(GOdata)
    all_res=GenTable(GOdata, KS=resultKS, orderBy="KS", topNodes=length(allGO))
    par(cex = 0.4)
    # Plotting results
    showSigOfNodes(GOdata, score(resultKS), firstSigNodes = 5, useInfo ='all')
    # red: significant, yellow = not sig
    printGraph(GOdata, resultKS, firstSigNodes = 5, fn.prefix = paste0("tGO_", study,"_parasite_", GeneOnt, collapse = ''), useInfo = "all", pdfSW = TRUE)
    
    # Get genes in a particular GO term:
    GenesForGOterm <- c()
    myterms = all_res$GO.ID
    mygenes <- genesInTerm(GOdata, myterms)
    for (i in 1:length(myterms))
    {
      myterm <- mygenes[myterms[i]][[1]]
      mygenesforterm <- myterm[which(myterm %in% para_in == TRUE)]
      mygenesforterm <- paste(mygenesforterm, collapse=',')
      GenesForGOterm[i] <- mygenesforterm
    }

    all_res$GenesForGOterm <- GenesForGOterm
    write.table(all_res, paste0("Pb_topGO_", GeneOnt,"_", study, "_para_result.txt", collapse = ''), sep = '\t', row.names = F)

    results.ks <- runTest(GOdata, algorithm="weight01", statistic="ks")
    goEnrichment <- GenTable(GOdata, KS=results.ks, orderBy="KS", topNodes=20)
    goEnrichment <- goEnrichment[goEnrichment$KS<0.05,]
    goEnrichment <- goEnrichment[,c("GO.ID","Term","KS", "Significant")]
    goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
    goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
    goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
    goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))
    goEnrichment$KS <- as.numeric(goEnrichment$KS)
    # 
    # b <- read.delim("~/topGO/Mmus_topGO_BP_ERP106451_SRP118996_concat_bipartite_host_result.txt",
    #                 stringsAsFactors=FALSE)
    # b <- b[-which(b$Significant == 0),]
    # b <- b[-which(b$Significant == 1),]
    # b <- b[b$KS<0.05,]
    # b <- b[,c("GO.ID","Term","KS", "Significant")]
    # b$Term <- gsub(" [a-z]*\\.\\.\\.$", "", b$Term)
    # b$Term <- gsub("\\.\\.\\.$", "", b$Term)
    # b$Term <- paste(b$GO.ID, b$Term, sep=", ")
    # b$Term <- factor(b$Term, levels=rev(b$Term))
    # b$KS <- as.numeric(b$KS)
    
    ggplot(goEnrichment, aes(x=Term, y=Significant, fill = KS)) +
      stat_summary(geom = "bar", width = 0.3, fun = mean, position = "dodge") +
      xlab(GeneOnt) +
      ylab("Number of parasite genes in GO term") +
      ggtitle(study) +
      #scale_y_continuous(breaks = round(seq(0, max(-log10(goEnrichment$KS)), by = 2), 1)) +
      theme_bw(base_size=14) +
      theme(
        #legend.position='none',
        #legend.background=element_rect(),
        plot.title=element_text(angle=0, size=10, face="bold", vjust=1),
        axis.text.x=element_text(angle=0, size=10, face="bold", hjust=1.10),
        axis.text.y=element_text(angle=0, size=10, face="bold", vjust=0.5),
        axis.title=element_text(size=10, face="bold"),
        #legend.key=element_blank(),     #removes the border
        #legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
        #legend.text=element_text(size=8),  #Text size
        title=element_text(size=8)) +
      #guides(colour=guide_legend(override.aes=list(size=2.5))) +
      coord_flip()
      ggsave(paste0(study, "_", GeneOnt, "_parasite_Enrichment.png"), height = 20, width = 20, units = "cm")
    # Host
    # 1. gene-to-GO annotation file
    
    # 2. Background genes
    ############################################## good code ################
    host_genes <- read.csv(paste0("~/Documents/Data/", study ,"_host_genes.txt",
                                  collapse = ""), sep="",
                           stringsAsFactors=FALSE)
    host_orthogroups <- read.delim("~/Downloads/host_orthogroups.txt", stringsAsFactors=FALSE)
    colnames(host_genes)[1] <- "Orthogroup"
    host_in <- inner_join(host_genes, host_orthogroups)
    host_in <- host_in[,c(1,3)]
    host_in <- unique(as.character(host_in[,2]))
    
    host_bg <- host_orthogroups[,3]
    ########################################################################
    
    # getBM() wasn't returning results
    # ensembl <- useMart("ensembl")
    # datasets <- listDatasets(ensembl)
    # ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
    # go_ids = getBM(attributes=c('go_id', 'external_gene_name', 'namespace_1003'), 
    #               filters='external_gene_name', values=host_bg, mart=ensembl)
    
    # downloaded dataset from biomart webpage
    # Mmus_annot <- read.delim("~/Downloads/Mmus_annot.txt", stringsAsFactors=FALSE)
    # Mmus_annot <- Mmus_annot[-which(Mmus_annot[,3]==""),c(1,3)]
    # Mmus_annot <- aggregate(Mmus_annot$GO.term.accession, Mmus_annot['Gene.stable.ID'],paste,collapse=', ')
    # write.table(Mmus_annot, "Mmus_agg_annot.txt", row.names = F, sep = "\t", quote = F)
    
    ##################### good code ##################################
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
    
    sg <- sigGenes(hGOdata)
    resultKS=runTest(hGOdata, algorithm='weight01', statistic='KS')
    allGO=usedGO(hGOdata)
    all_res=GenTable(hGOdata, KS=resultKS, orderBy="KS", topNodes=length(allGO))
    ################################################################
    
    # Plotting results
    #par(cex = 0.4)
    #showSigOfNodes(hGOdata, score(resultKS), firstSigNodes = 5, useInfo ='all')
    # # # red: significant, yellow = not sig
    #printGraph(hGOdata, resultKS, firstSigNodes = 5, fn.prefix = paste0("tGO_", study,"_host_", GeneOnt, collapse = ''), useInfo = "all", pdfSW = TRUE)
    # 
    # # Get genes in a particular GO term:
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
    write.table(all_res, paste0("Mmus_topGO_", GeneOnt,"_", study, "_host_result.txt", collapse = ''), sep = '\t', row.names = F)
    
    results.ks <- runTest(hGOdata, algorithm="weight01", statistic="ks")
    goEnrichment <- GenTable(hGOdata, KS=results.ks, orderBy="KS", topNodes=20)
    goEnrichment <- goEnrichment[goEnrichment$KS<0.05,]
    goEnrichment <- goEnrichment[,c("GO.ID","Term","KS", "Significant")]
    goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
    goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
    goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
    goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))
    goEnrichment$KS <- as.numeric(goEnrichment$KS)
    
    ggplot(goEnrichment, aes(x=Term, y=Significant, fill = KS)) +
      stat_summary(geom = "bar", width = 0.3, fun = mean, position = "dodge") +
      xlab(GeneOnt) +
      ylab("Number of host genes in GO term") +
      ggtitle(study) +
      #scale_y_continuous(breaks = round(seq(0, max(-log10(goEnrichment$KS)), by = 2), 1)) +
      theme_bw(base_size=14) +
      theme(
        #legend.position='none',
        #legend.background=element_rect(),
        plot.title=element_text(angle=0, size=10, face="bold", vjust=1),
        axis.text.x=element_text(angle=0, size=10, face="bold", hjust=1.10),
        axis.text.y=element_text(angle=0, size=10, face="bold", vjust=0.5),
        axis.title=element_text(size=10, face="bold"),
        #legend.key=element_blank(),     #removes the border
        #legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
        #legend.text=element_text(size=8),  #Text size
        title=element_text(size=8)) +
      #guides(colour=guide_legend(override.aes=list(size=2.5))) +
      coord_flip()
    ggsave(paste0(study, "_", GeneOnt, "_host_Enrichment.png"), height = 20, width = 20, units = "cm")
    
  }
}
