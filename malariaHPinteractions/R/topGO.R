#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
  
#BiocManager::install("topGO")
library(topGO)
library(dplyr)
library(Rgraphviz)
library(biomaRt)
library(org.Mm.eg.db)
library(ggplot2)

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

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
Pb_Gdb <- read.delim("~/topGO/Pberghei.gaf", header=FALSE, 
                     stringsAsFactors=FALSE)
Pb_bm <- read.delim("~/topGO/Pb_mart_export.txt", header=T,
                   stringsAsFactors=FALSE)

Pf_Gdb <- read.delim("~/topGO/Pfalciparum.gaf", header=FALSE, 
                     stringsAsFactors=FALSE)
Pf_bm <- read.delim("~/topGO/Pf_mart_export.txt", header=T,
                    stringsAsFactors=FALSE)

Pv_Gdb <- read.delim("~/topGO/PvivaxP01.gaf", header=FALSE, 
                     stringsAsFactors=FALSE)
Pv_bm <- read.delim("~/topGO/Pv_mart_export.txt", header=T,
                    stringsAsFactors=FALSE)

pOG <- read.delim("~/Documents/Data/parasite_orthogroups.txt", 
                                   stringsAsFactors=FALSE)

# treating Pf tables
# remove rows with empty GO accession from both
Pf_bm <- Pf_bm[which(Pf_bm$GO.term.accession != ""),]
Pf_Gdb <- Pf_Gdb[which(Pf_Gdb$V5 != ""),]

# remove .1 from Pf_Gdb
Pf_Gdb$V2 <- substr(Pf_Gdb$V2, 1, 13)

# collect rows wih gene IDs and GO term
Pf_bm_red <- Pf_bm[,c("Gene.stable.ID", "GO.term.accession")]
Pf_Gdb_red <- Pf_Gdb[,c("V2", "V5")]
colnames(Pf_Gdb_red) <- colnames(Pf_bm_red)

# rbind bm and Gdb
Pf <- rbind(Pf_bm_red, Pf_Gdb_red)
Pf <- aggregate(Pf$GO.term.accession,Pf['Gene.stable.ID'],paste,collapse=',')
 
# keep unique GOterms
Pf$GO <- sapply(Pf$x, function(y) paste(unique(strsplit(y, split = ",")[[1]]), collapse = ","))
Pf <- Pf[, -c(2)]

# treating Pb tables
# remove rows with empty GO accession from both
Pb_bm <- Pb_bm[which(Pb_bm$GO.term.accession != ""),]
Pb_Gdb <- Pb_Gdb[which(Pb_Gdb$V5 != ""),]

# remove .1 from Pb_Gdb
Pb_Gdb$V2 <- substr(Pb_Gdb$V2, 1, 13)

# collect rows wih gene IDs and GO term
Pb_bm_red <- Pb_bm[,c("Gene.stable.ID", "GO.term.accession")]
Pb_Gdb_red <- Pb_Gdb[,c("V2", "V5")]
colnames(Pb_Gdb_red) <- colnames(Pb_bm_red)

# rbind bm and Gdb
Pb <- rbind(Pb_bm_red, Pb_Gdb_red)
Pb <- aggregate(Pb$GO.term.accession,Pb['Gene.stable.ID'],paste,collapse=',')

# keep unique GOterms
Pb$GO <- sapply(Pb$x, function(y) paste(unique(strsplit(y, split = ",")[[1]]), collapse = ","))
Pb <- Pb[, -c(2)]

# treating Pv tables
# remove rows with empty GO accession from both
Pv_bm <- Pv_bm[which(Pv_bm$GO.term.accession != ""),]
Pv_Gdb <- Pv_Gdb[which(Pv_Gdb$V5 != ""),]

# remove .1 from Pb_Gdb
Pv_Gdb$V2 <- substr(Pv_Gdb$V2, 1, 13)

# collect rows wih gene IDs and GO term
Pv_bm_red <- Pv_bm[,c("Gene.stable.ID", "GO.term.accession")]
Pv_Gdb_red <- Pv_Gdb[,c("V2", "V5")]
colnames(Pv_Gdb_red) <- colnames(Pv_bm_red)

# rbind bm and Gdb
Pv <- rbind(Pv_bm_red, Pv_Gdb_red)
Pv <- aggregate(Pv$GO.term.accession,Pv['Gene.stable.ID'],paste,collapse=',')

# keep unique GOterms
Pv$GO <- sapply(Pv$x, function(y) paste(unique(strsplit(y, split = ",")[[1]]), collapse = ","))
Pv <- Pv[, -c(2)]

save(Pf, file = "Pf_annot.RData")
save(Pb, file = "Pb_annot.RData")
save(Pv, file = "Pv_annot.RData")

# merge with orthogroups
Pb_OG <- merge(Pb, pOG[,c("Orthogroup", "Pb_g")], by.x = "Gene.stable.ID", by.y = "Pb_g")
Pf_OG <- merge(Pf, pOG[,c("Orthogroup", "Pf_g")], by.x = "Gene.stable.ID", by.y = "Pf_g")
Pv_OG <- merge(Pv, pOG[,c("Orthogroup", "Pv_g")], by.x = "Gene.stable.ID", by.y = "Pv_g")

Pb_Pf <- full_join(Pb_OG, Pf_OG, by = "Orthogroup")
Pb_Pf_Pv <- full_join(Pb_Pf, Pv_OG, by = "Orthogroup")

# keeping only orthogroups
GO <- data.frame(paste(as.character(Pb_Pf_Pv$GO.x), as.character(Pb_Pf_Pv$GO.y), as.character(Pb_Pf_Pv$GO), sep = ","))
GO_ <- data.frame()
for(i in 1:nrow(GO))
{
  t <- unique(strsplit(as.character(GO[i,]), split = ",")[[1]])
  t <- t[!t %in% "NA"]
  GO_[i,1] <- paste(t, collapse = ",")
}

Pb_Pf_Pv_OG <- data.frame(Orthogroup = Pb_Pf_Pv$Orthogroup, GOterm = GO_)
colnames(Pb_Pf_Pv_OG)[2] <- "GO"
save(Pb_Pf_Pv_OG, file = "p_OG_GOterms.RData")
colnames(Pb_Pf_Pv_OG) <- NULL

write.table(Pb_Pf_Pv_OG, "p_OG_GOterms.txt", sep = '\t', row.names = F, quote = F)

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

DRP000987 <- loadRData("/SAN/Plasmo_compare/SRAdb/Output/DRP000987/cor/DRP000987_str_bipartite.RData")
DRP000987_para_genes <- unique(as.character(DRP000987[,2]))
DRP000987_host_genes <- unique(as.character(DRP000987[,1]))
write.table(DRP000987_para_genes, "DRP000987_para_genes.txt", quote = F, row.names = F)
write.table(DRP000987_host_genes, "DRP000987_host_genes.txt", quote = F, row.names = F)

ERP110375 <- loadRData("/SAN/Plasmo_compare/SRAdb/Output/ERP110375/cor/ERP110375_int_bipartite.RData")
ERP110375_para_genes <- unique(as.character(ERP110375[,2]))
ERP110375_host_genes <- unique(as.character(ERP110375[,1]))
write.table(ERP110375_para_genes, "ERP110375_para_genes.txt", quote = F, row.names = F)
write.table(ERP110375_host_genes, "ERP110375_host_genes.txt", quote = F, row.names = F)

SRP118827 <- loadRData("/SAN/Plasmo_compare/SRAdb/Output/SRP118827/cor/SRP118827_int_bipartite.RData")
SRP118827_para_genes <- unique(as.character(SRP118827[,2]))
SRP118827_host_genes <- unique(as.character(SRP118827[,1]))
write.table(SRP118827_para_genes, "SRP118827_para_genes.txt", quote = F, row.names = F)
write.table(SRP118827_host_genes, "SRP118827_host_genes.txt", quote = F, row.names = F)

SRP116793 <- loadRData("/SAN/Plasmo_compare/SRAdb/Output/SRP116793/cor/SRP116793_all_bipartite.RData")
SRP116793_para_genes <- unique(as.character(SRP116793[,2]))
SRP116793_host_genes <- unique(as.character(SRP116793[,1]))
write.table(SRP116793_para_genes, "SRP116793_para_genes.txt", quote = F, row.names = F)
write.table(SRP116793_host_genes, "SRP116793_host_genes.txt", quote = F, row.names = F)

ERP106451 <- loadRData("/SAN/Plasmo_compare/SRAdb/Output/ERP106451/cor/ERP106451_int_bipartite.RData")
ERP106451_para_genes <- unique(as.character(ERP106451[,2]))
ERP106451_host_genes <- unique(as.character(ERP106451[,1]))
write.table(ERP106451_para_genes, "ERP106451_para_genes.txt", quote = F, row.names = F)
write.table(ERP106451_host_genes, "ERP106451_host_genes.txt", quote = F, row.names = F)

ERP004598 <- loadRData("/SAN/Plasmo_compare/SRAdb/Output/ERP004598/cor/ERP004598_all_bipartite.RData")
ERP004598_para_genes <- unique(as.character(ERP004598[,2]))
ERP004598_host_genes <- unique(as.character(ERP004598[,1]))
write.table(ERP004598_para_genes, "ERP004598_para_genes.txt", quote = F, row.names = F)
write.table(ERP004598_host_genes, "ERP004598_host_genes.txt", quote = F, row.names = F)

SRP250329 <- loadRData("/SAN/Plasmo_compare/SRAdb/Output/SRP250329/cor/SRP250329_int_bipartite.RData")
SRP250329_para_genes <- unique(as.character(SRP250329[,2]))
SRP250329_host_genes <- unique(as.character(SRP250329[,1]))
write.table(SRP250329_para_genes, "SRP250329_para_genes.txt", quote = F, row.names = F)
write.table(SRP250329_host_genes, "SRP250329_host_genes.txt", quote = F, row.names = F)

ERP105548 <- loadRData("/SAN/Plasmo_compare/SRAdb/Output/ERP105548/cor/ERP105548_int_bipartite.RData")
ERP105548_para_genes <- unique(as.character(ERP105548[,2]))
ERP105548_host_genes <- unique(as.character(ERP105548[,1]))
write.table(ERP105548_para_genes, "ERP105548_para_genes.txt", quote = F, row.names = F)
write.table(ERP105548_host_genes, "ERP105548_host_genes.txt", quote = F, row.names = F)

SRP110282 <- loadRData("/SAN/Plasmo_compare/SRAdb/Output/SRP110282/cor/SRP110282_int_bipartite.RData")
SRP110282_para_genes <- unique(as.character(SRP110282[,2]))
SRP110282_host_genes <- unique(as.character(SRP110282[,1]))
write.table(SRP110282_para_genes, "SRP110282_para_genes.txt", quote = F, row.names = F)
write.table(SRP110282_host_genes, "SRP110282_host_genes.txt", quote = F, row.names = F)

SRP096160 <- loadRData("/SAN/Plasmo_compare/SRAdb/Output/SRP096160/cor/SRP096160_int_bipartite.RData")
SRP096160_para_genes <- unique(as.character(SRP096160[,2]))
SRP096160_host_genes <- unique(as.character(SRP096160[,1]))
write.table(SRP096160_para_genes, "SRP096160_para_genes.txt", quote = F, row.names = F)
write.table(SRP096160_host_genes, "SRP096160_host_genes.txt", quote = F, row.names = F)

liver.int.overall <- loadRData("/SAN/Plasmo_compare/SRAdb/Output/liver.int.overall/cor/liver.int.overall_bipartite.RData")
liver.int.overall_para_genes <- unique(as.character(liver.int.overall[,2]))
liver.int.overall_host_genes <- unique(as.character(liver.int.overall[,1]))
write.table(liver.int.overall_para_genes, "liver.int.overall_para_genes.txt", quote = F, row.names = F)
write.table(liver.int.overall_host_genes, "liver.int.overall_host_genes.txt", quote = F, row.names = F)



#bip_studies <- c("SRP250329", "ERP105548", 
 #                "SRP110282", "SRP096160",
  #               "liver.int.overall")
bip_studies <- c("DRP000987","ERP106451",
                 "ERP110375", "ERP004598",
                 "SRP118827", "SRP116793",
                 "df_concat_allhosts")
geneont <- c("BP", "CC", "MF")

for(m in 1:length(bip_studies))
{
  print(bip_studies[m])
  for(n in 1:length(geneont))
  {
    print(geneont[n])
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
    
    # parasite_orthogroups <- read.delim("~/Documents/Data/parasite_orthogroups.txt",
    #                                    stringsAsFactors=FALSE) %>%
    #   dplyr::as_tibble() %>%
    #   dplyr::select(Orthogroup, Pb_g) %>%
    #   mutate(Pb_g = paste(Pb_g, "0", sep = ""))
    
    geneID2GO <- readMappings(file = "topGO/p_OG_GOterms.txt") # 3659 genes
    
    # 2. Make list of interesting genes. They don't need p.values
    # Pvalues would be required only to categorise interesting genes from the universe
    # Here we achieve this by 0/1 using %in%
    para_genes <- read.delim(paste0("~/Documents/Data/", study ,"_para_genes.txt",
                                  collapse = ""), stringsAsFactors=FALSE, header = T)
    #para_in <- inner_join(para_genes, parasite_orthogroups)
    #para_in <- para_in[,c(1,4)]
    para_in <- unique(as.character(para_genes[,1]))
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
    para_bg <- names(geneID2GO)
    #colnames(para_bg) <- "p_OG"
    #para_annot <- parasite_orthogroups[,2]
    #para_bg <- unlist(c(inner_join(para_bg, para_annot)))
    
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
    resultKS=runTest(GOdata, algorithm='weight01', statistic='Fisher') 
    #resultFisher=runTest(GOdata, algorithm='weight01', statistic='Fisher') 
    allGO=usedGO(GOdata)
    all_res=GenTable(GOdata, Fisher=resultKS, orderBy="Fisher", topNodes=length(allGO))
    #par(cex = 0.4)
    # Plotting results
    #showSigOfNodes(GOdata, score(resultKS), firstSigNodes = 5, useInfo ='all')
    # red: significant, yellow = not sig
    #printGraph(GOdata, resultKS, firstSigNodes = 5, fn.prefix = paste0("tGO_", study,"_parasite_", GeneOnt, collapse = ''), useInfo = "all", pdfSW = TRUE)
    
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
    write.table(all_res, paste0("p_OG_topGO_", GeneOnt,"_", study, "_para_result_Fisher.txt", collapse = ''), sep = '\t', row.names = F)

    results.ks <- runTest(GOdata, algorithm="weight01", statistic="Fisher")
    goEnrichment <- GenTable(GOdata, Fisher=results.ks, orderBy="Fisher", topNodes=20)
    goEnrichment <- goEnrichment[goEnrichment$Fisher<0.05,]
    goEnrichment <- goEnrichment[,c("GO.ID","Term","Fisher", "Significant")]
    goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
    goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
    goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
    goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))
    goEnrichment$Fisher <- as.numeric(goEnrichment$Fisher)
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
    
    ggplot(goEnrichment, aes(x=Term, y=log10(Significant), fill = Fisher)) +
      stat_summary(geom = "bar", width = 0.5, fun = mean, position = "dodge") +
      xlab("Molecular Function") +
      ylab("Number of parasite genes in GO term") +
      ggtitle("Significant Plasmodium GOterms") +
      #scale_y_continuous(breaks = round(seq(0, max(-log10(goEnrichment$KS)), by = 2), 1)) +
      theme_bw(base_size=20) +
      theme(
        #legend.position='none',
        #legend.background=element_rect(),
        plot.title=element_text(angle=0, size=20, face="bold", vjust=1),
        axis.text.x=element_text(angle=0, size=20, face="bold", hjust=1.10),
        axis.text.y=element_text(angle=0, size=20, face="bold", vjust=0.5),
        axis.title=element_text(size=20, face="bold"),
        #legend.key=element_blank(),     #removes the border
        #legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
        #legend.text=element_text(size=8),  #Text size
        title=element_text(size=18)) +
      #guides(colour=guide_legend(override.aes=list(size=2.5))) +
      coord_flip()
      ggsave(paste0(study, "_", GeneOnt, "_p_OG_Enrichment_Fisher.png"), height = 30, width = 45, units = "cm")
      # Host
    # 1. gene-to-GO annotation file
    
    # 2. Background genes
    ############################################## good code ################
    host_genes <- read.delim(paste0("~/Documents/Data/", study ,"_host_genes.txt",
                                  collapse = ""), header =T,
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
    resultKS=runTest(hGOdata, algorithm='weight01', statistic='Fisher')
    allGO=usedGO(hGOdata)
    all_res=GenTable(hGOdata, Fisher=resultKS, orderBy="Fisher", topNodes=length(allGO))
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
    write.table(all_res, paste0("Mmus_topGO_", GeneOnt,"_", study, "_host_result_Fisher.txt", collapse = ''), sep = '\t', row.names = F)
    
    results.ks <- runTest(hGOdata, algorithm="weight01", statistic="Fisher")
    goEnrichment <- GenTable(hGOdata, Fisher=results.ks, orderBy="Fisher", topNodes=20)
    goEnrichment <- goEnrichment[goEnrichment$Fisher<0.05,]
    goEnrichment <- goEnrichment[,c("GO.ID","Term","Fisher", "Significant")]
    goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
    goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
    goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
    goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))
    goEnrichment$Fisher <- as.numeric(goEnrichment$Fisher)
    
    ggplot(goEnrichment, aes(x=Term, y=Significant, fill = Fisher)) +
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
    ggsave(paste0(study, "_", GeneOnt, "_host_Enrichment_Fisher.png"), height = 20, width = 20, units = "cm")
    
  }
}

