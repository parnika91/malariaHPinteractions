if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
  
#BiocManager::install("topGO")
library(topGO)
library(dplyr)
library(Rgraphviz)
#library(org.Pf.plasmo.db)
#In the first step a convenientRobject of class topGOdata is created 
#containing all the information requiredfor  the  remaining  two  steps.

# Parasite
#paraGO <- readGAF("/home/parnika/Downloads/Pberghei.gaf")
# para <- read.delim("~/Downloads/Pberghei_GO.txt", header=FALSE, 
#                    stringsAsFactors=FALSE)
# parasite_orthogroups <- read.delim("~/Documents/Data/parasite_orthogroups.txt", stringsAsFactors=FALSE)
# 
# para_annot <- para[,c(2,5)]
# para_annot[,1] <- sapply(para_annot[,1], function(x)
#   strsplit(x, split = "\\.")[[1]][[1]])
# colnames(para_annot) <- c("Gene", "GO")
# para_annot <- aggregate(para_annot$GO,para_annot['Gene'],paste,collapse=', ')
# colnames(para_annot) <- NULL
# write.table(para_annot, "Pb_annot.txt", row.names = F, sep = "\t", quote = F)

# 1. Backgound genes and their annotations
# Found via GeneDB redirecting to
# Index of ftp://ftp.sanger.ac.uk/pub/genedb/releases/latest/Pberghei/
geneID2GO <- readMappings(file = "Pb_annot.txt")

# 2. MAke list of interesting genes. They don't need p.values
# Pvalues would be required only to categorise interesting genes from the universe
# Here we achieve this by 0/1 using %in%
parasite_orthogroups <- read.delim("~/Documents/Data/parasite_orthogroups.txt", stringsAsFactors=FALSE)
ERP004598_all_bipartite_paragenes <- read.csv("~/Documents/Data/ERP004598_all_bipartite_paragenes.txt", sep="", stringsAsFactors=FALSE)
colnames(ERP004598_all_bipartite_paragenes)[1] <- "Orthogroup"
para_in <- inner_join(ERP004598_all_bipartite_paragenes, parasite_orthogroups)
para_in <- para_in[,c(1,4)]
para_in <- unique(as.character(para_in[,2]))
para_in <- paste(para_in, "0", sep = "")

# 3. Interesting genes
# para_in <- inner_join(ERP004598_all_bipartite_paragenes, parasite_orthogroups)
# para_in <- para_in[,c(1,4)]
# para_in <- unique(as.character(para_in[,2]))
# para_in <- as.data.frame(para_in)
# para_in$p.value <- rep(0.000001, times = nrow(para_in))
# para_int <- as.character(para_in[,1])
# para_in <- para_in$p.value
# para_names <- para_int
# names(para_in) <- para_in$para_names
# Make topGOdata

#3. To know which genes are interesting in the universe, we do %in% with background genes
para_bg <- names(geneID2GO)
para_bg <- gsub("\"", "", para_bg, fixed = T)

geneList=as.integer(para_bg %in% para_in)
names(geneList)= para_bg

GOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = geneList,
              annot = annFUN.gene2GO,
              geneSelectionFun = function(x)x==1,
              gene2GO = geneID2GO,
              nodeSize = 1)

sg <- sigGenes(GOdata)
resultFisher=runTest(GOdata, algorithm='weight01', statistic='fisher') 
allGO=usedGO(GOdata)
all_res=GenTable(GOdata, weightFisher=resultFisher, 
                 orderBy='weightFisher', topNodes=length(allGO))
write.table(all_res, "Pf_topGO_BP_ERP004598_para_result.txt", sep = '\t', row.names = F)

# Plotting results
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 5, useInfo ='all')

# Host
# 1. gene-to-GO annotation file

# 2. Background genes
ERP004598_all_bipartite_hostgenes <- read.csv("~/Documents/Data/ERP004598_all_bipartite_hostgenes.txt", 
                                              sep="", stringsAsFactors=FALSE)
host_orthogroups <- read.delim("~/Downloads/host_orthogroups.txt", stringsAsFactors=FALSE)
colnames(ERP004598_all_bipartite_hostgenes)[1] <- "Orthogroup"
host_in <- inner_join(ERP004598_all_bipartite_hostgenes, host_orthogroups)
#para_in <- parasite_orthogroups[,c(1,5)]
host_in <- host_in[,c(1,2)]
host_in <- unique(as.character(host_in[,2]))
Hs_geneList <- rep(0.000001, length(host_in))
names(Hs_geneList) <- host_in

x <- org.Hs.egENSEMBL
# Get the entrez gene IDs that are mapped to an Ensembl ID
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])

BiocManager::install("biomaRt")
library(biomaRt)
db= useMart('ENSEMBL_MART_ENSEMBL',dataset='hsapiens_gene_ensembl', host="www.ensembl.org")
go_ids= getBM(attributes=c('go_id', 'external_gene_name', 'namespace_1003'), filters='external_gene_name', values=host_orthogroups[,2], mart=db)

# 3. Interesting genes

