########### orthologous parasite genes ########

# get one-to-one orthologs for tree

# max.copies <- function(common)
# {
#   copies.max <- c()
#   for(i in 1:nrow(common))
#   {
#     common_singlecopy <- length(strsplit(as.character(common[i,2]), split = ",")[[1]])
#     group1_singlecopy <- length(strsplit(as.character(common[i,3]), split = ",")[[1]])
#     group2_singlecopy <- length(strsplit(as.character(common[i,4]), split = ",")[[1]])
#     copies.max[i] <- max(common_singlecopy, group1_singlecopy, group2_singlecopy)
#   }
#   return (copies.max)
# }
# 
# # first branch = Pv, Pcy and Pco
# 
# Pv_vs_Pco <-
#   read.csv2("../OrthoFinder/Parasites/Results/Orthologues_Apr10/Orthologues/Orthologues_Pvivax/Pvivax__v__Pcoatneyi.csv", sep = '\t', header = T)
# 
# Pcy_vs_Pco <- read.csv2("../OrthoFinder/Parasites/Results/Orthologues_Apr10/Orthologues/Orthologues_Pcynomolgi/Pcynomolgi__v__Pcoatneyi.csv", sep = '\t', header = T)
# 
# Pco_vs_PvPcy <- merge(Pv_vs_Pco, Pcy_vs_Pco)
#   #Pv_vs_Pco %>%
#   #filter(Orthogroup %in% Pcy_vs_Pco$Orthogroup)#
#   
# 
# Pco_vs_PvPcy$maximum.copies <- max.copies(Pco_vs_PvPcy)
# Pco_vs_PvPcy_1to1 <- Pco_vs_PvPcy[Pco_vs_PvPcy$maximum.copies==1,]
# unique_Pco_vs_PvPcy <- length(unique(Pco_vs_PvPcy_1to1[,1]))
# 
# # second branch = Pb, Py and Pch
# 
# Pb_vs_Pch <-
#   read.csv2("../OrthoFinder/Parasites/Results/Orthologues_Apr10/Orthologues/Orthologues_Pberghei/Pberghei__v__Pchabaudi.csv", sep = '\t', header = T)
# 
# Py_vs_Pch <- read.csv2("../OrthoFinder/Parasites/Results/Orthologues_Apr10/Orthologues/Orthologues_Pyoelii/Pyoelii__v__Pchabaudi.csv", sep = '\t', header = T)
# 
# Pch_vs_PbPy <- merge(Pb_vs_Pch, Py_vs_Pch)
# 
# Pch_vs_PbPy$maximum.copies <- max.copies(Pch_vs_PbPy)
# Pch_vs_PbPy_1to1 <- Pch_vs_PbPy[Pch_vs_PbPy$maximum.copies==1,]
# unique_Pch_vs_PbPy <- length(unique(Pch_vs_PbPy_1to1[,1]))
# 
# 
# all_except_Pf_1to1 <- merge(Pco_vs_PvPcy_1to1, Pch_vs_PbPy_1to1)
# all_except_Pf_1to1 <- all_except_Pf_1to1[,-which(names(all_except_Pf_1to1)%in%"maximum.copies")]
# unique_all_except_Pf <- length(unique(all_except_Pf_1to1[,1]))
# 
# Pf_vs_Pch <- read.csv2("../OrthoFinder/Parasites/Results/Orthologues_Apr10/Orthologues/Orthologues_Pfalciparum/Pfalciparum__v__Pchabaudi.csv", sep = '\t', header = T)
# Pf_vs_all <- merge(all_except_Pf_1to1, Pf_vs_Pch)
# 
# maximum.copies <- c()
# for(i in 1:nrow(Pf_vs_all))
# {
#   maximum.copies[i] <- length(strsplit(as.character(Pf_vs_all[i,8]), split = ",")[[1]])
# }
# Pf_vs_all$maximum.copies <- maximum.copies
# 
# Pf_vs_all_1to1 <- Pf_vs_all[Pf_vs_all$maximum.copies==1,]
# unique_all <- length(unique(Pf_vs_all_1to1[,1]))
# 
# all <- Pf_vs_all_1to1 %>% 
#   aggregate(list(Pf_vs_all_1to1$Orthogroup), FUN = paste) %>%
#   as.data.frame()
# 
# all_1to1 <- data.frame()
# a = 1
# for(i in 1:nrow(all))
# {
#   if(length(all[i,2][[1]])==1)
#     all_1to1[a,1:ncol(all)] <- all[i,]
#   a = a+1
# }
# 
# all_1to1 <- na.omit(all_1to1) # This does it!
# all_1to1 <- all_1to1[,-which(names(all_1to1)%in%"Group.1")]

########################## from orthogroups.csv ###########################

ortho <- readLines("../OrthoFinder/Parasites/Results/Orthogroups.csv")
para_species <- ortho[1] # extract the parasite species names
para_species <- strsplit(para_species, split = "\t")[[1]][-1]
ortho <- ortho[2:length(ortho)] # remove species names from dataframe

orthogroup <- c()
ortho_1to1 <- c()
a = 1

for(i in 1:length(ortho))
{
  # split the lines by tab
  line.split <- strsplit(ortho[i], split = "\t")
  
  # get the orthogroup
  # orthogroup[i] <- strsplit(ortho[i], split = "\t")[[1]][1] 
  
  # check the contents of the rest of the split part
  # check if there are in fact all species represented (there are 7 elements)
  # and check if the length of these elements is exactly one
  
  # get the number of genes/transcripts/proteins named for each species in that line
  split.orthologs <- sapply(sapply(line.split, function(x) strsplit(x, split = ",")), function(x) length(x))
  
  if(length(unique(split.orthologs))==1 & all(split.orthologs==1))
  {
    ortho_1to1[a] <- paste0(unlist(line.split), collapse = "\t")
    a = a+1
  }
}

# change ortho_1to1 into a data frame
ortho_1to1_df <- data.frame()

ortho_1to1_df[1:length(ortho_1to1),1] <- sapply(ortho_1to1, function(x) strsplit(x, split = "\t")[[1]][1]) # orthogroup
ortho_1to1_df[1:length(ortho_1to1),2] <- sapply(ortho_1to1, function(x) strsplit(x, split = "\t")[[1]][2])
ortho_1to1_df[1:length(ortho_1to1),3] <- sapply(ortho_1to1, function(x) strsplit(x, split = "\t")[[1]][3])
ortho_1to1_df[1:length(ortho_1to1),4] <- sapply(ortho_1to1, function(x) strsplit(x, split = "\t")[[1]][4])
ortho_1to1_df[1:length(ortho_1to1),5] <- sapply(ortho_1to1, function(x) strsplit(x, split = "\t")[[1]][5])
ortho_1to1_df[1:length(ortho_1to1),6] <- sapply(ortho_1to1, function(x) strsplit(x, split = "\t")[[1]][6])
ortho_1to1_df[1:length(ortho_1to1),7] <- sapply(ortho_1to1, function(x) strsplit(x, split = "\t")[[1]][7])
ortho_1to1_df[1:length(ortho_1to1),8] <- sapply(ortho_1to1, function(x) strsplit(x, split = "\t")[[1]][8])
  
colnames(ortho_1to1_df) <- c("Orthogroup", para_species)

######### change gene names #########
### Pf, Pv, Pb and Pch need changing from transcript ID to gene ID ###
### Pco and Pcy needs splitting at -t and Py needs splitting .1 ###

library(rtracklayer)
Pch <- as.data.frame(import("/SAN/Plasmo_compare/Genomes/annotation/Pchabaudi.gtf", format = "gtf"))
Pch <- data.frame(Pch_t = Pch$transcript_id, Pch_g = Pch$gene_id)
Pch <- unique(na.omit(Pch))

Pch_ortho <- data.frame(Orthogroup = ortho_1to1_df$Orthogroup, Pch_t = ortho_1to1_df$Pchabaudi)

ID_change <- merge(as.data.frame(Pch), as.data.frame(Pch_ortho), by = "Pch_t")

Pv <- as.data.frame(import("/SAN/Plasmo_compare/Genomes/annotation/Pvivax.gtf", format = "gtf"))
Pv <- data.frame(Pv_t = Pv$transcript_id, Pv_g = Pv$gene_id)
Pv <- unique(na.omit(Pv))

Pv_ortho <- data.frame(Orthogroup = ortho_1to1_df$Orthogroup, Pv_t = ortho_1to1_df$Pvivax)
Pv_ortho_gene <- merge(as.data.frame(Pv), as.data.frame(Pv_ortho), by = "Pv_t")
ID_change <- merge(ID_change, as.data.frame(Pv_ortho_gene), by = "Orthogroup")

Pb <- as.data.frame(import("/SAN/Plasmo_compare/Genomes/annotation/Pberghei.gtf", format = "gtf"))
Pb <- data.frame(Pb_t = Pb$transcript_id, Pb_g = Pb$gene_id)
Pb <- unique(na.omit(Pb))

Pb_ortho <- data.frame(Orthogroup = ortho_1to1_df$Orthogroup, Pb_t = ortho_1to1_df$Pberghei)
Pb_ortho_gene <- merge(as.data.frame(Pb), as.data.frame(Pb_ortho), by = "Pb_t")
ID_change <- merge(ID_change, as.data.frame(Pb_ortho_gene), by = "Orthogroup")

Pf <- as.data.frame(import("/SAN/Plasmo_compare/Genomes/annotation/Pfalciparum.gtf", format = "gtf"))
Pf <- data.frame(Pf_t = Pf$transcript_id, Pf_g = Pf$gene_id)
Pf <- unique(na.omit(Pf))

Pf_ortho <- data.frame(Orthogroup = ortho_1to1_df$Orthogroup, Pf_t = ortho_1to1_df$Pfalciparum)
Pf_ortho_gene <- merge(as.data.frame(Pf), as.data.frame(Pf_ortho), by = "Pf_t")
ID_change <- merge(ID_change, as.data.frame(Pf_ortho_gene), by = "Orthogroup")

ID_change <- ID_change[, -grep("t$", colnames(ID_change))]

###

ortho_1to1_df$Pco_g <- sapply(ortho_1to1_df$Pcoatneyi, function(x) strsplit(unlist(x), split = "-")[[1]])[1,]
ortho_1to1_df$Pcy_g <- sapply(ortho_1to1_df$Pcynomolgi, function(x) strsplit(unlist(x), split = "-")[[1]])[1,]
ortho_1to1_df$Py_g <- sapply(ortho_1to1_df$Pyoelii, function(x) strsplit(unlist(x), split = "\\.")[[1]][1])
Pco_Pcy_Py <- data.frame(Orthogroup = ortho_1to1_df$Orthogroup, Pco_g = ortho_1to1_df$Pco_g, Pcy_g = ortho_1to1_df$Pcy_g, Py_g = ortho_1to1_df$Py_g)

###

orthogroups <- merge(ID_change, Pco_Pcy_Py, by = "Orthogroup")
orthogroups$Orthogroup <- paste("p", orthogroups$Orthogroup, sep = "_")

write.table(orthogroups, "/SAN/Plasmo_compare/OrthoFinder/parasite_orthogroups.txt", sep = '\t', row.names = F)
save(orthogroups, file = "/SAN/Plasmo_compare/OrthoFinder/parasite_orthogroups.RData")

########### orthologous host genes ########

# library(seqinr)
# 
# human <- read.fasta("/SAN/Plasmo_compare/OrthoFinder/Hosts/human.fasta", seqtype = "AA", as.string = T)
# # df of pep name and gene name
# human_gene_pep <- data.frame(pep = attr(human, "name"), gene = sapply(strsplit(sapply(human, function(x) attr(x, "Annot")), split = " "), function(x) x[4]))
# human_gene_pep$h_g <- sapply(strsplit(sapply(strsplit(sapply(strsplit(sapply(human, function(x) attr(x, "Annot")), split = " "), function(x) x[4]), split = ":"), function(x) x[2]), split = "\\."), function(x) x[1])
# human_gene_pep <- data.frame(h_pep = human_gene_pep$pep, h_g = human_gene_pep$h_g)
# 
# mouse <- read.fasta("/SAN/Plasmo_compare/OrthoFinder/Hosts/mouse.fasta", seqtype = "AA", as.string = T)
# mouse_gene_pep <- data.frame(pep = attr(mouse, "name"), gene = sapply(strsplit(sapply(mouse, function(x) attr(x, "Annot")), split = " "), function(x) x[4]))
# mouse_gene_pep$m_g <- sapply(strsplit(sapply(strsplit(sapply(strsplit(sapply(mouse, function(x) attr(x, "Annot")), split = " "), function(x) x[4]), split = ":"), function(x) x[2]), split = "\\."), function(x) x[1])
# mouse_gene_pep <- data.frame(m_pep = mouse_gene_pep$pep, m_g = mouse_gene_pep$m_g)
# 
# monkey <- read.fasta("/SAN/Plasmo_compare/OrthoFinder/Hosts/monkey.fasta", seqtype = "AA", as.string = T)
# monkey_gene_pep <- data.frame(pep = attr(monkey, "name"), gene = sapply(strsplit(sapply(monkey, function(x) attr(x, "Annot")), split = " "), function(x) x[4]))
# monkey_gene_pep$mo_g <- sapply(strsplit(sapply(strsplit(sapply(strsplit(sapply(monkey, function(x) attr(x, "Annot")), split = " "), function(x) x[4]), split = ":"), function(x) x[2]), split = "\\."), function(x) x[1])
# monkey_gene_pep <- data.frame(mo_pep = monkey_gene_pep$pep, mo_g = monkey_gene_pep$mo_g)

# host_ortho <- read.delim("/SAN/Plasmo_compare/OrthoFinder/Hosts/Results_Apr10/Orthogroups.csv", header = T)
# colnames(host_ortho)[1] <- "Orthogroups"
# 
# h_ortho <- data.frame(Orthogroups = host_ortho$Orthogroups, h_pep = host_ortho$human)
# h_ortho <- merge(h_ortho, human_gene_pep, by = "h_pep")
# 
# m_ortho <- data.frame(Orthogroups = host_ortho$Orthogroups, m_pep = host_ortho$mouse)
# m_ortho <- merge(m_ortho, mouse_gene_pep, by = "m_pep")
# 
# mo_ortho <- data.frame(Orthogroups = host_ortho$Orthogroups, mo_pep = host_ortho$monkey)
# mo_ortho <- merge(mo_ortho, monkey_gene_pep, by = "mo_pep")
# 
# host_ortho <- merge(h_ortho, m_ortho, by = "Orthogroups")
# host_ortho <- merge(host_ortho, mo_ortho, by = "Orthogroups")
# host_ortho$Orthogroups <- paste("h_", host_ortho$Orthogroups, sep = "_")
# host_ortho <- host_ortho[, -grep("pep$", colnames(host_ortho))] # remove pep names
# only around 4000 genes

# write.table(host_ortho, "/SAN/Plasmo_compare/OrthoFinder/host_orthogroups.txt", sep = '\t', row.names = F)
# save(host_ortho, file = "/SAN/Plasmo_compare/OrthoFinder/host_orthogroups.RData")
# 
# ## get 1:1 orthologs
# 
# max.copies <- c()
# for(i in 1:nrow(host_ortho)) #get max.copies
# {
#   human.copies <- length(strsplit(as.character(host_ortho$human[i]), split = ",")[[1]])
#   mouse.copies <- length(strsplit(as.character(host_ortho$mouse[i]), split = ",")[[1]])
#   monkey.copies <- length(strsplit(as.character(host_ortho$monkey[i]), split = ",")[[1]])
#   max.copies[i] <- max(human.copies, mouse.copies, monkey.copies)
# }
# 
# host_ortho$max.copies <- max.copies
# host_ortho <- host_ortho[host_ortho$max.copies==1,]
# 
# h_vs_m <- 
#   read.csv2("../OrthoFinder/Hosts/Results_Apr10/Orthologues_Apr10/Orthologues/Orthologues_human/human__v__mouse.csv", sep = '\t', header = T)
# 
# # first expand all peptides separated by comma to get gene names for each peptide
# # separate human and mouse and then put them back together again
# # use dplyr and tidyr
# h_vs_m %>%
#   select(-mouse) %>% # remove mouse row
#   mutate(h_pep = as.character(str_split(as.character(human), ","))) %>% # split contents of each peptide row
#   unnest(h_pep) %>% # spread out the peptides
#   select(-human) %>% # remove nested human peptide column
#   as.data.frame() %>%
#   inner_join(., human_gene_pep, by = "h_pep") -> h_vs_m_h # get gene names for peptides
# 
# # same for mouse
# 
# h_vs_m %>%
#   select(-human) %>% # remove mouse row
#   mutate(m_pep = as.character(str_split(as.character(mouse), ","))) %>% # split contents of each peptide row
#   unnest(m_pep) %>% # spread out the peptides
#   select(-mouse) %>% # remove nested human peptide column
#   as.data.frame() %>%
#   inner_join(., mouse_gene_pep, by = "m_pep") -> h_vs_m_m # get gene names for peptides
# 
# # merge mouse and human based on orthologous groups
# h_vs_m <- merge(h_vs_m_h, h_vs_m_m, by = "Orthogroup") # why is the number reducing here?
# h_vs_m <- h_vs_m[,-grep("_pep$", colnames(h_vs_m))]
# 
# h_vs_m %>%
#   group_by(Orthogroup) %>%
#   summarise(h_g = unique(toString(as.character(h_g))), m_g = unique(toString(as.character(m_g)))) -> h_vs_mo
# ##
# h_vs_mo <- 
#   read.csv2("../OrthoFinder/Hosts/Results_Apr10/Orthologues_Apr10/Orthologues/Orthologues_human/human__v__monkey.csv", sep = '\t', header = T)
# 
# h_vs_mo %>%
#   select(-monkey) %>% # remove mouse row
#   mutate(h_pep = as.character(str_split(as.character(human), ","))) %>% # split contents of each peptide row
#   unnest(h_pep) %>% # spread out the peptides
#   select(-human) %>% # remove nested human peptide column
#   as.data.frame() %>%
#   inner_join(., human_gene_pep, by = "h_pep") -> h_vs_mo_h # get gene names for peptides
# 
# # same for monkey
# 
# h_vs_mo %>%
#   select(-human) %>% # remove mouse row
#   mutate(mo_pep = as.character(str_split(as.character(monkey), ","))) %>% # split contents of each peptide row
#   unnest(mo_pep) %>% # spread out the peptides
#   select(-monkey) %>% # remove nested human peptide column
#   as.data.frame() %>%
#   inner_join(., monkey_gene_pep, by = "mo_pep") -> h_vs_mo_mo # get gene names for peptides
# 
# # merge monkey and human based on orthologous groups
# h_vs_mo <- merge(h_vs_mo_h, h_vs_mo_mo, by = "Orthogroup") # why is the number reducing here?
# h_vs_mo <- h_vs_mo[,-grep("_pep$", colnames(h_vs_mo))]
# 
# h_vs_mo %>%
#   group_by(Orthogroup) %>%
#   summarise(h_g = unique(toString(as.character(h_g))), mo_g = unique(toString(as.character(mo_g)))) -> h_vs_mo
# 
# 
# # get max copies for mouse
# max.copies <- c()
# for(i in 1:nrow(h_vs_m)) #get max.copies
# {
#   human.copies <- length(unique(strsplit(as.character(h_vs_m$h_g[i]), split = ",")[[1]]))
#   mouse.copies <- length(unique(strsplit(as.character(h_vs_m$m_g[i]), split = ",")[[1]]))
#   #monkey.copies <- length(unique(strsplit(as.character(h_vs_m$mo_g[i]), split = ",")[[1]]))
#   max.copies[i] <- max(human.copies, mouse.copies)
# }
# h_vs_m$max.copies <- max.copies
# h_vs_m <- h_vs_m[h_vs_m$max.copies==1,]
# 
# # get max copies for monkey
# 
# max.copies <- c()
# for(i in 1:nrow(h_vs_mo)) #get max.copies
# {
#   human.copies <- length(unique(strsplit(as.character(h_vs_mo$h_g[i]), split = ",")[[1]]))
#   #mouse.copies <- length(unique(strsplit(as.character(h_vs_mo$m_g[i]), split = ",")[[1]]))
#   monkey.copies <- length(unique(strsplit(as.character(h_vs_mo$mo_g[i]), split = ",")[[1]]))
#   max.copies[i] <- max(human.copies, monkey.copies)
# }
# 
# h_vs_mo$max.copies <- max.copies
# h_vs_mo <- h_vs_mo[h_vs_mo$ma x.copies==1,]
# 
# h_m_mo <- merge(h_vs_m, h_vs_mo, by = "Orthogroup") # really low number -> check laptop dataset. Similar numbers from files.

# checking with orthogroups.genecount
# genecount <- read.csv2("/SAN/Plasmo_compare/OrthoFinder/Hosts/Results_Apr10/Orthogroups.GeneCount.csv", sep = '\t', header = T)[,-5] %>%
#   rename(orthogroup=X) %>%
#   filter(human==1) %>%
#   filter(mouse==1) %>%
#   filter(monkey==1) # gives only 4565 orthogroups


###################################### from ensembl biomart ############################################

# read mart and get only one to one orthologs
mart <- read.delim("../OrthoFinder/mart_export.txt") %>% #dated April 11, 2019
  filter(Mouse.homology.type=="ortholog_one2one") %>%
  #filter(Macaque.homology.type=="ortholog_one2one") %>%
  # get rid of transcript IDs
  select(-Transcript.stable.ID) %>%
  unique.data.frame() %>%
  mutate(Orthogroup = row_number())
mart_ <- stringr::str_pad(mart$Orthogroup, width = 7, pad = "0", side = "left")
mart$Orthogroup <- sub("^", "h_OG", mart_)
mart <-data.frame(Orthogroup = mart$Orthogroup, h_g = mart$Gene.stable.ID, m_g = mart$Mouse.gene.stable.ID, mo_g = mart$Macaque.gene.stable.ID)

write.table(mart, "../OrthoFinder/host_orthogroups.txt", sep = '\t', row.names = F)
save(mart, file = "../OrthoFinder/host_orthogroups.RData")

################## find how many genes are missing from proteome fasta files compared to gtf files ################

  annotation_differences <- list(list())
  parasites <- c("Pfalciparum", "Pberghei", "Pvivax", "Pcoatneyi", "Pcynomolgi", "Pchabaudi", "Pyoelii")
  for(i in parasites) # loop does not work, data is not consistent
  {
    mylist <- list()
    # load fasta files
    
    library(seqinr)
    fasta <- read.fasta(paste0("../OrthoFinder/Parasites/", i, ".fasta", collapse = ""), seqtype = "AA", as.string = T)
    # get all names and gene names from fasta file
    names <- attr(fasta, "name")
    genes <- sapply(sapply(sapply(sapply(fasta, function(x) strsplit(attr(x, "Annot"), split = " | ")), function(x) x[5]), function(x) strsplit(x, split = "=")), function(x) x[2])
    
    require(rtracklayer)
    gtf <- as.data.frame(import(paste0("../Genomes/annotation/",i,".gtf", collapse = ""), format = "gtf"))
    gtf <- gtf[gtf$gene_biotype=="protein_coding",]
    gtf_genes <- unique(gtf$gene_id)
    
    diff_fa_gtf <- setdiff(genes, gtf_genes)
    diff_gtf_fa <- setdiff(gtf_genes, genes)
    
    # if(length(diff_fa_gtf) > 0)
    # {annotation_differences[i,"diff_fa_gtf"] <- I(list(diff_fa_gtf))}else
    #   annotation_differences[i,"diff_fa_gtf"] <- "none"
    # 
    # if(length(diff_gtf_fa) > 0)
    # {annotation_differences[i,"diff_gtf_fa"] <- I(list(diff_gtf_fa))}else
    #   annotation_differences[i,"diff_gtf_fa"] <- "none"
    
    if(length(diff_fa_gtf) > 0)
      {mylist[[1]] <- diff_fa_gtf}else
      mylist[[1]] <- "none"
    
    if(length(diff_gtf_fa) > 0)
    {mylist[[2]] <- diff_gtf_fa}else
      mylist[[2]] <- "none"
    
    annotation_differences[[i]] <- mylist
  }
  
annotation_differences_plasmodb <- annotation_differences
annotation_differences_ensembl <- annotation_differences
  
save(annotation_differences_ensembl, file = "annotation_differences_ensembl.RData")

 