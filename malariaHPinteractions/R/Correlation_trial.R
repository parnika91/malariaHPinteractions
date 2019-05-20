# Script to test correlation time requirement and GPU time requirement

# Script to test correlation with existing code

# make the dataset

# Step 1: For each study, get 1:1 orthologous genes and get their orthogroup IDs
p_o <- read.delim("/SAN/Plasmo_compare/OrthoFinder/parasite_orthogroups.txt")
h_o <- read.delim("/SAN/Plasmo_compare/OrthoFinder/host_orthogroups.txt")
allHPexp <- read.delim("/SAN/Plasmo_compare/SRAdb/allHPexp.txt", sep = ',')
hp.species <- c(unique(allHPexp$Host), unique(allHPexp$Parasite))
studies <- as.character(unique(allHPexp$Study))

# use disp files for this
for(i in studies)
{
  # get <study>_disp.txt file
  data <- read.delim(paste0("/SAN/Plasmo_compare/SRAdb/Output/", i, "/", i, "_disp.txt", collapse = ""))
  data <- data[,-which(names(data)=="study_disp")] %>%
    tibble::rownames_to_column("gene")
  
  # get the host and parasite
  host = unique(as.character(allHPexp[allHPexp$Study==i,"Host"]))
  para = unique(as.character(allHPexp[allHPexp$Study==i,"Parasite"]))
  
  study = as.data.frame(data)
  study[str_detect(study[,"gene"], "^P+"),(ncol(study)+1)] <- "parasite"
  study[str_detect(study[,"gene"], "^E+"),ncol(study)] <- "host"
  study[str_detect(study[,1], "^b+"),ncol(study)] <- "parasite"
  study[str_detect(study[,1], "^t+"),ncol(study)] <- "parasite"
  colnames(study)[ncol(study)] <- "Species"
  host_genes <- study[study$Species=="host",]
  para_genes <- study[study$Species=="parasite",]
  
  if(host == "human")
    check_h <- "h_g"
  if(host == "mouse")
    check_h <- "m_g"
  if(host == "monkey")
    check_h <- "mo_g"
  
  if(para == "Pfalciparum")
    check_p <- "Pf_g"
  if(para == "Pberghei")
    check_p <- "Pb_g"
  if(para == "Pvivax")
    check_p <- "Pv_g"
  if(para == "Pyoelii")
    check_p <- "Py_g"
  if(para == "Pchabaudi")
    check_p <- "Pch_g"
  if(para == "Pcoatneyi")
    check_p <- "Pco_g"
  if(para == "Pcynomolgi")
    check_p <- "Pcy_g"
  
  h_o1 <- h_o[,c("Orthogroup", check_h)]
  colnames(h_o1)[2] <- "gene"
  p_o1 <- p_o[,c("Orthogroup", check_p)]
  colnames(p_o1)[2] <- "gene"
  
  
  ortho_host <- right_join(host_genes, h_o1)
  ortho_para <- right_join(para_genes, p_o1)
  
  ortho <- rbind(ortho_host, ortho_para)
  write.table(ortho, paste0("/SAN/Plasmo_compare/SRAdb/Output/", i, "/", i, "_ortho.txt", collapse = ''), row.names = F, sep = '\t')
}


# Step 2: concat all runs from all studies based on orthogroup IDs
ortho_data <- data.frame(Orthogroup = c(as.character(h_o$Orthogroup), as.character(p_o$Orthogroup)))
for(j in studies)
{
  print(j)
  
  data <- read.delim(paste0("/SAN/Plasmo_compare/SRAdb/Output/", j, "/", j, "_ortho.txt", collapse = ""))
  data <- data[,-which(names(data)=="Species")]
  data <- data[,-which(names(data)=="gene")]
  
  #allruns sudo apt get update
  ortho_data <- merge(ortho_data, data, by = "Orthogroup")
}
write.table(ortho_data, "ortho_data.txt", sep = '\t', row.names = F)

# BCV all
ortho_data <- ortho_data %>%
  tibble::rownames_to_column("Orthogroup")
row.names(ortho_data) <- NULL
host_ortho_data <- ortho_data[grep(ortho_data$Orthogroup, pattern = "h"),] %>%
  tibble::column_to_rownames("Orthogroup")

para_ortho_data <- ortho_data[grep(ortho_data$Orthogroup, pattern = "p"),] %>%
  tibble::column_to_rownames("Orthogroup")

host_disp <- estimateDisp(host_ortho_data, prior.n = 10)$tagwise.dispersion
para_disp <- estimateDisp(para_ortho_data, prior.n = 10)$tagwise.dispersion


all_w_disp <- cbind(ortho_data, all_disp)


# correlation trial

library(doParallel)
library(doMC)
library(foreach)
library(data.table)
library(WGCNA)
library(qgraph)

study <- read.csv2("ortho_data.txt", sep = "\t", header = T) %>%
  tibble::column_to_rownames("Orthogroup")
system.time(ori_cor <- WGCNA::cor(t(study), use = 'pairwise.complete.obs'))

  reps <- 2
  
  registerDoParallel(cores = 5)
  
  PermAsso <- function()
  {
    res <- list()
    system.time(res <- foreach(i=1:reps, .combine = '+') %dopar%
                  {
                    rand.col <- study[,sample(ncol(study), (ncol(study)))] %>%
                      setNames(., colnames(study))
                    
                    (abs(WGCNA::cor(t(study), t(rand.col), use = 'pairwise.complete.obs') >= abs(ori_cor) +0))
                  })
    system.time(red <- Reduce('+', res))
    return(res)
  }
  
  PermAssoCompiled <- compiler::cmpfun(PermAsso)
  
  outer_list <- list()
  outer_reps <- 1
  
  ptm <- proc.time() 
  for(i in 1:outer_reps)
  {
    outer_list[[i]] <- PermAssoCompiled()
  }
  proc.time() - ptm
  stopImplicitCluster()
  
  ptm1 <- proc.time()
  outer <- Reduce('+', outer_list)
  proc.time() - ptm1

pval <- outer/(reps*outer_reps)
