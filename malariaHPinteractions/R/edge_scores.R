## Script to assign scores to all edges from the overall and outside

## first load all studies - blood, liver, blood-liver
studies <- c("blood_all",
             # human studies
             "DRP000987_str", "ERP106451_all", "SRP032775_str", "SRP233153_int", "ERP023982_int", "hpv_bl_int", #hpv_bl is the collection of human - P.vivax blood runs
             # mouse studies
             "ERP110375_str", "ERP004598_int", "m_bl_int", # m_bl is a collection of mouse - Plasmodium studies
             # monkey studies
             "SRP118827_all", "SRP116793_all", "SRP118996_all", "SRP116593_int", "SRP108356_str", "SRP118503_int",
             
             #liver
             "liver_all",
             # human 
             "SRP034011_all", "SRP071199_all", "SRP126641_all", "SRP250329_str",
             # mouse
             "ERP105548_int", "SRP261098_all", "SRP110282_int", 
             "SRP131855_all", "SRP018945_all", "SRP150689_all", "SRP171171_all",
             # monkey
             "SRP096160_str", "ERP020067_all",
             
             "blood_liver_overall"
             )

for(i in studies)
  readRDS(paste0(substr(i, 1, 9), "/cor/", i, "_bipartite.rds", collapse = "")) %>%
  unite("hp", ends_with("gene1"):ends_with("gene2"), remove = F) %>%
  select(hp) %>%
  mutate("{substr(i, 1, 9)}_c" := rep(1, length(.))) -> study_list[[i]]


## join all of the edges together
edge_counter <- list(study_list) %>% 
  reduce(full_join, by = "hp")

# adds up the number of datasets an interaction is found in
edge_rowsums <- rowSums(edge_counter[,2:ncol(edge_counter)], na.rm = T)
edge_counter$edge_rowsums <- edge_rowsums

## make a list/df of which organism, organ a study belongs to; overall doesn't have an organism

organism_list <- list(
  human = c("DRP000987_str", "ERP106451_all", "SRP032775_str", "SRP233153_int", "ERP023982_int", "hpv_bl_int", #hpv_bl is the collection of human - P.vivax blood runs
            "SRP034011_all", "SRP071199_all", "SRP126641_all", "SRP250329_str"),
  mouse = c("ERP110375_str", "ERP004598_int", "m_bl_int",
            "ERP105548_int", "SRP261098_all", "SRP110282_int", 
            "SRP131855_all", "SRP018945_all", "SRP150689_all", "SRP171171_all"),
  monkey = c("SRP118827_all", "SRP116793_all", "SRP118996_all", "SRP116593_int", "SRP108356_str", "SRP118503_int",
             "SRP096160_str", "ERP020067_all")
)

## find the proportion of studies that an edge is found in each host
## the proportion for each host is the coefficient of the that host on that axis

liv <- readRDS("liver_edge_counter_all_studies.rds")

overall_indices <- which(liv$liver_all_c == 1 & liv$edge_rowsums == 1)
overall_interactions <- liv[which(liv$liver_all_c == 1 & liv$edge_rowsums == 1),]

liv <- liv[-overall_indices, -2]

prop_df <- data.frame()
# for each edge, find the props

prop_score <- function(liv)
{
  prop_df <- data.frame()
  organism_list <- list(
    human = c("DRP000987_str", "ERP106451_all", "SRP032775_str", "SRP233153_int", "ERP023982_int", "hpv_bl_int", #hpv_bl is the collection of human - P.vivax blood runs
              "SRP034011_all", "SRP071199_all", "SRP126641_all", "SRP250329_str"),
    mouse = c("ERP110375_str", "ERP004598_int", "m_bl_int",
              "ERP105548_int", "SRP261098_all", "SRP110282_int", 
              "SRP131855_all", "SRP018945_all", "SRP150689_all", "SRP171171_all"),
    monkey = c("SRP118827_all", "SRP116793_all", "SRP118996_all", "SRP116593_int", "SRP108356_str", "SRP118503_int",
               "SRP096160_str", "ERP020067_all")
  )
  
    for(i in 1:nrow(liv))
  {
    print(i)
    interaction <- liv[i,"hp"]
    presence <- liv[i,c(2:(ncol(liv)-1))]
    
    # find which studies the edge is present in
    studies <- colnames(presence[which(!is.na(presence) == T)])
    
    h_score = 0
    mo_score = 0
    m_score = 0
    
    for(j in 1:length(studies))
    {
      if(any(grepl(pattern = substr(studies[j], 1, 9), organism_list[["human"]])))
        h_score <- h_score+1
      if(any(grepl(pattern = substr(studies[j], 1, 9), organism_list[["mouse"]])))
        m_score <- m_score+1
      if(any(grepl(pattern = substr(studies[j], 1, 9), organism_list[["monkey"]])))
        mo_score <- mo_score+1
    }
    
    h_prop = h_score/length(organism_list[["human"]])
    mo_prop = mo_score/length(organism_list[["monkey"]])
    m_prop = m_score/length(organism_list[["mouse"]])
    
    prop_df[i,1] <- interaction
    prop_df[i,2] <- h_prop
    prop_df[i,3] <- mo_prop
    prop_df[i,4] <- m_prop
    }
  return(prop_df)
}

colnames(prop_df) <- c("interaction", "h_prop", "mo_prop", "m_prop")

n = detectCores()  #number of cluster
library(foreach)
library(doParallel)
cl = makeCluster(n)
registerDoParallel(cl)

df <- liv
z = nrow(df)
y = floor(z/n) 
x = nrow(df)%%n

ris = foreach(i = split(df[1:(z-x),],rep(1:n,each=y)), .combine = rbind) %dopar% prop_score(i)

stopCluster(cl)
colnames(ris) <- c("interaction", "h_prop", "mo_prop", "m_prop")

#### blood interactions scores

blood <- readRDS("blood_edge_counter_all_studies.rds")

prop_score <- function(liv)
{
  prop_df <- data.frame()
  organism_list <- list(
    human = c("DRP000987_str", "ERP106451_all", "SRP032775_str", "SRP233153_int", "ERP023982_int", "hpv_bl_int", #hpv_bl is the collection of human - P.vivax blood runs
              "SRP034011_all", "SRP071199_all", "SRP126641_all", "SRP250329_str"),
    mouse = c("ERP110375_str", "ERP004598_int", "mbl_int",
              "ERP105548_int", "SRP261098_all", "SRP110282_int", 
              "SRP131855_all", "SRP018945_all", "SRP150689_all", "SRP171171_all"),
    monkey = c("SRP118827_all", "SRP116793_all", "SRP118996_all", "SRP116593_int", "SRP108356_str", "SRP118503_int",
               "SRP096160_str", "ERP020067_all")
  )
  
  for(i in 1:nrow(liv))
  {
    print(i)
    interaction <- liv[i,"hp"]
    presence <- liv[i,c(2:(ncol(liv)-1))]
    
    # find which studies the edge is present in
    studies <- colnames(presence[which(!is.na(presence) == T)])
    
    h_score = 0
    mo_score = 0
    m_score = 0
    
    for(j in 1:length(studies))
    {
      if(any(grepl(pattern = substr(studies[j], 1, 9), organism_list[["human"]])))
        h_score <- h_score+1
      if(any(grepl(pattern = substr(studies[j], 1, 9), organism_list[["mouse"]])))
        m_score <- m_score+1
      if(any(grepl(pattern = substr(studies[j], 1, 9), organism_list[["monkey"]])))
        mo_score <- mo_score+1
    }
    
    h_prop = h_score/length(organism_list[["human"]])
    mo_prop = mo_score/length(organism_list[["monkey"]])
    m_prop = m_score/length(organism_list[["mouse"]])
    
    prop_df[i,1] <- interaction
    prop_df[i,2] <- h_prop
    prop_df[i,3] <- mo_prop
    prop_df[i,4] <- m_prop
  }
  return(prop_df)
}

n = 20  #number of cluster
library(foreach)
library(doParallel)
cl = makeCluster(n)
registerDoParallel(cl)

df <- blood
z = nrow(df)
y = floor(z/n) 
x = nrow(df)%%n

ris = foreach(i = split(df[1:(z-x),],rep(1:n,each=y)), .combine = rbind) %dopar% prop_score(i)

stopCluster(cl)
colnames(ris) <- c("interaction", "h_prop", "mo_prop", "m_prop")


############ blood-liver #############

blood_liver_score <- function(df)
{
  df = data.frame(interaction = df)
  for(i in 1:nrow(df))
  {
    print(i)
    interaction <- df[i,"interaction"]
    
    if(any(grepl(pattern = interaction, blood$interaction)) & any(grepl(pattern = interaction, liver$interaction)))
      score = (blood[grep(pattern = interaction, blood$interaction), "prop_sum"])*(liver[grep(pattern = interaction, liver$interaction), "prop_sum"])
    if(any(grepl(pattern = interaction, blood$interaction)) & any(!grepl(pattern = interaction, liver$interaction)))
      score = 5 + (blood[grep(pattern = interaction, blood$interaction), "prop_sum"])
    if(any(!grepl(pattern = interaction, blood$interaction)) & any(grepl(pattern = interaction, liver$interaction)))
      score = -5 - (liver[grep(pattern = interaction, liver$interaction), "prop_sum"])
    
    df[i,2] <- score
  }
  return(df)
}

n = 20  #number of clusters
library(foreach)
library(doParallel)
cl = makeCluster(n)
registerDoParallel(cl)

df <- blood_liver_edge_score
z = nrow(df)
y = floor(z/n) 
x = nrow(df)%%n

ris = foreach(i = split(df[1:(z-x),1],rep(1:n,each=y)), .combine = rbind) %dopar% blood_liver_score(i)

stopCluster(cl)
colnames(ris) <- c("interaction", "prop_sum")


############## blood para matrix style ######

organism_list <- list(
  human = c("DRP000987_str", "ERP106451_all", "SRP032775_str", "SRP233153_int", "ERP023982_int", "hpv_bl_int", #hpv_bl is the collection of human - P.vivax blood runs
            "SRP034011_all", "SRP071199_all", "SRP126641_all", "SRP250329_str"),
  mouse = c("ERP110375_str", "ERP004598_int", "mbl_int",
            "ERP105548_int", "SRP261098_all", "SRP110282_int", 
            "SRP131855_all", "SRP018945_all", "SRP150689_all", "SRP171171_all"),
  monkey = c("SRP118827_all", "SRP116793_all", "SRP118996_all", "SRP116593_int", "SRP108356_str", "SRP118503_int",
             "SRP096160_str", "ERP020067_all")
)
# join all possible para para pairs with those from the study
# all possible para pairs
pOG <- read.delim("/SAN/Plasmo_compare/OrthoFinder/parasite_orthogroups.txt")
para_genes <- pOG$Orthogroup
para_pairs <- expand.grid(para_genes, para_genes)
colnames(para_pairs) <- c("gene1", "gene2")

human_summary_mat <- matrix(rep(0, length(para_genes)*length(para_genes)), nrow = length(para_genes))
colnames(human_summary_mat) = para_genes; rownames(human_summary_mat) = para_genes

for(i in 1:length(organism_list[["human"]]))
{
  print(organism_list[["human"]][i])
  if(organism_list[["human"]][i] == "hpv_bl_int")
  {  
    study <- loadRData("/SAN/Plasmo_compare/SRAdb/Output/hpv_bl/cor/hpv_bl_int_para.RData")%>%
      mutate(c = rep(1, nrow(.)))
  } else if(organism_list[["human"]][i] == "mbl_int") {
    study <- loadRData("/SAN/Plasmo_compare/SRAdb/Output/m_bl/cor/m_bl_int_para.RData")%>%
    mutate(c = rep(1, nrow(.)))
  } else {
  study <- loadRData(paste0("/SAN/Plasmo_compare/SRAdb/Output/", substr(organism_list[["human"]][i], 1, 9), "/cor/", organism_list[["human"]][i], "_para.RData", collapse = ''))%>%
    mutate(c = rep(1, nrow(.)))
  }
  colnames(study)[1] <- c("gene1")
  colnames(study)[2] <- c("gene2")
  study <- study[,c("gene1","gene2", "c")]
  
  standard_df <- left_join(para_pairs, study, by = c("gene1", "gene2"))
  # dcast all host genes vs all parasite genes (in this case, only para by para)
  
  study_mat <- standard_df %>% 
    dcast(gene1 ~ gene2, value.var = "c") %>%
    tibble::column_to_rownames("gene1")
  
  study_mat[is.na(study_mat)] <- 0
  # add all human studies matrices, then monkey and then mouse matrices, get proportion for each
  
  human_summary_mat <- human_summary_mat + study_mat
  
}

human_summary_mat_prop <- human_summary_mat/length(organism_list[["human"]])

# monkey

monkey_summary_mat <- matrix(rep(0, length(para_genes)*length(para_genes)), nrow = length(para_genes))
colnames(monkey_summary_mat) = para_genes; rownames(monkey_summary_mat) = para_genes

for(i in 1:length(organism_list[["monkey"]]))
{
  print(organism_list[["monkey"]][i])
  if(organism_list[["monkey"]][i] == "hpv_bl_int")
  {  
    study <- loadRData("/SAN/Plasmo_compare/SRAdb/Output/hpv_bl/cor/hpv_bl_int_para.RData")%>%
      mutate(c = rep(1, nrow(.)))
  } else if(organism_list[["monkey"]][i] == "mbl_int") {
    study <- loadRData("/SAN/Plasmo_compare/SRAdb/Output/m_bl/cor/m_bl_int_para.RData")%>%
      mutate(c = rep(1, nrow(.)))
  } else {
    study <- loadRData(paste0("/SAN/Plasmo_compare/SRAdb/Output/", substr(organism_list[["monkey"]][i], 1, 9), "/cor/", organism_list[["monkey"]][i], "_para.RData", collapse = ''))%>%
      mutate(c = rep(1, nrow(.)))
  }
  colnames(study)[1] <- c("gene1")
  colnames(study)[2] <- c("gene2")
  study <- study[,c("gene1","gene2", "c")]
  
  standard_df <- left_join(para_pairs, study, by = c("gene1", "gene2"))
  # dcast all host genes vs all parasite genes (in this case, only para by para)
  
  study_mat <- standard_df %>% 
    dcast(gene1 ~ gene2, value.var = "c") %>%
    tibble::column_to_rownames("gene1")
  
  study_mat[is.na(study_mat)] <- 0
  # add all monkey studies matrices, then monkey and then mouse matrices, get proportion for each
  
  monkey_summary_mat <- monkey_summary_mat + study_mat
}
monkey_summary_mat_prop <- monkey_summary_mat/length(organism_list[["monkey"]])

# mouse

mouse_summary_mat <- matrix(rep(0, length(para_genes)*length(para_genes)), nrow = length(para_genes))
colnames(mouse_summary_mat) = para_genes; rownames(mouse_summary_mat) = para_genes

for(i in 1:length(organism_list[["mouse"]]))
{
  print(organism_list[["mouse"]][i])
  if(organism_list[["mouse"]][i] == "hpv_bl_int")
  {  
    study <- loadRData("/SAN/Plasmo_compare/SRAdb/Output/hpv_bl/cor/hpv_bl_int_para.RData")%>%
      mutate(c = rep(1, nrow(.)))
  } else if(organism_list[["mouse"]][i] == "mbl_int") {
    study <- loadRData("/SAN/Plasmo_compare/SRAdb/Output/m_bl/cor/m_bl_int_para.RData")%>%
      mutate(c = rep(1, nrow(.)))
  } else {
    study <- loadRData(paste0("/SAN/Plasmo_compare/SRAdb/Output/", substr(organism_list[["mouse"]][i], 1, 9), "/cor/", organism_list[["mouse"]][i], "_para.RData", collapse = ''))%>%
      mutate(c = rep(1, nrow(.)))
  }
  colnames(study)[1] <- c("gene1")
  colnames(study)[2] <- c("gene2")
  study <- study[,c("gene1","gene2", "c")]
  
  standard_df <- left_join(para_pairs, study, by = c("gene1", "gene2"))
  # dcast all host genes vs all parasite genes (in this case, only para by para)
  
  study_mat <- standard_df %>% 
    dcast(gene1 ~ gene2, value.var = "c") %>%
    tibble::column_to_rownames("gene1")
  
  study_mat[is.na(study_mat)] <- 0
  # add all mouse studies matrices, then mouse and then mouse matrices, get proportion for each
  
  mouse_summary_mat <- mouse_summary_mat + study_mat
}
mouse_summary_mat_prop <- mouse_summary_mat/length(organism_list[["mouse"]])

M %>% tibble::rownames_to_column("gene1") -> M
melt(M, variable.name = "gene2",value.name ="prop", id="gene1")


######
library(tidyverse)
library(reshape2)
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# join all possible para para pairs with those from the study
# all possible para pairs
type = "para" # or "para" for essentiality
organ = "liver" # or "liver" or "both"

edge_scoring <- function(type, organ)
{
  hOG <- read.delim("/SAN/Plasmo_compare/OrthoFinder/host_orthogroups.txt")
  host_genes <- hOG$Orthogroup
  pOG <- read.delim("/SAN/Plasmo_compare/OrthoFinder/parasite_orthogroups.txt")
  para_genes <- pOG$Orthogroup
  
  # host_pairs <- expand.grid(host_genes, host_genes)
  # colnames(host_pairs) <- c("gene1", "gene2")
  para_pairs <- expand.grid(para_genes, para_genes)
  colnames(para_pairs) <- c("gene1", "gene2")
  
  host_para <- expand.grid(host_genes, para_genes)
  colnames(host_para) <- c("gene1", "gene2")
  
  both_list <- list(
    human = c("DRP000987_str", "ERP106451_all", "SRP032775_str", "SRP233153_int", "ERP023982_int", "hpv_bl_int",
              "SRP034011_all", "ERP105548_int", "SRP071199_all", "SRP126641_all", "SRP250329_str"),
    mouse = c("ERP110375_str", "ERP004598_int", "mbl_int",
              "SRP261098_all", "SRP110282_int", "SRP131855_all", "SRP018945_all", "SRP150689_all", "SRP171171_all"),
    monkey = c("SRP118827_all", "SRP116793_all", "SRP118996_all", "SRP116593_int", "SRP108356_str", "SRP118503_int",
               "SRP096160_str", "ERP020067_all")
  )
  
  blood_list <- list(
    human = c("DRP000987_str", "ERP106451_all", "SRP032775_str", "SRP233153_int", "ERP023982_int", "hpv_bl_int"),
    mouse = c("ERP110375_str", "ERP004598_int", "mbl_int"),
    monkey = c("SRP118827_all", "SRP116793_all", "SRP118996_all", "SRP116593_int", "SRP108356_str", "SRP118503_int")
  )
  
  liver_list <- list(
    human = c("SRP034011_all", "ERP105548_int", "SRP071199_all", "SRP126641_all", "SRP250329_str"),
    mouse = c("SRP261098_all", "SRP110282_int",
              "SRP131855_all", "SRP018945_all", "SRP150689_all", "SRP171171_all"),
    monkey = c("SRP096160_str", "ERP020067_all")
  )
  
  if(organ == "blood"){
    organism_list = blood_list
  } else if(organ == "liver"){
    organism_list = liver_list
  } else {
    organism_list = both_list
  }
  
  if(type == "bipartite"){
    genes = c(host_genes, para_genes)
    pairs = expand.grid(host_genes, para_genes)
  } else if(type == "para"){
    genes = para_genes
    pairs = expand.grid(para_genes, para_genes)
  }
  colnames(pairs) = c("gene1", "gene2")
  
  human_summary_mat <- matrix(rep(0, length(genes)*length(genes)), nrow = length(genes))
  colnames(human_summary_mat) = genes; rownames(human_summary_mat) = genes
  
  monkey_summary_mat <- matrix(rep(0, length(genes)*length(genes)), nrow = length(genes))
  colnames(monkey_summary_mat) = genes; rownames(monkey_summary_mat) = genes
  
  mouse_summary_mat <- matrix(rep(0, length(genes)*length(genes)), nrow = length(genes))
  colnames(mouse_summary_mat) = genes; rownames(mouse_summary_mat) = genes
  
  for(j in c("human", "monkey", "mouse"))
  {
    for(i in 1:length(organism_list[[j]]))
    {
      print(organism_list[[j]][i])
      if(organism_list[[j]][i] == "hpv_bl_int")
      {  
        study <- loadRData("/SAN/Plasmo_compare/SRAdb/Output/hpv_bl/cor/hpv_bl_int_para.RData")%>%
          mutate(c = rep(1, nrow(.)))
      } else if(organism_list[[j]][i] == "mbl_int") {
        study <- loadRData("/SAN/Plasmo_compare/SRAdb/Output/m_bl/cor/m_bl_int_para.RData")%>%
          mutate(c = rep(1, nrow(.)))
      } else {
        study <- loadRData(paste0("/SAN/Plasmo_compare/SRAdb/Output/", substr(organism_list[[j]][i], 1, 9), "/cor/", organism_list[[j]][i], "_", type, ".RData", collapse = ''))%>%
          mutate(c = rep(1, nrow(.)))
      }
      colnames(study)[1] <- c("gene1")
      colnames(study)[2] <- c("gene2")
      study <- study[,c("gene1","gene2", "c")]
      
      standard_df <- left_join(pairs, study, by = c("gene1", "gene2"))
      # dcast all host genes vs all parasite genes (in this case, only para by para)
      
      study_mat <- standard_df %>% 
        dcast(gene1 ~ gene2, value.var = "c") %>%
        tibble::column_to_rownames("gene1")
      
      study_mat[is.na(study_mat)] <- 0
      # add all mouse studies matrices, then mouse and then mouse matrices, get proportion for each
      
      if(j == "human")
        human_summary_mat <- human_summary_mat + study_mat
      if(j == "mouse")
        mouse_summary_mat <- mouse_summary_mat + study_mat
      if(j == "monkey")
        monkey_summary_mat <- monkey_summary_mat + study_mat
    }
  
  human_summary_mat_prop <- human_summary_mat/length(organism_list[["human"]])
  mouse_summary_mat_prop <- mouse_summary_mat/length(organism_list[["mouse"]])
  monkey_summary_mat_prop <- monkey_summary_mat/length(organism_list[["monkey"]])
  
  human_summary_mat_prop %>% tibble::rownames_to_column("gene1") -> human_summary_mat_prop
  human_melt <- melt(human_summary_mat_prop, variable.name = "gene2",value.name ="h_prop", id="gene1"); human_melt$gene2 = as.character(human_melt$gene2)
  mouse_summary_mat_prop %>% tibble::rownames_to_column("gene1") -> mouse_summary_mat_prop
  mouse_melt <- melt(mouse_summary_mat_prop, variable.name = "gene2",value.name ="m_prop", id="gene1"); mouse_melt$gene2 = as.character(mouse_melt$gene2)
  monkey_summary_mat_prop %>% tibble::rownames_to_column("gene1") -> monkey_summary_mat_prop
  monkey_melt <- melt(monkey_summary_mat_prop, variable.name = "gene2",value.name ="mo_prop", id="gene1"); monkey_melt$gene2 = as.character(monkey_melt$gene2)
  
  all_org_df <- plyr::join_all(list(human_melt, mouse_melt, monkey_melt), type = "full", by = c("gene1", "gene2"))
  }
  
return(all_org_df)
}

#exp = edge_scoring("para", "blood")
exp = all_org_df[which(rowSums(all_org_df[,c(3,4,5)]) > 0),]
exp = distinct(exp)

sum_prop = rowSums(exp[,c(3:5)])
exp$sum_prop = sum_prop
saveRDS(exp, file = "blood_liver_para_prop.rds")

mul_prop = apply(exp[,3:5], 1, FUN = function(x) (1-x[1])*(1-x[2])*(1-x[3]))
exp$mul_prop = mul_prop
saveRDS(exp, file = "blood_liver_para_prop.rds")

### plot the drop in scores ###

# y axis: proportion; x-axis: sorted interactions, based on score
df <- readRDS("blood_liver_para_prop.rds")
x = df[order(df$sum_prop),]
png("blood_liver_para_sumprop_plot.png")
plot(x = c(1:nrow(x)), y = x$sum_prop,
     xlim = c(1,nrow(x)),
     ylim = c(min(df$sum_prop), 
              max(df$sum_prop)),
     xaxt = 'n')
dev.off()

x = df[order(df$mul_prop),]
png("blood_liver_para_mulprop_plot.png")
plot(x = c(1:nrow(x)), y = x$mul_prop,
     xlim = c(1,nrow(x)),
     ylim = c(min(df$mul_prop), 
              max(df$mul_prop)),
     xaxt = 'n')
dev.off()


### MDS ###

library(tidyverse)
library(reshape2)
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# join all possible para para pairs with those from the study
# all possible para pairs
type = "bipartite" # or "para" for essentiality
organ = "both" # or "liver" or "both"

edge_scoring <- function(type, organ)
{
  hOG <- read.delim("/SAN/Plasmo_compare/OrthoFinder/host_orthogroups.txt")
  host_genes <- hOG$Orthogroup
  pOG <- read.delim("/SAN/Plasmo_compare/OrthoFinder/parasite_orthogroups.txt")
  para_genes <- pOG$Orthogroup
  
  # host_pairs <- expand.grid(host_genes, host_genes)
  # colnames(host_pairs) <- c("gene1", "gene2")
  para_pairs <- expand.grid(para_genes, para_genes)
  colnames(para_pairs) <- c("gene1", "gene2")
  
  host_para <- expand.grid(host_genes, para_genes)
  colnames(host_para) <- c("gene1", "gene2")
  
  both_list <- list(
    human = c("DRP000987_str", "ERP106451_all", "SRP032775_str", "SRP233153_int", "ERP023982_int", "hpv_bl_int",
              "SRP034011_all", "ERP105548_int", "SRP071199_all", "SRP126641_all", "SRP250329_str"),
    mouse = c("ERP110375_str", "ERP004598_int", "mbl_int",
              "SRP261098_all", "SRP110282_int", "SRP131855_all", "SRP018945_all", "SRP150689_all", "SRP171171_all"),
    monkey = c("SRP118827_all", "SRP116793_all", "SRP118996_all", "SRP116593_int", "SRP108356_str", "SRP118503_int",
               "SRP096160_str", "ERP020067_all")
  )
  
  blood_list <- list(
    human = c("DRP000987_str", "ERP106451_all", "SRP032775_str", "SRP233153_int", "ERP023982_int", "hpv_bl_int"),
    mouse = c("ERP110375_str", "ERP004598_int", "mbl_int"),
    monkey = c("SRP118827_all", "SRP116793_all", "SRP118996_all", "SRP116593_int", "SRP108356_str", "SRP118503_int")
  )
  
  liver_list <- list(
    human = c("SRP034011_all", "ERP105548_int", "SRP071199_all", "SRP126641_all", "SRP250329_str", "SRP110282_int"),
    mouse = c("SRP261098_all", "SRP131855_all", "SRP018945_all", "SRP150689_all", "SRP171171_all"),
    monkey = c("SRP096160_str", "ERP020067_all")
  )
  
  if(organ == "blood"){
    organism_list = blood_list
  } else if(organ == "liver"){
    organism_list = liver_list
  } else {
    organism_list = both_list
  }
  
  if(type == "bipartite"){
    genes = c(host_genes, para_genes)
    pairs = expand.grid(host_genes, para_genes)
  } else if(type == "para"){
    genes = para_genes
    pairs = expand.grid(para_genes, para_genes)
  }
  colnames(pairs) = c("gene1", "gene2")
  
  human_summary_mat <- list()
  
  monkey_summary_mat <- list()
  
  mouse_summary_mat <- list()
  
  for(j in c("human", "monkey", "mouse"))
  {
    for(i in 1:length(organism_list[[j]]))
    {
      print(organism_list[[j]][i])
      if(organism_list[[j]][i] == "hpv_bl_int")
      {  
        study <- loadRData(paste0("/SAN/Plasmo_compare/SRAdb/Output/hpv_bl/cor/hpv_bl_int_", type,".RData", collapse = ''))%>%
          mutate(hpv =  rep(1, nrow(.)))
      } else if(organism_list[[j]][i] == "mbl_int") {
        study <- loadRData(paste0("/SAN/Plasmo_compare/SRAdb/Output/m_bl/cor/m_bl_int_", type, ".RData", collapse = ''))%>%
          mutate(mbl = rep(1, nrow(.)))
      } else {
        study <- loadRData(paste0("/SAN/Plasmo_compare/SRAdb/Output/", substr(organism_list[[j]][i], 1, 9), "/cor/", organism_list[[j]][i], "_", type, ".RData", collapse = ''))%>%
          mutate(!! substr(organism_list[[j]][i], 1, 9) := !! rep(1, nrow(.)))
      }
      colnames(study)[1] <- c("gene1")
      colnames(study)[2] <- c("gene2")
      study <- study[,c("gene1","gene2", colnames(study[ncol(study)]))]
      
      standard_df <- left_join(pairs, study, by = c("gene1", "gene2"))
      # dcast all host genes vs all parasite genes (in this case, only para by para)
      
      #study_mat <- standard_df %>% 
      #  dcast(gene1 ~ gene2, value.var = "c") %>%
      #  tibble::column_to_rownames("gene1")
      
      standard_df[is.na(standard_df[,3]),3] <- 0
      # add all mouse studies matrices, then mouse and then mouse matrices, get proportion for each
      
      if(j == "human")
        human_summary_mat[[i]] <- standard_df
      if(j == "mouse")
        mouse_summary_mat[[i]] <- standard_df
      if(j == "monkey")
        monkey_summary_mat[[i]] <- standard_df
    }
    human_df <- plyr::join_all(human_summary_mat, type = "full", by = c("gene1", "gene2"))
    mouse_df <- plyr::join_all(mouse_summary_mat, type = "full", by = c("gene1", "gene2"))
    monkey_df <- plyr::join_all(monkey_summary_mat, type = "full", by = c("gene1", "gene2"))
    
    all_org_df <- plyr::join_all(list(human_df, monkey_df, mouse_df), type = "full", by = c("gene1", "gene2"))
    all_org_df <- all_org_df[which(rowSums(all_org_df[,c(3:ncol(all_org_df))]) > 0),]
}
  
  return(all_org_df)
}
rownames(all_org_df) <- NULL
df <- all_org_df %>%
  unite("hp", ends_with("gene1"):ends_with("gene2"), remove = T) %>%
  tibble::column_to_rownames("hp")

mds <- t.df %>%
  dist() %>%          
  cmdscale() %>%
  as_tibble()
colnames(mds) <- c("Dim.1", "Dim.2")

ggscatter(mds, x = "Dim.1", y = "Dim.2", 
          size = 1,
          repel = TRUE)

library(factoextra)
t.df <- as.data.frame(t(df))
res.pca <- prcomp(t.df[,1:(ncol(t.df)-1)], scale = F, center = F)
t.df$org <- c("human", "human", "human", "human", "human", "monkey", 
              "monkey", "mouse", "human", "mouse", "mouse", "mouse", 
              "mouse")
group = as.factor(t.df$org)
fviz_eig(res.pca)
fviz_pca_ind(res.pca_nozero,
             col.ind = t.df_nozero$tissue, # Color by the quality of representation
             palette = c("#0073C2FF", "#EFC000FF", "#868686FF"),          
             repel = TRUE,     # Avoid text overlapping
             geom.ind = "point", pointshape = 21,
             pointsize = 2,
             legend.title = "Organism") +
  ggtitle("PCA for blood and liver studies") +
  theme(plot.title = element_text(hjust = 0.5))

xd = data.frame(PC1 = res.pca$x[,"PC1"], PC2 = res.pca$x[,"PC2"], org = t.df$org, tissue = t.df$tissue)
ggplot(xd, aes(x = PC1, y = PC2, fill = org, color = org, shape = tissue)) + 
  geom_point(alpha = 0.4, size = 3) + 
  theme(text = element_text(size=10),axis.text = element_text(size = 10)) + 
  theme_bw()
          

fviz_pca_ind(res.pca.t,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

library(ggfortify)
autoplot(res.pca.t, data = t.df, colour = as.factor(t.df$tissue), label = T, label.size = 2)

t.df_reduced$org = c("human", "human", "human", "human", "human", "human", 
             "human", "human", "human", "human", "human", 
               "monkey", "monkey", "monkey", "monkey", 
             "monkey", "mouse", "mouse", "mouse", "mouse", 
             "mouse", "mouse", "mouse", "mouse")
t.df_reduced$tissue = c("blood", "blood", "blood", "blood", "blood", "blood",
                "liver", "liver", "liver", "liver", "liver", "blood",
                  "blood", "blood", "blood", "liver", "blood", "blood",
               "blood", "liver", "liver", "liver", "liver", "liver")

t.df_nozero <- t.df_reduced[,which(colSums(t.df_reduced[,1:(ncol(t.df_reduced)-2)]) > 0)]


############### logisticPCA ############

logsvd_model = logisticSVD(t.df_nozero[,1:(ncol(t.df_nozero)-2)], k = 2)

cols = t.df_nozero$tissue
plot(logsvd_model, type = "scores") + geom_point(aes(colour = cols)) + 
  ggtitle("Exponential Family PCA") + scale_colour_manual(values = c("blue", "red"))

## hierarchical clustering ######

distances = dist(t.df_nozero[,1:(ncol(t.df_nozero)-2)], method = "euclidean")
clusters = hclust(distances, method = "ward.D2")
dend = as.dendrogram(clusters)
dendata = dendro_data(dend)
dendata$labels$tissue = c("blood", "liver", "liver", "blood",
                            "liver", "blood", "liver", "blood",
                            "blood", "blood", "blood", "blood",
                            "blood", "liver", "liver", "liver",
                            "liver", "blood", "liver", "liver",
                            "blood", "blood", "liver", "blood")
dendata$labels$org = c("mouse", "human", "mouse", "human",
                      "mouse", "monkey", "mouse", "monkey",
                      "human", "human", "mouse", "human",
                      "human", "monkey", "mouse", "human",
                      "human", "human", "human", "mouse",
                       "monkey", "mouse", "human", "monkey")
labs <- label(dendata)
p <- ggplot(segment(dendata)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +coord_flip()+scale_y_reverse()
p <- p + geom_text(data=labs,
                   aes(label=label, x=x, y=-0.5, colour=tissue)) +
  theme(text = element_text(size=2),axis.text = element_text(size = 10)) + 
  theme_dendro()
  scale_colour_manual(values=c("blue", "red", "darkgreen"))

   plot_ggdendro(dendata, direction   = "tb",
  scale.color = cols,
  label.size  = 2.5,
  branch.size = 0.5,
  expand.y    = 0.2) + expand_limits(x = c(-1, 32)) +
    geom_text(data=labs,
              aes(label=label, x=x, y=-0.5, colour=tissue)) +
    theme(text = element_text(size=2),axis.text = element_text(size = 10)) + 
    theme_bw()
   
   
