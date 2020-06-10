setwd("/home/parnika/Documents/")
library(dplyr)
allHPexp <- read.delim("Data/allHPexp.txt", sep = ',')
colnames(allHPexp)
ortho_data <- read.delim("Data/ortho_data.txt")# %>%
#tibble::column_to_rownames("Orthogroup")

#SRP108356 <- allHPexp[which(allHPexp$Study=="SRP108356"),]
#SRP108632 <- allHPexp[which(allHPexp$Study=="SRP108632"),]
ERP106451 <- allHPexp[which(allHPexp$Study=="ERP106451"),]
#SRP032775 <- allHPexp[which(allHPexp$Study=="SRP032775"),]
SRP118996 <- allHPexp[which(allHPexp$Study=="SRP118996"),]
DRP000987 <- allHPexp[which(allHPexp$Study=="DRP000987"),]
#SRP118503 <- allHPexp[which(allHPexp$Study=="SRP118503"),] # 6 in intermediate and 3 in stringent
SRP118827 <- allHPexp[which(allHPexp$Study=="SRP118827"),]
SRP116793 <- allHPexp[which(allHPexp$Study=="SRP116793"),]
SRP116593 <- allHPexp[which(allHPexp$Study=="SRP116593"),] # 13 rins in intermed and 8 runs in stringent
# DRP001953 <- allHPexp[which(allHPexp$Study=="DRP001953"),]
ERP023982 <- allHPexp[which(allHPexp$Study=="ERP023982"),]
# SRP108632 <- allHPexp[which(allHPexp$Study=="SRP108632"),]
ERP004598 <- allHPexp[which(allHPexp$Study=="ERP004598"),]
ERP110375 <- allHPexp[which(allHPexp$Study=="ERP110375"),]
ERP002273 <- allHPexp[which(allHPexp$Study=="ERP002273"),]
SRP096160 <- allHPexp[which(allHPexp$Study=="SRP096160"),]
SRP110282 <- allHPexp[which(allHPexp$Study=="SRP110282"),]


liver <- allHPexp[which(allHPexp$Tissue=="liver"),]
table(liver$Host)
table(liver$HostParasite)
table(liver$Study)
liver.screen <- liver[which(liver$ProteinCodHost >= 1e6 & liver$ProteinCodPara >= 1e5),]
liver.screen <- liver.screen[which(liver.screen$NumberProtCodGenesHost >= 10000 &
                                     liver.screen$NumberProtCodGenesPara > 3000),]
liver.screen <- liver.screen[which(liver.screen$MapPercent >= 70),]
# SRP108356 #

# # Reads for protein-coding genes: Host - >= 10^6, Parasite - >= 10^5
# SRP108356.screen <- SRP108356[which(SRP108356$ProteinCodHost >= 1e6 & SRP108356$ProteinCodPara >= 1e5),]
# # Protein-coding genes: Host - >= 10^4, Parasite - >= 3000
# SRP108356.screen <- SRP108356.screen[which(SRP108356.screen$NumberProtCodGenesHost >= 10000 & SRP108356.screen$NumberProtCodGenesPara > 3000),]
# # Unique map percent: >= 70%
# SRP108356.screen <- SRP108356.screen[which(SRP108356.screen$MapPercent >= 70),]
# # I have 12 samples from SRP108356
# SRP108356.screen.runs <- as.character(SRP108356.screen[,"RunID"])

# SRP032775 #

# # Reads for protein-coding genes: Host - >= 10^6, Parasite - >= 10^5
# SRP032775.screen <- SRP032775[which(SRP032775$ProteinCodHost >= 1e6 & SRP032775$ProteinCodPara >= 1e5),]
# # Protein-coding genes: Host - >= 10^4, Parasite - >= 3000
# SRP032775.screen <- SRP032775.screen[which(SRP032775.screen$NumberProtCodGenesHost >= 10000 & SRP032775.screen$NumberProtCodGenesPara > 3000),]
# # Unique map percent: >= 70%
# SRP032775.screen <- SRP032775.screen[which(SRP032775.screen$MapPercent >= 70),]
# # parasite proportion > 5%
# SRP032775.screen <- SRP032775.screen[which(SRP032775.screen$Parasite_percent >= 5),]
# # I have 35 samples from SRP032775
# SRP032775.screen.runs <- as.character(SRP032775.screen[,"RunID"])

# ERP106451 #

# Reads for protein-coding genes: Host - >= 10^6, Parasite - >= 10^5
ERP106451.screen <- ERP106451[which(ERP106451$ProteinCodHost >= 1e7 & ERP106451$ProteinCodPara >= 1e6),]
# Protein-coding genes: Host - >= 10^4, Parasite - >= 3000
ERP106451.screen <- ERP106451.screen[which(ERP106451.screen$NumberProtCodGenesHost >= 10000 & ERP106451.screen$NumberProtCodGenesPara > 3000),]
# Unique map percent: >= 70%
ERP106451.screen <- ERP106451.screen[which(ERP106451.screen$MapPercent >= 70),]
# # parasite proportion > 5%
# ERP106451.screen <- ERP106451.screen[which(ERP106451.screen$Parasite_percent >= 5),]
# I have 48 samples from ERP106451
ERP106451.screen.runs <- as.character(ERP106451.screen[,"RunID"])

# SRP118996 #

# Reads for protein-coding genes: Host - >= 10^6, Parasite - >= 10^5
SRP118996.screen <- SRP118996[which(SRP118996$ProteinCodHost >= 1e7 & SRP118996$ProteinCodPara >= 1e6),]
# Protein-coding genes: Host - >= 10^4, Parasite - >= 3000
SRP118996.screen <- SRP118996.screen[which(SRP118996.screen$NumberProtCodGenesHost >= 10000 & SRP118996.screen$NumberProtCodGenesPara > 3000),]
# Unique map percent: >= 70%
SRP118996.screen <- SRP118996.screen[which(SRP118996.screen$MapPercent >= 70),]
# # host parasite proportion: keep parasite percent == 5
# SRP118996.screen <- SRP118996.screen[which(SRP118996.screen$Parasite_percent >= 5),]
# I have 47 samples from SRP118996 
SRP118996.screen.runs <- as.character(SRP118996.screen[,"RunID"])

# DRP000987 #

# Int: Reads for protein-coding genes: Host - >= 10^6, Parasite - >= 10^5
DRP000987.screen <- DRP000987[which(DRP000987$ProteinCodHost >= 1e7 & DRP000987$ProteinCodPara >= 1e6),]
# # Protein-coding genes: Host - >= 10^4, Parasite - >= 3000
DRP000987.screen <- DRP000987.screen[which(DRP000987.screen$NumberProtCodGenesHost >= 10000 & DRP000987.screen$NumberProtCodGenesPara > 3000),]
# # Unique map percent: >= 70%
DRP000987.screen <- DRP000987.screen[which(DRP000987.screen$MapPercent >= 70),]
# # # parasite proportion > 5%
# DRP000987.screen <- DRP000987.screen[which(DRP000987.screen$Parasite_percent >= 5),]
# # I have 48 samples from DRP000987
DRP000987.screen.runs <- as.character(DRP000987.screen[,"RunID"])

# DRP001953 #

# Int: Reads for protein-coding genes: Host - >= 10^6, Parasite - >= 10^5
DRP001953.screen <- DRP001953[which(DRP001953$ProteinCodHost >= 1e6 & DRP001953$ProteinCodPara >= 1e5),]
# # Protein-coding genes: Host - >= 10^4, Parasite - >= 3000
DRP001953.screen <- DRP001953.screen[which(DRP001953.screen$NumberProtCodGenesHost >= 10000 & DRP001953.screen$NumberProtCodGenesPara > 3000),]
# # Unique map percent: >= 70%
DRP001953.screen <- DRP001953.screen[which(DRP001953.screen$MapPercent >= 70),]
# # # parasite proportion > 5%
# DRP001953.screen <- DRP001953.screen[which(DRP001953.screen$Parasite_percent >= 5),]
# # I have 48 samples from DRP001953
DRP001953.screen.runs <- as.character(DRP001953.screen[,"RunID"])

# ERP023982 #

# Int: Reads for protein-coding genes: Host - >= 10^6, Parasite - >= 10^5
ERP023982.screen <- ERP023982[which(ERP023982$ProteinCodHost >= 1e7 & ERP023982$ProteinCodPara >= 1e6),]
# # Protein-coding genes: Host - >= 10^4, Parasite - >= 3000
ERP023982.screen <- ERP023982.screen[which(ERP023982.screen$NumberProtCodGenesHost >= 10000 & ERP023982.screen$NumberProtCodGenesPara > 3000),]
# # Unique map percent: >= 70%
ERP023982.screen <- ERP023982.screen[which(ERP023982.screen$MapPercent >= 70),]
# # # parasite proportion > 5%
# ERP023982.screen <- ERP023982.screen[which(ERP023982.screen$Parasite_percent >= 5),]
# # I have 48 samples from ERP023982
ERP023982.screen.runs <- as.character(ERP023982.screen[,"RunID"])

# SRP118503 #

# Reads for protein-coding genes: Host - >= 10^6, Parasite - >= 10^5
# SRP118503.screen <- SRP118503[which(SRP118503$ProteinCodHost >= 1e7 & SRP118503$ProteinCodPara >= 1e6),]
# # Protein-coding genes: Host - >= 10^4, Parasite - >= 3000
# SRP118503.screen <- SRP118503.screen[which(SRP118503.screen$NumberProtCodGenesHost >= 10000 & SRP118503.screen$NumberProtCodGenesPara > 3000),]
# # Unique map percent: >= 70%
# SRP118503.screen <- SRP118503.screen[which(SRP118503.screen$MapPercent >= 70),]
# # # parasite proportion > 5%
# # SRP118503.screen <- SRP118503.screen[which(SRP118503.screen$Parasite_percent >= 5),]
# # I have 48 samples from SRP118503
# SRP118503.screen.runs <- as.character(SRP118503.screen[,"RunID"])

# SRP118827 #

# Reads for protein-coding genes: Host - >= 10^6, Parasite - >= 10^5
SRP118827.screen <- SRP118827[which(SRP118827$ProteinCodHost >= 1e6 & SRP118827$ProteinCodPara >= 1e5),]
# Protein-coding genes: Host - >= 10^4, Parasite - >= 3000
SRP118827.screen <- SRP118827.screen[which(SRP118827.screen$NumberProtCodGenesHost >= 10000 & SRP118827.screen$NumberProtCodGenesPara > 3000),]
# Unique map percent: >= 70%
SRP118827.screen <- SRP118827.screen[which(SRP118827.screen$MapPercent >= 70),]
# # parasite proportion > 5%
# SRP118827.screen <- SRP118827.screen[which(SRP118827.screen$Parasite_percent >= 5),]
# I have 48 samples from SRP118827
SRP118827.screen.runs <- as.character(SRP118827.screen[,"RunID"])

# SRP116793 #

# Reads for protein-coding genes: Host - >= 10^6, Parasite - >= 10^5
SRP116793.screen <- SRP116793[which(SRP116793$ProteinCodHost >= 1e6 & SRP116793$ProteinCodPara >= 1e5),]
# Protein-coding genes: Host - >= 10^4, Parasite - >= 3000
SRP116793.screen <- SRP116793.screen[which(SRP116793.screen$NumberProtCodGenesHost >= 10000 & SRP116793.screen$NumberProtCodGenesPara > 3000),]
# Unique map percent: >= 70%
SRP116793.screen <- SRP116793.screen[which(SRP116793.screen$MapPercent >= 70),]
# # parasite proportion > 5%
# SRP116793.screen <- SRP116793.screen[which(SRP116793.screen$Parasite_percent >= 5),]
# I have 48 samples from SRP116793
SRP116793.screen.runs <- as.character(SRP116793.screen[,"RunID"])

# SRP116593 #

# Reads for protein-coding genes: Host - >= 10^6, Parasite - >= 10^5
SRP116593.screen <- SRP116593[which(SRP116593$ProteinCodHost >= 1e6 & SRP116593$ProteinCodPara >= 1e5),]
# Protein-coding genes: Host - >= 10^4, Parasite - >= 3000
SRP116593.screen <- SRP116593.screen[which(SRP116593.screen$NumberProtCodGenesHost >= 10000 & SRP116593.screen$NumberProtCodGenesPara > 3000),]
# Unique map percent: >= 70%
SRP116593.screen <- SRP116593.screen[which(SRP116593.screen$MapPercent >= 70),]
# # parasite proportion > 5%
# SRP116593.screen <- SRP116593.screen[which(SRP116593.screen$Parasite_percent >= 5),]
# I have 48 samples from SRP116593
SRP116593.screen.runs <- as.character(SRP116593.screen[,"RunID"])

# ERP004598 #

# Reads for protein-coding genes: Host - >= 10^6, Parasite - >= 10^5
ERP004598.screen <- ERP004598[which(ERP004598$ProteinCodHost >= 1e7 & ERP004598$ProteinCodPara >= 1e6),]
# Protein-coding genes: Host - >= 10^4, Parasite - >= 3000
ERP004598.screen <- ERP004598.screen[which(ERP004598.screen$NumberProtCodGenesHost >= 10000 & ERP004598.screen$NumberProtCodGenesPara > 3000),]
# Unique map percent: >= 70%
ERP004598.screen <- ERP004598.screen[which(ERP004598.screen$MapPercent >= 70),]
# # parasite proportion > 5%
# ERP004598.screen <- ERP004598.screen[which(ERP004598.screen$Parasite_percent >= 5),]
# I have 48 samples from ERP004598
ERP004598.screen.runs <- as.character(ERP004598.screen[,"RunID"])

# ERP110375 #

# Reads for protein-coding genes: Host - >= 10^6, Parasite - >= 10^5
ERP110375.screen <- ERP110375[which(ERP110375$ProteinCodHost >= 1e7 & ERP110375$ProteinCodPara >= 1e6),]
# Protein-coding genes: Host - >= 10^4, Parasite - >= 3000
ERP110375.screen <- ERP110375.screen[which(ERP110375.screen$NumberProtCodGenesHost >= 10000 & ERP110375.screen$NumberProtCodGenesPara > 3000),]
# Unique map percent: >= 70%
ERP110375.screen <- ERP110375.screen[which(ERP110375.screen$MapPercent >= 70),]
# # parasite proportion > 5%
# ERP110375.screen <- ERP110375.screen[which(ERP110375.screen$Parasite_percent >= 5),]
# I have 48 samples from ERP110375
ERP110375.screen.runs <- as.character(ERP110375.screen[,"RunID"])

# ERP002273 #

# Reads for protein-coding genes: Host - >= 10^6, Parasite - >= 10^5
ERP002273.screen <- ERP002273[which(ERP002273$ProteinCodHost >= 1e6 & ERP002273$ProteinCodPara >= 1e5),]
# Protein-coding genes: Host - >= 10^4, Parasite - >= 3000
ERP002273.screen <- ERP002273.screen[which(ERP002273.screen$NumberProtCodGenesHost >= 10000 & ERP002273.screen$NumberProtCodGenesPara > 3000),]
# Unique map percent: >= 70%
ERP002273.screen <- ERP002273.screen[which(ERP002273.screen$MapPercent >= 70),]
# # parasite proportion > 5%
# ERP002273.screen <- ERP002273.screen[which(ERP002273.screen$Parasite_percent >= 5),]
# I have 48 samples from ERP002273
ERP002273.screen.runs <- as.character(ERP002273.screen[,"RunID"])

# col.num <- c()
# for(i in 1:length(SRP032775.screen.runs))
# {
#   id <- grep(pattern = SRP032775.screen.runs[i], colnames(ortho_data))
#   if(length(id)==1)
#     col.num[i] <- id
# }
# 
# SRP032775.ortho.data2 <- ortho_data[,col.num]
# save(SRP032775.ortho.data2, file = "Data/SRP032775.ortho.data2.RData")

# col.num <- c()
# for(i in 1:length(SRP108356.screen.runs))
# {
#   id <- grep(pattern = SRP108356.screen.runs[i], colnames(ortho_data))
#   if(length(id)==1)
#     col.num[i] <- id
# }
# 
# SRP108356.ortho.data <- ortho_data[,col.num]
# save(SRP108356.ortho.data, file = "Data/SRP108356.ortho.data.RData")

# Reads for protein-coding genes: Host - >= 10^6, Parasite - >= 10^5
SRP096160.screen <- SRP096160[which(SRP096160$ProteinCodHost >= 1e6 & SRP096160$ProteinCodPara >= 1e5),]
# Protein-coding genes: Host - >= 10^4, Parasite - >= 3000
SRP096160.screen <- SRP096160.screen[which(SRP096160.screen$NumberProtCodGenesHost >= 10000 & SRP096160.screen$NumberProtCodGenesPara > 3000),]
# Unique map percent: >= 70%
SRP096160.screen <- SRP096160.screen[which(SRP096160.screen$MapPercent >= 70),]
# # parasite proportion > 5%
# SRP096160.screen <- SRP096160.screen[which(SRP096160.screen$Parasite_percent >= 5),]
# I have 48 samples from SRP096160
SRP096160.screen.runs <- as.character(SRP096160.screen[,"RunID"])

# SRP110282 #

# Reads for protein-coding genes: Host - >= 10^6, Parasite - >= 10^5
SRP110282.screen <- SRP110282[which(SRP110282$ProteinCodHost >= 1e6 & SRP110282$ProteinCodPara >= 1e5),]
# Protein-coding genes: Host - >= 10^4, Parasite - >= 3000
SRP110282.screen <- SRP110282.screen[which(SRP110282.screen$NumberProtCodGenesHost >= 10000 & SRP110282.screen$NumberProtCodGenesPara > 3000),]
# Unique map percent: >= 70%
SRP110282.screen <- SRP110282.screen[which(SRP110282.screen$MapPercent >= 70),]
# # parasite proportion > 5%
# SRP110282.screen <- SRP110282.screen[which(SRP110282.screen$Parasite_percent >= 5),]
# I have 48 samples from SRP110282
SRP110282.screen.runs <- as.character(SRP110282.screen[,"RunID"])


col.num <- c()
for(i in 1:length(ERP106451.screen.runs))
{
  id <- grep(pattern = ERP106451.screen.runs[i], colnames(ortho_data))
  if(length(id)==1)
    col.num[i] <- id
}

ERP106451.ortho.data.str <- ortho_data[,col.num]
save(ERP106451.ortho.data.str, file = "Data/ERP106451.ortho.data.str.RData")
    
col.num <- c()
for(i in 1:length(SRP118996.screen.runs))
{
  id <- grep(pattern = SRP118996.screen.runs[i], colnames(ortho_data))
  if(length(id)==1)
    col.num[i] <- id
}

SRP118996.ortho.data.stringent22 <- ortho_data[,col.num]
save(SRP118996.ortho.data.stringent22, file = "Data/SRP118996.ortho.data.stringent22.RData")
  
## DRP000987 allruns ##
col.num <- c()
for(i in 1:length(DRP000987.screen.runs))
{
  id <- grep(pattern = DRP000987.screen.runs[i], colnames(ortho_data))
  if(length(id)==1)
    col.num[i] <- id
}

DRP000987.ortho.data.str <- ortho_data[,col.num]
save(DRP000987.ortho.data.str, file = "Data/DRP000987.ortho.data.str.RData")
##

## SRP118503  ##
# col.num <- c()
# for(i in 1:length(SRP118503.screen.runs))
# {
#   id <- grep(pattern = SRP118503.screen.runs[i], colnames(ortho_data))
#   if(length(id)==1)
#     col.num[i] <- id
# }
# 
# SRP118503.ortho.data.allruns <- ortho_data[,col.num]
# save(SRP118503.ortho.data.allruns, file = "Data/SRP118503.ortho.data.allruns.RData")
##

## SRP118827  ##
col.num <- c()
for(i in 1:length(SRP118827.screen.runs))
{
  id <- grep(pattern = SRP118827.screen.runs[i], colnames(ortho_data))
  if(length(id)==1)
    col.num[i] <- id
}

SRP118827.ortho.data.str <- ortho_data[,col.num]
save(SRP118827.ortho.data.str, file = "Data/SRP118827.ortho.data.str.RData")

## SRP116793  ##
col.num <- c()
for(i in 1:length(SRP116793.screen.runs))
{
  id <- grep(pattern = SRP116793.screen.runs[i], colnames(ortho_data))
  if(length(id)==1)
    col.num[i] <- id
}

SRP116793.ortho.data.all <- ortho_data[,col.num]
save(SRP116793.ortho.data.all, file = "Data/SRP116793.ortho.data.all.RData")

## SRP116593  ##
col.num <- c()
for(i in 1:length(SRP116593.screen.runs))
{
  id <- grep(pattern = SRP116593.screen.runs[i], colnames(ortho_data))
  if(length(id)==1)
    col.num[i] <- id
}

SRP116593.ortho.data.all <- ortho_data[,col.num]
save(SRP116593.ortho.data.all, file = "Data/SRP116593.ortho.data.all.RData")

## ERP023982  ##
col.num <- c()
for(i in 1:length(ERP023982.screen.runs))
{
  id <- grep(pattern = ERP023982.screen.runs[i], colnames(ortho_data))
  if(length(id)==1)
    col.num[i] <- id
}

ERP023982.ortho.data.str <- ortho_data[,col.num]
save(ERP023982.ortho.data.str, file = "Data/ERP023982.ortho.data.str.RData")

## ERP004598  ##
col.num <- c()
for(i in 1:length(ERP004598.screen.runs))
{
  id <- grep(pattern = ERP004598.screen.runs[i], colnames(ortho_data))
  if(length(id)==1)
    col.num[i] <- id
}

ERP004598.ortho.data.str <- ortho_data[,col.num]
save(ERP004598.ortho.data.str, file = "Data/ERP004598.ortho.data.str.RData")

## ERP110375  ##
col.num <- c()
for(i in 1:length(ERP110375.screen.runs))
{
  id <- grep(pattern = ERP110375.screen.runs[i], colnames(ortho_data))
  if(length(id)==1)
    col.num[i] <- id
}

ERP110375.ortho.data.str <- ortho_data[,col.num]
save(ERP110375.ortho.data.str, file = "Data/ERP110375.ortho.data.str.RData")

## ERP002273  ##
col.num <- c()
for(i in 1:length(ERP002273.screen.runs))
{
  id <- grep(pattern = ERP002273.screen.runs[i], colnames(ortho_data))
  if(length(id)==1)
    col.num[i] <- id
}

ERP002273.ortho.data.all <- ortho_data[,col.num]
save(ERP002273.ortho.data.all, file = "Data/ERP002273.ortho.data.all.RData")

## RP096160  ##
col.num <- c()
for(i in 1:length(SRP096160.screen.runs))
{
  id <- grep(pattern = SRP096160.screen.runs[i], colnames(ortho_data))
  if(length(id)==1)
    col.num[i] <- id
}

SRP096160.ortho.data.int <- ortho_data[,col.num]
save(SRP096160.ortho.data.int, file = "Data/SRP096160.ortho.data.int.RData")

## RP110282  ##
col.num <- c()
for(i in 1:length(SRP110282.screen.runs))
{
  id <- grep(pattern = SRP110282.screen.runs[i], colnames(ortho_data))
  if(length(id)==1)
    col.num[i] <- id
}

SRP110282.ortho.data.int <- ortho_data[,col.num]
save(SRP110282.ortho.data.int, file = "Data/SRP110282.ortho.data.int.RData")


####### Make the big dataset #######

df <- DRP000987.ortho.data.str
df <- cbind(df, ERP106451.ortho.data, ERP110375.ortho.data.int, ERP004598.ortho.data.all, SRP118827.ortho.data.int, SRP116793.ortho.data.all)
df_concat_allhosts <- df
save(df_concat_allhosts2, file = "df_concat_allhosts2.RData")
write.table(df_concat_allhosts2, "df_concat_allhosts2.txt", sep = '\t')

############## Prep Kai's data #############

Mph <- read.delim("~/Documents/Data/Kai/Macrophage_RNAseqReadCount.txt", 
                  stringsAsFactors=FALSE) %>%
  tibble::column_to_rownames("ID")
# # divide into host and parasite part
Mph.h <- Mph[grep(pattern = "ENS", rownames(Mph)),] %>% tibble::rownames_to_column("ID")
Mph.p <- Mph[-grep(pattern = "ENS", rownames(Mph)),] %>% tibble::rownames_to_column("ID")
# # get runs that have 1e6 host reads and 1e5 para reads
geneSums.h <- colSums(Mph.h)
geneSums.p <- colSums(Mph.p)
runs.h <- colnames(Mph.h)[which(geneSums.h >= 1e6)]
runs.p <- colnames(Mph.p)[which(geneSums.p >= 1e5)]
# # Protein-coding genes: Host - >= 10^4, Parasite - >= 3000
geneCounts.h <- c()
for(i in 1:ncol(Mph.h))
{
  counts = 0
  for(j in 1:nrow(Mph.h))
  {
    if(Mph.h[j,i]>0)
      counts = counts +1
  }
  geneCounts.h[i] <- counts
}

geneCounts.p <- c()
for(i in 1:ncol(Mph.p))
{
  counts = 0
  for(j in 1:nrow(Mph.p))
  {
    if(Mph.p[j,i]>0)
      counts = counts +1
  }
  geneCounts.p[i] <- counts
}

names(geneCounts.h) <- colnames(Mph.h)
names(geneCounts.p) <- colnames(Mph.p)

runs.g.h <- colnames(Mph.h)[which(geneCounts.h >= 10000)]
runs.g.p <- colnames(Mph.p)[which(geneCounts.p >= 3000)]

# # Unique map percent: >= 70%
# holds good for all runs

# to find the runs that fit all of the criteria:
Reduce(intersect, list(runs.p, runs.h, runs.g.p, runs.g.h))

##  or use Mph_all
# map data to orthogroups

parasite_orthogroups <- read.delim("~/Downloads/parasite_orthogroups.txt", stringsAsFactors=FALSE) 
parasite_orthogroups <- parasite_orthogroups[,c(1,4)]
colnames(parasite_orthogroups)[2] <- "ID"
host_orthogroups <- read.delim("~/Downloads/host_orthogroups.txt", stringsAsFactors=FALSE)
host_orthogroups <- host_orthogroups[,c(1,3)]
colnames(host_orthogroups)[2] <- "ID"

h <- inner_join(host_orthogroups, Mph.h) %>% 
  tibble %>%
  dplyr::select(-c(2)) %>% 
  tibble::column_to_rownames("Orthogroup")

p <- inner_join(parasite_orthogroups, Mph.p) %>% 
  tibble %>%
  dplyr::select(-c(2)) %>% 
  tibble::column_to_rownames("Orthogroup")

Mph.ortho.data <- rbind(h, p)
save(Mph.ortho.data, file = "Mph.ortho.data.RData")
write.table(Mph.ortho.data, "Mph.ortho.data", sep = '\t', row.names = T)
