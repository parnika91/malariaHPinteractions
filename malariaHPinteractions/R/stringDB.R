# libraries
library(tidyverse)
library(seqinr)


# Get orthogroups
pOG <- read.delim("~/Documents/Data/parasite_orthogroups.txt")

### Get Pb PPIs
Pb <- read.csv("~/Downloads/5821.protein.links.full.v11.0.txt", sep="", stringsAsFactors=FALSE)
Mm <- read.csv("~/Downloads/10090.protein.links.full.v11.0.txt", sep="", stringsAsFactors=FALSE)

### Annotated protein sequences downloaded of version 49, 2020-10-27 from PlasmoDB
## load protein sequence fasta file

# Pb_seq <- read.fasta("Downloads/Pberghei_proteins.fasta", seqtype = "AA", as.string = T)
# ## remove ".1-p1" from the names
# names(Pb_seq)<- substr(names(Pb_seq), 1, nchar(names(Pb_seq))-5)
# Pb_entries_names <- sapply(pOG$Pberghei, function(x) grep(pattern = x, names(Pb_seq)))
# Pb_entries_seq <- Pb_seq[unlist(Pb_entries_names)]
# 
# write.fasta(Pb_entries_seq[1:2000], names = names(Pb_entries_seq)[1:2000], file.out = "Pberghei_Pfam_entries_first.fasta")
# write.fasta(Pb_entries_seq[2001:4006], names = names(Pb_entries_seq)[2001:4006], file.out = "Pberghei_Pfam_entries_second.fasta")

# get domains and remove superfamily rows
d1 <- read.csv("Downloads/PbergheiPFAM_1_2_hitdata.txt", stringsAsFactors=TRUE) %>%
  filter(!(Hit.type == "superfamily"))
d2 <- read.csv("Downloads/PbergheiPFAM_2_2_hitdata.txt", stringsAsFactors=TRUE) %>%
  filter(!(Hit.type == "superfamily"))
d <- rbind(d1, d2)

### for mouse
dm <- data.frame()
for(i in 1:17)
{
  d <- read.delim(paste0("Downloads/mouse", i, "_hitdata.txt", collapse = ""), stringsAsFactors=TRUE) %>%
    filter(grepl(pattern = "pfam", Accession))
  dm <- rbind(dm, d)
}

### intra-species DDIs function
intraspecies_DDIs <- function(PPI, OG, D){
  
  # remove taxonID from PPIs
  PPI[,1] <- sapply(PPI[,1], function(x) strsplit(x, split = "\\.")[[1]][2])
  PPI[,2] <- sapply(PPI[,2], function(x) strsplit(x, split = "\\.")[[1]][2])
  
  # Put Orthogroup names next to protein names
  P1 <- as.data.frame(PPI[,1]); colnames(P1) <- "Pberghei"
  P2 <- as.data.frame(PPI[,2]); colnames(P2) <- "Pberghei"
  OG1 <- inner_join(P1, OG[,c("Orthogroup", "Pberghei")])
  OG2 <- inner_join(P2, OG[,c("Orthogroup", "Pberghei")])
  
  ppi_og <- cbind(OG1, OG2); colnames(ppi_og) <- c("P1", "OG1", "P2", "OG2")
  
  # add their domain and expand the rows
  ## remove "Q - >" part of Query column of domains to retain only the protein name
  D[,1] <- sapply(D[,1], function(x) strsplit(as.character(x), split = ">")[[1]][2])
  D[,1] <- substr(D[,"Query"], 1, 13) # conversion to previous ID
  # add domains to proteins -> could be multiple rows for one protein
  # P1 OG1 D1
  # P1 OG1 D2
  ppi_og_d <- rbind(OG1, OG2) %>%
    distinct() %>%
    left_join(., D[,c("Query", "Accession")], by = c("Pberghei" = "Query"))
  
  # put this info into the ppi_og dataframe
  dom1 = left_join(ppi_og, ppi_og_d, by = c("P1" = "Pberghei")) %>%
    rename(Dom1 = Accession)
  
  dom2 = left_join(dom1, ppi_og_d, by = c("P2" = "Pberghei")) %>%
    select(-c(5,7)) %>%
    na.omit() %>%
    rename(Dom2 = Accession)
  
  ddi <- dom2[,c("Dom1", "Dom2")]
  ddi_mat <- as.matrix(table(as.character(ddi$Dom1), as.character(ddi$Dom2)))
  ddi_df <- as.data.frame(table(as.character(ddi$Dom1), as.character(ddi$Dom2)))
  ddi_df[,1:2] <- lapply(ddi_df[,1:2], as.character)
  
  ddi_mat <- spread(ddi_df, Var2, Freq) %>% tibble::column_to_rownames("Var1")
 
  # universe of domain interactions
  domains <- na.omit(unique(as.character(ppi_og_d$Accession)))
  d_uni <- expand.grid(domains, domains)
  
  # takes a very long time
  d_uni <- d_uni[!duplicated(t(apply(d_uni, 1, sort))), ]
  d_uni[] <- lapply(d_uni[], as.character)
  
  #y <- ddi[!duplicated(t(apply(ddi, 1, sort))), ]
  
  # ddpair <- matrix(nrow = 10, ncol = 10)
  # 
  # ptm <- proc.time()
  # for(i in 1:10)
  #   for(j in 1:10)
  #   {
  #     ddpair[i,j] = 
  #       chisq.test(table(Dom1 = ddi$Dom1 %in% domains[i], 
  #                        Dom2 = ddi$Dom2 %in% domains[j]) + 
  #                  table(Dom2 = ddi$Dom2 %in% domains[i], 
  #                        Dom1 = ddi$Dom1 %in% domains[j]))$p.value
  #   }
  # proc.time() - ptm
  # rownames(ddpair) <- colnames(ddpair) <- domains
  
  # if(any(grepl(pattern = d_uni$Var1[i], ddi_df$Var1)) & 
  #    any(grepl(pattern = d_uni$Var2[i], ddi_df$Var2)) & 
  #    any(grepl(pattern = d_uni$Var2[i], ddi_df$Var1)) &
  #    any(grepl(pattern = d_uni$Var1[i], ddi_df$Var2)))
  # {
  # tt = sum(ddi_df[which(ddi_df$Var1==d_uni$Var1[i] & ddi_df$Var2==d_uni$Var2[i]), "Freq"],
  #          ddi_df[which(ddi_df$Var1==d_uni$Var2[i] & ddi_df$Var2==d_uni$Var1[i]), "Freq"])
  # 
  # tf = sum(ddi_df[which(ddi_df$Var1 == d_uni$Var1[i]), "Freq"]) + 
  #   sum(ddi_df[which(ddi_df$Var2 == d_uni$Var1[i]), "Freq"]) - 
  #   tt
  # 
  # ft = sum(ddi_df[which(ddi_df$Var1 == d_uni$Var2[i]), "Freq"]) + 
  #   sum(ddi_df[which(ddi_df$Var2 == d_uni$Var2[i]), "Freq"]) - 
  #   tt
  # 
  # ff = sum(ddi_df$Freq) - 
  #   (sum(ddi_df[which(ddi_df$Var2==d_uni$Var1[i]), "Freq"]) + sum(ddi_df[which(ddi_df$Var1==d_uni$Var1[i]), "Freq"]) + 
  #      sum(ddi_df[which(ddi_df$Var2==d_uni$Var2[i]), "Freq"]) + sum(ddi_df[which(ddi_df$Var1==d_uni$Var2[i]), "Freq"])) + 
  #   sum(ddi_df[which(ddi_df$Var1==d_uni$Var1[i] & ddi_df$Var2==d_uni$Var2[i]), "Freq"], 
  #       ddi_df[which(ddi_df$Var1==d_uni$Var2[i] & ddi_df$Var2==d_uni$Var1[i]), "Freq"], 
  #       ddi_df[which(ddi_df$Var1==d_uni$Var1[i] & ddi_df$Var2==d_uni$Var1[i]), "Freq"],
  #       ddi_df[which(ddi_df$Var1==d_uni$Var2[i] & ddi_df$Var2==d_uni$Var2[i]), "Freq"])
  # 
  # ___________________________
  
 enriched_DDI <- function(d_uni)
 {
     ddipair <- data.frame()
     for(i in 1:nrow(d_uni))
   {
     ddipair[i,1] <- d_uni$Var1[i]
     ddipair[i,2] <- d_uni$Var2[i]
     
     if(any(grepl(pattern = d_uni$Var1[i], rownames(ddi_mat))) &
        any(grepl(pattern = d_uni$Var2[i], colnames(ddi_mat))) & 
        any(grepl(pattern = d_uni$Var2[i], rownames(ddi_mat))) &
        any(grepl(pattern = d_uni$Var1[i], colnames(ddi_mat))))
     {
       tt = sum(ddi_mat[which(rownames(ddi_mat)==d_uni$Var1[i]), which(colnames(ddi_mat)==d_uni$Var2[i])],
                ddi_mat[which(rownames(ddi_mat)==d_uni$Var2[i]), which(colnames(ddi_mat)==d_uni$Var1[i])])
       
       tf = sum(ddi_mat[which(rownames(ddi_mat) == d_uni$Var1[i]),]) + 
         sum(ddi_mat[,which(colnames(ddi_mat) == d_uni$Var1[i])]) - 
         tt
       
       ft = sum(ddi_mat[which(rownames(ddi_mat) == d_uni$Var2[i]),]) + 
         sum(ddi_mat[,which(colnames(ddi_mat) == d_uni$Var2[i])]) - 
         tt
       
       ff = sum(colSums(ddi_mat)) - 
         (sum(ddi_mat[,which(colnames(ddi_mat) == d_uni$Var1[i])]) + sum(ddi_mat[which(rownames(ddi_mat) == d_uni$Var1[i]),]) + 
            sum(ddi_mat[,which(colnames(ddi_mat) == d_uni$Var2[i])]) + sum(ddi_mat[which(rownames(ddi_mat) == d_uni$Var2[i]),])) + 
         sum(ddi_mat[which(rownames(ddi_mat) == d_uni$Var1[i]), which(colnames(ddi_mat) == d_uni$Var2[i])], 
             ddi_mat[which(rownames(ddi_mat) == d_uni$Var2[i]), which(colnames(ddi_mat) == d_uni$Var1[i])], 
             ddi_mat[which(rownames(ddi_mat) == d_uni$Var1[i]), which(colnames(ddi_mat) == d_uni$Var1[i])],
             ddi_mat[which(rownames(ddi_mat) == d_uni$Var2[i]), which(colnames(ddi_mat) == d_uni$Var2[i])])
       
       ddipair[i,3] <- chisq.test(matrix(c(ff, tf, ft, tt), nrow = 2, byrow = T))$p.value
     } else {
       ddipair[i,3] <- NA
     }
   }
   return(ddipair)
 }
}
 
 system.time(f <- enriched_DDI(d_uni[1:10,]))
 d_uni_split <- split(d_uni, cut(1:nrow(d_uni), 80)) # ncores = 80
 system.time(Pres <- mclapply(d_uni_split, enriched_DDI, mc.cores = 80))
 

#   ddipair <- data.frame()
#   a = 1
#   
#   system.time(
#     for(i in 1)
#       for(j in 1:ncol(ddi_mat))
#       {
#         domain1 <- rownames(ddi_mat)[i]
#         domain2 <- colnames(ddi_mat)[j]
#         
#         ddipair[a,1] <- domain1
#         ddipair[a,2] <- domain2
#         
#        tt = sum(ddi_mat[which(rownames(ddi_mat)==domain1), which(colnames(ddi_mat)==domain2)],
#                    ddi_mat[which(rownames(ddi_mat)==domain2), which(colnames(ddi_mat)==domain1)])
#           
#         tf = sum(ddi_mat[which(rownames(ddi_mat) == domain1),]) + 
#           sum(ddi_mat[,which(colnames(ddi_mat) == domain1)]) - 
#           tt
#         
#         ft = sum(ddi_mat[which(rownames(ddi_mat) == domain2),]) + 
#           sum(ddi_mat[,which(colnames(ddi_mat) == domain2)]) - 
#           tt
#         
#         ff = sum(colSums(ddi_mat)) - 
#           (sum(ddi_mat[,which(colnames(ddi_mat) == domain1)]) + sum(ddi_mat[which(rownames(ddi_mat) == domain1),]) + 
#              sum(ddi_mat[,which(colnames(ddi_mat) == domain2)]) + sum(ddi_mat[which(rownames(ddi_mat) == domain2),])) + 
#           sum(ddi_mat[which(rownames(ddi_mat) == domain1), which(colnames(ddi_mat) == domain2)], 
#               ddi_mat[which(rownames(ddi_mat) == domain2), which(colnames(ddi_mat) == domain1)], 
#               ddi_mat[which(rownames(ddi_mat) == domain1), which(colnames(ddi_mat) == domain1)],
#               ddi_mat[which(rownames(ddi_mat) == domain2), which(colnames(ddi_mat) == domain2)])
#         
#         ddipair[a,3] <- chisq.test(matrix(c(ff, tf, ft, tt), nrow = 2, byrow = T))$p.value
#         
#         a = a+1
#       })
#     
#  }
# 
# chisq.test(table(Pfam00004 = ddi$Dom1 %in% "pfam00004", Pfam00026 = ddi$Dom2 %in% "pfam00026") + table(Pfam00004 = ddi$Dom2 %in% "pfam00004", Pfam00026 = ddi$Dom1 %in% "pfam00026"))
# 
# table(Pfam00271 = ddi$Dom1 %in% "pfam00271", pfam00270 = ddi$Dom2 %in% "pfam00270") + table(pfam00271 = ddi$Dom2 %in% "pfam00271", pfam00270 = ddi$Dom1 %in% "pfam00270")
# 
# table(as.character(dom2$Dom1[1:10000]), as.character(dom2$Dom2[1:10000]))[1:5, 1:5]
# 
# df <- data.frame(A=c("j","K","Z"), B=c("i","P","Z"), C=c(100,101,102), ntimes=c(2,4,1))
# df <- as.data.frame(lapply(df, rep, df$ntimes))
# mat2 <- as.data.frame(lapply(mat2, rep, mat2$value))

 
 ################## mouse
 

 fa <- read.fasta("~/Documents/Data/Mus_musculus.GRCm39.pep.all.fa.txt", seqtype = "AA", as.string = T)
 
 protein <- sapply(fa, function(x) substr(strsplit(attr(x, "Annot"), split = " ")[[1]][1], 
                                          2, 
                                          nchar(strsplit(attr(x, "Annot"), split = " ")[[1]][1]))) %>%
   as.data.frame() %>%
   separate(".", c("A", NA))
 
 gene <- sapply(fa, function(x) substr(strsplit(attr(x, "Annot"), split = " ")[[1]][4], 
                                       6, 
                                       nchar(strsplit(attr(x, "Annot"), split = " ")[[1]][4]))) %>%
   as.data.frame() %>%
   separate(".", c("B", NA))
 
 
 pro_gene <- data.frame(Protein = protein, Gene = gene)