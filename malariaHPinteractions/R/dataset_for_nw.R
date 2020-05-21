library(dplyr)

load("Output/df_concat_allhosts/cor/df_concat_allhosts_all_bipartite.RData")
po<-read.delim("../OrthoFinder/parasite_orthogroups.txt", stringsAsFactors=F)
ho<-read.delim("../OrthoFinder/host_orthogroups.txt", stringsAsFactors=F)
ho_m <- as.data.frame(ho[,c(1,3)])
po_m <- as.data.frame(po[,c(1,4)])
bh <- as.data.frame(as.character(bipartite[,1]))
colnames(bh) <- "Orthogroup"
ph <- as.data.frame(as.character(bipartite[,2]))
colnames(ph) <- "Orthogroup"
hj <- left_join(bh, ho_m, by = "Orthogroup")
colnames(hj)[1] <- "h_OG"
pj <- left_join(ph, po_m, by = "Orthogroup")
colnames(pj)[1] <- "p_OG"

RON12 <- new_bip[which(new_bip$Pb_g == "PBANKA_050140"),]
RALP1 <- new_bip[which(new_bip$Pb_g == "PBANKA_061970"),]
RON5 <- new_bip[which(new_bip$Pb_g == "PBANKA_071310"),]
AMP <- new_bip[which(new_bip$Pb_g == "PBANKA_072180"),]
ASP <- new_bip[which(new_bip$Pb_g == "PBANKA_100360"),]
RON11 <- new_bip[which(new_bip$Pb_g == "PBANKA_132710"),]
rhoptry <- rbind(RON6, RON12, RALP1, RON5, AMP, ASP, RON11)

load("Documents/Data/rhoptry.RData")
