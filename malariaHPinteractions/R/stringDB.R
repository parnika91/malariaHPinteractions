# libraries
library(tidyverse)
library(seqinr)

# load string db interactions
Pb <- read.csv("~/Downloads/5821.protein.links.full.v11.0.txt", sep="")
Mm <- read.csv("~/Downloads/10090.protein.links.full.v11.0.txt", sep="")

# load orthogroups
pOG <- read.delim("~/Documents/Data/parasite_orthogroups.txt")
hOG <- read.delim("~/Documents/Data/host_orthogroups.txt")

# remove taxon ID from protein names from stringdb files
Pb[,1] <- substr(Pb[,1], 6, nchar(Pb[,1]))
Pb[,2] <- substr(Pb[,2], 6, nchar(Pb[,2]))

Mm[,1] <- substr(Mm[,1], 7, nchar(Mm[,1]))            
Mm[,2] <- substr(Mm[,2], 7, nchar(Mm[,2]))

# Make all possible combinations with the genes in orthogroups
Pb_comb <- expand.grid(pOG$Pberghei, pOG$Pberghei); colnames(Pb_comb) <- c("protein1", "protein2")
Pb_stringdb <- inner_join(Pb[,1:2], Pb_comb)

Mm_comb <- expand.grid(hOG$mouse, hOG$mouse); colnames(Mm_comb) <- c("protein1", "protein2")
Mm_stringdb <- inner_join(Mm[,1:2], Mm_comb)

# get a list of proteins so that we can obtain their PFAM domains
Pb_list <- unique(c(Pb_stringdb[,1], Pb_stringdb[,2]))
Pb_list <- as.data.frame(Pb_list)
write.csv(Pb_list, "Pb_list_Pfram_entries.csv", row.names = F)

### Annotated protein sequences downloaded of version 49, 2020-10-27 from PlasmoDB
# load protein sequence fasta file
Pb_seq <- read.fasta("Downloads/Pberghei_proteins.fasta", seqtype = "AA", as.string = T)
# remove ".1-p1" from the names
names(Pb_seq)<- substr(names(Pb_seq), 1, nchar(names(Pb_seq))-5)
Pb_entries_names <- sapply(pOG$Pberghei, function(x) grep(pattern = x, names(Pb_seq)))
Pb_entries_seq <- Pb_seq[unlist(Pb_entries_names)]

write.fasta(Pb_entries_seq[1:2000], names = names(Pb_entries_seq)[1:2000], file.out = "Pberghei_Pfam_entries_first.fasta")
write.fasta(Pb_entries_seq[2001:4006], names = names(Pb_entries_seq)[2001:4006], file.out = "Pberghei_Pfam_entries_second.fasta")
 