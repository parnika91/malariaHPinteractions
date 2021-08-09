# libraries
library(tidyverse)
library(seqinr)

# load string db interactions
Pb <- read.csv("Downloads/5821.protein.links.full.v11.0.txt", sep="")
Mm <- read.csv("10090.protein.links.full.v11.0.txt", sep="")

# load orthogroups
pOG <- read.delim("Downloads//parasite_orthogroups.txt")
hOG <- read.delim("Downloads//host_orthogroups.txt")

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
 
  # universe of domain interactions
  domains <- na.omit(unique(as.character(ppi_og_d$Accession)))
  d_uni <- expand.grid(domains, domains)
  
  # takes a very long time
  x <- d_uni[!duplicated(t(apply(d_uni, 1, sort))), ]
  y <- ddi[!duplicated(t(apply(ddi, 1, sort))), ]
  
 }

 