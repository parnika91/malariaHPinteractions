library(dplyr)
library(rtracklayer)

allhosts_bipartite_0.95 <- read.delim("~/Documents/Data/allhosts_bipartite_0.95.txt", stringsAsFactors=FALSE)
View(allhosts_bipartite_0.95)
ho <- read.delim("~/Downloads/host_orthogroups.txt", stringsAsFactors=FALSE)
View(ho)
po <- read.delim("~/Downloads/parasite_orthogroups.txt", stringsAsFactors=FALSE)
View(po)
ho <- ho[,c(1,3)]
po <- po[,c(1,4)]
View(ho)
View(po)

host_part <- as.data.frame(allhosts_bipartite_0.95[,1])
colnames(host_part) <- c("Orthogroup")
para_part <- as.data.frame(allhosts_bipartite_0.95[,2])
colnames(para_part) <- c("Orthogroup")

h <- merge(host_part, ho)
p <- merge(para_part, po)

hp <- cbind(h,p)
write.table(hp, "hp_0.95.txt", sep = '\t', row.names = F)
