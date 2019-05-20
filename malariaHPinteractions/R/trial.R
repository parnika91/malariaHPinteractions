library(doParallel)
registerDoParallel(cores=9)

x <- iris[which(iris[,5] != "setosa"), c(1,5)]
trials <- 10000
ptime <- system.time({
r <- foreach(icount(trials), .combine=cbind) %dopar% {
     ind <- sample(100, 100, replace=TRUE)
     result1 <- glm(x[ind,2]~x[ind,1], family=binomial(logit))
     coefficients(result1)
   }
})[3]
ptime

library(dplyr)
library(multidplyr)
library(nycflights13)
planes %>% partition() %>% group_by(type) %>% summarize(n())

concatfile <- read.delim("Output/SRP034011/ConcatRunsToStudy_SRP034011.txt")

parasite_row <- concatfile[grep("P+", concatfile[,1]),]
host_row <- setdiff(concatfile, parasite_row)

host_row$seqnames <- paste0("Hs_chr", host_row$seqnames)

concat <- rbind(host_row, parasite_row)

write.table(concat, "Output/SRP034011/ConcatRunsToStudy_SRP034011.txt", sep = '\t', row.names = F)
