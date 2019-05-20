options(echo=TRUE)
args <- commandArgs(TRUE)
study <- args[1]

positive_experiments <- read.table("/SAN/Plasmo_compare/SRAdb/Input/positive_experiments.txt", header = F, sep = ',')
host <- positive_experiments[grep(study,positive_experiments[,1]),2]
para <- positive_experiments[grep(study,positive_experiments[,1]),3]
concated <- read.csv2(paste0("/SAN/Plasmo_compare/SRAdb/Output/", study,"/ConcatRunsToStudy_", study,".txt", collapse = ''), sep = '\t', header = T)
anno <- concated[,1:5]

current_runs <- read.csv2(paste0("/SAN/Plasmo_compare/SRAdb/Output/",study,"/runs_",study,".txt", collapse=''), header = FALSE, sep = ',')
runs_in_study <- current_runs[grep(study, current_runs[,2]),1]
runs <- runs_in_study

for(i in 1:length(runs))
{
  if(file.exists(paste0("/SAN/Plasmo_compare/SRAdb/Output/",study,"/count_",runs[i],".txt")))
  {
    print(i)
    count_df <- read.csv2(file=paste0("/SAN/Plasmo_compare/SRAdb/Output/",study,"/count_",runs[i],".txt"), sep='\t', header=F)
    
    parasite_rows <- grep(paste0(substr(para,1,3),"_chr",collapse=''), anno[,1])
    all_parasite <- anno[parasite_rows,]
    all_host <- anno[-parasite_rows,]
    
    all_host1 <- sapply(all_host[,1], function(x) paste0("Hs_chr",x, collapse = ''))
    all_host[,1] <- all_host1
    host_para <- rbind(all_host, all_parasite)
    run_file <- cbind(host_para, count_df)
    colnames(run_file)[6] <- c("counts")
    
    write.table(run_file, paste0("/SAN/Plasmo_compare/SRAdb/Output/",study,"/countWithGFF3_",runs[i],".txt",collapse=''), row.names = F, sep = '\t')
  }
}
