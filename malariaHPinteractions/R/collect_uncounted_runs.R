# get all runs of study

ERP004042_all <- read.csv2(paste0("/SAN/Plasmo_compare/SRAdb/Output/",study,"/runs_",study,".txt", collapse=''), header = FALSE, sep = ',')

# get runs that were counted

l <- c()
for(i in 1:nrow(ERP004042_all))
  l <- c(l, list.files(path="/SAN/Plasmo_compare/SRAdb/Output/ERP004042/", pattern = paste0("countWithGFF3_",ERP004042_all[i,1],".txt",collapse='')))

counted_runs.txt <- sapply(strsplit(l,"_"), "[[", 2)
counted_runs <- unlist(strsplit(counted_runs.txt,".txt"))

# get uncounted runs

uncounted_runs <- ERP004042_all[!(ERP004042_all[,1] %in% counted_runs),1]

  
# remove them from blacklist