positive_experiments <- read.table("/SAN/Plasmo_compare/SRAdb/Input/positive_experiments.txt", header = F, sep = ',')
studyIDs <- positive_experiments[,1]
studyIDs <- studyIDs[-32]

study_inputReads <- data.frame()

for(j in 1:length(studyIDs))
{
  input_reads_sum <- 0
  runIDs <- read.table(paste0("/SAN/Plasmo_compare/SRAdb/Output/",studyIDs[j],"/runs_",studyIDs[j],".txt",collapse=''), header = F, sep = ',')
  number_of_runs <- nrow(runIDs)
  
  for(i in 1:number_of_runs)
  {
    # get runID
    runID <- runIDs[i,1]
    
    inputReadLine <- readLines(paste0("/SAN/Plasmo_compare/SRAdb/Output/",studyIDs[j],"/",runID,"_",studyIDs[j],".final.out",collapse=''), n=6)[6]
    inputReads <- as.numeric(strsplit(inputReadLine, split="\t")[[1]][2])
    
    input_reads_sum <- input_reads_sum + inputReads
  }
  
  study_inputReads[j,1] <- studyIDs[j]
  study_inputReads[j,2] <- input_reads_sum
  
}

colnames(study_inputReads) <- c("Study", "InputReads")
write.table(study_inputReads, "study_inputReads.txt", sep = '\t', row.names=F)