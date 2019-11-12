# "Only when using tagwise dispersion will genes that are
# consistent between replicates be ranked more highly than genes that are not"

studyIDs <- as.character(read.table("/SAN/Plasmo_compare/SRAdb/Input/positive_experiments.txt", sep = '\t', header = F)[,1])

for(i in 5:length(studyIDs))
{
  print(studyIDs[i])
  study <- read.delim(paste0("Output/", studyIDs[i], "/", studyIDs[i], "_coding_genes.txt", collapse = ''))
  
  #study <- t(study)
  study_disp <- estimateDisp(study, prior.n = 10)$tagwise.dispersion

  study_w_disp <- cbind(study, study_disp)
  disp_range <- range(study_w_disp[,"study_disp"])
  
  # order according to dispersion values
  disp_ordered <- study_w_disp[order(-study_w_disp[,"study_disp"]),]
  
  # get number of rows to be included: 0.95 * number of rows (from ordered matrix)
  disp_95 <- as.integer(0.95 * nrow(disp_ordered))
  
  # obtain the top 95% rows (== top 95% dispersion)
  disp_include <- disp_ordered[1:disp_95,]
  
  # write table
  write.table(disp_include, paste0("Output/",studyIDs[i],"/",studyIDs[i],"_filteredDisp.txt", collapse=''), sep = '\t', row.names = T)

  write.table(study_w_disp, paste0("/SAN/Plasmo_compare/SRAdb/Output/",studyIDs[i],"/",studyIDs[i],"_disp.txt",collapse = ''), sep = '\t', row.names = T)
}

# I push changes! Tra la la la la
