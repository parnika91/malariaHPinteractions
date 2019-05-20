# get positive experiments
positive_experiments <- read.table("/SAN/Plasmo_compare/SRAdb/Input/positive_experiments.txt", header = F, sep = ',')
studyIDs <- positive_experiments[,1]

# for each study, go to every countWithGFF_<run>.txt and add gene ids

for(i in 1:length(studyIDs))
{
  runs_in_study <- read.table(paste0("/SAN/Plasmo_compare/SRAdb/Output/",studyIDs[i],"/runs_",studyIDs[i],".txt",collapse=''), header = F, sep = ',')
  runs_in_study <- runs_in_study[,1]
  
  host <- positive_experiments[grep(studyIDs[i],positive_experiments[,1]),2]
  para <- positive_experiments[grep(studyIDs[i],positive_experiments[,1]),3]
  
  hp.gff <- makeTxDbFromGFF(paste0("/SAN/Plasmo_compare/SRAdb/Genomes/annotation/",host,para,".gff3",collapse=''), format = "gff3")
  hp.genes <- genes(hp.gff)
  
  for(j in 1:length(runs_in_study))
  {
    countFile <- read.table(paste0("/SAN/Plasmo_compare/SRAdb/Output/",studyIDs[i],"/countWithGFF3_",runs_in_study[j],".txt", collapse=''), sep='\t', header = T)
    hp.genes.df <- as.data.frame(hp.genes)
    countFilegeneID <- merge(hp.genes.df, countFile, by = c("seqnames","start", "end", "width", "strand"))
    write.table(countFilegeneID, paste0("/SAN/Plasmo_comapre/SRAdb/Output/",studyIDs[i],"/countWithGFFgeneID_",runs_in_study[j],".txt", collapse=''), row.names=F, sep='\t')
  }
}

# hp.gff <- makeTxDbFromGFF(paste0("/SAN/Plasmo_compare/SRAdb/Genomes/annotation/",host,para,".gff3",collapse=''), format = "gff3")
# # Warning message:
# #  In .extract_exons_from_GRanges(exon_IDX, gr, ID, Name, Parent, feature = "exon",  :
# #                                   The following orphan exon were dropped (showing only the 6 first):
# #                                   seqid    start      end strand   ID                        Parent               Name
# #                                 1 Mmul_chr1  8251952  8252121      + <NA> transcript:ENSMMUT00000053157 ENSMMUE00000321551
# #                                 2 Mmul_chr1 12513775 12513920      - <NA> transcript:ENSMMUT00000052526 ENSMMUE00000320920
# #                                 3 Mmul_chr1 14422713 14422852      + <NA> transcript:ENSMMUT00000051348 ENSMMUE00000319742
# #                                4 Mmul_chr1 24808763 24808886      + <NA> transcript:ENSMMUT00000051782 ENSMMUE00000320176
# #                                 5 Mmul_chr1 24808972 24809054      + <NA> transcript:ENSMMUT00000051421 ENSMMUE00000319815
# #                                6 Mmul_chr1 26715847 26716012      + <NA> transcript:ENSMMUT00000052545 ENSMMUE00000320939
# hp.genes <- genes(hp.gff) # gives gene id and gene range but not ensembl gene ids
# write.table(hp.genes, "hp.genes.txt", sep ='\t', row.names = F)
# 
# temp.countERR1744068 <- read.table("/SAN/Plasmo_compare/SRAdb/Output/ERP020067/countWithGFF3_ERR1744068.txt", sep='\t', header = T)
# hp.genes.df <- as.data.frame(hp.genes)
# temp.countERR1744068.geneID <- merge(hp.genes.df, temp.countERR1744068, by = c("seqnames","start", "end", "width", "strand"))
# write.table(temp.countERR1744068.geneID, "temp.count.txt", row.names=F, sep='\t')
