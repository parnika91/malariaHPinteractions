# suppressMessages(source("https://bioconductor.org/biocLite.R")) # for countOverlaps
# suppressMessages(biocLite("GenomicRanges"))
#suppressMessages(biocLite("GenomicFeatures"))
#suppressMessages(biocLite("GenomicAlignments"))
#suppressMessages(biocLite("rtracklayer"))

#suppressMessages(library(GenomicRanges))
#suppressMessages(library(GenomicFeatures))
#suppressMessages(library(GenomicAlignments))
#suppressMessages(library(rtracklayer))

 require(GenomicRanges)
 require(GenomicFeatures)
 require(GenomicAlignments)
 require(rtracklayer)

# check.packages <- function(pkg){
#   new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
#   if (length(new.pkg)) 
#     install.packages(new.pkg, dependencies = TRUE)
#   sapply(pkg, library, character.only = TRUE)
# }
# 
# # Usage example
# packages<-c("GenomicRanges", "GenomicFeatures", "GenomicAlignments", "rtracklayer")
# check.packages(packages)


# take host and parasite names as arguements from cmd
options(echo=TRUE)
args <- commandArgs(TRUE)
host <- args[1]
para <- args[2]
run <- args[3]
study <- args[4]

# # read gene annotation files for host and parasite
# host <- readGFF(file=paste0(host, ".gff3", collapse=''))
# # host[,1] <- paste0("Mm_chr", host[,1])
# # mouse_ens.gff3 <- host
# # save(mouse_ens.gff3, file="mouse_ens.gff3.RData")
# 
# 
# # parasite <- readGFF(file=paste0(parasite, ".gff3", collapse=''))
# # parasite[,1] <- paste0("Pyo_chr", parasite[,1])
# # Pyoelii_ens.gff3 <- parasite
# # save(Pyoelii_ens.gff3, file="Pyoelii_ens.gff3.RData")
# parasite <- readRDS(paste0(para, "_ens.gff3.rds"))
# 
# # host_GR <- GenomicRanges::makeGRangesFromDataFrame(host)
# #
# # add 'HS_chr' to all seqid
# # host$seqid <- paste("chr", host$seqid, sep="")
# # parasite$seqid <- paste("chr", parasite$seqid, sep="")
# 
# # extract only exons from host and parasite and convert into df
# exon_host <- host[host$type%in%"exon",]
# exon_parasite <- parasite[parasite$type%in%"exon",]
# 
# # exon_host_trimmed <- exon_host[, c("seqid", "source", "type", "start", "end","score", "strand",
#                                    # "phase","ID","Alias","Name","biotype","description","gene_id","logic_name",
#                                    # "Parent","transcript_id","constitutive","ensembl_end_phase",
#                                    # "ensembl_phase","exon_id","rank", "version" , "protein_id","external_name")]
# exon_host_trimmed <- exon_host[,c(colnames(parasite))]
# # levels(factor(exon_parasite[,1]))
# 
# # colnames(exon_parasite_trimmed)[6] <- "exon_id"
# # colnames(exon_host_trimmed)[7] <- "Parent"
# # exon_host_trimmed[,1] <-rep("Contig_HS")
# # exon_parasite_trimmed[,1] <-rep("Contig_PF")
# 
# # bind host df and parasite df
# exon_hp <- rbind(exon_host_trimmed, exon_parasite)
# # exon_hp_GR <- GenomicRanges::makeGRangesFromDataFrame(exon_hp)
# export(exon_hp, "exon_hp.gff3", format = "GFF3")

# concatGFF3 <- function(h, p)
# {
#   host <- readGFF(file=paste0(h, ".gff3", collapse=''))
#   parasite <- readGFF(file=paste0(p, ".gff3", collapse=''))
#   
#   exon_host <- host[host$type%in%"exon",]
#   exon_parasite <- parasite[parasite$type%in%"exon",]
#   
#   exon_host_trimmed <- exon_host[,c(colnames(parasite))]
#   exon_hp <- rbind(exon_host_trimmed, exon_parasite)
#   export(exon_hp, con=paste0(h,p,".gff3", collapse=''), format = "GFF3")
# }
# go to study folder
# setwd(paste0("/SAN/Plasmo_compare/SRAdb/",study))
# get list of bam files in folder; I am already inside a study folder,
# so I don't have to worry about going into a study folder

# make counts for every bam file 

#setwd("/SAN/Plasmo_compare/SRAdb/Genomes/annotation")

gff3 <- import(paste0("/SAN/Plasmo_compare/Genomes/annotation/", host, para, ".gtf",collapse=''), format = "gtf")
gff3<- gff3[gff3$type%in%"exon",]
gff3 <- GenomicRanges::makeGRangesFromDataFrame(gff3)

#setwd("/SAN/Plasmo_compare/fastq_download_tmp/")
#bamFile <- import(paste0("/SAN/Plasmo_compare/fastq_download_tmp/",run,"Aligned.out.bam", collapse=''))

setMethod("seqinfo", "BamFile", function (x) {
 h <- scanBamHeader(x, what = "targets")[["targets"]]
 h <- h[!duplicated(names(h))]
 Seqinfo(names(h), h)
})

bamFile <- readGAlignments(paste0("/SAN/Plasmo_compare/fastq_download_tmp/",run,"Aligned.sortedByCoord.out.bam", collapse=''))
counts <- countOverlaps(gff3, bamFile)
counts_df <- data.frame(gff3, counts)

#studyIDs <- "Macrophage_Kai"

#try(for(study in studyIDs)
#{
#  run <- read.csv2(paste0("Output/",study,"/runs_",study,".txt", collapse = ''), sep = ',', header = F)[,1]
#  for(i in 1:length(run))
#  {
#    print(run[i])
#    if(file.exists(paste0("/SAN/Plasmo_compare/fastq_download_tmp/",run[i],"Aligned.out.bam", collapse='')) & !file.exists(paste0("/SAN/Plasmo_compare/SRAdb/Output/",study,"/countWithGFF3_",run[i],".txt", collapse='')))
#      {
#        bamFile <- readGAlignments(paste0("/SAN/Plasmo_compare/fastq_download_tmp/",run[i],"Aligned.out.bam", collapse=''))
#        counts <- countOverlaps(gff3, bamFile)
#        counts_df <- data.frame(gff3, counts)
#          
#          #setwd(paste0("/SAN/Plasmo_compare/SRAdb/Output/",study, collapse=''))
#        write.table(counts_df, file=paste0("/SAN/Plasmo_compare/SRAdb/Output/",study,"/countWithGFF3_",run[i],".txt"), sep='\t', row.names = FALSE)
#        write.table(counts, paste0("/SAN/Plasmo_compare/SRAdb/Output/",study,"/count_",run[i],".txt"), sep='\t', row.names = FALSE, col.names = FALSE)
#    }
#  }
#}
#)

write.table(counts_df, file=paste0("/SAN/Plasmo_compare/SRAdb/Output/", study,"/countWithGFF3_",run,".txt"), sep='\t', row.names = FALSE)

#file.copy(from=paste0("/SAN/Plasmo_compare/SRAdb/countWithGFF3_",run,".txt",collapse=''), to=paste0("/SAN/Plasmo_compare/SRAdb/Output/",study,"/countWithGFF3_",run,".txt",collapse=''))
#file.copy(from=paste0("/SAN/Plasmo_compare/SRAdb/count_",run,".txt",collapse=''), to=paste0("/SAN/Plasmo_compare/SRAdb/Output/",study,"/count_",run,".txt",collapse=''))

# for(i in 1:length(bamList))
# {
#   print(i)
#   bamFile_1 <- import(bamList[i])
#   # levels(factor(seqnames(bamFile)))
#   # Pf_bam <- bamFile[grep("Pf", seqnames(bamFile)),]
#   counts <- countOverlaps(exon_hp_GR, bamFile_1)
#   counts_df <- data.frame(exon_hp_GR, counts)
#   
#   # sum of all counts in host genes and in parasite genes, then divide by the total
#   # host_percent = (total counts in host * 100) / (total counts in host + total counts in parasite )
#   parasite_rows <- grep(substr(para,1,3), counts_df[,1])
#   all_parasite <- counts_df[parasite_rows,]
#   all_host <- counts_df[-parasite_rows,]
#   parasite_count_percent <- (sum(all_parasite[,6])*100)/ (sum(all_parasite[,6]) + sum(all_host[,6]))
#   host_count_percent <- (sum(all_host[,6])*100) / (sum(all_parasite[,6]) + sum(all_host[,6]))
#   
#   study_map_percent[i,1] <- bamList[i]
#   study_map_percent[i,2] <- host_count_percent
#   study_map_percent[i,3] <- parasite_count_percent
#   write.table(counts_df, paste0("count_df_",bamList[i],".txt"), sep='\t', row.names = FALSE)
# }

# colnames(study_map_percent) <- c("Run", "Host_percent", "Parasite_percent")
# write.table(study_map_percent, "study_map_percent.txt", sep='\t', row.names = FALSE)


# setwd("..")
