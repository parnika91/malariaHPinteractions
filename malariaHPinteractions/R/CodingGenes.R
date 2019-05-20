# extract coding genes
# biomaRt for the hosts and gtf files for the parasites

################ hosts ###############

# install BiomaRt
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt", version = "3.8")
library("biomaRt")

# have a look at all the marts availablr. We want the protists mart
listMarts()

ensembl = useEnsembl(biomart="ensembl")
human_ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
human_gb <- getBM(attributes=c("ensembl_gene_id","external_gene_name",
                         "description"), 
            filters='biotype', 
            values=c('protein_coding'), 
            mart=human_ensembl)
write.table(human_gb, "human_protein_coding_genes.txt", sep = '\t', row.names = F)


mouse_ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
mouse_gb <- getBM(attributes=c("ensembl_gene_id","external_gene_name",
                               "description"), 
                  filters='biotype', 
                  values=c('protein_coding'), 
                  mart=mouse_ensembl)
write.table(mouse_gb, "mouse_protein_coding_genes.txt", sep = '\t', row.names = F)


monkey_ensembl = useDataset("mmulatta_gene_ensembl",mart=ensembl)
monkey_gb <- getBM(attributes=c("ensembl_gene_id","external_gene_name",
                               "description"), 
                  filters='biotype', 
                  values=c('protein_coding'), 
                  mart=monkey_ensembl)
write.table(monkey_gb, "monkey_protein_coding_genes.txt", sep = '\t', row.names = F)


################# parasites ################
library(rtracklayer)

Pfal_coding <- as.data.frame(import("/SAN/Plasmo_compare/Genomes/annotation/Pfalciparum.gtf")) %>%
  filter(type%in%"exon") %>%
  filter(gene_biotype%in%"protein_coding") %>%
  distinct(gene_id)

Pber_coding <- as.data.frame(import("/SAN/Plasmo_compare/Genomes/annotation/Pberghei.gtf")) %>%
  filter(type%in%"exon") %>%
  filter(gene_biotype%in%"protein_coding") %>%
  distinct(gene_id)

Pviv_coding <- as.data.frame(import("/SAN/Plasmo_compare/Genomes/annotation/Pvivax.gtf")) %>%
  filter(type%in%"exon") %>%
  filter(gene_biotype%in%"protein_coding") %>%
  distinct(gene_id)

Pyoe_coding <- as.data.frame(import("/SAN/Plasmo_compare/Genomes/annotation/Pyoelii.gtf")) %>%
  filter(type%in%"exon") %>%
  filter(gene_biotype%in%"protein_coding") %>%
  distinct(gene_id)

Pcha_coding <- as.data.frame(import("/SAN/Plasmo_compare/Genomes/annotation/Pchabaudi.gtf")) %>%
  filter(type%in%"exon") %>%
  filter(gene_biotype%in%"protein_coding") %>%
  distinct(gene_id)

Pcyn_coding <- as.data.frame(import("/SAN/Plasmo_compare/Genomes/annotation/Pcynomolgi.gtf")) %>%
  filter(type%in%"exon") %>%
  filter(gene_biotype%in%"protein_coding") %>%
  distinct(gene_id)

Pcoa_coding <- as.data.frame(import("/SAN/Plasmo_compare/Genomes/annotation/Pcoatneyi.gtf")) %>%
  filter(type%in%"exon") %>%
  filter(gene_biotype%in%"protein_coding") %>%
  distinct(gene_id)

################ host-parasite pairs ###############

library(rtracklayer)

hPfal_coding <- as.data.frame(import("/SAN/Plasmo_compare/Genomes/annotation/humanPfalciparum.gtf")) %>%
  filter(type%in%"exon") %>%
  filter(gene_biotype%in%"protein_coding") %>%
  distinct(gene_id)

hPber_coding <- as.data.frame(import("/SAN/Plasmo_compare/Genomes/annotation/humanPberghei.gtf")) %>%
  filter(type%in%"exon") %>%
  filter(gene_biotype%in%"protein_coding") %>%
  distinct(gene_id)

hPviv_coding <- as.data.frame(import("/SAN/Plasmo_compare/Genomes/annotation/humanPvivax.gtf")) %>%
  filter(type%in%"exon") %>%
  filter(gene_biotype%in%"protein_coding") %>%
  distinct(gene_id)

mPyoe_coding <- as.data.frame(import("/SAN/Plasmo_compare/Genomes/annotation/mousePyoelii.gtf")) %>%
  filter(type%in%"exon") %>%
  filter(gene_biotype%in%"protein_coding") %>%
  distinct(gene_id)

mPcha_coding <- as.data.frame(import("/SAN/Plasmo_compare/Genomes/annotation/mousePchabaudi.gtf")) %>%
  filter(type%in%"exon") %>%
  filter(gene_biotype%in%"protein_coding") %>%
  distinct(gene_id)

mPber_coding <- as.data.frame(import("/SAN/Plasmo_compare/Genomes/annotation/mousePberghei.gtf")) %>%
  filter(type%in%"exon") %>%
  filter(gene_biotype%in%"protein_coding") %>%
  distinct(gene_id)

moPcoa_coding <- as.data.frame(import("/SAN/Plasmo_compare/Genomes/annotation/monkeyPcoatneyi.gtf")) %>%
  filter(type%in%"exon") %>%
  filter(gene_biotype%in%"protein_coding") %>%
  distinct(gene_id)

moPcyn_coding <- as.data.frame(import("/SAN/Plasmo_compare/Genomes/annotation/monkeyPcynomolgi.gtf")) %>%
  filter(type%in%"exon") %>%
  filter(gene_biotype%in%"protein_coding") %>%
  distinct(gene_id)

########## do it for every study #########

studyIDs <- as.character(read.table("/SAN/Plasmo_compare/SRAdb/Input/positive_experiments.txt", sep = '\t', header = F)[,1])
allHPexp <- read.csv2("/SAN/Plasmo_compare/SRAdb/allHPexp.txt", sep = ',', header = T)

for(i in 1:length(studyIDs))
{
  print(i)
  # get study.txt including all runs
  study <- as.data.frame(t(read.delim(paste0("Output/", studyIDs[i], "/", studyIDs[i], ".txt", collapse = ''), sep = '\t', header = T)))
  # get host parasite from allHPexp
  hp <- as.character(unique(allHPexp[allHPexp$Study==studyIDs[i],"HostParasite"]))
  
  if(hp == "humanPfalciparum"){ coding = hPfal_coding }
  if(hp == "humanPberghei"){ coding = hPber_coding }
  if(hp == "humanPvivax"){ coding = hPviv_coding }
  if(hp == "mousePberghei"){ coding = mPber_coding }
  if(hp == "mousePyoelii"){ coding = mPyoe_coding }
  if(hp == "mousePchabaudi"){ coding = mPcha_coding }
  if(hp == "monkeyPcoatneyi"){ coding = moPcoa_coding }
  if(hp == "monkeyPcynomolgi"){ coding = moPcyn_coding }
  
  # filter study to keep only protein-coding genes
  study_coding_genes <- study %>%
    tibble::rownames_to_column('gene') %>%
    filter(rownames(study)%in%coding$gene_id) %>%
    tibble::column_to_rownames('gene')
   
  # write the table out
  write.table(study_coding_genes, paste0("/SAN/Plasmo_compare/SRAdb/Output/", studyIDs[i],"/", studyIDs[i], "_coding_genes.txt", collapse = ''), sep ='\t', row.names = T)
}