
#library(vegan)
library(dplyr)
library(ggplot2)
library(ggpubr)

####### get allHPexp and ortho_data #######
allHPexp <- read.delim("Output/allHPexp.txt", sep = ',')
ortho_data <- as.data.frame(read.delim("Output/ortho_data.txt")) #%>%
 # tibble::column_to_rownames("Orthogroup")

####### get runs from category #######
# first category: all studies of one host #
# human
#human.runs <- as.character(allHPexp[which(allHPexp$Host=="human"), "RunID"])
#human <- ortho_data[,sapply(human.runs, function(x) grep(pattern = x, colnames(ortho_data)))]
#mouse.runs <- as.character(allHPexp[which(allHPexp$Host=="mouse"), "RunID"])
#monkey.runs <- as.character(allHPexp[which(allHPexp$Host=="monkey"), "RunID"
incl.runs <- allHPexp[which(allHPexp$ProteinCodHost >=1e6 & allHPexp$ProteinCodPara >= 1e5 & allHPexp$NumberProtCodGenesHost >= 10000 & allHPexp$NumberProtCodGenesPara >= 3000 & allHPexp$MapPercent >= 70),]
blood.runs <- as.character(incl.runs[which(incl.runs$Tissue=="blood"), "RunID"])
col.num <- c()
for(i in 1:length(blood.runs))
{
  id <- grep(pattern = blood.runs[i], colnames(ortho_data))
  if(length(id)==1)
    col.num[i] <- id
}

blood <- ortho_data[,col.num]

# monkey
# mo.col.num <- c()
# for(i in 1:length(monkey.runs))
# {
#   id <- grep(pattern = monkey.runs[i], colnames(ortho_data))
#   if(length(id)==1)
#     mo.col.num[i] <- id
# }
# monkey <- ortho_data[,mo.col.num]
# 
# # mouse
# m.col.num <- c()
# for(i in 1:length(mouse.runs))
# {
#   id <- grep(pattern = mouse.runs[i], colnames(ortho_data))
#   if(length(id)==1)
#     m.col.num[i] <- id
# }
# mouse <- ortho_data[,m.col.num]

## nMDS ##
# human_NMDS <- metaMDS(human[1:1000,], k = 2)
# stressplot(human_NMDS)
# 
# plot(human_NMDS)
# ordiplot(human_NMDS,type="n")
# orditorp(human_NMDS,display="species",col="red",air=0.01)
# orditorp(human_NMDS,display="sites",cex=1.25,air=0.01)
# 
# treat=c(rep("Treatment1",5),rep("Treatment2",5))
# ordiplot(example_NMDS,type="n")
# ordihull(example_NMDS,groups=treat,draw="polygon",col="grey90",label=F)
# orditorp(example_NMDS,display="species",col="red",air=0.01)
# orditorp(example_NMDS,display="sites",col=c(rep("green",5),rep("blue",5)),
#          air=0.01,cex=1.25)

library(edgeR)
#install.packages("circlize")
library(circlize)
run <- sapply(colnames(blood), function(x) strsplit(x, "_")[[1]][1])
#human <- human[1:500,560:760]
#y <- DGEList(counts = human)
#colour.vector <- c()
#cols <- rand_color(55, hue = NULL, luminosity = "random", transparency = 0)
#studies <- unique(as.character(allHPexp$Study))

#for(i in 1:length(cols))
#{
#  for(j in 1:ncol(human))
#  {
#    study <- strsplit(colnames(human)[j], "_")[[1]][2]
#    if(study == studies[i])
#      colour.vector[j] <- cols[i]
#  }
#}

#y$samples$group <- colour.vector
#y.tmm <- calcNormFactors(y, na.rm = T)
#y.cpm <- cpm(y.tmm)
#pdf("MDS_human.pdf", width = 30, height = 30)
#plotMDS(y.tmm, col=colour.vector, pch = 1, cex= 0.95, main = "Human studies MDS", cex.lab = 0.75)
#dev.off()

# mouse
#ortho_data <- ortho_data[1:500,600:800]
y <- DGEList(counts = blood, group = run)
colour.vector <- c()
cols <- rand_color(55, hue = NULL, luminosity = "random", transparency = 0)
#studies <- unique(as.character(allHPexp$Study))
study.vector <- c()
host <- unique(as.character(allHPexp$Host))
(host)
host.vector <- c()

#for(i in 1:length(cols))
#{
#  for(j in 1:ncol(blood))        
#  {
#    study <- strsplit(colnames(blood)[j], "_")[[1]][2]       
#    if(study == studies[i])
#      colour.vector[j] <- cols[i]
#  }
#}

#host.cols <- rand_color(length(host), hue = NULL, luminosity = "random", transparency = 0)
#(host.cols)
for(k in 1:ncol(blood))
{
   study <- strsplit(colnames(blood)[k], "_")[[1]][2]
   h <- unique(as.character(allHPexp[which(allHPexp$Study==study),"Host"]))
   if(h == "human")
     host.vector[k] <- "human"
   if(h == "monkey")                              
     host.vector[k] <- "monkey"
   if(h == "mouse")                              
     host.vector[k] <- "mouse"
   study.vector[k] <- study
}




#y$samples$group <- colour.vector
#y.tmm <- calcNormFactors(y, na.rm = T)
#y.cpm <- cpm(y.tmm)
pdf("MDS_allstudies_wBCV_blood_anno_threshold_dots.pdf")
plotMDS(y, method = "bcv", col=host.vector,pch = 1, main = "All blood studies after threshold MDS")
dev.off()

# monkey

# mo.y <- DGEList(counts = monkey)
# colour.vector <- c()
# cols <- rand_color(55, hue = NULL, luminosity = "random", transparency = 0)
# studies <- unique(as.character(allHPexp$Study))            
# 
# for(i in 1:length(cols))
# {
#   for(j in 1:ncol(monkey))        
#   {
#     study <- strsplit(colnames(monkey)[j], "_")[[1]][2]       
#     if(study == studies[i])
#       colour.vector[j] <- cols[i]
#   }
# }       
# 
# mo.y$samples$group <- colour.vector
# mo.y.tmm <- calcNormFactors(mo.y, na.rm = T)
# mo.y.cpm <- cpm(mo.y.tmm)
# pdf("MDS_monkey.pdf", width = 30, height = 30)
# plotMDS(mo.y.tmm, col=colour.vector, pch = 1, cex= 0.95, main = "Monkey studies MDS", cex.lab = 1.2)
# dev.off()

res.cor <- cor(blood, method = "p")
mds.cor <- (1-res.cor) %>%
  cmdscale() %>%
  as_tibble() %>%
  mutate(group = host.vector) %>%
  mutate(study = study.vector)

colnames(mds.cor) <- c("Dim.1", "Dim.2", "Group", "Study")
ggscatter(mds.cor, x = "Dim.1", y = "Dim.2", 
          size = 1,
          col = "Group",
          shape = "Study",
          repel = TRUE,
          ) + scale_shape_manual(values=1:nlevels(mds.cor$Study))
ggsave("mds.blood.png")

ggplot() +
  geom_polygon(data=mds.cor,aes(x=Dim.1,y=Dim.2,fill=Group,group=Group),alpha=0.30) + # add the convex hulls
  #geom_text(data=mds.cor,aes(x=Dim.1,y=Dim.2,label=Study),alpha=0.5) +  # add the species labels
  geom_point(data=mds.cor,aes(x=Dim.1,y=Dim.2,shape=Study,colour=Group),size=1) + # add the point markers
  scale_shape_manual(values=1:nlevels(as.factor(mds.cor$Study))) +
  coord_equal() +
  theme_bw() + 
  theme(axis.title.x = element_text(size=10), # remove x-axis labels
        axis.title.y = element_text(size=10), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())
ggsave("mds.blood.png")

