#suppressMessages(install.packages("ggpubr"))
suppressMessages(library(ggpubr))
#suppressMessages(install.packages("reshape2"))
suppressMessages(library(reshape2))
#suppressMessages(install.packages("ggplot2"))
suppressMessages(library(ggplot2))
#suppressMessages(install.packages("plyr"))
suppressMessages(library(plyr))
#suppressMessages(install.packages("gridExtra"))
suppressMessages(library(gridExtra))
# 
# options(echo=TRUE)
# args <- commandArgs(TRUE)
# host <- args[1]
# para <- args[2]
# study <- args[3]



# check for runs in the study that already have count_run.txt files


HostParasitePercent <- function(runs)
{
  content <- data.frame()
  for( m in 1:length(runs))
  {
    print(paste0(m,"/",length(runs)))
    count_df <- read.csv2(paste0("/SAN/Plasmo_compare/SRAdb/Output/",study[i],"/countWithGFF3_",runs[m],".txt",collapse=''), sep='\t', header=TRUE)
    
    parasite_rows <- grep(paste0(substr(para,1,3),"_chr",collapse=''), count_df[,1])
    all_parasite <- count_df[parasite_rows,]
    all_host <- count_df[-parasite_rows,]
    parasite_count_percent <- (sum(all_parasite[,6])*100)/ (sum(all_parasite[,6]) + sum(all_host[,6]))
    host_count_percent <- (sum(all_host[,6])*100) / (sum(all_parasite[,6]) + sum(all_host[,6]))
    
    
    content[m,1] <- runs[m]
    content[m,2] <- host_count_percent
    content[m,3] <- parasite_count_percent
    
    colnames(content) <- c("Run", "Host_percent", "Parasite_percent")
  }
  return(content)
}

positive_experiments <- read.delim("/SAN/Plasmo_compare/SRAdb/Input/positive_experiments.txt", header = F, sep = ',')
study <- c("SRP250329", "SRP171171", "ERP109432", "SRP261098")


for(i in 1:length(study))
{
  print(study[i])
  host <- positive_experiments[grep(study[i],positive_experiments[,1]),2]
  para <- positive_experiments[grep(study[i],positive_experiments[,1]),3]
  
  current_runs <- read.csv2(paste0("/SAN/Plasmo_compare/SRAdb/Output/",study[i],"/runs_",study[i],".txt", collapse=''), header = FALSE, sep = ',')
  runs_in_study <- as.character(current_runs[grep(study[i], current_runs[,2]),1])
  
  runs_counted_number <- c()
  runs_counted_ID <- c()
  for(j in 1:length(runs_in_study))
  { 
    if(file.exists(paste0("/SAN/Plasmo_compare/SRAdb/Output/",study[i],"/countWithGFF3_",runs_in_study[j], ".txt", collapse = '')))
      runs_counted_ID <- c(runs_counted_ID, as.character(runs_in_study[j]))

  }

  runs_counted_number <- length(runs_counted_ID)

  hp_percent <- data.frame()
  hp_percent <- HostParasitePercent(runs_counted_ID)

  if(nrow(hp_percent) > 10)
  {
    l <- ceiling(runs_counted_number/10)
    Label <- c()
    a=0
    for (k in 1:l)
    {
      a=a+1
      Label <- c(Label,rep(k,10))
      a=a+9
    }
    Label <- as.character(Label[1:nrow(hp_percent)])
    hp_percent[,4] <- Label
  } else
  {
    l <- runs_counted_number
    Label <- rep(1, nrow(hp_percent))
    hp_percent[,4] <- Label
  }
  
  colnames(hp_percent) <- c("Run", "Host_percent", "Parasite_percent", "Label")
  hp_percent <- na.omit(hp_percent)
  
  ## add map%
  runs <- gsub("^\\_","",runs)
  distribution <- hp_percent[,c(1:3)]
  for(m in 1:nrow(distribution))
  {
    runID <- distribution[m,1]
    if(grep(pattern = paste0(runs[m],"Log.final.out",collapse=''), list.files()))
    {
      file <- list.files()[grep(pattern = paste0(runs[m],"Log.final.out",collapse=''), list.files())]
      line <- readLines(paste0(runs[m],"Log.final.out",collapse=''), n=10)[10]

      map_percent <- as.numeric(strsplit(strsplit(line, split="%")[[1]][2], split="\t")[[1]][2])
      distribution[m,4] <- map_percent
    }else
      distribution[m,4] <- "NA"
  }
  colnames(distribution)[4] <- c("Map_percent")
  
  ## add host and para percent taking into account map %
  
  for(n in 1:nrow(distribution))
  {
    host_overall <- (as.numeric(distribution[n,2])/100)*(as.numeric(distribution[n,4]))
    para_overall <- (as.numeric(distribution[n,3])/100)*(as.numeric(distribution[n,4]))

    distribution[n,5] <- host_overall
    distribution[n,6] <- para_overall
  }
  colnames(distribution)[5] <- c("Host_overall")
  colnames(distribution)[6] <- c("Parasite_overall")
  
  # examining hp difference brackets
  hp_difference <- abs(hp_percent[,2] - hp_percent[,3])
  
  a=0;b=0;c=0;d=0
  
  a=length(which(hp_difference <= 40))
  b=length(which(hp_difference > 40 & hp_difference <= 80))
  c=length(which(hp_difference > 80 & hp_difference <= 98))
  d=length(which(hp_difference > 98 & hp_difference <= 100))
  #
  table_output <- data.frame()
  table_output[1,1] <- study[i]
  table_output[1,2] <- paste0(nrow(hp_percent), "/", length(runs_in_study), collapse='')
  table_output[1,3] <- a
  table_output[1,4] <- b
  table_output[1,5] <- c
  table_output[1,6] <- d
  #
  colnames(table_output) <- c("Study", "Runs_finished/Runs_total", "At most 70%-30%", "At most 90%-10%", "At most 99%-1%", "Worse than 99%-1%")
  table_output <- t(table_output)
  

  h<-as.data.frame(hp_percent[,2])
  colnames(h) <- "Percentage"
  histogram_h_ggplot <- ggplot(h, aes(Percentage))+geom_histogram(binwidth = 0.5)+ scale_x_continuous(breaks=seq(0,100,10), labels=seq(0,100,10), position="bottom")+ labs(x="Host gene expression", y = "Count", title = "Histogram: Host gene expression") + theme_classic() + scale_y_continuous(breaks=seq(0,runs_counted_number,ceiling(runs_counted_number/10)), labels=seq(0,runs_counted_number,ceiling(runs_counted_number/10)),position="left") + theme(plot.margin = unit(c(1,1,1,1), "cm")) + theme(axis.title.y=element_text(size=12)) + theme(axis.title.x = element_text(size=12)) +  theme(plot.title=element_text(size=15)) + theme(axis.text.x = element_text(face="bold",size=10),axis.text.y = element_text(face="bold",size=10)) + coord_cartesian(xlim=c(0,100), ylim=c(0,nrow(hp_percent)))
  #
  ggsave(paste0("/SAN/Plasmo_compare/SRAdb/Output/",study,"/Histogram_host_",study,".png",collapse=''))
  #
  #
  p<- as.data.frame(hp_percent[,3])
  colnames(p) <- "Percentage"
  histogram_p_ggplot <- ggplot(p, aes(Percentage))+geom_histogram(binwidth = 0.5)+ scale_x_continuous(breaks=seq(0,100,10), labels=seq(0,100,10), position="bottom")+ labs(x="Parasite gene expression", y = "Count", title = "Histogram: Parasite gene expression") + theme_classic() + scale_y_continuous(breaks=seq(0,runs_counted_number,ceiling(runs_counted_number/10)), labels=seq(0,runs_counted_number,ceiling(runs_counted_number/10)),position="left") + theme(plot.margin = unit(c(1,1,1,1), "cm")) + theme(axis.title.y=element_text(size=12)) + theme(axis.title.x = element_text(size=12)) +  theme(plot.title=element_text(size=15)) + theme(axis.text.x = element_text(face="bold",size=10),axis.text.y = element_text(face="bold",size=10)) + coord_cartesian(xlim=c(0,100), ylim=c(0,nrow(hp_percent)))
  #
  ggsave(paste0("/SAN/Plasmo_compare/SRAdb/Output/",study,"/Histogram_parasite_",study,".png",collapse=''))
  #
  
  require(reshape2)
  
  hp_percent_melt <- melt(hp_percent, id=c("Run", "Label"))
  col=c("Host_percent"="salmon", "Parasite_percent"="skyblue")

  hp_bar_gg <- ggplot(hp_percent_melt, aes(x = Run, y=value, fill=variable))

  p =  hp_bar_gg + geom_bar(stat='identity', width=0.3) + scale_fill_manual(values=col) + xlab("Runs in the study") + ylab("Host% and Parasite% distribution") + ggtitle("Distribution of host% and parasite%")  + theme_classic() + theme(legend.title=element_blank(), legend.key.size = unit(0.2, "in"), legend.text = element_text(size=4),legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.2, "cm")) + theme(strip.text.x = element_text(margin = margin(.05, 0, .05, 0, "cm"))) + theme(plot.margin = unit(c(1,0.35,1,0.35), "cm"))  + theme(axis.title.y=element_text(size=6)) + theme(axis.title.x = element_text(size=6))  +  theme(plot.title=element_text(size=7))  + theme(axis.text.x = element_text(size=3),axis.text.y = element_text(size=5))  + facet_wrap(~Label,scale="free_x", nrow=ceiling(l^0.5))
  
  plots = dlply(hp_percent_melt , "Label", `%+%`, e1 = p)
  ml = do.call(marrangeGrob, list(plots, nrow=2, ncol=2))
  ggsave(paste0("/SAN/Plasmo_compare/SRAdb/Output/",study,"/Distribution_",study,".pdf", collapse=''), ml)
  
  write.table(table_output, paste0("/SAN/Plasmo_compare/SRAdb/Output/",study,"/table_",study,".txt", collapse=''), row.names=TRUE, sep="\t")
  write.table(distribution, paste0("/SAN/Plasmo_compare/SRAdb/Output/",study,"/hp_percent_",study,".txt", collapse=''), row.names=FALSE, sep="\t")
  
  
}
