### Script to visualise highly represented edges starting from the most represented and following that
## This will be a nested list data structure
## First get the edge with the highest representation
## the input data frame has "hp" as the first column, then all the studies, and the last column is edge_rowsums

edge_counter <- readRDS("liver_edge_counter.rds")

max.edge <- which.max(edge_counter$edge_rowsums)

# Level 0 of the nested list is the data frame of hp interactions with max edge score
pre_df0 <- edge_counter[max.edge,c("hp", "edge_rowsums")]
# separate host and parasite genes
hp <- strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\1 \\2', pre_df0$hp), ' ')
h <- sapply(hp, function(x) x[[1]])
p <- sapply(hp, function(x) x[[2]])

df0 <- data.frame(h = h, p = p, weight = pre_df0$edge_rowsums)

# separate host and parasite genes in edge_counter
ec_hp <- strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\1 \\2', edge_counter$hp), ' ')
h <- sapply(ec_hp, function(x) x[[1]])
p <- sapply(ec_hp, function(x) x[[2]])

edge_counter <- data.frame(h = h, p = p, weight = edge_counter$edge_rowsums)

## initiate list
rooted_network <- list(list())
rooted_network[[1]] <- df0

rooted_edgelist <- function(df)
{
  host_interactors <- c()
  para_interactors <- c()
  edges_host <- data.frame()
  edges_para <- data.frame()
  
  host_unique <- unique(df$h)
  para_unique <- unique(df$p)
  
  for(i in 1:length(host_unique))
  {
    print(paste0(i, "/", length(host_unique), collapse = ''))
    # for all host genes in the df data frame, find the interactors and same for all parasite genes 
    host_gene <- host_unique[i]
    
    # find interactors
    para_interactors <- c(para_interactors, edge_counter[grep(pattern = host_gene, edge_counter$h), "p"])
    
    host_side_interactions <- data.frame(h = rep(host_gene, length(para_interactors)),
                                         p = para_interactors)
     
  
    
    edges_host <- rbind(edges_host, host_side_interactions)
  }
  
  for(j in 1:length(para_unique))
  {
    print(paste0(j, "/", length(para_unique), collapse = ''))
    # for all host genes in the df data frame, find the interactors and same for all parasite genes 
    para_gene <- para_unique[j]
    
    # find interactors
    host_interactors <- c(host_interactors, edge_counter[grep(pattern = para_gene, edge_counter$p), "h"])
    
    para_side_interactions <- data.frame(h = host_interactors,
                                         p = rep(para_gene, length(host_interactors)))
    
    
    
    edges_para <- rbind(edges_para, para_side_interactions)
  }
  
  host_para <- rbind(edges_host, edges_para)
  
  df1 <- merge(host_para, edge_counter, by = c("h", "p"))
  return(df1)
}

## call the function to make the edgelist for visualisation

# I want to run it 5 times and see how big the list gets
edgelist <- df0
df2 <- data.frame()
for(k in 1:5)
{
  df2 <- rbind(df2, edgelist)
  edgelist <- rooted_edgelist(edgelist)
}

E(ig)$weight <- df1$weight
length(unique(df1$h))
#[1] 1869
length(unique(df1$p))
#[1] 3036
V(ig)$color = c(rep("grey", length(unique(df1$h))), rep("tomato", length(unique(df1$p))))

E(ig)$color[E(ig)$weight == 1] <- 'yellow'
E(ig)$color[E(ig)$weight == 2] <- 'salmon'
E(ig)$color[E(ig)$weight == 3] <- 'darkmagenta'
E(ig)$color[E(ig)$weight == 4] <- 'skyblue'
E(ig)$color[E(ig)$weight == 5] <- 'green'
E(ig)$color[E(ig)$weight == 6] <- 'navy'

png("liver_rooted_network.png" )
plot(ig, vertex.color = V(ig)$color,
     edge.color = E(ig)$color,
     vertex.size = 2, vertex.label=NA,
     edge.width = 2,
     edge.curved=F)
dev.off()


###################### second type of algorithm ##################

# read edge list and wights as data frame
edge_counter <- readRDS("liver_edge_counter.rds")
hp <- strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\1 \\2', edge_counter$hp), ' ')
h <- sapply(hp, function(x) x[[1]])
p <- sapply(hp, function(x) x[[2]])

edge_counter <- data.frame(h = h, p = p, weight = edge_counter$edge_rowsums)

# convert the edge list data frame into a matrix of weights.
library(reshape2)
edge_matrix <- dcast(edge_counter, h~p, value.var="weight")
edge_matrix[is.na(edge_matrix)] <- 0
rownames(edge_matrix) <- edge_matrix[,1]; edge_matrix <- edge_matrix[,2:ncol(edge_matrix)]

max.cell <- max(edge_matrix)
m1 = as.matrix(edge_matrix)
ind <- which(m1==max.cell, arr.ind=TRUE)
host <- sapply(ind, function(x) rownames(edge_matrix[x[,1],]))
para <- sapply(ind, function(x) colnames(edge_matrix[x[,2]]))

df0 = data.frame(h = host, p = para)
df0 <- merge(edge_counter, df0, by = c("h", "p"))

rooted_edgelist <- function(df)
{
  host <- unique(df$h)
  para <- unique(df$p)
  
  interactors_of_host <- c()
  interactors_of_para <- c()
  host_df<- data.frame()
  para_df<- data.frame()
  
  #if(length(host) <= 100)
  #  n = length(host)
  #if(length(host) > 100)
  #  n = 100
  
  for(i in 1:n)
  {
    print(paste0(i, "/", length(host), collapse = ''))
    interactors_h <- c()
    interactors_h <- as.vector(edge_matrix[which(rownames(edge_matrix) == host[i]), ])
    interactors_h <- interactors_h[which(interactors_h != 0)]
    
    interactors_of_host <- c(interactors_of_host, unlist(names(interactors_h)))
    
    host_df <- rbind(host_df, data.frame(h = rep(host[i], length(interactors_of_host)), p = interactors_of_host))
  }
  
  if(length(para) <= 100)
    m = length(para)
  if(length(para) > 100)
    m = 100
  
  for(j in 1:m)
  {
    print(paste0(j, "/", length(para), collapse = ''))
    
    interactors_p <- c()
    interactors_p <- rownames(edge_matrix[which(edge_matrix[which(colnames(edge_matrix) == para[j])] != 0),])

    interactors_of_para <- c(interactors_of_para, interactors_p)
    
    para_df <- rbind(para_df, data.frame(h = interactors_of_para, p = rep(para[j], length(interactors_of_para))))
  }
  
  host_para <- rbind(host_df, para_df)
  host_para<-merge(edge_counter, host_para, by = c("h", "p"))
  return(host_para)
}

edgelist <- df0
df2 <- data.frame()
for(k in 1:5)
{
  df2 <- rbind(df2, edgelist)%>%
    distinct()
  edgelist <- rooted_edgelist(edgelist) %>%
    distinct()
}
