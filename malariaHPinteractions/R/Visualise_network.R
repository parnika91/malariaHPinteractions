library(tibble)
library(dplyr)
library(visNetwork)
library(networkD3)
SRP118996_bipartite <- read.delim("~/Documents/SRP118996_bipartite.txt", stringsAsFactors=FALSE)
nodes <- as_tibble(c(as.character(SRP118996_bipartite[,1])),
                   c(as.character(SRP118996_bipartite[,2])))# %>% 
  #rowid_to_column("id")
edges <- data.frame(from = SRP118996_bipartite[,1], 
                    to = SRP118996_bipartite[,2], 
                    weight = SRP118996_bipartite[,3])
nw <- network(edges, vertex.attr = nodes, matrix.type = "edgelist", ignore.eval = FALSE)
visNetwork(nodes, edges)


load("Documents/Data/rhoptry.RData")
nodes <- as_tibble(c(as.character(rhoptry[,2])),
                   c(as.character(rhoptry[,4]))) #%>% 
  #rowid_to_column("id")
edges <- data.frame(from = rhoptry[,4], 
                    to = rhoptry[,2], 
                    weight = rhoptry[,5])
nw <- network(edges, vertex.attr = nodes, matrix.type = "edgelist", ignore.eval = FALSE)
visNetwork(nodes, edges)
