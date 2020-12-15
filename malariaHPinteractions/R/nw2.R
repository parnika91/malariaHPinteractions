library(igraph)
# Assume we examine (fictive) train connections of 4 countries: Switzerland, Italy, France, Spain
# in the Swiss cities "Genf" and "Lugano" there are different languages/ethnicities

#construct the graph
g <- graph (c( "Zurich","Bern","Zurich","Bern", "Genf","Bern","Lugano","Zurich",
               "Genf","Zurich","Lugano","Bern",
               "Rome","Venice","Rome","Milano","Venice","Milano",
               "Marseille","Lyon","Marseille","Toulouse","Lyon","Toulouse",
               "Barcelona","Saragosa","Barcelona","Valencia","Saragosa","Valencia",
               "Milano","Lugano","Genf","Lyon","Milano","Marseille","Marseille","Barcelona"
))

#set major language/ethnicities
V(g)$etnic <- c("Swiss", "Swiss","French","Italian",  #for Genf and Lugano respectively!
                "Italian","Italian","Italian",
                "French","French","French",
                "Spanish","Spanish","Spanish")

V(g)$color <- ifelse(V(g)$etnic == "Italian", "#61D04F", ifelse(V(g)$etnic =="French", "#2297E6", ifelse(V(g)$etnic == "Spanish","#F5C710","red")))

#when we simply plot this graph, everything looks good
plot(g, vertex.label.color="black", vertex.label.dist=1.8, edge.arrow.size=.5,
     vertex.color = V(g)$color) 

# now let's see, whether the clustering finds the four countries
clust <- edge.betweenness.community(g)

#but when we plot this, the clustered graph loses the color of the vertices
plot(clust, g, edge.arrow.size=.15, edge.curved=0, vertex.frame.color="black",
     vertex.label=V(g)$city, vertex.label.color="black",
     vertex.label.cex=.8, layout=layout_with_dh(g))

# use the mark.groups argument
plot(g, mark.groups=communities(clust),  
     edge.arrow.size=.15, edge.curved=0, vertex.frame.color="black",
     vertex.label=V(g)$city, vertex.label.color="black",
     vertex.label.cex=.8, layout=layout_with_dh(g))
# also check out the other arguments for the grouping:
# mark.shape, mark.border, mark.col and mark.expand
net <- read.csv("blood_core_edges_human_Pb.csv", stringsAsFactors = F)
features <- read.csv("blood_core_fast_greedy_40comm.csv", stringsAsFactors = F)

ig <- graph_from_data_frame(net, directed = F, vertices = features)
pdf("blood_core_fast_greedy_40comm.pdf")
plot(ig, mark.groups=V(ig)$Membership, vertex.size = 2,
  edge.arrow.size=.15, edge.curved=0.2, 
  vertex.frame.color="black",vertex.label=NA, 
  vertex.label.color="black",vertex.label.cex=.8,
  layout=layout_with_graphopt(ig)
   )
dev.off()

clust <- cluster_fast_greedy(ig, merges = TRUE, modularity = TRUE,
  membership = TRUE, weights = abs(E(ig)$weight))

pdf("blood_core_fast_greedy_40comm.pdf")
plot(clust, ig)
dev.off()