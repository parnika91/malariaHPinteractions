library(igraph)
library(grid)

study <- "ERP106451"
load(paste0(study,"_bipartite.RData", collapse =''))
df1 <- ERP106451_bipartite[,c(1,2,3)]
vnames = c(unique(as.character(df1$gene1)), unique(as.character(df1$gene2)))
hnum = length(unique(as.character(df1$gene1)))
pnum = length(unique(as.character(df1$gene2)))

net <- graph_from_data_frame(d = df1, vertices=vnames, directed=F)
V(net)$color <- c(rep("gray50", hnum), rep("tomato", pnum))

png(paste0(study,"_bipartite.png")) # , width = 10, height = 10, unit = "cm", res = 300
plot(net, edge.arrow.size=.4,vertex.label=NA, vertex.size=4, vertex.frame.width=0.2, layout = layout_nicely(net))
dev.off()
