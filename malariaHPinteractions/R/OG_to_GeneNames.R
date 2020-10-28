host <- sapply(liver_core_edged[,1], function(x) host_orthogroups[grep(pattern = x, host_orthogroups$Orthogroup),"mouse"])
para <- sapply(liver_core_edged[,2], function(x) pOG[grep(pattern = x, pOG$Orthogroup),"Pb_g"])
liver_core_edge_genenames <- data.frame(host = host, para = para, weight = liver_core_edged$weight)
liver_core_nodes <- c(as.character(unique(liver_core_edge_genenames$host)), as.character(unique(liver_core_edge_genenames$para)))
liver_core_nodes <- data.frame(name = liver_core_nodes)
liver_core_nodes$colour = rep(1, nrow(liver_core_nodes))
hrows = grep(pattern = "ENS", liver_core_nodes[,1])
liver_core_nodes$colour[hrows] <- 2
write.table(liver_core_edge_genenames, "liver_core_edge_genenames.csv", sep = '\t', quote = F, row.names = F)
write.table(liver_core_nodes, "liver_core_nodes.csv", sep = '\t', quote = F, row.names = F)


# blood

host <- sapply(blood_core_edges[,1], function(x) host_orthogroups[grep(pattern = x, host_orthogroups$Orthogroup),"mouse"])
para <- sapply(blood_core_edges[,2], function(x) pOG[grep(pattern = x, pOG$Orthogroup),"Pb_g"])
blood_core_edge_genenames <- data.frame(host = host, para = para, weight = blood_core_edges$weight)
blood_core_nodes <- c(as.character(unique(blood_core_edge_genenames$host)), as.character(unique(blood_core_edge_genenames$para)))
blood_core_nodes <- data.frame(name = blood_core_nodes)
blood_core_nodes$colour = rep(1, nrow(blood_core_nodes))
hrows = grep(pattern = "ENS", blood_core_nodes[,1])
blood_core_nodes$colour[hrows] <- 2
write.table(blood_core_edge_genenames, "blood_core_edge_genenames.csv", sep = '\t', quote = F, row.names = F)
write.table(blood_core_nodes, "blood_core_nodes.csv", sep = '\t', quote = F, row.names = F)
