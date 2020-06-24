library(dplyr)
library(igraph)
library(betareg)
library(ggeffects)
library(ggplot2)
library(reshape2)

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

study <- "ERP106451"
type <- "int" # str, int or all
nw <- "bp" # all, pp or bp

if(nw == "all")
  net = "all"
if(nw == "bp")
  net = "bipartite"
if(nw == "pp")
  net = "para"

if(grep(pattern = "df_concat", study))
  org = "ov"
if(grep(pattern = "ERP004598", study))
  org = "pb"
if(grep(pattern = "DRP000987", study))
  org = "pf"
if(grep(pattern = "ERP106451", study))
  org = "pfE"

data <- loadRData(paste0(study,"/","cor", "/", study, "_", type, "_", net, ".RData", collapse = ""))
ig <- graph_from_data_frame(data, directed = F)

dg <- degree(ig, v = V(ig))
bw <- betweenness(ig, v = V(ig), directed = FALSE)
cl <- closeness(ig, vids = V(ig))
ec <- eigen_centrality(ig, directed = FALSE)
kc <- coreness(ig, mode = c("all"))

dg_p <- dg[grep(pattern = "p_OG", names(dg))]
dg_df <- as.data.frame(dg_p) %>%
  tibble::rownames_to_column("Orthogroup")

bw_p <- bw[grep(pattern = "p_OG", names(bw))]
bw_df <- as.data.frame(bw_p) %>%
  tibble::rownames_to_column("Orthogroup")

cl_p <- cl[grep(pattern = "p_OG", names(cl))]
cl_df <- as.data.frame(cl_p) %>%
  tibble::rownames_to_column("Orthogroup")

ec_p <- ec$vector[grep(pattern = "p_OG", names(ec$vector))]
ec_df <- as.data.frame(ec_p) %>%
  tibble::rownames_to_column("Orthogroup")

kc_p <- kc[grep(pattern = "p_OG", names(kc))]
kc_df <- as.data.frame(kc_p) %>%
  tibble::rownames_to_column("Orthogroup")

join_df <-plyr::join_all(list(dg_df, bw_df, cl_df,
                              ec_df, kc_df), by = "Orthogroup", type = "full")
join_df[is.na(join_df)] <- 0

colnames(join_df) <- c(paste0(org, "_", net, "_dg", collapse = ""),
                       paste0(org, "_", net, "_bw", collapse = ""),
                       paste0(org, "_", net, "_cl", collapse = ""),
                       paste0(org, "_", net, "_ec", collapse = ""),
                       paste0(org, "_", net, "_kc", collapse = ""))

save(join_df, file = paste0(org, "_", net, "_nwstats.RData", collapse = ""))


load("Unweighted_RGR_MIS_cleaned.RData")
RGR_MIS_cleaned <- inner_join(join_df, RGR_MIS_cleaned)
save(RGR_MIS_cleaned, file = "Unwebp_ighted_RGR_MIS_cleaned.RData")