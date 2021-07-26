library(dplyr)
library(igraph)
library(betareg)
library(ggeffects)
library(ggplot2)
library(reshape2)
library(stringr)

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

models_list <- list()

study <- "blood"
type <- "all" # str, int or all
nw <- "bp" # pp_bp, pp or bp

if(nw == "pp_bp")
  net = "pp_bp"
if(nw == "bp")
  net = "bipartite"
if(nw == "pp")
  net = "para"

if(grepl(pattern = "df_concat", study))
  org = "ov"
if(grepl(pattern = "ERP004598", study))
  org = "pfE"
if(grepl(pattern = "DRP000987", study))
  org = "pf"
if(grepl(pattern = "ERP106451", study))
  org = "pfE"
if(grepl(pattern = "hpv_bl", study))
  org = "pv"
if(grepl(pattern = "overall_addblood", study))
  org = "ov_ad"

if(net=="pp_bp")
  {
    data_pp <- loadRData(paste0(study,"/","cor", "/", study, "_", type, "_para.RData", collapse = "")); colnames(data_pp)[1] <- "h"; colnames(data_pp)[2] <- "p"
    data_bp <- loadRData(paste0(study,"/","cor", "/", study, "_", type, "_bipartite.RData", collapse = "")); colnames(data_bp)[1] <- "h"; colnames(data_bp)[2] <- "p"
    data <- rbind(data_pp[,c(1,2)], data_bp[,c(1,2)])
  } else {
    data <- loadRData(paste0(study,"/","cor", "/", study, "_", type, "_", net, ".RData", collapse = ""))
  }
data_<- loadRData("blood_all_para.RData")
d <- data.frame(h = as.character(data[,1]), p = as.character(data[,2]))
ig <- graph_from_data_frame(d, directed = F)

dg <- degree(ig, v = V(ig), loops = F, normalized = F)
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

colnames(join_df) <- c("Orthogroup",paste0(org, "_", nw, "_dg", collapse = ""),
                       paste0(org, "_", nw, "_bw", collapse = ""),
                       paste0(org, "_", nw, "_cl", collapse = ""),
                       paste0(org, "_", nw, "_ec", collapse = ""),
                       paste0(org, "_", nw, "_kc", collapse = ""))

save(join_df, file = paste0(org, "_", net, "_",nw,"_nwstats.RData", collapse = ""))


load("Unweighted_RGR_MIS_cleaned.RData")
RGR_MIS_cleaned <- left_join(RGR_MIS_cleaned, join_df)
save(RGR_MIS_cleaned, file = "Unweighted_RGR_MIS_cleaned.RData")

### models ###

eff <- c("dg_m", "ec_m", "bw_m", "dg_bw_m", "dg_ec_m", "m")

dg_rgr <- ec_rgr <- bw_rgr <- dg_bw_rgr <- dg_ec_rgr <- dg_bw_ec_rgr <- dg_mis <- ec_mis <- bw_mis <- dg_bw_mis <- dg_ec_mis <- dg_bw_ec_mis <- 0

dg_rgr <- betareg(data = RGR_MIS_cleaned, RGR ~ RGR_MIS_cleaned[,paste0(org, "_", nw, "_dg", collapse = "")])
ec_rgr <- betareg(data = RGR_MIS_cleaned, RGR ~ RGR_MIS_cleaned[,paste0(org, "_", nw, "_ec", collapse = "")])
bw_rgr <- betareg(data = RGR_MIS_cleaned, RGR ~ RGR_MIS_cleaned[,paste0(org, "_", nw, "_bw", collapse = "")])
dg_bw_rgr <- betareg(data = RGR_MIS_cleaned, RGR ~ RGR_MIS_cleaned[,paste0(org, "_", nw, "_dg", collapse = "")] + RGR_MIS_cleaned[,paste0(org, "_", nw, "_bw", collapse = "")])
dg_ec_rgr <- betareg(data = RGR_MIS_cleaned, RGR ~ RGR_MIS_cleaned[,paste0(org, "_", nw, "_dg", collapse = "")] + RGR_MIS_cleaned[,paste0(org, "_", nw, "_ec", collapse = "")])
dg_bw_ec_rgr <- betareg(data = RGR_MIS_cleaned, RGR ~ RGR_MIS_cleaned[,paste0(org, "_", nw, "_dg", collapse = "")] + RGR_MIS_cleaned[,paste0(org, "_", nw, "_bw", collapse = "")] +
    RGR_MIS_cleaned[,paste0(org, "_", net, "_ec", collapse = "")])

dg_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ RGR_MIS_cleaned[,paste0(org, "_", nw, "_dg", collapse = "")])
ec_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ RGR_MIS_cleaned[,paste0(org, "_", nw, "_ec", collapse = "")])
bw_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ RGR_MIS_cleaned[,paste0(org, "_", nw, "_bw", collapse = "")])
dg_bw_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ RGR_MIS_cleaned[,paste0(org, "_", nw, "_dg", collapse = "")] + RGR_MIS_cleaned[,paste0(org, "_", nw, "_bw", collapse = "")])
dg_ec_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ RGR_MIS_cleaned[,paste0(org, "_", nw, "_dg", collapse = "")] + RGR_MIS_cleaned[,paste0(org, "_", nw, "_ec", collapse = "")])
dg_bw_ec_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ RGR_MIS_cleaned[,paste0(org, "_", nw, "_dg", collapse = "")] + RGR_MIS_cleaned[,paste0(org, "_", nw, "_bw", collapse = "")] +
    RGR_MIS_cleaned[,paste0(org, "_", net, "_ec", collapse = "")])

models_list <- list(summary(dg_rgr), summary(ec_rgr), summary(bw_rgr), summary(dg_bw_rgr), summary(dg_ec_rgr), summary(dg_bw_ec_rgr),
  summary(dg_mis), summary(ec_mis), summary(bw_mis), summary(dg_bw_mis), summary(dg_ec_mis), summary(dg_bw_ec_mis))
for(i in 1:length(models_list))
{
  if(i<=6)
    names(models_list)[i] <- paste0(org, "_", nw, "_", eff[i], "_rgr")
  if(i>6)
    names(models_list)[i] <- paste0(org, "_", nw, "_", eff[i-6], "_mis")
}
models <- list()
s = 1
for(u in 1:length(models_list))
{
  if(length(models_list[u][[1]]) == 24)
    {models[s] <- models_list[u]
      names(models)[s] <- names(models_list[u])
      s = s+1}

}

### intersection of overall network and P. berghei network

ov_ad_pp <- loadRData("overall_addblood/cor/overall_addblood_all_para.RData"); colnames(ov_ad_pp)[1] <- "host"; colnames(ov_ad_pp)[2] <- "para"
ov_ad_bp <- loadRData("overall_addblood/cor/overall_addblood_all_bipartite.RData"); colnames(ov_ad_bp)[1] <- "host"; colnames(ov_ad_bp)[2] <- "para"
ov_ad_pp_bp <- rbind(ov_ad_pp[,c(1,2)], ov_ad_bp[,c(1,2)])

pb_pp <- loadRData("ERP004598/cor/ERP004598_all_para.RData"); colnames(pb_pp)[1] <- "host"; colnames(pb_pp)[2] <- "para"
pf_pp <- loadRData("DRP000987/cor/DRP000987_str_para.RData"); colnames(pf_pp)[1] <- "host"; colnames(pf_pp)[2] <- "para"
pfE_pp <- loadRData("ERP106451/cor/ERP106451_int_para.RData"); colnames(pfE_pp)[1] <- "host"; colnames(pfE_pp)[2] <- "para"
pv_pp <- loadRData("hpv_bl/cor/hpv_bl_int_para.RData"); colnames(pv_pp)[1] <- "host"; colnames(pv_pp)[2] <- "para"

pb_bp <- loadRData("ERP004598/cor/ERP004598_all_bipartite.RData"); colnames(pb_bp)[1] <- "host"; colnames(pb_bp)[2] <- "para"
pf_bp <- loadRData("DRP000987/cor/DRP000987_str_bipartite.RData"); colnames(pf_bp)[1] <- "host"; colnames(pf_bp)[2] <- "para"
pfE_bp <- loadRData("ERP106451/cor/ERP106451_int_bipartite.RData"); colnames(pfE_bp)[1] <- "host"; colnames(pfE_bp)[2] <- "para"
pv_bp <- loadRData("hpv_bl/cor/hpv_bl_int_bipartite.RData"); colnames(pv_bp)[1] <- "host"; colnames(pv_bp)[2] <- "para"

pb_pp_bp <- rbind(pb_pp[,c(1,2)], pb_bp[,c(1,2)])
pf_pp_bp <- rbind(pf_pp[,c(1,2)], pf_bp[,c(1,2)])
pfE_pp_bp <- rbind(pfE_pp[,c(1,2)], pfE_bp[,c(1,2)])
pv_pp_bp <- rbind(pv_pp[,c(1,2)], pv_bp[,c(1,2)])

ov_ad_pb_pp_merge <- merge(ov_ad_pp, pb_pp, by=c("host","para"))
ov_ad_pf_pp_merge <- merge(ov_ad_pp, pf_pp, by=c("host","para"))
ov_ad_pfE_pp_merge <- merge(ov_ad_pp, pfE_pp, by=c("host","para"))
ov_ad_pv_pp_merge <- merge(ov_ad_pp, pv_pp, by=c("host","para"))

ov_ad_pb_pp_bp_merge <- merge(ov_ad_pp_bp, pb_pp_bp, by=c("host","para"))
ov_ad_pf_pp_bp_merge <- merge(ov_ad_pp_bp, pf_pp_bp, by=c("host","para"))
ov_ad_pfE_pp_bp_merge <- merge(ov_ad_pp_bp, pfE_pp_bp, by=c("host","para"))
ov_ad_pv_pp_bp_merge <- merge(ov_ad_pp_bp, pv_pp_bp, by=c("host","para"))

ov_ad_pb_bp_merge <- merge(ov_ad_bp, pb_bp, by=c("host","para"))
ov_ad_pf_bp_merge <- merge(ov_ad_bp, pf_bp, by=c("host","para"))
ov_ad_pfE_bp_merge <- merge(ov_ad_bp, pfE_bp, by=c("host","para"))
ov_ad_pv_bp_merge <- merge(ov_ad_bp, pv_bp, by=c("host","para"))

intersect_edges <- list(ov_ad_pb_pp_merge = ov_ad_pb_pp_merge,
  ov_ad_pf_pp_merge = ov_ad_pf_pp_merge, 
  ov_ad_pfE_pp_merge = ov_ad_pfE_pp_merge,
  ov_ad_pv_pp_merge = ov_ad_pv_pp_merge,

  ov_ad_pb_pp_bp_merge = ov_ad_pb_pp_bp_merge,
  ov_ad_pf_pp_bp_merge = ov_ad_pf_pp_bp_merge, 
  ov_ad_pfE_pp_bp_merge = ov_ad_pfE_pp_bp_merge,
  ov_ad_pv_pp_bp_merge = ov_ad_pv_pp_bp_merge,

  ov_ad_pb_bp_merge = ov_ad_pb_bp_merge,
  ov_ad_pf_bp_merge = ov_ad_pf_bp_merge,
  ov_ad_pfE_bp_merge = ov_ad_pfE_bp_merge,
  ov_ad_pv_bp_merge = ov_ad_pv_bp_merge)

property_list <- list()
a = 1

for(i in 1:length(intersect_edges))
{
  print(i)
  n <- stringr::str_sub(string = names(intersect_edges[i]),start = 1, end = -7)
  d <- data.frame(h = as.character(intersect_edges[[i]][,1]), p = as.character(intersect_edges[[i]][,2]))
  ig <- graph_from_data_frame(d, directed = F)

  dg <- degree(ig, v = V(ig), loops = F, normalized = F)


  dg_p <- dg[grep(pattern = "p_OG", names(dg))]
  dg_df <- as.data.frame(dg_p) %>%
    tibble::rownames_to_column("Orthogroup")
  colnames(dg_df)[2] <- paste0(n, "_dg")
## betweenness
  bw <- betweenness(ig, v = V(ig), directed = FALSE)
  bw_p <- bw[grep(pattern = "p_OG", names(bw))]
  bw_df <- as.data.frame(bw_p) %>%
    tibble::rownames_to_column("Orthogroup")
  colnames(bw_df)[2] <- paste0(n, "_bw")

## eigencentrality
  ec <- eigen_centrality(ig, directed = FALSE)
  ec_p <- ec$vector[grep(pattern = "p_OG", names(ec$vector))]
  ec_df <- as.data.frame(ec_p) %>%
    tibble::rownames_to_column("Orthogroup")
  colnames(ec_df)[2] <- paste0(n, "_ec")

property_list[[a]] <- assign(paste0(n, "_dg"),get("dg_df"))
names(property_list)[a] <- paste0(n, "_dg")
a = a+1
property_list[[a]] <- assign(paste0(n, "_bw"),get("bw_df"))
names(property_list)[a] <- paste0(n, "_bw")
a = a+1
property_list[[a]] <- assign(paste0(n, "_ec"),get("ec_df"))
names(property_list)[a] <- paste0(n, "_ec")
a = a+1
  
}

join_df <-plyr::join_all(property_list, by = "Orthogroup", type = "full")

join_df[is.na(join_df)] <- 0

load("Unweighted_RGR_MIS_cleaned.RData")
RGR_MIS_cleaned <- left_join(RGR_MIS_cleaned, join_df)
save(RGR_MIS_cleaned, file = "Unweighted_RGR_MIS_cleaned.RData")


#### blood core edges

# bipartite
blood_core_edges <- loadRData("blood_core_edges.RData")
#blood_core_host <- blood_core_edges[,"host"]
#blood_core_para <- blood_core_edges[,"para"]
d  = blood_core_edges
ig <- graph_from_data_frame(d, directed = F)

dg <- degree(ig, v = V(ig), loops = F, normalized = F)


  dg_p <- dg[grep(pattern = "p_OG", names(dg))]
  dg_df <- as.data.frame(dg_p) %>%
    tibble::rownames_to_column("Orthogroup")
colnames(dg_df)[2] <- "bl_core_bp_dg"
## betweenness
bw <- betweenness(ig, v = V(ig), directed = FALSE)
  bw_p <- bw[grep(pattern = "p_OG", names(bw))]
  bw_df <- as.data.frame(bw_p) %>%
    tibble::rownames_to_column("Orthogroup")
colnames(bw_df)[2] <- "bl_core_bp_bw"

## eigencentrality
ec <- eigen_centrality(ig, directed = FALSE)
  ec_p <- ec$vector[grep(pattern = "p_OG", names(ec$vector))]
  ec_df <- as.data.frame(ec_p) %>%
    tibble::rownames_to_column("Orthogroup")
colnames(ec_df)[2] <- "bl_core_bp_ec"

property_list <- list(bl_core_bp_dg = dg_df, bl_core_bp_bw = bw_df, bl_core_bp_ec = ec_df)

join_df <-plyr::join_all(property_list, by = "Orthogroup", type = "full")

join_df[is.na(join_df)] <- 0

load("Unweighted_RGR_MIS_cleaned.RData")
RGR_MIS_cleaned <- left_join(RGR_MIS_cleaned, join_df)
save(RGR_MIS_cleaned, file = "Unweighted_RGR_MIS_cleaned.RData")

# para
blood_core_para <- loadRData("blood_core_para_edges.RData")
#blood_core_host <- blood_core_edges[,"host"]
#blood_core_para <- blood_core_edges[,"para"]
d  = blood_core_para
ig <- graph_from_data_frame(d, directed = F)

dg <- degree(ig, v = V(ig), loops = F, normalized = F)


dg_p <- dg[grep(pattern = "p_OG", names(dg))]
  dg_df <- as.data.frame(dg_p) %>%
    tibble::rownames_to_column("Orthogroup")
colnames(dg_df)[2] <- "bl_core_pp_dg"
## betweenness
bw <- betweenness(ig, v = V(ig), directed = FALSE)
  bw_p <- bw[grep(pattern = "p_OG", names(bw))]
  bw_df <- as.data.frame(bw_p) %>%
    tibble::rownames_to_column("Orthogroup")
colnames(bw_df)[2] <- "bl_core_pp_bw"

## eigencentrality
ec <- eigen_centrality(ig, directed = FALSE)
  ec_p <- ec$vector[grep(pattern = "p_OG", names(ec$vector))]
  ec_df <- as.data.frame(ec_p) %>%
    tibble::rownames_to_column("Orthogroup")
colnames(ec_df)[2] <- "bl_core_pp_ec"
}

property_list <- list(bl_core_pp_dg = dg_df, bl_core_pp_bw = bw_df, bl_core_pp_ec = ec_df)

join_df <-plyr::join_all(property_list, by = "Orthogroup", type = "full")

join_df[is.na(join_df)] <- 0

load("Unweighted_RGR_MIS_cleaned.RData")
RGR_MIS_cleaned <- left_join(RGR_MIS_cleaned, join_df)
save(RGR_MIS_cleaned, file = "Unweighted_RGR_MIS_cleaned.RData")

# pp_bp
blood_core <- rbind(blood_core_para, blood_core_edges)
#blood_core_host <- blood_core_edges[,"host"]
#blood_core_para <- blood_core_edges[,"para"]
d  = blood_core
ig <- graph_from_data_frame(d, directed = F)

dg <- degree(ig, v = V(ig), loops = F, normalized = F)


dg_p <- dg[grep(pattern = "p_OG", names(dg))]
  dg_df <- as.data.frame(dg_p) %>%
    tibble::rownames_to_column("Orthogroup")
colnames(dg_df)[2] <- "bl_core_pp_bp_dg"
## betweenness
bw <- betweenness(ig, v = V(ig), directed = FALSE)
  bw_p <- bw[grep(pattern = "p_OG", names(bw))]
  bw_df <- as.data.frame(bw_p) %>%
    tibble::rownames_to_column("Orthogroup")
colnames(bw_df)[2] <- "bl_core_pp_bp_bw"

## eigencentrality
ec <- eigen_centrality(ig, directed = FALSE)
  ec_p <- ec$vector[grep(pattern = "p_OG", names(ec$vector))]
  ec_df <- as.data.frame(ec_p) %>%
    tibble::rownames_to_column("Orthogroup")
colnames(ec_df)[2] <- "bl_core_pp_bp_ec"

property_list <- list(bl_core_pp_bp_dg = dg_df, bl_core_pp_bp_bw = bw_df, bl_core_pp_bp_ec = ec_df)

join_df <-plyr::join_all(property_list, by = "Orthogroup", type = "full")

join_df[is.na(join_df)] <- 0

load("Unweighted_RGR_MIS_cleaned.RData")
RGR_MIS_cleaned <- left_join(RGR_MIS_cleaned, join_df)
save(RGR_MIS_cleaned, file = "Unweighted_RGR_MIS_cleaned.RData")

pb_pp_dg_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pb_pp_dg)
pb_pp_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pb_pp_ec)

# ov models
ov_pp_dg_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_pp_dg)
ov_pp_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_pp_ec)
ov_pp_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_pp_bw)

ov_pp_dg_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_pp_dg + ov_pp_bw)
ov_pp_dg_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_pp_dg + ov_pp_ec)

ov_pp_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_pp_dg + ov_pp_bw + ov_pp_ec)



ov_bp_dg_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_bp_dg)
ov_bp_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_bp_ec)
ov_bp_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_bp_bw)

ov_bp_dg_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_bp_dg + ov_bp_bw)
ov_bp_dg_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_bp_dg + ov_bp_ec)

ov_bp_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_bp_dg + ov_bp_bw + ov_bp_ec)

ov_pp_bp_dg_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_pp_bp_dg)
ov_pp_bp_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_pp_bp_ec)
ov_pp_bp_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_pp_bp_bw)

ov_pp_bp_dg_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_pp_bp_dg + ov_pp_bp_bw)
ov_pp_bp_dg_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_pp_bp_dg + ov_pp_bp_ec)

ov_pp_bp_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_pp_bp_dg + ov_pp_bp_bw + ov_pp_bp_ec)


ov_all_dg_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_all_dg)
ov_all_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_all_ec)
ov_all_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_all_bw)

ov_all_dg_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_all_dg + ov_all_bw)
ov_all_dg_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_all_dg + ov_all_ec)

ov_all_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_all_dg + ov_all_bw + ov_all_ec)

# pfE models

pb_pp_dg_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pb_pp_dg)


pb_pp_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pb_pp_ec)
pb_pp_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pb_pp_bw)

pb_pp_dg_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pb_pp_dg + pb_pp_bw)
pb_pp_dg_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pb_pp_dg + pb_pp_ec)

pb_pp_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pb_pp_dg + pb_pp_bw + pb_pp_ec)

pb_bp_dg_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pb_bp_dg)
pb_bp_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pb_bp_ec)
pb_bp_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pb_bp_bw)

pb_bp_dg_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pb_bp_dg + pb_bp_bw)
pb_bp_dg_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pb_bp_dg + pb_bp_ec)

pb_bp_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pb_bp_dg + pb_bp_bw + pb_bp_ec)

pb_pp_bp_dg_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pb_pp_bp_dg)
pb_pp_bp_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pb_pp_bp_ec)
pb_pp_bp_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pb_pp_bp_bw)

pb_pp_bp_dg_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pb_pp_bp_dg + pb_pp_bp_bw)
pb_pp_bp_dg_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pb_pp_bp_dg + pb_pp_bp_ec)

pb_pp_bp_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pb_pp_bp_dg + pb_pp_bp_bw + pb_pp_bp_ec)


pb_all_dg_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pb_all_dg)
pb_all_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pb_all_ec)
pb_all_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pb_all_bw)

pb_all_dg_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pb_all_dg + pb_all_bw)
pb_all_dg_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pb_all_dg + pb_all_ec)

pb_all_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pb_all_dg + pb_all_bw + pb_all_ec)

# pf models

pf_pp_dg_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pf_pp_dg)
pf_pp_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pf_pp_ec)
pf_pp_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pf_pp_bw)

pf_pp_dg_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pf_pp_dg + pf_pp_bw)
pf_pp_dg_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pf_pp_dg + pf_pp_ec)

pf_pp_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pf_pp_dg + pf_pp_bw + pf_pp_ec)

pf_bp_dg_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pf_bp_dg)
pf_bp_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pf_bp_ec)
pf_bp_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pf_bp_bw)

pf_bp_dg_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pf_bp_dg + pf_bp_bw)
pf_bp_dg_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pf_bp_dg + pf_bp_ec)

pf_bp_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pf_bp_dg + pf_bp_bw + pf_bp_ec)

pf_pp_bp_dg_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pf_pp_bp_dg)
pf_pp_bp_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pf_pp_bp_ec)
pf_pp_bp_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pf_pp_bp_bw)

pf_pp_bp_dg_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pf_pp_bp_dg + pf_pp_bp_bw)
pf_pp_bp_dg_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pf_pp_bp_dg + pf_pp_bp_ec)

pf_pp_bp_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pf_pp_bp_dg + pf_pp_bp_bw + pf_pp_bp_ec)


pf_all_dg_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pf_all_dg)
pf_all_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pf_all_ec)
pf_all_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pf_all_bw)

pf_all_dg_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pf_all_dg + pf_all_bw)
pf_all_dg_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pf_all_dg + pf_all_ec)

pf_all_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pf_all_dg + pf_all_bw + pf_all_ec)

# pfE models

pfE_pp_dg_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pfE_pp_dg)
pfE_pp_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pfE_pp_ec)
pfE_pp_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pfE_pp_bw)

pfE_pp_dg_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pfE_pp_dg + pfE_pp_bw)
pfE_pp_dg_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pfE_pp_dg + pfE_pp_ec)

pfE_pp_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pfE_pp_dg + pfE_pp_bw + pfE_pp_ec)

pfE_bp_dg_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pfE_bp_dg)
pfE_bp_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pfE_bp_ec)
pfE_bp_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pfE_bp_bw)

pfE_bp_dg_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pfE_bp_dg + pfE_bp_bw)
pfE_bp_dg_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pfE_bp_dg + pfE_bp_ec)

pfE_bp_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pfE_bp_dg + pfE_bp_bw + pfE_bp_ec)

pfE_pp_bp_dg_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pfE_pp_bp_dg)
pfE_pp_bp_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pfE_pp_bp_ec)
pfE_pp_bp_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pfE_pp_bp_bw)

pfE_pp_bp_dg_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pfE_pp_bp_dg + pfE_pp_bp_bw)
pfE_pp_bp_dg_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pfE_pp_bp_dg + pfE_pp_bp_ec)

pfE_pp_bp_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pfE_pp_bp_dg + pfE_pp_bp_bw + pfE_pp_bp_ec)


pfE_all_dg_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pfE_all_dg)
pfE_all_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pfE_all_ec)
pfE_all_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pfE_all_bw)

pfE_all_dg_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pfE_all_dg + pfE_all_bw)
pfE_all_dg_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pfE_all_dg + pfE_all_ec)

pfE_all_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pfE_all_dg + pfE_all_bw + pfE_all_ec)

# pv models
pv_pp_dg_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pv_pp_dg)
pv_pp_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pv_pp_ec)
pv_pp_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pv_pp_bw)

pv_pp_dg_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pv_pp_dg + pv_pp_bw)
pv_pp_dg_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pv_pp_dg + pv_pp_ec)

pv_pp_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pv_pp_dg + pv_pp_bw + pv_pp_ec)

pv_bp_dg_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pv_bp_dg)
pv_bp_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pv_bp_ec)
pv_bp_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pv_bp_bw)

pv_bp_dg_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pv_bp_dg + pv_bp_bw)
pv_bp_dg_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pv_bp_dg + pv_bp_ec)

pv_bp_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pv_bp_dg + pv_bp_bw + pv_bp_ec)

pv_pp_bp_dg_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pv_pp_bp_dg)
pv_pp_bp_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pv_pp_bp_ec)
pv_pp_bp_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pv_pp_bp_bw)

pv_pp_bp_dg_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pv_pp_bp_dg + pv_pp_bp_bw)
pv_pp_bp_dg_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pv_pp_bp_dg + pv_pp_bp_ec)

pv_pp_bp_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pv_pp_bp_dg + pv_pp_bp_bw + pv_pp_bp_ec)

pv_all_dg_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pv_all_dg)
pv_all_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pv_all_ec)
pv_all_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pv_all_bw)

pv_all_dg_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pv_all_dg + pv_all_bw)
pv_all_dg_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pv_all_dg + pv_all_ec)

pv_all_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ pv_all_dg + pv_all_bw + pv_all_ec)


# ov_ad models

ov_ad_pp_dg_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_ad_pp_dg)
ov_ad_pp_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_ad_pp_ec)
ov_ad_pp_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_ad_pp_bw)

ov_ad_pp_dg_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_ad_pp_dg + ov_ad_pp_bw)
ov_ad_pp_dg_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_ad_pp_dg + ov_ad_pp_ec)

ov_ad_pp_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_ad_pp_dg + ov_ad_pp_bw + ov_ad_pp_ec)

ov_ad_bp_dg_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_ad_bp_dg)
ov_ad_bp_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_ad_bp_ec)
ov_ad_bp_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_ad_bp_bw)

ov_ad_bp_dg_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_ad_bp_dg + ov_ad_bp_bw)
ov_ad_bp_dg_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_ad_bp_dg + ov_ad_bp_ec)

ov_ad_bp_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_ad_bp_dg + ov_ad_bp_bw + ov_ad_bp_ec)

ov_ad_pp_bp_dg_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_ad_pp_bp_dg)
ov_ad_pp_bp_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_ad_pp_bp_ec)
ov_ad_pp_bp_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_ad_pp_bp_bw)

ov_ad_pp_bp_dg_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_ad_pp_bp_dg + ov_ad_pp_bp_bw)
ov_ad_pp_bp_dg_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_ad_pp_bp_dg + ov_ad_pp_bp_ec)

ov_ad_pp_bp_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_ad_pp_bp_dg + ov_ad_pp_bp_bw + ov_ad_pp_bp_ec)

ov_ad_all_dg_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_ad_all_dg)
ov_ad_all_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_ad_all_ec)
ov_ad_all_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_ad_all_bw)

ov_ad_all_dg_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_ad_all_dg + ov_ad_all_bw)
ov_ad_all_dg_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_ad_all_dg + ov_ad_all_ec)

ov_ad_all_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ ov_ad_all_dg + ov_ad_all_bw + ov_ad_all_ec)

# blood core

bl_core_pp_dg_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ core_pp_dg)
bl_core_pp_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ core_pp_ec)
bl_core_pp_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ core_pp_bw)

bl_core_pp_dg_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ core_pp_dg + core_pp_bw)
bl_core_pp_dg_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ core_pp_dg + core_pp_ec)

bl_core_pp_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ core_pp_dg + core_pp_bw + core_pp_ec)

bl_core_bp_dg_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ bl_core_bp_dg)
bl_core_bp_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ bl_core_bp_ec)
bl_core_bp_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ bl_core_bp_bw)

bl_core_bp_dg_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ bl_core_bp_dg + bl_core_bp_bw)
bl_core_bp_dg_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ bl_core_bp_dg + bl_core_bp_ec)

bl_core_bp_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ bl_core_bp_dg + bl_core_bp_bw + bl_core_bp_ec)

bl_core_pp_bp_dg_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ bl_core_pp_bp_dg)
bl_core_pp_bp_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ bl_core_pp_bp_ec)
bl_core_pp_bp_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ bl_core_pp_bp_bw)

bl_core_pp_bp_dg_bw_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ bl_core_pp_bp_dg + bl_core_pp_bp_bw)
bl_core_pp_bp_dg_ec_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ bl_core_pp_bp_dg + bl_core_pp_bp_ec)

bl_core_pp_bp_m <- betareg(data = rgr_mis_allgenes_clean, RGR ~ bl_core_pp_bp_dg + bl_core_pp_bp_bw + bl_core_pp_bp_ec)

# MIS models

# pb models 
# pb models
pb_pp_dg_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pb_pp_dg)
pb_pp_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pb_pp_ec)
pb_pp_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pb_pp_bw)

pb_pp_dg_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pb_pp_dg + pb_pp_bw)
pb_pp_dg_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pb_pp_dg + pb_pp_ec)

pb_pp_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pb_pp_dg + pb_pp_bw + pb_pp_ec)

pb_bp_dg_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pb_bp_dg)
pb_bp_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pb_bp_ec)
pb_bp_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pb_bp_bw)

pb_bp_dg_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pb_bp_dg + pb_bp_bw)
pb_bp_dg_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pb_bp_dg + pb_bp_ec)

pb_bp_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pb_bp_dg + pb_bp_bw + pb_bp_ec)

pb_pp_bp_dg_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pb_pp_bp_dg)
pb_pp_bp_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pb_pp_bp_ec)
pb_pp_bp_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pb_pp_bp_bw)

pb_pp_bp_dg_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pb_pp_bp_dg + pb_pp_bp_bw)
pb_pp_bp_dg_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pb_pp_bp_dg + pb_pp_bp_ec)

pb_pp_bp_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pb_pp_bp_dg + pb_pp_bp_bw + pb_pp_bp_ec)


pb_all_dg_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pb_all_dg)
pb_all_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pb_all_ec)
pb_all_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pb_all_bw)

pb_all_dg_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pb_all_dg + pb_all_bw)
pb_all_dg_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pb_all_dg + pb_all_ec)

pb_all_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pb_all_dg + pb_all_bw + pb_all_ec)


# ov models
ov_pp_dg_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_pp_dg)
ov_pp_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_pp_ec)
ov_pp_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_pp_bw)

ov_pp_dg_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_pp_dg + ov_pp_bw)
ov_pp_dg_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_pp_dg + ov_pp_ec)

ov_pp_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_pp_dg + ov_pp_bw + ov_pp_ec)

ov_bp_dg_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_bp_dg)
ov_bp_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_bp_ec)
ov_bp_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_bp_bw)

ov_bp_dg_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_bp_dg + ov_bp_bw)
ov_bp_dg_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_bp_dg + ov_bp_ec)

ov_bp_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_bp_dg + ov_bp_bw + ov_bp_ec)

ov_pp_bp_dg_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_pp_bp_dg)
ov_pp_bp_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_pp_bp_ec)
ov_pp_bp_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_pp_bp_bw)

ov_pp_bp_dg_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_pp_bp_dg + ov_pp_bp_bw)
ov_pp_bp_dg_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_pp_bp_dg + ov_pp_bp_ec)

ov_pp_bp_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_pp_bp_dg + ov_pp_bp_bw + ov_pp_bp_ec)


ov_all_dg_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_all_dg)
ov_all_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_all_ec)
ov_all_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_all_bw)

ov_all_dg_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_all_dg + ov_all_bw)
ov_all_dg_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_all_dg + ov_all_ec)

ov_all_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_all_dg + ov_all_bw + ov_all_ec)

# pfE models

pfE_pp_dg_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_pp_dg)
pfE_pp_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_pp_ec)
pfE_pp_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_pp_bw)

pfE_pp_dg_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_pp_dg + pfE_pp_bw)
pfE_pp_dg_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_pp_dg + pfE_pp_ec)

pfE_pp_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_pp_dg + pfE_pp_bw + pfE_pp_ec)

pfE_bp_dg_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_bp_dg)
pfE_bp_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_bp_ec)
pfE_bp_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_bp_bw)

pfE_bp_dg_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_bp_dg + pfE_bp_bw)
pfE_bp_dg_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_bp_dg + pfE_bp_ec)

pfE_bp_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_bp_dg + pfE_bp_bw + pfE_bp_ec)

pfE_pp_bp_dg_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_pp_bp_dg)
pfE_pp_bp_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_pp_bp_ec)
pfE_pp_bp_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_pp_bp_bw)

pfE_pp_bp_dg_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_pp_bp_dg + pfE_pp_bp_bw)
pfE_pp_bp_dg_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_pp_bp_dg + pfE_pp_bp_ec)

pfE_pp_bp_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_pp_bp_dg + pfE_pp_bp_bw + pfE_pp_bp_ec)


pfE_all_dg_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_all_dg)
pfE_all_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_all_ec)
pfE_all_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_all_bw)

pfE_all_dg_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_all_dg + pfE_all_bw)
pfE_all_dg_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_all_dg + pfE_all_ec)

pfE_all_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_all_dg + pfE_all_bw + pfE_all_ec)

# pf models

pf_pp_dg_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pf_pp_dg)
pf_pp_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pf_pp_ec)
pf_pp_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pf_pp_bw)

pf_pp_dg_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pf_pp_dg + pf_pp_bw)
pf_pp_dg_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pf_pp_dg + pf_pp_ec)

pf_pp_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pf_pp_dg + pf_pp_bw + pf_pp_ec)

pf_bp_dg_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pf_bp_dg)
pf_bp_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pf_bp_ec)
pf_bp_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pf_bp_bw)

pf_bp_dg_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pf_bp_dg + pf_bp_bw)
pf_bp_dg_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pf_bp_dg + pf_bp_ec)

pf_bp_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pf_bp_dg + pf_bp_bw + pf_bp_ec)

pf_pp_bp_dg_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pf_pp_bp_dg)
pf_pp_bp_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pf_pp_bp_ec)
pf_pp_bp_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pf_pp_bp_bw)

pf_pp_bp_dg_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pf_pp_bp_dg + pf_pp_bp_bw)
pf_pp_bp_dg_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pf_pp_bp_dg + pf_pp_bp_ec)

pf_pp_bp_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pf_pp_bp_dg + pf_pp_bp_bw + pf_pp_bp_ec)


pf_all_dg_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pf_all_dg)
pf_all_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pf_all_ec)
pf_all_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pf_all_bw)

pf_all_dg_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pf_all_dg + pf_all_bw)
pf_all_dg_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pf_all_dg + pf_all_ec)

pf_all_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pf_all_dg + pf_all_bw + pf_all_ec)

# pfE models

pfE_pp_dg_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_pp_dg)
pfE_pp_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_pp_ec)
pfE_pp_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_pp_bw)

pfE_pp_dg_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_pp_dg + pfE_pp_bw)
pfE_pp_dg_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_pp_dg + pfE_pp_ec)

pfE_pp_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_pp_dg + pfE_pp_bw + pfE_pp_ec)

pfE_bp_dg_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_bp_dg)
pfE_bp_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_bp_ec)
pfE_bp_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_bp_bw)

pfE_bp_dg_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_bp_dg + pfE_bp_bw)
pfE_bp_dg_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_bp_dg + pfE_bp_ec)

pfE_bp_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_bp_dg + pfE_bp_bw + pfE_bp_ec)

pfE_pp_bp_dg_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_pp_bp_dg)
pfE_pp_bp_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_pp_bp_ec)
pfE_pp_bp_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_pp_bp_bw)

pfE_pp_bp_dg_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_pp_bp_dg + pfE_pp_bp_bw)
pfE_pp_bp_dg_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_pp_bp_dg + pfE_pp_bp_ec)

pfE_pp_bp_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_pp_bp_dg + pfE_pp_bp_bw + pfE_pp_bp_ec)


pfE_all_dg_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_all_dg)
pfE_all_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_all_ec)
pfE_all_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_all_bw)

pfE_all_dg_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_all_dg + pfE_all_bw)
pfE_all_dg_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_all_dg + pfE_all_ec)

pfE_all_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pfE_all_dg + pfE_all_bw + pfE_all_ec)

# pv
pv_pp_dg_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pv_pp_dg)
pv_pp_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pv_pp_ec)
pv_pp_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pv_pp_bw)

pv_pp_dg_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pv_pp_dg + pv_pp_bw)
pv_pp_dg_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pv_pp_dg + pv_pp_ec)

pv_pp_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pv_pp_dg + pv_pp_bw + pv_pp_ec)

pv_bp_dg_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pv_bp_dg)
pv_bp_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pv_bp_ec)
pv_bp_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pv_bp_bw)

pv_bp_dg_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pv_bp_dg + pv_bp_bw)
pv_bp_dg_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pv_bp_dg + pv_bp_ec)

pv_bp_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pv_bp_dg + pv_bp_bw + pv_bp_ec)

pv_pp_bp_dg_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pv_pp_bp_dg)
pv_pp_bp_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pv_pp_bp_ec)
pv_pp_bp_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pv_pp_bp_bw)

pv_pp_bp_dg_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pv_pp_bp_dg + pv_pp_bp_bw)
pv_pp_bp_dg_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pv_pp_bp_dg + pv_pp_bp_ec)

pv_pp_bp_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pv_pp_bp_dg + pv_pp_bp_bw + pv_pp_bp_ec)


pv_all_dg_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pv_all_dg)
pv_all_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pv_all_ec)
pv_all_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pv_all_bw)

pv_all_dg_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pv_all_dg + pv_all_bw)
pv_all_dg_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pv_all_dg + pv_all_ec)

pv_all_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ pv_all_dg + pv_all_bw + pv_all_ec)

# ov_ad models

ov_ad_pp_dg_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_ad_pp_dg)
ov_ad_pp_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_ad_pp_ec)
ov_ad_pp_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_ad_pp_bw)

ov_ad_pp_dg_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_ad_pp_dg + ov_ad_pp_bw)
ov_ad_pp_dg_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_ad_pp_dg + ov_ad_pp_ec)

ov_ad_pp_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_ad_pp_dg + ov_ad_pp_bw + ov_ad_pp_ec)

ov_ad_bp_dg_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_ad_bp_dg)
ov_ad_bp_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_ad_bp_ec)
ov_ad_bp_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_ad_bp_bw)

ov_ad_bp_dg_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_ad_bp_dg + ov_ad_bp_bw)
ov_ad_bp_dg_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_ad_bp_dg + ov_ad_bp_ec)

ov_ad_bp_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_ad_bp_dg + ov_ad_bp_bw + ov_ad_bp_ec)

ov_ad_pp_bp_dg_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_ad_pp_bp_dg)
ov_ad_pp_bp_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_ad_pp_bp_ec)
ov_ad_pp_bp_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_ad_pp_bp_bw)

ov_ad_pp_bp_dg_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_ad_pp_bp_dg + ov_ad_pp_bp_bw)
ov_ad_pp_bp_dg_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_ad_pp_bp_dg + ov_ad_pp_bp_ec)

ov_ad_pp_bp_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_ad_pp_bp_dg + ov_ad_pp_bp_bw + ov_ad_pp_bp_ec)


ov_ad_all_dg_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_ad_all_dg)
ov_ad_all_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_ad_all_ec)
ov_ad_all_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_ad_all_bw)

ov_ad_all_dg_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_ad_all_dg + ov_ad_all_bw)
ov_ad_all_dg_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_ad_all_dg + ov_ad_all_ec)

ov_ad_all_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ ov_ad_all_dg + ov_ad_all_bw + ov_ad_all_ec)

## blood core

bl_core_pp_dg_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ core_pp_dg)
bl_core_pp_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ core_pp_ec)
bl_core_pp_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ core_pp_bw)

bl_core_pp_dg_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ core_pp_dg + core_pp_bw)
bl_core_pp_dg_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ core_pp_dg + core_pp_ec)

bl_core_pp_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ core_pp_dg + core_pp_bw + core_pp_ec)

bl_core_bp_dg_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ bl_core_bp_dg)
bl_core_bp_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ bl_core_bp_ec)
bl_core_bp_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ bl_core_bp_bw)

bl_core_bp_dg_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ bl_core_bp_dg + bl_core_bp_bw)
bl_core_bp_dg_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ bl_core_bp_dg + bl_core_bp_ec)

bl_core_bp_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ bl_core_bp_dg + bl_core_bp_bw + bl_core_bp_ec)

bl_core_pp_bp_dg_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ bl_core_pp_bp_dg)
bl_core_pp_bp_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ bl_core_pp_bp_ec)
bl_core_pp_bp_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ bl_core_pp_bp_bw)

bl_core_pp_bp_dg_bw_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ bl_core_pp_bp_dg + bl_core_pp_bp_bw)
bl_core_pp_bp_dg_ec_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ bl_core_pp_bp_dg + bl_core_pp_bp_ec)

bl_core_pp_bp_mis <- betareg(data = rgr_mis_allgenes_clean, MIS ~ bl_core_pp_bp_dg + bl_core_pp_bp_bw + bl_core_pp_bp_ec)


# trying to make the table automated
# input: list of models to enter into the table in the desired order
# | model name | dg | bw | cl | ec | kc | prsq |

#ls()[which(class(eval(parse(text = ls()))) == "betareg")]

summary(bl_core_pp_dg_m)$coef

f <-c()
for(i in ls())
{
  if(class(eval(parse(text = i))) == "betareg")
    f <- append(f, i)
}
models <- f


model_table <- function(models)
{
  remove <- c(0)
  for(x in 1:length(models))
  {
    if(!exists(models[x]))
      remove <- append(remove, x)
  }

  if(length(remove) > 1)
    models <- models[-remove]

  tab <- data.frame(model_name = rep("mn", 150),
    dg = rep("dg", 150),
    bw = rep("bw", 150),
    cl = rep("cl", 150),
    ec = rep("ec", 150),
    kc = rep("kc", 150),
    psrq = rep("psrq", 150))
  tab[] <- lapply(tab, as.character)

  for(i in 1:length(models))
  {
    tab[i,1] <- names(models[i])
    #as.character(substitute(models[i]))

    m <- eval(parse(text = models[i]))

    coefs <- as.data.frame(m$call$coefficients$mean)
    coef_names <- names(m$call$coefficients$mean)

    for(j in 2:length(coef_names))
    {
        #pred <- strsplit(coef_names[j], "_")[[1]][3]
        pred <- str_sub(coef_names[j], 42, -19)
        
        effect_size <- format(exp(coefs[j, 1]), nsmall = 3, scientific = T, digits = 5)
        #pval <- format(coefs[j, 4], nsmall = 3, scientific = T, digits = 5)
        #pval <- format(models[i][[1]]$coefficients$precision[4], nsmall = 3, scientific = T, digits = 5) # this is wrong, this is the pval for precision model, not coefs

        eff_pv <- paste(effect_size, pval, sep = "; ")

        tab[i,grep(pattern = pred, colnames(tab))] = eff_pv
    }

    tab[i,7] <- m[["pseudo.r.squared"]]
  }
  return(tab)
}

model_table(models_list)

write.table(tab, paste0(org,"_", nw, "_","models_table_exp.txt", collapse = ""), sep = '\t', row.names = F)



## RGR plots ##
ggplot(data=RGR_MIS_cleaned, aes(x=pfE_pp_ec, y=RGR, colour = phenotype.y))+
  geom_point(alpha = 0.8) +
  #scale_colour_manual(values=unique(as.character(RGR_MIS_cleaned$color.codes))) +
  theme_bw() + 
  ggtitle("Relative growth rate vs eigen centrality of P. berghei dataset") +
  xlab("Eigenvector centrality") + 
  ylab("Relative Growth Rate")
ggsave("RGR_vs_pfE_pp_ec.png")

ggplot(data=RGR_MIS_cleaned, aes(x=pfE_pp_dg, y=RGR, colour = phenotype.y))+
  geom_point(alpha = 0.8) +
  #scale_colour_manual(values=unique(as.character(RGR_MIS_cleaned$color.codes))) +
  theme_bw() + 
  ggtitle("Relative growth rate vs node degree of P. berghei dataset") +
  xlab("Node degree") + 
  ylab("Relative Growth Rate")
ggsave("RGR_vs_pfE_pp_dg.png")

ggplot(data=RGR_MIS_cleaned, aes(x=ov_pp_ec, y=RGR, colour = phenotype.y))+
  geom_point(alpha = 0.8) +
  #scale_colour_manual(values=unique(as.character(RGR_MIS_cleaned$color.codes))) +
  theme_bw() + 
  ggtitle("Relative growth rate vs eigen centrality of Overall dataset") +
  xlab("Eigenvector centrality") + 
  ylab("Relative Growth Rate")
ggsave("RGR_vs_ov_pp_ec.png")

ggplot(data=RGR_MIS_cleaned, aes(x=ov_pp_dg, y=RGR, colour = phenotype.y))+
  geom_point(alpha = 0.8) +
  #scale_colour_manual(values=unique(as.character(RGR_MIS_cleaned$color.codes))) +
  theme_bw() + 
  ggtitle("Relative growth rate vs node degree of Overall dataset") +
  xlab("Node degree") + 
  ylab("Relative Growth Rate")
ggsave("RGR_vs_ov_pp_dg.png")

ggplot(data=RGR_MIS_cleaned, aes(x=pf_pp_ec, y=RGR, colour = phenotype.y))+
  geom_point(alpha = 0.8) +
  #scale_colour_manual(values=unique(as.character(RGR_MIS_cleaned$color.codes))) +
  theme_bw() + 
  ggtitle("Relative growth rate vs eigen centrality of P. falci (Indonesia) dataset") +
  xlab("Eigenvector centrality") + 
  ylab("Relative Growth Rate")
ggsave("RGR_vs_pf_pp_ec.png")

ggplot(data=RGR_MIS_cleaned, aes(x=pf_pp_dg, y=RGR, colour = phenotype.y))+
  geom_point(alpha = 0.8) +
  #scale_colour_manual(values=unique(as.character(RGR_MIS_cleaned$color.codes))) +
  theme_bw() + 
  ggtitle("Relative growth rate vs node degree of P. falci (Indonesia) dataset") +
  xlab("Node degree") + 
  ylab("Relative Growth Rate")
ggsave("RGR_vs_pf_pp_dg.png")

ggplot(data=RGR_MIS_cleaned, aes(x=pfE_pp_ec, y=RGR, colour = phenotype.y))+
  geom_point(alpha = 0.8) +
  #scale_colour_manual(values=unique(as.character(RGR_MIS_cleaned$color.codes))) +
  theme_bw() + 
  ggtitle("Relative growth rate vs eigen centrality of P. falci (Gambia) dataset") +
  xlab("Eigenvector centrality") + 
  ylab("Relative Growth Rate")
ggsave("RGR_vs_pfE_pp_ec.png")

ggplot(data=RGR_MIS_cleaned, aes(x=pfE_pp_dg, y=RGR, colour = phenotype.y))+
  geom_point(alpha = 0.8) +
  #scale_colour_manual(values=unique(as.character(RGR_MIS_cleaned$color.codes))) +
  theme_bw() + 
  ggtitle("Relative growth rate vs node degree of P. falci (Gambia) dataset") +
  xlab("Node degree") + 
  ylab("Relative Growth Rate")
ggsave("RGR_vs_pfE_pp_dg.png")

## MIS plots

ggplot(data=RGR_MIS_cleaned, aes(x=pfE_pp_ec, y=MIS, colour = phenotype.y))+
  geom_point(alpha = 0.8) +
  #scale_colour_manual(values=unique(as.character(RGR_MIS_cleaned$color.codes))) +
  theme_bw() + 
  ggtitle("Mutagenesis Index Score vs eigen centrality of P. berghei dataset") +
  xlab("Eigenvector centrality") + 
  ylab("Mutagenesis Index Score")
ggsave("MIS_vs_pfE_pp_ec.png")

ggplot(data=RGR_MIS_cleaned, aes(x=pfE_pp_dg, y=MIS, colour = phenotype.y))+
  geom_point(alpha = 0.8) +
  #scale_colour_manual(values=unique(as.character(RGR_MIS_cleaned$color.codes))) +
  theme_bw() + 
  ggtitle("Mutagenesis Index Score vs node degree of P. berghei dataset") +
  xlab("Node degree") + 
  ylab("Mutagenesis Index Score")
ggsave("MIS_vs_pfE_pp_dg.png")

ggplot(data=RGR_MIS_cleaned, aes(x=pf_pp_ec, y=MIS, colour = phenotype.y))+
  geom_point(alpha = 0.8) +
  #scale_colour_manual(values=unique(as.character(RGR_MIS_cleaned$color.codes))) +
  theme_bw() + 
  ggtitle("Mutagenesis Index Score vs eigen centrality of P. falciparum (Indonesia) dataset") +
  xlab("Eigenvector centrality") + 
  ylab("Mutagenesis Index Score")
ggsave("MIS_vs_pf_pp_ec.png")

ggplot(data=RGR_MIS_cleaned, aes(x=pf_pp_dg, y=MIS, colour = phenotype.y))+
  geom_point(alpha = 0.8) +
  #scale_colour_manual(values=unique(as.character(RGR_MIS_cleaned$color.codes))) +
  theme_bw() + 
  ggtitle("Mutagenesis Index Score vs node degree of P. falciparum (Indonesia) dataset") +
  xlab("Node degree") + 
  ylab("Mutagenesis Index Score")
ggsave("MIS_vs_pf_pp_dg.png")

ggplot(data=RGR_MIS_cleaned, aes(x=ov_pp_ec, y=MIS, colour = phenotype.y))+
  geom_point(alpha = 0.8) +
  #scale_colour_manual(values=unique(as.character(RGR_MIS_cleaned$color.codes))) +
  theme_bw() + 
  ggtitle("Mutagenesis Index Score vs eigen centrality of Overall dataset") +
  xlab("Eigenvector centrality") + 
  ylab("Mutagenesis Index Score")
ggsave("MIS_vs_ov_pp_ec.png")

ggplot(data=RGR_MIS_cleaned, aes(x=ov_pp_dg, y=MIS, colour = phenotype.y))+
  geom_point(alpha = 0.8) +
  #scale_colour_manual(values=unique(as.character(RGR_MIS_cleaned$color.codes))) +
  theme_bw() + 
  ggtitle("Mutagenesis Index Score vs node degree of Overall dataset") +
  xlab("Node degree") + 
  ylab("Mutagenesis Index Score")
ggsave("MIS_vs_ov_pp_dg.png")

ggplot(data=RGR_MIS_cleaned, aes(x=pfE_pp_ec, y=MIS, colour = phenotype.y))+
  geom_point(alpha = 0.8) +
  #scale_colour_manual(values=unique(as.character(RGR_MIS_cleaned$color.codes))) +
  theme_bw() + 
  ggtitle("Mutagenesis Index Score vs eigen centrality of P. falciparum (Gambia) dataset") +
  xlab("Eigenvector centrality") + 
  ylab("Mutagenesis Index Score")
ggsave("MIS_vs_pfE_pp_ec.png")

ggplot(data=RGR_MIS_cleaned, aes(x=pfE_pp_dg, y=MIS, colour = phenotype.y))+
  geom_point(alpha = 0.8) +
  #scale_colour_manual(values=unique(as.character(RGR_MIS_cleaned$color.codes))) +
  theme_bw() + 
  ggtitle("Mutagenesis Index Score vs node degree of P. falciparum (Gambia) dataset") +
  xlab("Node degree") + 
  ylab("Mutagenesis Index Score")
ggsave("MIS_vs_pfE_pp_dg.png")
##

n_ess <- df[which(df$phenotype == "Essential"), "n"]
n_dis <- df[which(df$phenotype == "Dispensable"), "n"]

png("pfE_PlasmoGEM.png")
ggplot(data=RGR_MIS_cleaned, aes(x=pfE_pp_dg, y=RGR, colour = phenotype, label=gene_id))+
  geom_point(alpha = 0.4) +
  # geom_text(aes(label=df$gene_id),hjust=0, vjust=0, size = 2) +
  ggtitle("Plasmodium berghei gene interactors in correlation >= 0.9 (PlasmoGEM)") +
  xlab("node degree")
dev.off()

library(ggeffects)

## RGR

# pb
# bl_core_pp_ec_ef <- ggeffect(bl_core_pp_ec_m, terms = c("bl_core_pp_ec"))
# ggplot(pfE_pp_dg_ef, aes(x, predicted)) +
#   geom_line() +
#   geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
#   xlab("P. berghei dataset node degree")
# ggsave("pfE_pp_dg_ef.png")
# 
# plot(bl_core_pp_ec_ef, rawdata = T)#, xlab = "Eigenvector centrality of blood core network", main = "RGR prediction")

#################
############
#############
### golden code
############
##############
#############

library(ggplot2)
p <- ggplot(rgr_mis_allgenes_clean, aes(x = pb_pp_ec, y = RGR, fill = phenotype.y, colour = phenotype.y)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_line(aes(y = predict(pb_pp_ec_m, rgr_mis_allgenes_clean))) +
  theme_bw() + 
  ggtitle("Relative Growth Rate prediction") + 
  xlab("P. berghei network eigenvector centrality") + 
  ylab("Relative Growth Rate") +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 18), 
        plot.title = element_text(size = 18))
ggsave("RGR_vs_Pberg.png", width = 30, height = 20, units = "cm")

p <- ggplot(rgr_mis_allgenes_clean, aes(x = ov_pp_ec, y = RGR, fill = phenotype.y, colour = phenotype.y)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_line(aes(y = predict(ov_pp_ec_m, rgr_mis_allgenes_clean))) +
  theme_bw() + 
  ggtitle("Relative Growth Rate prediction") + 
  xlab("Overall network centrality") + 
  ylab("Relative Growth Rate") +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 18), 
        plot.title = element_text(size = 18))
ggsave("RGR_vs_overall.png", width = 30, height = 20, units = "cm")

p <- ggplot(rgr_mis_allgenes_clean, aes(x = bl_core_pp_ec, y = RGR, fill = phenotype.y, colour = phenotype.y)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_line(aes(y = predict(bl_core_pp_ec_m, rgr_mis_allgenes_clean))) +
  theme_bw() + 
  ggtitle("Relative Growth Rate prediction") + 
  xlab("Core network centrality") + 
  ylab("Relative Growth Rate") +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 18), 
        plot.title = element_text(size = 18))
ggsave("RGR_vs_core.png", width = 30, height = 20, units = "cm")

p <- ggplot(rgr_mis_allgenes_clean, aes(x = pb_pp_ec, y = MIS, fill = phenotype.y, colour = phenotype.y)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_line(aes(y = predict(pb_pp_ec_mis, rgr_mis_allgenes_clean))) +
  theme_bw() + 
  ggtitle("Mutagenesis Index Score prediction") + 
  xlab("P. berghei network eigenvector centrality") + 
  ylab("Mutagenesis Index Score") +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 18), 
        plot.title = element_text(size = 18))
ggsave("MIS_vs_Pberg.png", width = 30, height = 20, units = "cm")

p <- ggplot(rgr_mis_allgenes_clean, aes(x = ov_pp_ec, y = MIS, fill = phenotype.y, colour = phenotype.y)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_line(aes(y = predict(ov_pp_ec_mis, rgr_mis_allgenes_clean))) +
  theme_bw() + 
  ggtitle("Mutagenesis Index Score prediction") + 
  xlab("Overall network centrality") + 
  ylab("Mutagenesis Index Score") +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 18), 
        plot.title = element_text(size = 18))
ggsave("MIS_vs_overall.png", width = 30, height = 20, units = "cm")

p <- ggplot(rgr_mis_allgenes_clean, aes(x = bl_core_pp_ec, y = MIS, fill = phenotype.y, colour = phenotype.y)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_line(aes(y = predict(bl_core_pp_ec_mis, rgr_mis_allgenes_clean))) +
  theme_bw() + 
  ggtitle("Mutagenesis Index Score prediction") + 
  xlab("Core network centrality") + 
  ylab("Mutagenesis Index Score") +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 18), 
        plot.title = element_text(size = 18))
ggsave("MIS_vs_core.png", width = 30, height = 20, units = "cm")

p <- ggplot(rgr_mis_allgenes_clean, aes(x = pfE_pp_ec, y = MIS, fill = phenotype.y, colour = phenotype.y)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_line(aes(y = predict(pfE_pp_ec_mis, rgr_mis_allgenes_clean))) +
  theme_bw() + 
  ggtitle("Mutagenesis Index Score prediction") + 
  xlab("P. falciparum (E) network centrality") + 
  ylab("Mutagenesis Index Score") +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 18), 
        plot.title = element_text(size = 18))
ggsave("MIS_vs_pfE.png", width = 30, height = 20, units = "cm")

##pfE
pfE_pp_dg_ef <- ggeffect(pfE_pp_dg_m, terms = c("pfE_pp_dg"))

ggplot(pfE_pp_dg_ef, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
  xlab("P. berghei dataset node degree")
ggsave("pfE_pp_dg_ef.png")

pfE_pp_ec_ef <- ggeffect(pfE_pp_ec_m, terms = c("pfE_pp_ec"))

ggplot(pfE_pp_ec_ef, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  xlab("P. berghei dataset eigen centrality")
ggsave("pfE_pp_ec_ef.png")

##Ov
ov_pp_dg_ef <- ggeffect(ov_pp_dg_m, terms = c("ov_pp_dg"))

ggplot(ov_pp_dg_ef, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
  xlab("Overall dataset node degree")
ggsave("ov_pp_dg_ef.png")

ov_pp_ec_ef <- ggeffect(ov_pp_ec_m, terms = c("ov_pp_ec"))

ggplot(ov_pp_ec_ef, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  xlab("Overall dataset eigen centrality")
ggsave("ov_pp_ec_ef.png")

#Pf

pf_pp_dg_ef <- ggeffect(pf_pp_dg_m, terms = c("pf_pp_dg"))

ggplot(pf_pp_dg_ef, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
  xlab("P. falciparum (Indonesia) dataset node degree")
ggsave("pf_pp_dg_ef.png")

pf_pp_ec_ef <- ggeffect(pf_pp_ec_m, terms = c("pf_pp_ec"))

ggplot(pf_pp_ec_ef, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  xlab("P. falciparum (Indonesia) dataset eigen centrality")
ggsave("pf_pp_ec_ef.png")

pfE_pp_dg_ef <- ggeffect(pfE_pp_dg_m, terms = c("pfE_pp_dg"))

ggplot(pfE_pp_dg_ef, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
  xlab("P. falciparum (Gambia) dataset node degree")
ggsave("pfE_pp_dg_ef.png")

pfE_pp_ec_ef <- ggeffect(pfE_pp_ec_m, terms = c("pfE_pp_ec"))

ggplot(pfE_pp_ec_ef, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  xlab("P. falciparum (Gambia) dataset eigen centrality")
ggsave("pfE_pp_ec_ef.png")

## MIS

##pfE
pfE_pp_dg_mis_ef <- ggeffect(pfE_pp_dg_mis, terms = c("pfE_pp_dg"))

ggplot(pfE_pp_dg_mis_ef, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
  xlab("P. berghei dataset node degree (MIS)")
ggsave("pfE_pp_dg_mis_ef.png")

pfE_pp_ec_mis_ef <- ggeffect(pfE_pp_ec_mis, terms = c("pfE_pp_ec"))

ggplot(pfE_pp_ec_mis_ef, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  xlab("P. berghei dataset eigen centrality (MIS)")
ggsave("pfE_pp_ec_mis_ef.png")

##Ov
ov_pp_dg_mis_ef <- ggeffect(ov_pp_dg_mis, terms = c("ov_pp_dg"))

ggplot(ov_pp_dg_mis_ef, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
  xlab("Overall dataset node degree (MIS)")
ggsave("ov_pp_dg_mis_ef.png")

ov_pp_ec_mis_ef <- ggeffect(ov_pp_ec_mis, terms = c("ov_pp_ec"))

ggplot(ov_pp_ec_mis_ef, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  xlab("Overall dataset eigen centrality (MIS)")
ggsave("ov_pp_ec_mis_ef.png")

##Pf

pf_pp_dg_mis_ef <- ggeffect(pf_pp_dg_mis, terms = c("pf_pp_dg"))

ggplot(pf_pp_dg_mis_ef, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
  xlab("P.falciparum (Indonesia) dataset node degree (MIS)")
ggsave("pf_pp_dg_mis_ef.png")

pf_pp_ec_mis_ef <- ggeffect(pf_pp_ec_mis, terms = c("pf_pp_ec"))

ggplot(pf_pp_ec_mis_ef, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  xlab("P.falciparum (Indonesia) eigen centrality (MIS)")
ggsave("pf_pp_ec_mis_ef.png")


pfE_pp_dg_mis_ef <- ggeffect(pfE_pp_dg_mis, terms = c("pfE_pp_dg"))

ggplot(pfE_pp_dg_mis_ef, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
  xlab("P.falciparum (Gambia) dataset node degree (MIS)")
ggsave("pfE_pp_dg_mis_ef.png")

pfE_pp_ec_mis_ef <- ggeffect(pfE_pp_ec_mis, terms = c("pfE_pp_ec"))

ggplot(pfE_pp_ec_mis_ef, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  xlab("P.falciparum (Gambia) eigen centrality (MIS)")
ggsave("pfE_pp_ec_mis_ef.png")

stargazer(pb_pp_ec_m, ov_pp_ec_m, bl_core_pp_ec_m, pf_pp_ec_m, pfE_pp_ec_m,
          pb_pp_ec_mis, ov_pp_ec_mis, bl_core_pp_ec_mis, pf_pp_ec_mis, pfE_pp_ec_mis,
          title="Results", align=TRUE, type = "text")


### liver

# ov
ov_pp_dg_m <- betareg(data = liver_RGR_MIS, RGR ~ liver_ov_dg)
ov_pp_ec_m <- betareg(data = liver_RGR_MIS, RGR ~ liver_ov_ec)
ov_pp_bw_m <- betareg(data = liver_RGR_MIS, RGR ~ liver_ov_bw)

ov_pp_dg_bw_m <- betareg(data = liver_RGR_MIS, RGR ~ liver_ov_dg + liver_ov_bw)
ov_pp_dg_ec_m <- betareg(data = liver_RGR_MIS, RGR ~ liver_ov_dg + liver_ov_ec)

ov_pp_m <- betareg(data = liver_RGR_MIS, RGR ~ liver_ov_dg + liver_ov_bw + liver_ov_ec)

ov_pp_dg_mis <- betareg(data = liver_RGR_MIS, MIS ~ liver_ov_dg)
ov_pp_ec_mis <- betareg(data = liver_RGR_MIS, MIS ~ liver_ov_ec)
ov_pp_bw_mis <- betareg(data = liver_RGR_MIS, MIS ~ liver_ov_bw)

ov_pp_dg_bw_mis <- betareg(data = liver_RGR_MIS, MIS ~ liver_ov_dg + liver_ov_bw)
ov_pp_dg_ec_mis <- betareg(data = liver_RGR_MIS, MIS ~ liver_ov_dg + liver_ov_ec)

ov_pp_mis <- betareg(data = liver_RGR_MIS, MIS ~ liver_ov_dg + liver_ov_bw + liver_ov_ec)
### 

# SRP250329
SRP250329_pp_dg_m <- betareg(data = liver_RGR_MIS, RGR ~ h_SRP250329_dg)
SRP250329_pp_ec_m <- betareg(data = liver_RGR_MIS, RGR ~ h_SRP250329_ec)
SRP250329_pp_bw_m <- betareg(data = liver_RGR_MIS, RGR ~ h_SRP250329_bw)

SRP250329_pp_dg_bw_m <- betareg(data = liver_RGR_MIS, RGR ~ h_SRP250329_dg + h_SRP250329_bw)
SRP250329_pp_dg_ec_m <- betareg(data = liver_RGR_MIS, RGR ~ h_SRP250329_dg + h_SRP250329_ec)

SRP250329_pp_m <- betareg(data = liver_RGR_MIS, RGR ~ h_SRP250329_dg + h_SRP250329_bw + h_SRP250329_ec)

SRP250329_pp_dg_mis <- betareg(data = liver_RGR_MIS, MIS ~ h_SRP250329_dg)
SRP250329_pp_ec_mis <- betareg(data = liver_RGR_MIS, MIS ~ h_SRP250329_ec)
SRP250329_pp_bw_mis <- betareg(data = liver_RGR_MIS, MIS ~ h_SRP250329_bw)

SRP250329_pp_dg_bw_mis <- betareg(data = liver_RGR_MIS, MIS ~ h_SRP250329_dg + h_SRP250329_bw)
SRP250329_pp_dg_ec_mis <- betareg(data = liver_RGR_MIS, MIS ~ h_SRP250329_dg + h_SRP250329_ec)

SRP250329_pp_mis <- betareg(data = liver_RGR_MIS, MIS ~ h_SRP250329_dg + h_SRP250329_bw + h_SRP250329_ec)
### 

# ERP105548
liver_RGR_MIS <- left_join(liver_RGR_MIS, join_df)
ERP105548_pp_dg_m <- betareg(data = liver_RGR_MIS, RGR ~ h_ERP105548_dg)
ERP105548_pp_ec_m <- betareg(data = liver_RGR_MIS, RGR ~ h_ERP105548_ec)
ERP105548_pp_bw_m <- betareg(data = liver_RGR_MIS, RGR ~ h_ERP105548_bw)

ERP105548_pp_dg_bw_m <- betareg(data = liver_RGR_MIS, RGR ~ h_ERP105548_dg + h_ERP105548_bw)
ERP105548_pp_dg_ec_m <- betareg(data = liver_RGR_MIS, RGR ~ h_ERP105548_dg + h_ERP105548_ec)

ERP105548_pp_m <- betareg(data = liver_RGR_MIS, RGR ~ h_ERP105548_dg + h_ERP105548_bw + h_ERP105548_ec)

ERP105548_pp_dg_mis <- betareg(data = liver_RGR_MIS, MIS ~ h_ERP105548_dg)
ERP105548_pp_ec_mis <- betareg(data = liver_RGR_MIS, MIS ~ h_ERP105548_ec)
ERP105548_pp_bw_mis <- betareg(data = liver_RGR_MIS, MIS ~ h_ERP105548_bw)

ERP105548_pp_dg_bw_mis <- betareg(data = liver_RGR_MIS, MIS ~ h_ERP105548_dg + h_ERP105548_bw)
ERP105548_pp_dg_ec_mis <- betareg(data = liver_RGR_MIS, MIS ~ h_ERP105548_dg + h_ERP105548_ec)

ERP105548_pp_mis <- betareg(data = liver_RGR_MIS, MIS ~ h_ERP105548_dg + h_ERP105548_bw + h_ERP105548_ec)
### 

# SRP110282
liver_RGR_MIS <- left_join(liver_RGR_MIS, join_df)
SRP110282_pp_dg_m <- betareg(data = liver_RGR_MIS, RGR ~ h_SRP110282_dg)
SRP110282_pp_ec_m <- betareg(data = liver_RGR_MIS, RGR ~ h_SRP110282_ec)
SRP110282_pp_bw_m <- betareg(data = liver_RGR_MIS, RGR ~ h_SRP110282_bw)

SRP110282_pp_dg_bw_m <- betareg(data = liver_RGR_MIS, RGR ~ h_SRP110282_dg + h_SRP110282_bw)
SRP110282_pp_dg_ec_m <- betareg(data = liver_RGR_MIS, RGR ~ h_SRP110282_dg + h_SRP110282_ec)

SRP110282_pp_m <- betareg(data = liver_RGR_MIS, RGR ~ h_SRP110282_dg + h_SRP110282_bw + h_SRP110282_ec)

SRP110282_pp_dg_mis <- betareg(data = liver_RGR_MIS, MIS ~ h_SRP110282_dg)
SRP110282_pp_ec_mis <- betareg(data = liver_RGR_MIS, MIS ~ h_SRP110282_ec)
SRP110282_pp_bw_mis <- betareg(data = liver_RGR_MIS, MIS ~ h_SRP110282_bw)

SRP110282_pp_dg_bw_mis <- betareg(data = liver_RGR_MIS, MIS ~ h_SRP110282_dg + h_SRP110282_bw)
SRP110282_pp_dg_ec_mis <- betareg(data = liver_RGR_MIS, MIS ~ h_SRP110282_dg + h_SRP110282_ec)

SRP110282_pp_mis <- betareg(data = liver_RGR_MIS, MIS ~ h_SRP110282_dg + h_SRP110282_bw + h_SRP110282_ec)
### 

# SRP034011
liver_RGR_MIS <- left_join(liver_RGR_MIS, join_df)
SRP034011_pp_dg_m <- betareg(data = liver_RGR_MIS, RGR ~ h_SRP034011_dg)
SRP034011_pp_ec_m <- betareg(data = liver_RGR_MIS, RGR ~ h_SRP034011_ec)
SRP034011_pp_bw_m <- betareg(data = liver_RGR_MIS, RGR ~ h_SRP034011_bw)

SRP034011_pp_dg_bw_m <- betareg(data = liver_RGR_MIS, RGR ~ h_SRP034011_dg + h_SRP034011_bw)
SRP034011_pp_dg_ec_m <- betareg(data = liver_RGR_MIS, RGR ~ h_SRP034011_dg + h_SRP034011_ec)

SRP034011_pp_m <- betareg(data = liver_RGR_MIS, RGR ~ h_SRP034011_dg + h_SRP034011_bw + h_SRP034011_ec)

SRP034011_pp_dg_mis <- betareg(data = liver_RGR_MIS, MIS ~ h_SRP034011_dg)
SRP034011_pp_ec_mis <- betareg(data = liver_RGR_MIS, MIS ~ h_SRP034011_ec)
SRP034011_pp_bw_mis <- betareg(data = liver_RGR_MIS, MIS ~ h_SRP034011_bw)

SRP034011_pp_dg_bw_mis <- betareg(data = liver_RGR_MIS, MIS ~ h_SRP034011_dg + h_SRP034011_bw)
SRP034011_pp_dg_ec_mis <- betareg(data = liver_RGR_MIS, MIS ~ h_SRP034011_dg + h_SRP034011_ec)

SRP034011_pp_mis <- betareg(data = liver_RGR_MIS, MIS ~ h_SRP034011_dg + h_SRP034011_bw + h_SRP034011_ec)
### 

# SRP071199
liver_RGR_MIS <- left_join(liver_RGR_MIS, join_df)
SRP071199_pp_dg_m <- betareg(data = liver_RGR_MIS, RGR ~ h_SRP071199_dg)
SRP071199_pp_ec_m <- betareg(data = liver_RGR_MIS, RGR ~ h_SRP071199_ec)
SRP071199_pp_bw_m <- betareg(data = liver_RGR_MIS, RGR ~ h_SRP071199_bw)

SRP071199_pp_dg_bw_m <- betareg(data = liver_RGR_MIS, RGR ~ h_SRP071199_dg + h_SRP071199_bw)
SRP071199_pp_dg_ec_m <- betareg(data = liver_RGR_MIS, RGR ~ h_SRP071199_dg + h_SRP071199_ec)

SRP071199_pp_m <- betareg(data = liver_RGR_MIS, RGR ~ h_SRP071199_dg + h_SRP071199_bw + h_SRP071199_ec)

SRP071199_pp_dg_mis <- betareg(data = liver_RGR_MIS, MIS ~ h_SRP071199_dg)
SRP071199_pp_ec_mis <- betareg(data = liver_RGR_MIS, MIS ~ h_SRP071199_ec)
SRP071199_pp_bw_mis <- betareg(data = liver_RGR_MIS, MIS ~ h_SRP071199_bw)

SRP071199_pp_dg_bw_mis <- betareg(data = liver_RGR_MIS, MIS ~ h_SRP071199_dg + h_SRP071199_bw)
SRP071199_pp_dg_ec_mis <- betareg(data = liver_RGR_MIS, MIS ~ h_SRP071199_dg + h_SRP071199_ec)

SRP071199_pp_mis <- betareg(data = liver_RGR_MIS, MIS ~ h_SRP071199_dg + h_SRP071199_bw + h_SRP071199_ec)
###

# SRP126641
liver_RGR_MIS <- left_join(liver_RGR_MIS, join_df)
SRP126641_pp_dg_m <- betareg(data = liver_RGR_MIS, RGR ~ h_SRP126641_dg)
SRP126641_pp_ec_m <- betareg(data = liver_RGR_MIS, RGR ~ h_SRP126641_ec)
SRP126641_pp_bw_m <- betareg(data = liver_RGR_MIS, RGR ~ h_SRP126641_bw)

SRP126641_pp_dg_bw_m <- betareg(data = liver_RGR_MIS, RGR ~ h_SRP126641_dg + h_SRP126641_bw)
SRP126641_pp_dg_ec_m <- betareg(data = liver_RGR_MIS, RGR ~ h_SRP126641_dg + h_SRP126641_ec)

SRP126641_pp_m <- betareg(data = liver_RGR_MIS, RGR ~ h_SRP126641_dg + h_SRP126641_bw + h_SRP126641_ec)

SRP126641_pp_dg_mis <- betareg(data = liver_RGR_MIS, MIS ~ h_SRP126641_dg)
SRP126641_pp_ec_mis <- betareg(data = liver_RGR_MIS, MIS ~ h_SRP126641_ec)
SRP126641_pp_bw_mis <- betareg(data = liver_RGR_MIS, MIS ~ h_SRP126641_bw)

SRP126641_pp_dg_bw_mis <- betareg(data = liver_RGR_MIS, MIS ~ h_SRP126641_dg + h_SRP126641_bw)
SRP126641_pp_dg_ec_mis <- betareg(data = liver_RGR_MIS, MIS ~ h_SRP126641_dg + h_SRP126641_ec)

SRP126641_pp_mis <- betareg(data = liver_RGR_MIS, MIS ~ h_SRP126641_dg + h_SRP126641_bw + h_SRP126641_ec)
### 

# SRP096160
liver_RGR_MIS <- left_join(liver_RGR_MIS, join_df)
SRP096160_pp_dg_m <- betareg(data = liver_RGR_MIS, RGR ~ mo_SRP096160_dg)
SRP096160_pp_ec_m <- betareg(data = liver_RGR_MIS, RGR ~ mo_SRP096160_ec)
SRP096160_pp_bw_m <- betareg(data = liver_RGR_MIS, RGR ~ mo_SRP096160_bw)

SRP096160_pp_dg_bw_m <- betareg(data = liver_RGR_MIS, RGR ~ mo_SRP096160_dg + mo_SRP096160_bw)
SRP096160_pp_dg_ec_m <- betareg(data = liver_RGR_MIS, RGR ~ mo_SRP096160_dg + mo_SRP096160_ec)

SRP096160_pp_m <- betareg(data = liver_RGR_MIS, RGR ~ mo_SRP096160_dg + mo_SRP096160_bw + mo_SRP096160_ec)

SRP096160_pp_dg_mis <- betareg(data = liver_RGR_MIS, MIS ~ mo_SRP096160_dg)
SRP096160_pp_ec_mis <- betareg(data = liver_RGR_MIS, MIS ~ mo_SRP096160_ec)
SRP096160_pp_bw_mis <- betareg(data = liver_RGR_MIS, MIS ~ mo_SRP096160_bw)

SRP096160_pp_dg_bw_mis <- betareg(data = liver_RGR_MIS, MIS ~ mo_SRP096160_dg + mo_SRP096160_bw)
SRP096160_pp_dg_ec_mis <- betareg(data = liver_RGR_MIS, MIS ~ mo_SRP096160_dg + mo_SRP096160_ec)

SRP096160_pp_mis <- betareg(data = liver_RGR_MIS, MIS ~ mo_SRP096160_dg + mo_SRP096160_bw + mo_SRP096160_ec)
### 
#

# ERP020067
liver_RGR_MIS <- left_join(liver_RGR_MIS, join_df)
ERP020067_pp_dg_m <- betareg(data = liver_RGR_MIS, RGR ~ mo_ERP020067_dg)
ERP020067_pp_ec_m <- betareg(data = liver_RGR_MIS, RGR ~ mo_ERP020067_ec)
ERP020067_pp_bw_m <- betareg(data = liver_RGR_MIS, RGR ~ mo_ERP020067_bw)

ERP020067_pp_dg_bw_m <- betareg(data = liver_RGR_MIS, RGR ~ mo_ERP020067_dg + mo_ERP020067_bw)
ERP020067_pp_dg_ec_m <- betareg(data = liver_RGR_MIS, RGR ~ mo_ERP020067_dg + mo_ERP020067_ec)

ERP020067_pp_m <- betareg(data = liver_RGR_MIS, RGR ~ mo_ERP020067_dg + mo_ERP020067_bw + mo_ERP020067_ec)

ERP020067_pp_dg_mis <- betareg(data = liver_RGR_MIS, MIS ~ mo_ERP020067_dg)
ERP020067_pp_ec_mis <- betareg(data = liver_RGR_MIS, MIS ~ mo_ERP020067_ec)
ERP020067_pp_bw_mis <- betareg(data = liver_RGR_MIS, MIS ~ mo_ERP020067_bw)

ERP020067_pp_dg_bw_mis <- betareg(data = liver_RGR_MIS, MIS ~ mo_ERP020067_dg + mo_ERP020067_bw)
ERP020067_pp_dg_ec_mis <- betareg(data = liver_RGR_MIS, MIS ~ mo_ERP020067_dg + mo_ERP020067_ec)

ERP020067_pp_mis <- betareg(data = liver_RGR_MIS, MIS ~ mo_ERP020067_dg + mo_ERP020067_bw + mo_ERP020067_ec)
### 

# SRP018945
liver_RGR_MIS <- left_join(liver_RGR_MIS, join_df)
SRP018945_pp_dg_m <- betareg(data = liver_RGR_MIS, RGR ~ m_SRP018945_dg)
SRP018945_pp_ec_m <- betareg(data = liver_RGR_MIS, RGR ~ m_SRP018945_ec)
SRP018945_pp_bw_m <- betareg(data = liver_RGR_MIS, RGR ~ m_SRP018945_bw)

SRP018945_pp_dg_bw_m <- betareg(data = liver_RGR_MIS, RGR ~ m_SRP018945_dg + m_SRP018945_bw)
SRP018945_pp_dg_ec_m <- betareg(data = liver_RGR_MIS, RGR ~ m_SRP018945_dg + m_SRP018945_ec)

SRP018945_pp_m <- betareg(data = liver_RGR_MIS, RGR ~ m_SRP018945_dg + m_SRP018945_bw + m_SRP018945_ec)

SRP018945_pp_dg_mis <- betareg(data = liver_RGR_MIS, MIS ~ m_SRP018945_dg)
SRP018945_pp_ec_mis <- betareg(data = liver_RGR_MIS, MIS ~ m_SRP018945_ec)
SRP018945_pp_bw_mis <- betareg(data = liver_RGR_MIS, MIS ~ m_SRP018945_bw)

SRP018945_pp_dg_bw_mis <- betareg(data = liver_RGR_MIS, MIS ~ m_SRP018945_dg + m_SRP018945_bw)
SRP018945_pp_dg_ec_mis <- betareg(data = liver_RGR_MIS, MIS ~ m_SRP018945_dg + m_SRP018945_ec)

SRP018945_pp_mis <- betareg(data = liver_RGR_MIS, MIS ~ m_SRP018945_dg + m_SRP018945_bw + m_SRP018945_ec)
### 

# SRP131855
liver_RGR_MIS <- left_join(liver_RGR_MIS, join_df)
SRP131855_pp_dg_m <- betareg(data = liver_RGR_MIS, RGR ~ m_SRP131855_dg)
SRP131855_pp_ec_m <- betareg(data = liver_RGR_MIS, RGR ~ m_SRP131855_ec)
SRP131855_pp_bw_m <- betareg(data = liver_RGR_MIS, RGR ~ m_SRP131855_bw)

SRP131855_pp_dg_bw_m <- betareg(data = liver_RGR_MIS, RGR ~ m_SRP131855_dg + m_SRP131855_bw)
SRP131855_pp_dg_ec_m <- betareg(data = liver_RGR_MIS, RGR ~ m_SRP131855_dg + m_SRP131855_ec)

SRP131855_pp_m <- betareg(data = liver_RGR_MIS, RGR ~ m_SRP131855_dg + m_SRP131855_bw + m_SRP131855_ec)

SRP131855_pp_dg_mis <- betareg(data = liver_RGR_MIS, MIS ~ m_SRP131855_dg)
SRP131855_pp_ec_mis <- betareg(data = liver_RGR_MIS, MIS ~ m_SRP131855_ec)
SRP131855_pp_bw_mis <- betareg(data = liver_RGR_MIS, MIS ~ m_SRP131855_bw)

SRP131855_pp_dg_bw_mis <- betareg(data = liver_RGR_MIS, MIS ~ m_SRP131855_dg + m_SRP131855_bw)
SRP131855_pp_dg_ec_mis <- betareg(data = liver_RGR_MIS, MIS ~ m_SRP131855_dg + m_SRP131855_ec)

SRP131855_pp_mis <- betareg(data = liver_RGR_MIS, MIS ~ m_SRP131855_dg + m_SRP131855_bw + m_SRP131855_ec)
###

# SRP150689
liver_RGR_MIS <- left_join(liver_RGR_MIS, join_df)
SRP150689_pp_dg_m <- betareg(data = liver_RGR_MIS, RGR ~ m_SRP150689_dg)
SRP150689_pp_ec_m <- betareg(data = liver_RGR_MIS, RGR ~ m_SRP150689_ec)
SRP150689_pp_bw_m <- betareg(data = liver_RGR_MIS, RGR ~ m_SRP150689_bw)

SRP150689_pp_dg_bw_m <- betareg(data = liver_RGR_MIS, RGR ~ m_SRP150689_dg + m_SRP150689_bw)
SRP150689_pp_dg_ec_m <- betareg(data = liver_RGR_MIS, RGR ~ m_SRP150689_dg + m_SRP150689_ec)

SRP150689_pp_m <- betareg(data = liver_RGR_MIS, RGR ~ m_SRP150689_dg + m_SRP150689_bw + m_SRP150689_ec)

SRP150689_pp_dg_mis <- betareg(data = liver_RGR_MIS, MIS ~ m_SRP150689_dg)
SRP150689_pp_ec_mis <- betareg(data = liver_RGR_MIS, MIS ~ m_SRP150689_ec)
SRP150689_pp_bw_mis <- betareg(data = liver_RGR_MIS, MIS ~ m_SRP150689_bw)

SRP150689_pp_dg_bw_mis <- betareg(data = liver_RGR_MIS, MIS ~ m_SRP150689_dg + m_SRP150689_bw)
SRP150689_pp_dg_ec_mis <- betareg(data = liver_RGR_MIS, MIS ~ m_SRP150689_dg + m_SRP150689_ec)

SRP150689_pp_mis <- betareg(data = liver_RGR_MIS, MIS ~ m_SRP150689_dg + m_SRP150689_bw + m_SRP150689_ec)
###

# SRP171171
liver_RGR_MIS <- left_join(liver_RGR_MIS, join_df)
SRP171171_pp_dg_m <- betareg(data = liver_RGR_MIS, RGR ~ m_SRP171171_dg)
SRP171171_pp_ec_m <- betareg(data = liver_RGR_MIS, RGR ~ m_SRP171171_ec)
SRP171171_pp_bw_m <- betareg(data = liver_RGR_MIS, RGR ~ m_SRP171171_bw)

SRP171171_pp_dg_bw_m <- betareg(data = liver_RGR_MIS, RGR ~ m_SRP171171_dg + m_SRP171171_bw)
SRP171171_pp_dg_ec_m <- betareg(data = liver_RGR_MIS, RGR ~ m_SRP171171_dg + m_SRP171171_ec)

SRP171171_pp_m <- betareg(data = liver_RGR_MIS, RGR ~ m_SRP171171_dg + m_SRP171171_bw + m_SRP171171_ec)

SRP171171_pp_dg_mis <- betareg(data = liver_RGR_MIS, MIS ~ m_SRP171171_dg)
SRP171171_pp_ec_mis <- betareg(data = liver_RGR_MIS, MIS ~ m_SRP171171_ec)
SRP171171_pp_bw_mis <- betareg(data = liver_RGR_MIS, MIS ~ m_SRP171171_bw)

SRP171171_pp_dg_bw_mis <- betareg(data = liver_RGR_MIS, MIS ~ m_SRP171171_dg + m_SRP171171_bw)
SRP171171_pp_dg_ec_mis <- betareg(data = liver_RGR_MIS, MIS ~ m_SRP171171_dg + m_SRP171171_ec)

SRP171171_pp_mis <- betareg(data = liver_RGR_MIS, MIS ~ m_SRP171171_dg + m_SRP171171_bw + m_SRP171171_ec)
###

# SRP261098
liver_RGR_MIS <- left_join(liver_RGR_MIS, join_df)
SRP261098_pp_dg_m <- betareg(data = liver_RGR_MIS, RGR ~ m_SRP261098_dg)
SRP261098_pp_ec_m <- betareg(data = liver_RGR_MIS, RGR ~ m_SRP261098_ec)
SRP261098_pp_bw_m <- betareg(data = liver_RGR_MIS, RGR ~ m_SRP261098_bw)

SRP261098_pp_dg_bw_m <- betareg(data = liver_RGR_MIS, RGR ~ m_SRP261098_dg + m_SRP261098_bw)
SRP261098_pp_dg_ec_m <- betareg(data = liver_RGR_MIS, RGR ~ m_SRP261098_dg + m_SRP261098_ec)

SRP261098_pp_m <- betareg(data = liver_RGR_MIS, RGR ~ m_SRP261098_dg + m_SRP261098_bw + m_SRP261098_ec)

SRP261098_pp_dg_mis <- betareg(data = liver_RGR_MIS, MIS ~ m_SRP261098_dg)
SRP261098_pp_ec_mis <- betareg(data = liver_RGR_MIS, MIS ~ m_SRP261098_ec)
SRP261098_pp_bw_mis <- betareg(data = liver_RGR_MIS, MIS ~ m_SRP261098_bw)

SRP261098_pp_dg_bw_mis <- betareg(data = liver_RGR_MIS, MIS ~ m_SRP261098_dg + m_SRP261098_bw)
SRP261098_pp_dg_ec_mis <- betareg(data = liver_RGR_MIS, MIS ~ m_SRP261098_dg + m_SRP261098_ec)

SRP261098_pp_mis <- betareg(data = liver_RGR_MIS, MIS ~ m_SRP261098_dg + m_SRP261098_bw + m_SRP261098_ec)
###