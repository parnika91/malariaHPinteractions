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
save(RGR_MIS_cleaned, file = "Unweighted_RGR_MIS_cleaned.RData")

### models ###

# ov models
ov_pp_dg_m <- betareg(data = RGR_MIS_cleaned, RGR ~ ov_pp_dg)
ov_pp_ec_m <- betareg(data = RGR_MIS_cleaned, RGR ~ ov_pp_ec)
ov_pp_bw_m <- betareg(data = RGR_MIS_cleaned, RGR ~ ov_pp_bw)

ov_pp_dg_bw_m <- betareg(data = RGR_MIS_cleaned, RGR ~ ov_pp_dg + ov_pp_bw)
ov_pp_dg_ec_m <- betareg(data = RGR_MIS_cleaned, RGR ~ ov_pp_dg + ov_pp_ec)

ov_pp_m <- betareg(data = RGR_MIS_cleaned, RGR ~ ov_pp_dg + ov_pp_bw + ov_pp_ec)

ov_bp_dg_m <- betareg(data = RGR_MIS_cleaned, RGR ~ ov_bp_dg)
ov_bp_ec_m <- betareg(data = RGR_MIS_cleaned, RGR ~ ov_bp_ec)
ov_bp_bw_m <- betareg(data = RGR_MIS_cleaned, RGR ~ ov_bp_bw)

ov_bp_dg_bw_m <- betareg(data = RGR_MIS_cleaned, RGR ~ ov_bp_dg + ov_bp_bw)
ov_bp_dg_ec_m <- betareg(data = RGR_MIS_cleaned, RGR ~ ov_bp_dg + ov_bp_ec)

ov_bp_m <- betareg(data = RGR_MIS_cleaned, RGR ~ ov_bp_dg + ov_bp_bw + ov_bp_ec)

ov_all_dg_m <- betareg(data = RGR_MIS_cleaned, RGR ~ ov_all_dg)
ov_all_ec_m <- betareg(data = RGR_MIS_cleaned, RGR ~ ov_all_ec)
ov_all_bw_m <- betareg(data = RGR_MIS_cleaned, RGR ~ ov_all_bw)

ov_all_dg_bw_m <- betareg(data = RGR_MIS_cleaned, RGR ~ ov_all_dg + ov_all_bw)
ov_all_dg_ec_m <- betareg(data = RGR_MIS_cleaned, RGR ~ ov_all_dg + ov_all_ec)

ov_all_m <- betareg(data = RGR_MIS_cleaned, RGR ~ ov_all_dg + ov_all_bw + ov_all_ec)

# pb models

pb_pp_dg_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pb_pp_dg)
pb_pp_ec_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pb_pp_ec)
pb_pp_bw_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pb_pp_bw)

pb_pp_dg_bw_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pb_pp_dg + pb_pp_bw)
pb_pp_dg_ec_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pb_pp_dg + pb_pp_ec)

pb_pp_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pb_pp_dg + pb_pp_bw + pb_pp_ec)

pb_bp_dg_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pb_bp_dg)
pb_bp_ec_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pb_bp_ec)
pb_bp_bw_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pb_bp_bw)

pb_bp_dg_bw_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pb_bp_dg + pb_bp_bw)
pb_bp_dg_ec_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pb_bp_dg + pb_bp_ec)

pb_bp_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pb_bp_dg + pb_bp_bw + pb_bp_ec)

pb_all_dg_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pb_all_dg)
pb_all_ec_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pb_all_ec)
pb_all_bw_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pb_all_bw)

pb_all_dg_bw_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pb_all_dg + pb_all_bw)
pb_all_dg_ec_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pb_all_dg + pb_all_ec)

pb_all_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pb_all_dg + pb_all_bw + pb_all_ec)

# pf models

pf_pp_dg_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pf_pp_dg)
pf_pp_ec_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pf_pp_ec)
pf_pp_bw_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pf_pp_bw)

pf_pp_dg_bw_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pf_pp_dg + pf_pp_bw)
pf_pp_dg_ec_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pf_pp_dg + pf_pp_ec)

pf_pp_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pf_pp_dg + pf_pp_bw + pf_pp_ec)

pf_bp_dg_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pf_bp_dg)
pf_bp_ec_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pf_bp_ec)
pf_bp_bw_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pf_bp_bw)

pf_bp_dg_bw_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pf_bp_dg + pf_bp_bw)
pf_bp_dg_ec_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pf_bp_dg + pf_bp_ec)

pf_bp_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pf_bp_dg + pf_bp_bw + pf_bp_ec)

pf_all_dg_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pf_all_dg)
pf_all_ec_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pf_all_ec)
pf_all_bw_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pf_all_bw)

pf_all_dg_bw_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pf_all_dg + pf_all_bw)
pf_all_dg_ec_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pf_all_dg + pf_all_ec)

pf_all_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pf_all_dg + pf_all_bw + pf_all_ec)

# pfE models

pfE_pp_dg_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pfE_pp_dg)
pfE_pp_ec_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pfE_pp_ec)
pfE_pp_bw_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pfE_pp_bw)

pfE_pp_dg_bw_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pfE_pp_dg + pfE_pp_bw)
pfE_pp_dg_ec_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pfE_pp_dg + pfE_pp_ec)

pfE_pp_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pfE_pp_dg + pfE_pp_bw + pfE_pp_ec)

pfE_bp_dg_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pfE_bp_dg)
pfE_bp_ec_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pfE_bp_ec)
pfE_bp_bw_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pfE_bp_bw)

pfE_bp_dg_bw_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pfE_bp_dg + pfE_bp_bw)
pfE_bp_dg_ec_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pfE_bp_dg + pfE_bp_ec)

pfE_bp_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pfE_bp_dg + pfE_bp_bw + pfE_bp_ec)

pfE_all_dg_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pfE_all_dg)
pfE_all_ec_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pfE_all_ec)
pfE_all_bw_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pfE_all_bw)

pfE_all_dg_bw_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pfE_all_dg + pfE_all_bw)
pfE_all_dg_ec_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pfE_all_dg + pfE_all_ec)

pfE_all_m <- betareg(data = RGR_MIS_cleaned, RGR ~ pfE_all_dg + pfE_all_bw + pfE_all_ec)


# MIS models

# ov models
ov_pp_dg_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ ov_pp_dg)
ov_pp_ec_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ ov_pp_ec)
ov_pp_bw_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ ov_pp_bw)

ov_pp_dg_bw_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ ov_pp_dg + ov_pp_bw)
ov_pp_dg_ec_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ ov_pp_dg + ov_pp_ec)

ov_pp_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ ov_pp_dg + ov_pp_bw + ov_pp_ec)

ov_bp_dg_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ ov_bp_dg)
ov_bp_ec_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ ov_bp_ec)
ov_bp_bw_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ ov_bp_bw)

ov_bp_dg_bw_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ ov_bp_dg + ov_bp_bw)
ov_bp_dg_ec_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ ov_bp_dg + ov_bp_ec)

ov_bp_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ ov_bp_dg + ov_bp_bw + ov_bp_ec)

ov_all_dg_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ ov_all_dg)
ov_all_ec_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ ov_all_ec)
ov_all_bw_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ ov_all_bw)

ov_all_dg_bw_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ ov_all_dg + ov_all_bw)
ov_all_dg_ec_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ ov_all_dg + ov_all_ec)

ov_all_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ ov_all_dg + ov_all_bw + ov_all_ec)

# pb models

pb_pp_dg_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pb_pp_dg)
pb_pp_ec_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pb_pp_ec)
pb_pp_bw_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pb_pp_bw)

pb_pp_dg_bw_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pb_pp_dg + pb_pp_bw)
pb_pp_dg_ec_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pb_pp_dg + pb_pp_ec)

pb_pp_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pb_pp_dg + pb_pp_bw + pb_pp_ec)

pb_bp_dg_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pb_bp_dg)
pb_bp_ec_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pb_bp_ec)
pb_bp_bw_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pb_bp_bw)

pb_bp_dg_bw_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pb_bp_dg + pb_bp_bw)
pb_bp_dg_ec_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pb_bp_dg + pb_bp_ec)

pb_bp_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pb_bp_dg + pb_bp_bw + pb_bp_ec)

pb_all_dg_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pb_all_dg)
pb_all_ec_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pb_all_ec)
pb_all_bw_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pb_all_bw)

pb_all_dg_bw_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pb_all_dg + pb_all_bw)
pb_all_dg_ec_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pb_all_dg + pb_all_ec)

pb_all_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pb_all_dg + pb_all_bw + pb_all_ec)

# pf models

pf_pp_dg_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pf_pp_dg)
pf_pp_ec_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pf_pp_ec)
pf_pp_bw_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pf_pp_bw)

pf_pp_dg_bw_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pf_pp_dg + pf_pp_bw)
pf_pp_dg_ec_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pf_pp_dg + pf_pp_ec)

pf_pp_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pf_pp_dg + pf_pp_bw + pf_pp_ec)

pf_bp_dg_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pf_bp_dg)
pf_bp_ec_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pf_bp_ec)
pf_bp_bw_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pf_bp_bw)

pf_bp_dg_bw_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pf_bp_dg + pf_bp_bw)
pf_bp_dg_ec_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pf_bp_dg + pf_bp_ec)

pf_bp_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pf_bp_dg + pf_bp_bw + pf_bp_ec)

pf_all_dg_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pf_all_dg)
pf_all_ec_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pf_all_ec)
pf_all_bw_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pf_all_bw)

pf_all_dg_bw_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pf_all_dg + pf_all_bw)
pf_all_dg_ec_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pf_all_dg + pf_all_ec)

pf_all_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pf_all_dg + pf_all_bw + pf_all_ec)

# pfE models

pfE_pp_dg_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pfE_pp_dg)
pfE_pp_ec_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pfE_pp_ec)
pfE_pp_bw_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pfE_pp_bw)

pfE_pp_dg_bw_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pfE_pp_dg + pfE_pp_bw)
pfE_pp_dg_ec_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pfE_pp_dg + pfE_pp_ec)

pfE_pp_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pfE_pp_dg + pfE_pp_bw + pfE_pp_ec)

pfE_bp_dg_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pfE_bp_dg)
pfE_bp_ec_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pfE_bp_ec)
pfE_bp_bw_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pfE_bp_bw)

pfE_bp_dg_bw_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pfE_bp_dg + pfE_bp_bw)
pfE_bp_dg_ec_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pfE_bp_dg + pfE_bp_ec)

pfE_bp_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pfE_bp_dg + pfE_bp_bw + pfE_bp_ec)

pfE_all_dg_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pfE_all_dg)
pfE_all_ec_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pfE_all_ec)
pfE_all_bw_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pfE_all_bw)

pfE_all_dg_bw_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pfE_all_dg + pfE_all_bw)
pfE_all_dg_ec_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pfE_all_dg + pfE_all_ec)

pfE_all_mis <- betareg(data = RGR_MIS_cleaned, MIS ~ pfE_all_dg + pfE_all_bw + pfE_all_ec)
