library(sna)
library(bipartite)
library(dplyr)
library(reshape2)

load("df_concat_allhosts/cor/df_concat_allhosts_all_bipartite.RData")

# nrow(bipartite)
#[1] 8320187
# length(unique(bipartite[,1]))
#[1] 13217
# length(unique(bipartite[,2]))
#[1] 4005
# 13217*4005
#[1] 52934085

bipartite$interactions <- rep(1, nrow(bipartite))
bipartitedc <- bipartite %>%
  dcast(df_concat_allhosts_all_gene2 ~ df_concat_allhosts_all_gene1) %>%
  tibble::column_to_rownames("df_concat_allhosts_all_gene2")
# Uses interactions as value column: use value.var to override.
bipartitedc[is.na(bipartitedc)] <- 0

onemode <- bipartite::as.one.mode(bipartitedc, project = "higher")
# str(onemode)
# num [1:17222, 1:17222] 0 0 0 0 0 0 0 0 0 0 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:17222] "p_OG0000293" "p_OG0000294" "p_OG0000295" "p_OG0000296" ...
# ..$ : chr [1:17222] "p_OG0000293" "p_OG0000294" "p_OG0000295" "p_OG0000296" ...

pdf("visweb.pdf", onefile = T)
bipartite::visweb(onemode, type="diagonal")

degree <- sna::degree(onemode, g=1, 
                 gmode="graph", 
                 diag=T,tmaxdev=FALSE, cmode="freeman", 
                 rescale=FALSE, ignore.eval=FALSE)

cl <- sna::closeness(onemode, g=1, 
          gmode="graph", diag=F,
          tmaxdev=FALSE, cmode="undirected",
          geodist.precomp=NULL,rescale=FALSE, 
          ignore.eval=TRUE)
save(cl, file = "cl_sna.RData")

bw <- betweenness(onemode, g=1, 
                  gmode="graph", 
                  diag=FALSE,tmaxdev=FALSE, cmode="undirected", 
                  geodist.precomp=NULL,rescale=FALSE, ignore.eval=TRUE)
save(bw, file = "bw_sna.RData")