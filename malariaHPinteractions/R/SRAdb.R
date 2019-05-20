# libraries
source("https://bioconductor.org/biocLite.R")
biocLite("SRAdb")
source("https://bioconductor.org/biocLite.R")
biocLite("GEOquery")
install.packages("pkgconfig")

library(SRAdb)
library(GEOquery)
library(pkgconfig)

# SRAdb first step from manual 
sqlfile <- 'SRAmetadb.sqlite'

# first run the function from https://github.com/seandavi/SRAdb/blob/master/R/getSRAdbFile.R
# otherwise the function getSRAdbFile() may not work - There were problems with download.file(...)


getSRAdbFile <-
  function (destdir = getwd(), destfile = "SRAmetadb.sqlite.gz",
            method)
  {
    if (missing(method))
      method <- ifelse(!is.null(getOption("download.file.method")),
                       getOption("download.file.method"), "auto")
    localfile <- file.path(destdir, destfile)   
    options(warn=-1)
    
    url_sra_1 = 'https://s3.amazonaws.com/starbuck1/sradb/SRAmetadb.sqlite.gz'
    url_sra_2 = 'https://gbnci-abcc.ncifcrf.gov/backup/SRAmetadb.sqlite.gz'
    
    if(! inherits(try(url(url_sra_1, open='rb'), silent = TRUE), "try-error") ) {
      url_sra = url_sra_1
    } else {
      url_sra = url_sra_2
    }
    
    download.file(url_sra, destfile = localfile, mode = "wb", method=method)
    message("Unzipping...\n")
    gunzip(localfile, overwrite = TRUE)
    unzippedlocalfile <- gsub("[.]gz$", "", localfile)
    con <- dbConnect(SQLite(), unzippedlocalfile)
    dat <- dbGetQuery(con, "SELECT * FROM metaInfo")
    dbDisconnect(con)
    message("Metadata associate with downloaded file:\n")
    message(dat)
    return(unzippedlocalfile)
  }

# download SRAdb file - 32.3 GB
if(!file.exists('SRAmetadb.sqlite')) sqlfile <<- getSRAdbFile()

file.info('SRAmetadb.sqlite')
# size isdir mode               mtime               ctime               atime  uid
# SRAmetadb.sqlite 34697342976 FALSE  644 2017-10-02 16:42:09 2017-10-02 16:42:09 2017-10-02 16:42:09 1009
# gid   uname  grname
# SRAmetadb.sqlite 1011 parnika parnika

# SRAmetadb.sqlite 35606409216 FALSE  644 2017-11-30 11:30:19 2017-11-30 11:30:19 2017-11-30 11:30:21 1009 1011
# uname  grname
# SRAmetadb.sqlite parnika parnika

# SRAmetadb.sqlite 37538463744 FALSE  644 2018-08-07 11:02:20 2018-08-07 11:02:20 2018-08-07 11:02:21 1009 1011 
#                   grname
# SRAmetadb.sqlite parnika
# size isdir mode               mtime               ctime               atime  uid  gid   uname  grname
# SRAmetadb.sqlite 35623591936 FALSE  644 2019-01-22 11:51:24 2019-01-22 11:51:24 2019-01-22 11:54:26 1009 1011 parnika parnika

sra_con <- dbConnect(SQLite(),sqlfile)

(sra_tables <- dbListTables(sra_con))
# [1] "col_desc"        "experiment"      "fastq"           "metaInfo"        "run"            
# [6] "sample"          "sra"             "sra_ft"          "sra_ft_content"  "sra_ft_segdir"  
# [11] "sra_ft_segments" "study"           "submission" 


dbListFields(sra_con,"study")
# [1] "study_ID"             "study_alias"          "study_accession"      "study_title"         
# [5] "study_type"           "study_abstract"       "broker_name"          "center_name"         
# [9] "center_project_name"  "study_description"    "related_studies"      "primary_study"       
# [13] "sra_link"             "study_url_link"       "xref_link"            "study_entrez_link"   
# [17] "ddbj_link"            "ena_link"             "study_attribute"      "submission_accession"
# [21] "sradb_updated"  


getTableCounts <- function(tableName,conn) {
  sql <- sprintf("select count(*) from %s",tableName)
  return(dbGetQuery(conn,sql)[1,1])
}

sapply(sra_tables,getTableCounts, sra_con, simplify=FALSE)
# $col_desc
# [1] 129
# 
# $experiment
# [1] 3326631
# 
# $fastq
# [1] 3434946
# 
# $metaInfo
# [1] 2
# 
# $run
# [1] 3732848
# 
# $sample
# [1] 2962581
# 
# $sra
# [1] 3679805
# 
# $sra_ft
# [1] 3679805
# 
# $sra_ft_content
# [1] 3679805
# 
# $sra_ft_segdir
# [1] 34
# 
# $sra_ft_segments
# [1] 297704
# 
# $study
# [1] 116014
# 
# $submission
# [1] 748819


do.call(rbind,sapply(sra_tables,getTableCounts, sra_con, simplify=FALSE)) # all the databases in SRAdb
# [,1]
# col_desc            129
# experiment      3326631
# fastq           3434946
# metaInfo              2
# run             3732848
# sample          2962581
# sra             3679805
# sra_ft          3679805
# sra_ft_content  3679805
# sra_ft_segdir        34
# sra_ft_segments  297704
# study            116014
# submission       748819

#rs <- dbGetQuery(sra_con, paste( "SELECT library_strategy AS 'Library Strategy', count( * ) AS Runs FROM `experiment` GROUP BY library_strategy order by Runs DESC", sep=""))

rs <- dbGetQuery(sra_con, paste( "SELECT library_strategy AS 'Library Strategy', count( * ) AS Runs, study_title AS 'Study Title' FROM `sra` WHERE library_strategy LIKE '%RNA%'AND study_title LIKE '%Plasmodium%' GROUP BY library_strategy order by Runs DESC", sep=""))

rs <- dbGetQuery(sra_con, paste( "SELECT * FROM `study`", sep=""))
rs <- dbGetQuery(sra_con, paste( "SELECT * FROM `col_desc`", sep="")) # overall view of the data : column description - submission, study, etc
rs <- dbGetQuery(sra_con, paste( "SELECT * FROM `experiment`", sep=""))
rs <- dbGetQuery(sra_con, paste( "SELECT * FROM `fastq`", sep=""))
rs <- dbGetQuery(sra_con, paste( "SELECT * FROM `metaInfo`", sep=""))
rs <- dbGetQuery(sra_con, paste( "SELECT * FROM `run`", sep=""))
rs <- dbGetQuery(sra_con, paste( "SELECT * FROM `sample`", sep=""))
rs <- dbGetQuery(sra_con, paste( "SELECT * FROM `sra`", sep=""))
rs <- dbGetQuery(sra_con, paste( "SELECT * FROM `sra_ft`", sep=""))
rs <- dbGetQuery(sra_con, paste( "SELECT * FROM `sra_ft_content`", sep=""))
rs <- dbGetQuery(sra_con, paste( "SELECT * FROM `sra_ft_segdir`", sep=""))
rs <- dbGetQuery(sra_con, paste( "SELECT * FROM `sra_ft_segment`", sep=""))
rs <- dbGetQuery(sra_con, paste( "SELECT * FROM `submission`", sep=""))

# Visualization
#
library(SRAdb)
library(Rgraphviz)
g <- sraGraph('primary thyroid cell line'  , sra_con)
attrs <- getDefaultAttrs(list(node=list(fillcolor='lightblue', shape='ellipse')))
plot(g, attrs=attrs)
## similiar search as the above, returned much larger data.frame and graph is too clouded
g <- sraGraph('Ewing Sarcoma', sra_con)
plot(g)

# has 8 library strategies
rs <- dbGetQuery(sra_con, paste( "SELECT library_strategy AS 'Library Strategy', count( * ) AS Runs, study_title AS 'Study Title' FROM `sra` WHERE library_strategy LIKE '%RNA%'AND study_title LIKE '%Plasmodium%' OR study_title LIKE '%Malaria%' GROUP BY library_strategy order by Runs DESC", sep=""))

# has RNA-seq, miRNA-seq and ncRNA-seq
rs <- dbGetQuery(sra_con, paste( "SELECT library_strategy AS 'Library Strategy', count( * ) AS Runs, study_title AS 'Study Title' FROM `sra` WHERE library_strategy LIKE '%RNA%'AND study_title LIKE '%Plasmodium%' order by Runs DESC", sep=""))

# has 8 library strategies ======> all columns
rs <- dbGetQuery(sra_con, paste( "SELECT *, count( * ) AS Runs FROM `sra` WHERE library_strategy LIKE '%RNA%'AND study_title LIKE '%Plasmodium%' OR study_title LIKE '%Malaria%' GROUP BY library_strategy order by Runs DESC", sep=""))

# only RNA-seq experiments
rs <- dbGetQuery(sra_con, paste( "SELECT *, count( * ) AS Runs FROM `sra` WHERE library_strategy LIKE 'RNA%' AND study_title LIKE '%Plasmodium%' OR study_title LIKE '%Malaria%' order by Runs DESC", sep=""))

# 23 experiments, 1427 runs: falciparum/malaria 
rs_Pf <- dbGetQuery(sra_con, paste( "SELECT *, count(*) AS Runs FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND ((study_abstract LIKE '%falciparum%') AND ((study_title LIKE '%malaria%') OR study_title LIKE '%falciparum%')) GROUP BY study_name", sep=""))

# RNA-Seq experiments with the words "human" and "falciparum" in study_abstract, grouped by study_names (experiment names) and with number of runs in each experiment
# 10 different experiments (study_name == "NA" counted as 1) with 300 runs
rs_Pf_human <- dbGetQuery(sra_con, paste( "SELECT *, count(*) AS Runs FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND study_abstract LIKE '%human%' AND study_abstract LIKE '%falciparum%' GROUP BY study_name", sep=""))
rs_Pf_human <- dbGetQuery(sra_con, paste( "SELECT *, count(*) AS Runs FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND study_abstract LIKE '%human%' AND study_title NOT LIKE '%mosquito%' AND study_abstract LIKE '%falciparum%' AND ((((study_abstract LIKE '%blood%') OR study_abstract LIKE '%erythrocyt%') OR study_abstract LIKE '%liver%') OR study_abstract LIKE '%hepatocytes%') AND ((((study_abstract NOT LIKE '%Anopheles%') AND study_abstract NOT LIKE '%gambiae%') AND study_abstract NOT LIKE '%stephensi%') AND study_abstract NOT LIKE '%midgut%') GROUP BY study_name", sep=""))


Plasmodium_host_table <- data.frame()

# 7 expt in human blood falciparum
rs_Pf_human_blood <- dbGetQuery(sra_con, paste( "SELECT *, count(*) AS Runs FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND study_abstract LIKE '%human%' AND study_abstract LIKE '%falciparum%' AND study_abstract LIKE '%blood%' GROUP BY study_name", sep=""))

rs_Pf_human_blood <- dbGetQuery(sra_con, paste( "SELECT *, count(*) AS Runs FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND study_title LIKE '%falciparum%' AND study_title NOT LIKE '%mosquito%' AND ((study_abstract LIKE '%blood%') OR study_abstract LIKE '%erythrocyt%') AND study_abstract LIKE '%human%' AND study_abstract LIKE '%falciparum%' GROUP BY study_name", sep=""))


# Categorised as malaria, plasmodium, anopheles and mosquito
rs_title_malaria <- dbGetQuery(sra_con, paste( "SELECT *, count(*) AS Runs FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND study_title LIKE '%malaria%' GROUP BY study_name", sep=""))
rs_title_plasmodium <- dbGetQuery(sra_con, paste( "SELECT *, count(*) AS Runs FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND study_title LIKE '%plasmodium%' GROUP BY study_name", sep=""))
rs_title_anopheles <- dbGetQuery(sra_con, paste( "SELECT *, count(*) AS Runs FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND study_title LIKE '%anopheles%' GROUP BY study_name", sep=""))
rs_title_mosquito <- dbGetQuery(sra_con, paste( "SELECT *, count(*) AS Runs FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND study_title LIKE '%mosquito%' GROUP BY study_name", sep=""))
rs_abstract_mosquito <- dbGetQuery(sra_con, paste( "SELECT *, count(*) AS Runs FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND study_abstract LIKE '%mosquito%' GROUP BY study_name", sep=""))
rs_abstract_anopheles <- dbGetQuery(sra_con, paste( "SELECT *, count(*) AS Runs FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND study_abstract LIKE '%anopheles%' GROUP BY study_name", sep=""))
rs_abstract_plasmodium <- dbGetQuery(sra_con, paste( "SELECT *, count(*) AS Runs FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND study_abstract LIKE '%plasmodium%' GROUP BY study_name", sep=""))
rs_abstract_malaria <- dbGetQuery(sra_con, paste( "SELECT *, count(*) AS Runs FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND study_abstract LIKE '%malaria%' GROUP BY study_name", sep=""))
rs_exptitle_malaria <- dbGetQuery(sra_con, paste( "SELECT *, count(*) AS Runs FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND experiment_title LIKE '%malaria%' GROUP BY study_name", sep=""))
rs_exptitle_plasmodium <- dbGetQuery(sra_con, paste( "SELECT *, count(*) AS Runs FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND experiment_title LIKE '%plasmodium%' GROUP BY study_name", sep=""))
rs_exptitle_anopheles <- dbGetQuery(sra_con, paste( "SELECT *, count(*) AS Runs FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND experiment_title LIKE '%anopheles%' GROUP BY study_name", sep=""))
rs_exptitle_mosquito <- dbGetQuery(sra_con, paste( "SELECT *, count(*) AS Runs FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND experiment_title LIKE '%mosquito%' GROUP BY study_name", sep=""))
rs_sample_plasmodium <- dbGetQuery(sra_con, paste( "SELECT *, count(*) AS Runs FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND sample_attribute LIKE '%plasmodium%' GROUP BY study_name", sep=""))
rs_sample_malaria <- dbGetQuery(sra_con, paste( "SELECT *, count(*) AS Runs FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND sample_attribute LIKE '%malaria%' GROUP BY study_name", sep=""))
rs_sample_anopheles <- dbGetQuery(sra_con, paste( "SELECT *, count(*) AS Runs FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND sample_attribute LIKE '%anopheles%' GROUP BY study_name", sep=""))
rs_sample_mosquito <- dbGetQuery(sra_con, paste( "SELECT *, count(*) AS Runs FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND sample_attribute LIKE '%mosquito%' GROUP BY study_name", sep=""))

# Same as above, but grouped by study_alias. Same number of total runs, but difference in number of experiments

rs_alias_title_malaria <- dbGetQuery(sra_con, paste( "SELECT *, count(*) AS Runs FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND study_title LIKE '%malaria%' GROUP BY study_alias", sep=""))
rs_alias_title_plasmodium <- dbGetQuery(sra_con, paste( "SELECT *, count(*) AS Runs FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND study_title LIKE '%plasmodium%' GROUP BY study_alias", sep=""))
rs_alias_title_anopheles <- dbGetQuery(sra_con, paste( "SELECT *, count(*) AS Runs FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND study_title LIKE '%anopheles%' GROUP BY study_alias", sep=""))
rs_alias_title_mosquito <- dbGetQuery(sra_con, paste( "SELECT *, count(*) AS Runs FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND study_title LIKE '%mosquito%' GROUP BY study_alias", sep=""))
rs_alias_abstract_mosquito <- dbGetQuery(sra_con, paste( "SELECT *, count(*) AS Runs FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND study_abstract LIKE '%mosquito%' GROUP BY study_alias", sep=""))
rs_alias_abstract_anopheles <- dbGetQuery(sra_con, paste( "SELECT *, count(*) AS Runs FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND study_abstract LIKE '%anopheles%' GROUP BY study_alias", sep=""))
rs_alias_abstract_plasmodium <- dbGetQuery(sra_con, paste( "SELECT *, count(*) AS Runs FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND study_abstract LIKE '%plasmodium%' GROUP BY study_alias", sep=""))
rs_alias_abstract_malaria <- dbGetQuery(sra_con, paste( "SELECT *, count(*) AS Runs FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND study_abstract LIKE '%malaria%' GROUP BY study_alias", sep=""))
rs_alias_exptitle_malaria <- dbGetQuery(sra_con, paste( "SELECT *, count(*) AS Runs FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND experiment_title LIKE '%malaria%' GROUP BY study_alias", sep=""))
rs_alias_exptitle_plasmodium <- dbGetQuery(sra_con, paste( "SELECT *, count(*) AS Runs FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND experiment_title LIKE '%plasmodium%' GROUP BY study_alias", sep=""))
rs_alias_exptitle_anopheles <- dbGetQuery(sra_con, paste( "SELECT *, count(*) AS Runs FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND experiment_title LIKE '%anopheles%' GROUP BY study_alias", sep=""))
rs_alias_exptitle_mosquito <- dbGetQuery(sra_con, paste( "SELECT *, count(*) AS Runs FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND experiment_title LIKE '%mosquito%' GROUP BY study_alias", sep=""))
rs_alias_sample_plasmodium <- dbGetQuery(sra_con, paste( "SELECT *, count(*) AS Runs FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND sample_attribute LIKE '%plasmodium%' GROUP BY study_alias", sep=""))
rs_alias_sample_malaria <- dbGetQuery(sra_con, paste( "SELECT *, count(*) AS Runs FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND sample_attribute LIKE '%malaria%' GROUP BY study_alias", sep=""))
rs_alias_sample_anopheles <- dbGetQuery(sra_con, paste( "SELECT *, count(*) AS Runs FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND sample_attribute LIKE '%anopheles%' GROUP BY study_alias", sep=""))
rs_alias_sample_mosquito <- dbGetQuery(sra_con, paste( "SELECT *, count(*) AS Runs FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND sample_attribute LIKE '%mosquito%' GROUP BY study_alias", sep=""))

# 283 unique experiments  with 11,631 runs
alias_experiments <- rbind.data.frame(rs_alias_abstract_anopheles, rs_alias_abstract_malaria, rs_alias_abstract_mosquito, rs_abstract_plasmodium, rs_alias_exptitle_anopheles, rs_alias_exptitle_malaria, rs_alias_exptitle_mosquito, rs_alias_exptitle_plasmodium, rs_alias_title_anopheles, rs_alias_title_malaria, rs_alias_title_mosquito, rs_alias_title_plasmodium, rs_alias_sample_anopheles, rs_alias_sample_malaria, rs_alias_sample_mosquito, rs_alias_sample_plasmodium)
unique_alias_experments <- unique(alias_experiments)

# Ungrouped experiments

rs_ungrouped_title_malaria <- dbGetQuery(sra_con, paste( "SELECT *  FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND study_title LIKE '%malaria%'   ", sep=""))
rs_ungrouped_title_plasmodium <- dbGetQuery(sra_con, paste( "SELECT *  FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND study_title LIKE '%plasmodium%'   ", sep=""))
rs_ungrouped_title_anopheles <- dbGetQuery(sra_con, paste( "SELECT *  FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND study_title LIKE '%anopheles%'   ", sep=""))
rs_ungrouped_title_mosquito <- dbGetQuery(sra_con, paste( "SELECT *  FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND study_title LIKE '%mosquito%'   ", sep=""))
rs_ungrouped_abstract_mosquito <- dbGetQuery(sra_con, paste( "SELECT *  FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND study_abstract LIKE '%mosquito%'   ", sep=""))
rs_ungrouped_abstract_anopheles <- dbGetQuery(sra_con, paste( "SELECT *  FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND study_abstract LIKE '%anopheles%'   ", sep=""))
rs_ungrouped_abstract_plasmodium <- dbGetQuery(sra_con, paste( "SELECT *  FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND study_abstract LIKE '%plasmodium%'   ", sep=""))
rs_ungrouped_abstract_malaria <- dbGetQuery(sra_con, paste( "SELECT *  FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND study_abstract LIKE '%malaria%'   ", sep=""))
rs_ungrouped_exptitle_malaria <- dbGetQuery(sra_con, paste( "SELECT *  FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND experiment_title LIKE '%malaria%'   ", sep=""))
rs_ungrouped_exptitle_plasmodium <- dbGetQuery(sra_con, paste( "SELECT *  FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND experiment_title LIKE '%plasmodium%'   ", sep=""))
rs_ungrouped_exptitle_anopheles <- dbGetQuery(sra_con, paste( "SELECT *  FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND experiment_title LIKE '%anopheles%'   ", sep=""))
rs_ungrouped_exptitle_mosquito <- dbGetQuery(sra_con, paste( "SELECT *  FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND experiment_title LIKE '%mosquito%'   ", sep=""))
rs_ungrouped_sample_plasmodium <- dbGetQuery(sra_con, paste( "SELECT *  FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND sample_attribute LIKE '%plasmodium%'   ", sep=""))
rs_ungrouped_sample_malaria <- dbGetQuery(sra_con, paste( "SELECT *  FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND sample_attribute LIKE '%malaria%'   ", sep=""))
rs_ungrouped_sample_anopheles <- dbGetQuery(sra_con, paste( "SELECT *  FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND sample_attribute LIKE '%anopheles%'   ", sep=""))
rs_ungrouped_sample_mosquito <- dbGetQuery(sra_con, paste( "SELECT *  FROM `sra` WHERE library_strategy IS 'RNA-Seq' AND sample_attribute LIKE '%mosquito%'   ", sep=""))

# 22148 rows
ungrouped_experiments <- rbind.data.frame(rs_ungrouped_abstract_malaria, rs_ungrouped_abstract_plasmodium, rs_ungrouped_exptitle_malaria, rs_ungrouped_exptitle_plasmodium, rs_ungrouped_title_malaria,rs_ungrouped_title_plasmodium, rs_ungrouped_sample_malaria, rs_ungrouped_sample_plasmodium)
unique_ungrouped_experiments <- unique(ungrouped_experiments) #7912 rows

db <- dbConnect(SQLite(), dbname="unique_ungrouped_experiments.sqlite")
dbWriteTable(conn = db, name = "unique_ungrouped_experiments", unique_ungrouped_experiments, overwrite=T, row.names=FALSE)

# 245 experiments, runs ranging from 1 to 1598, 7912 total runs
unique_ungrouped_experiments_runs <- dbGetQuery(db, paste("SELECT *, count(*) AS Runs FROM `unique_ungrouped_experiments` GROUP BY study_alias", sep=""))

# Using getSRA()
#rs_mosquito <- getSRA( search_terms = "mosquito", out_types = c('sra'), sra_con )
rs_murine <- getSRA( search_terms = c("murine AND malaria OR  Plasmodium"), out_types = c('sra'), sra_con )
rs_plasmodium <- getSRA( search_terms = "Plasmodium", out_types = c('sra'), sra_con )
rs_malaria <- getSRA( search_terms = "malaria", out_types = c('sra'), sra_con )
rs_mouse <- getSRA( search_terms = "mouse AND malaria OR Plasmodium", out_types = c('sra'), sra_con )
rs_rodent <- getSRA( search_terms = "rodent AND malaria OR Plasmodium", out_types = c('sra'), sra_con)
rs_human <- getSRA( search_terms = "human AND malaria OR Plasmodium", out_types = c('sra'), sra_con)
rs_monkey <- getSRA( search_terms = "monkey AND malaria OR Plasmodium", out_types = c('sra'), sra_con)
rs_macaque <- getSRA( search_terms = "Macaca AND malaria OR Plasmodium", out_types = c('sra'), sra_con)
rs_macaque2 <- getSRA( search_terms = "macaque AND malaria OR Plasmodium", out_types = c('sra'), sra_con)

ungrouped_rs <- rbind.data.frame(rs_malaria, rs_plasmodium, rs_murine, rs_mouse, rs_rodent, rs_human, rs_monkey, rs_macaque, rs_macaque2) #1812703 rows
unique_ungrouped_rs <- unique(ungrouped_rs) #1616183 rows without filtering RNA-Seq experiments

db <- dbConnect(SQLite(), dbname="unique_ungrouped_rs.sqlite")
dbWriteTable(conn = db, name = "unique_ungrouped_rs.sqlite", unique_ungrouped_rs, overwrite=T, row.names=FALSE, verbose = T)
# db_rodent <- dbConnect(SQLite(), dbname="rodent.sqlite")
# dbWriteTable(conn = db_rodent, name = "rodent", rs_rodent, overwrite=T, row.names=FALSE)

# 3886 experiments, runs ranging from 1 to 6665, 203341 total runs, 4052 experiments
unique_ungrouped_rs_runs <- dbGetQuery(db, paste("SELECT *, count(*) AS Runs FROM `unique_ungrouped_rs.sqlite` WHERE library_strategy IS 'RNA-Seq' GROUP BY study", sep=""))
# unique_rodents <- dbGetQuery(db_rodent, paste("SELECT *, count(*) AS Runs FROM `rodent` WHERE library_strategy IS 'RNA-Seq' GROUP BY study_alias", sep="")) # 76 rows, 2097 runs

# OR, use R subsetting to extraxt RNA-Seq expts

# get rows that have "malaria" or "Plasmodium" in them
# malaria <- sapply(unique_ungrouped_rs_runs[1:nrow(unique_ungrouped_rs_runs), 1:ncol(unique_ungrouped_rs_runs)], function(x) grep("malaria", x))
# plasmodium <- sapply(unique_ungrouped_rs_runs[1:nrow(unique_ungrouped_rs_runs), 1:ncol(unique_ungrouped_rs_runs)], function(x) grep("Plasmodium", x))
# malaria_rod <- sapply(unique_rodents[1:76, 1:71], function(x) grep("malaria", x))
# plasmodium_rod <- sapply(unique_rodents[1:76, 1:71], function(x) grep("Plasmodium", x))

# merge the row numbers obtained from malaria and plasmodium
# plasmodium_melt <- melt(plasmodium)
# malaria_melt <- melt(malaria)
# vector_rows <- c(malaria_melt[,1], plasmodium_melt[,1])
# unique_rows <- unique(vector_rows)

# plasmodium_melt_rod <- melt(plasmodium_rod)
# malaria_melt_rod <- melt(malaria_rod)
# vector_rows_rod <- c(malaria_melt_rod[,1], plasmodium_melt_rod[,1])
# unique_rows_rod <- unique(vector_rows_rod)

# 143 experiments, 6468 runs, run range: 1 to 1598
# check_unique_ungrouped_rs_runs <- unique_ungrouped_rs_runs[unique_rows,] # Some experiments to be excluded
# check_unique_ungrouped_rs_runs <- rbind(check_unique_ungrouped_rs_runs, check_unique_rod) # 140 experiments, 6643 runs
# check_unique_rod <- unique_rodents[unique_rows_rod,]

write.table(unique_ungrouped_rs_runs, "unique_ungrouped_rs_runs.txt", sep = "\t", row.names = FALSE)

# 3886 - 134 rows
# new_unique_ungrouped_rs_runs <- unique_ungrouped_rs_runs[-unique_rows,] # has no experiments related to malaria

# exclusion of undesired experiments

# Black_list <- c("ERP009392","ERP010191","SRP044861", "SRP070962", "SRP075802","SRP072137", "ERP004042", "SRP065966","ERP017920", "ERP005133", "ERP011524", "SRP028164", "SRP071926", "SRP066121", "SRP067884", "SRP067891", "SRP026535", "SRP078950", "SRP031506", "SRP075558", "SRP099346", "SRP103209", "SRP098800", "ERP012913","SRP009381", "DRP003241", "SRP055084", "SRP064783", "SRP067399", "ERP012209", "ERP001997", "SRP012327", "SRP050131", "ERP013211", "ERP004232", "ERP006797", "SRP012496", "SRP003813", "ERP005571", "ERP005646", "ERP005856", "ERP006110", "ERP021698", "ERP001849", "SRP063489", "SRP045347", "SRP070748", "SRP044683", "ERP020038", "SRP003507", "SRP016856", "SRP019362", "SRP034011", "SRP043116", "SRP048710", "SRP048711", "SRP055417", "SRP058108", "SRP058195", "SRP067125", "SRP073610", "SRP073801", "SRP079357", "SRP098628", "SRP013839", "SRP018078", "SRP028885", "SRP038137", "SRP045243", "SRP062654", "SRP068605", "SRP077596", "SRP082548", "SRP106793", "ERP016804", "ERP001696", "ERP101091", "ERP010044", "ERP010806", "SRP033414", "SRP043322", "SRP069075", "SRP106064", "SRP099925", "SRP013741", "SRP008152", "SRP115504", "SRP106032", "SRP090611", "SRP014155", "ERP001455", "ERP004740", "SRP009370", "SRP100893", "SRP090342", "SRP048791", "SRP026367", "SRP021890", "SRP017623")
# # 99 in black list
# remove_black_list <- c()
# for(i in 1:length(Black_list))
# {
#   remove_black_list <- c(remove_black_list, grep(Black_list[i], as.vector(check_unique_ungrouped_rs_runs[,"study"])))
# }

dbDisconnect(db, sra_con)

White_list <- c("DRP000987", "ERP004042", "SRP066796", "SRP070748", "ERP002116", "SRP018945", "SRP083918", "SRP094565", "SRP102683", "ERP004598", "ERP014302", "SRP110609", "SRP116117", "SRP029990", "SRP032775", "SRP034011", "SRP073801", "SRP074580", "SRP075802", "SRP106638", "SRP106798", "SRP108356", "SRP116593", "SRP116793", "SRP118503", "SRP118827", "SRP118996", "ERP021024", "ERP020067", "ERP004868", "ERP002273", "SRP059851", "SRP110282", "ERP005730", "SRP056443")


select_positive_studies <- function(White_list, check_unique_ungrouped_rs_runs)
{
  a <- nrow(check_unique_ungrouped_rs_runs)
  b <- length(White_list)
  c <- c()

  for(i in 1:b)
  {
    c[i] <- grep(White_list[i], as.vector(check_unique_ungrouped_rs_runs[,"study"]))
  }

  positive_studies <- check_unique_ungrouped_rs_runs[c,]

  return (positive_studies)
}

positive_studies <- select_positive_studies(White_list, check_unique_ungrouped_rs_runs)

range(positive_studies[, "Runs"])
sum(positive_studies[, "Runs"])
# 35 experiments in white list -> 4376 runs: Range -> 1 - 1598 ----> not all are published
write.table(positive_studies, "positive_studies.txt", sep = "\t", row.names = FALSE)
save(positive_studies, file = "positive_studies.RData")

positive_experiments <- positive_studies[,"study"] # extract experiment IDs
l <- length(positive_experiments)

write.table(positive_experiments, "positive_experiments.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
save(positive_experiments, file ="positive_experiments.RData")

# To check status: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA432151&go=go