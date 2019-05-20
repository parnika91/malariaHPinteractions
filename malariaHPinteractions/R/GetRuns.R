# This file is to  retrieve the run IDs from each study ID. It needs positive_studies.txt as an input.
# It reads each entry from the study column and runs getSRAinfo on each study ID. It then 
# stores the run IDs in a separate data frame.

# positive_experiments <- read.table("positive_experiments.txt", header = TRUE, sep = c("\t"))

#!/bin/env/Rscript
options(echo=TRUE)
args <- commandArgs(TRUE)

(studies <- args[1])

# source("https://bioconductor.org/biocLite.R")
# biocLite("SRAdb")
# source("https://bioconductor.org/biocLite.R")
# biocLite("GEOquery")
# install.packages("pkgconfig", repos="https://ftp.gwdg.de/pub/misc/cran/")

suppressMessages(library(SRAdb, verbose = FALSE))
suppressMessages(library(GEOquery, verbose = FALSE))
suppressMessages(library(pkgconfig, verbose = FALSE))

sqlfile <- 'SRAmetadb.sqlite'
sra_con <- dbConnect(SQLite(),sqlfile)
sraType = "sra"

# getSRAinfo <-
#   function( in_acc, sra_con, sraType = "sra" ) {
#     ## note: 'litesra' has phased out 
#     sraFile <- listSRAfile(in_acc, sra_con=sra_con, fileType=sraType, srcType='ftp')
#     sraFileDir <- paste(na.omit(unique(dirname(sraFile$ftp))), '/',
#                         sep='')
#     
#     file_name=NULL
#     file_size=NULL
#     file_date=NULL
#     require(RCurl)
#     opts=curlOptions(header=TRUE)
#     for( sraFileDir_1 in sraFileDir ) {
#       print(sraFileDir_1)
#       x <- getURL(sraFileDir_1,.opts=opts)
#       x1 <- strsplit(x[1], "\n")[[1]]
#       x2 <- sub('^.*\\s+anonymous\\s+','', x1, perl=TRUE)
#       
#       if(length(x2) == 0)
#       {
#         file_name <- c(file_name, "NA")
#         file_date <- c(file_date, "NA")
#         file_size <- c(file_size, "NA")
#       }else{
#       
#         file_name <-
#           c(file_name,
#             paste(sraFileDir_1, sub('^.*\\s+','', x2, perl=TRUE),
#                   sep='') )
#         file_date1 =sub('^\\d+\\s{1}','', x2, perl=TRUE)
#         file_date <-
#           c(file_date,
#             substr(file_date1, 1,
#                    gregexpr("\\s", file_date1, perl=TRUE)[[1]][4]-1 ) )
#         file_size <-
#           c(file_size,
#             ceiling(as.numeric(sub('\\s+.*$','', x2, perl=TRUE)) / 1024) )
#       }
#       ## curlPerform sometimes gives 'Access denied: 530' error
#       Sys.sleep(5)
#     }
#     Sys.sleep(5)
#     
#     file_info <-
#       as.data.frame(cbind('file_name'=file_name,
#                           'size(KB)'=file_size, 'date'=file_date))
#     sraFileInfo <-
#       merge(sraFile, file_info, by.x="ftp", by.y="file_name",
#             all.x=TRUE)
# 
#     
#     return(sraFileInfo)
#   }

#rs = getSRAinfo(in_acc=study, sra_con=sra_con, sraType = "sra")

# sraConvert function definition:

# sraConvert <-
#   function (in_acc,
#             out_type=c('sra','submission','study','sample','experiment','run'),
#             sra_con)
#   {
#     out_type <- tolower(out_type);
#     out_type <- match.arg(out_type, several.ok = T)
#     if( is.element('sra',out_type) )
#       out_type = c('submission','study','sample','experiment','run')
#     
#     ## validate in_acc
#     valid_in_acc_type <-
#       c('SRA', 'ERA', 'DRA', 'SRP', 'ERP', 'DRP', 'SRS', 'ERS',
#         'DRS', 'SRX', 'ERX', 'DRX', 'SRR', 'ERR', 'DRR')
#     valid_in_type <-
#       c('SRA'='submission', 'ERA'='submission', 'DRA'='submission',
#         'SRP'='study', 'ERP'='study', 'DRP'='study', 'SRS'='sample',
#         'ERS'='sample', 'DRS'='sample', 'SRX'='experiment',
#         'ERX'='experiment', 'DRX'='experiment', 'SRR'='run',
#         'ERR'='run', 'DRR'='run')
#     
#     ## trim leading or tailing spaces
#     in_acc <- sub('^\\s+|\\s+$','', in_acc, perl=TRUE)
#     ## the first three should be letters, not special characters, and
#     ## followed by numbers
#     if(any(grep('\\^W{3}|\\D+$', in_acc, perl=TRUE)))
#       stop("invalid input SRA accession(s), right ones are like 'SRA003625' or 'SRP000403', or 'SRS001834', 'SRR013350', or 'SRX002512'")
#     
#     ## extract the leading letters, which should be valid 
#     in_acc_type = toupper(unique(sub('\\d+$', '', in_acc, perl= TRUE)))
#     ## they should be valid
#     if( !all(in_acc_type %in%  valid_in_acc_type) )
#       stop("Input type shuld be in '",
#            paste(valid_in_acc_type, collapse="' '"),
#            "'")
#     in_type <- unique(valid_in_type[in_acc_type])
#     ## in_type should be only one type
#     if(length(in_type) != 1 )
#       stop("Only one type of SRA accession(s) is allowed in an input accession vector, either 'submission','study','sample','experiment' or 'run'")
#     
#     ## Exclude the in_type in the out_type	
#     out_type <- out_type[out_type != in_type];
#     select_type <- c(in_type, out_type)	
#     
#     ##Remove self converion
#     #	if(length(out_type) == 0) {		
#     #		sra_acc <- as.data.frame(cbind(run = in_acc))
#     #		return(sra_acc)
#     #		## print("Not necessary to convert to input itself");
#     #	}
#     
#     in_acc_sql = paste("'", paste(in_acc, collapse = "','"),"'", sep="");
#     select_type_sql <- paste(paste(select_type, "_accession", sep=''),
#                              collapse = "," );
#     sql <- paste ("SELECT DISTINCT ", select_type_sql,
#                   " FROM sra WHERE ", in_type ,
#                   "_accession IN (", in_acc_sql, ")", sep = "");			 
#     sra_acc <- dbGetQuery(sra_con, sql);
#     names(sra_acc) <- sub('_accession', '', names(sra_acc))
#     return(sra_acc);
#     
#   }
# 
# 
# # listSRAfile function definition:
# 
listSRAfile <-
  function (in_acc, sra_con, fileType = "sra", srcType = "ftp")
{
  if (fileType == "fastq") {
    sra_acc = sraConvert(in_acc, out_type = c("run"), sra_con = sra_con)
    sraFiles = getFASTQinfo(sra_acc$run, sra_con, srcType)
  }
  else if (fileType == "sra") {
    sraExt <- ".sra"
    sra_acc <- sraConvert(in_acc, out_type = c("study", "sample",
                                               "experiment", "run"), sra_con = sra_con)
    if (srcType == "fasp") {
      srcMain = "anonftp@ftp-trace.ncbi.nlm.nih.gov:"
    }
    else if (srcType == "ftp") {
      srcMain = "ftp://ftp-trace.ncbi.nlm.nih.gov"
    }
    sraFiles_1 = NULL
    for (i in 1:nrow(sra_acc)) {
      sraFileDir <- paste(srcMain, "/sra/sra-instant/reads/ByRun/",
                          fileType, "/", substring(sra_acc$run[i], 1, 3),
                          "/", substring(sra_acc$run[i], 1, 6), "/", sra_acc$run[i],
                          "/", sep = "")
      if (is.na(sra_acc$run[i])) {
        sraFiles1 <- NA
      }
      else {
        sraFiles1 <- paste(sraFileDir, sra_acc$run[i],
                           sraExt, sep = "")
      }
      sraFiles_1 = c(sraFiles_1, sraFiles1)
    }
    sraFiles <- cbind(sra_acc, sraFiles_1, stringsAsFactors = FALSE)
    colnames(sraFiles) <- c(names(sra_acc), srcType)
  }
  return(sraFiles)
  }

# sraConvert:
function (in_acc, out_type = c("sra", "submission", "study", 
                               "sample", "experiment", "run"), sra_con) 
{
  out_type <- tolower(out_type)
  out_type <- match.arg(out_type, several.ok = T)
  if (is.element("sra", out_type)) 
    out_type = c("submission", "study", "sample", "experiment", 
                 "run")
  valid_in_acc_type <- c("SRA", "ERA", "DRA", "SRP", "ERP", 
                         "DRP", "SRS", "ERS", "DRS", "SRX", "ERX", "DRX", "SRR", 
                         "ERR", "DRR")
  valid_in_type <- c(SRA = "submission", ERA = "submission", 
                     DRA = "submission", SRP = "study", ERP = "study", DRP = "study", 
                     SRS = "sample", ERS = "sample", DRS = "sample", SRX = "experiment", 
                     ERX = "experiment", DRX = "experiment", SRR = "run", 
                     ERR = "run", DRR = "run")
  in_acc <- sub("^\\s+|\\s+$", "", in_acc, perl = TRUE)
  if (any(grep("\\^W{3}|\\D+$", in_acc, perl = TRUE))) 
    stop("invalid input SRA accession(s), right ones are like 'SRA003625' or 'SRP000403', or 'SRS001834', 'SRR013350', or 'SRX002512'")
  in_acc_type = toupper(unique(sub("\\d+$", "", in_acc, perl = TRUE)))
  if (!all(in_acc_type %in% valid_in_acc_type)) 
    stop("Input type shuld be in '", paste(valid_in_acc_type, 
                                           collapse = "' '"), "'")
  in_type <- unique(valid_in_type[in_acc_type])
  if (length(in_type) != 1) 
    stop("Only one type of SRA accession(s) is allowed in an input accession vector, either 'submission','study','sample','experiment' or 'run'")
  out_type <- out_type[out_type != in_type]
  select_type <- c(in_type, out_type)
  in_acc_sql = paste("'", paste(in_acc, collapse = "','"), 
                     "'", sep = "")
  select_type_sql <- paste(paste(select_type, "_accession", 
                                 sep = ""), collapse = ",")
  sql <- paste("SELECT DISTINCT ", select_type_sql, " FROM sra WHERE ", 
               in_type, "_accession IN (", in_acc_sql, ")", sep = "")
  sra_acc <- dbGetQuery(sra_con, sql)
  names(sra_acc) <- sub("_accession", "", names(sra_acc))
  return(sra_acc)
}

current_runs <- data.frame()
for(i in 1:length(studies))
{
	rs_runs <- listSRAfile(in_acc=studies[i], sra_con=sra_con, fileType=sraType, srcType='ftp')[,c(1,4)]
	# randomize runs, keeping study intact
	rs_runs <- data.frame(run=sample(rs_runs$run), study=rs_runs$study)
	# rs = rs[!is.na(rs[,6]),]
	# rs_ = as.numeric(as.character(rs[,6]))
	# runs_rs <- c(rs[,5]) # ---> include this later instead of the following if condition

	# runs_rs <- c()
	# if(nrow(rs) <= 5)
	 #  {runs_rs <- rs[,5]} else
	 #  {runs_rs <- rs[sample(nrow(rs), 5),5] }


	#runs_rs <- as.data.frame(runs_rs[!is.na(runs_rs)])
	# # #colnames(runs_rs) <- "run"
	save(rs_runs, file = "runs_rs.RData")
	write.table(rs_runs, paste0("Output/",studies[i],"/runs_",studies[i],".txt", sep = ",", collapse = ''),row.names = FALSE, col.names = FALSE, quote = FALSE)
	#write.table(rs_, "rs_.txt", sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)

	# shell command: bash --verbose runDownloadRuns.sh
	current_runs <- rbind(current_runs, rs_runs)
}
write.table(current_runs, "Input/current_runs.txt", sep = ',', row.names = FALSE, col.names = FALSE, quote = FALSE)


# for(i in 1:length(studies))
# {
#   read <- read.csv2(paste0("Output/",studies[i],"/runs_",studies[i],".txt", collapse = ''), sep = ',', header = F)
#   current_runs <- rbind(current_runs, read)
# }