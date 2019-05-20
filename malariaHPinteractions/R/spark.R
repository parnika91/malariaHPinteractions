install.packages("sparklyr")
install.packages("rsparkling")
options(rsparkling.sparklingwater.version = "2.4.8")
library(rsparkling)
library(sparklyr)
library(dplyr)
spark_install(version = "2.4.0")
sc <- spark_connect("local", verison = "2.1")
detach("package:rsparkling", unload = TRUE)
if ("package:h2o" %in% search()) { detach("package:h2o", unload = TRUE) }
if (isNamespaceLoaded("h2o")){ unloadNamespace("h2o") }
remove.packages("h2o")
install.packages("h2o", type = "source", repos = "https://h2o-release.s3.amazonaws.com/h2o/rel-tverberg/2/R")
library(h2o)


#spark_read_csv(sc, "erp", path = "Output/ERP106451/ERP106451_filteredDisp.txt",memory = FALSE, delimiter= '\t', header = T) -> tab1
tab1 <- read.table("Output/ERP106451/ERP106451_filteredDisp.txt", sep = '\t', header = T)
tab1 <- tab1[,-ncol(tab1)]
df_spark <- sparklyr::sdf_copy_to(sc, tab1, overwrite = T)

tab1 %>% sdf_sample(fraction = 0.9) %>% collect -> df1
df <- tab1[,-1]
df_spark <- sparklyr::sdf_copy_to(sc, df, overwrite = T)
library(corrplot)
M <- ml_corr(t(df_spark), method = "pearson") # also method="spearman" is possible
#Error: Unable to retrieve a spark_connection from object of class NULL

iris_tbl <- sdf_copy_to(sc, iris, name = "iris_tbl", overwrite = TRUE)
features <- c("Petal_Width", "Petal_Length", "Sepal_Length", "Sepal_Width")
ml_corr(iris_tbl, columns = features , method = "pearson")
# Error: ml_corr requires Spark 2.2.0 or higher.
spark_installed_versions()
# spark hadoop                                                   dir
# 2.4.0    2.7 /localstorage/parnika/spark/spark-2.4.0-bin-hadoop2.7


h2o.init()
p <- system.file("extdata", path = "Output/ERP106451/ERP106451_filteredDisp.txt",package = "h2o")
prostate <- h2o.importFile(path = "Output/ERP106451/ERP106451_filteredDisp.txt")
h2o_tab <- as.h2o(tab1)
tab1 %>% sdf_sample(fraction = 0.9) %>% collect -> df1
corre <- h2o.cor(t(h2o_tab), y = NULL, na.rm = F, use = "everything")
#ERROR: Unexpected HTTP Status code: 500 Server Error (url = http://localhost:54321/99/Rapids)

# Error in is.character(urlSuffix) : 
#   lexical error: invalid char in json text.
# <html> <head> <meta http-equiv=
#   (right here) ------^
#
df_spark <- sparklyr::copy_to(sc, df, overwrite = T)

hc <- rsparkling::h2o_context(sc)
hf <- rsparkling::as_h20_frame(sc, df_spark)

localH2O = h2o.init()
demo(h2o.kmeans)