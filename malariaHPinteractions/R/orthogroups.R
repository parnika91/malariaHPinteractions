orthogroups <- read.delim("Downloads/Orthogroups.csv")
Pv_Pcy <- orthogroups[,c("X","Pvivax", "Pcynomolgi")]
Pv_Pv_Pcy <- Pv_Pcy[,c("X","Pvivax")]
Pv_copies <- data.frame(Orthogroup = Pv_Pv_Pcy$X, copies = sapply(Pv_Pv_Pcy$Pvivax, function(x) length(strsplit(as.character(x), split = ",")[[1]])))
Pcy_Pv_Pcy <- Pv_Pcy[,c("X","Pcynomolgi")]
Pcy_copies <- data.frame(Orthogroup = Pcy_Pv_Pcy$X, copies = sapply(Pcy_Pv_Pcy$Pcynomolgi, function(x) length(strsplit(as.character(x), split = ",")[[1]])))

Pv_Pcy <- merge(Pv_copies, Pcy_copies, by = "Orthogroup")
Pv_Pcy <- Pv_Pcy[(Pv_Pcy$copies.x==1 & Pv_Pcy$copies.y==1),]
# Pv_Pcy = 4617


Pco <- orthogroups[,c("X","Pcoatneyi")]
Pco_copies <- data.frame(Orthogroup = Pco$X, copies = sapply(Pco$Pcoatneyi, function(x) length(strsplit(as.character(x), split = ",")[[1]])))
Pv_Pcy_Pco <- merge(Pv_Pcy, Pco_copies, by = "Orthogroup")
Pv_Pcy_Pco <- Pv_Pcy_Pco[(Pv_Pcy_Pco$copies.x==1 & Pv_Pcy_Pco$copies.y==1 & Pv_Pcy_Pco$copies==1),]
# Pv_Pcy_Pco = 4365

Pb_Pch <- orthogroups[,c("X","Pberghei", "Pchabaudi")]
Pb_Pb_Pch <- Pb_Pch[,c("X","Pberghei")]
Pb_copies <- data.frame(Orthogroup = Pb_Pb_Pch$X, copies = sapply(Pb_Pb_Pch$Pberghei, function(x) length(strsplit(as.character(x), split = ",")[[1]])))
Pch_Pb_Pch <- Pb_Pch[,c("X","Pchabaudi")]
Pch_copies <- data.frame(Orthogroup = Pch_Pb_Pch$X, copies = sapply(Pch_Pb_Pch$Pchabaudi, function(x) length(strsplit(as.character(x), split = ",")[[1]])))

Pb_Pch <- merge(Pb_copies, Pch_copies, by = "Orthogroup")
Pb_Pch <- Pb_Pch[(Pb_Pch$copies.x==1 & Pb_Pch$copies.y==1),]
# Pb_Pch = 4611

Py <- orthogroups[,c("X","Pyoelii")]
Py_copies <- data.frame(Orthogroup = Py$X, copies = sapply(Py$Pyoelii, function(x) length(strsplit(as.character(x), split = ",")[[1]])))
Pb_Pch_Py <- merge(Pb_Pch, Py_copies, by = "Orthogroup")
Pb_Pch_Py <- Pb_Pch_Py[(Pb_Pch_Py$copies.x==1 & Pb_Pch_Py$copies.y==1 & Pb_Pch_Py$copies==1),]
# Pb_Pch_Py = 4605

Pv_Pcy_Pco_Pb_Pch_Py <- merge(Pb_Pch_Py, Pv_Pcy_Pco, by = "Orthogroup")
# Pv_Pcy_Pco_Pb_Pch_Py = 4113

Pf <- orthogroups[,c("X", "Pfalciparum")]
Pf_copies <- data.frame(Orthogroup = Pf$X, copies = sapply(Pf$Pfalciparum, function(x) length(strsplit(as.character(x), split = ",")[[1]])))
Pf <- Pf_copies[(Pf_copies$copies==1),]
Pf_Pv_Pcy_Pco_Pb_Pch_Py <- merge(Pf, Pv_Pcy_Pco_Pb_Pch_Py, by = "Orthogroup")
# Pf_Pv_Pcy_Pco_Pb_Pch_Py = 4010