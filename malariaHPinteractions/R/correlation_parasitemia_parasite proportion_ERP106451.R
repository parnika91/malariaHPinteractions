parasitemia_ERP106451 <- read.table("/SAN/Plasmo_compare/SRAdb/Output/ERP106451/parasitemia_ERP106451.txt", sep = '\t', header = T)
parasiteProportion_ERP106451 <- read.table("/SAN/Plasmo_compare/SRAdb/Output/ERP106451/hp_percent_ERP106451.txt", sep = '\t', header = T)

parasiteProportion_parasitemia <- merge(parasitemia_ERP106451, parasiteProportion_ERP106451, by = "Run")

parasiteProportion_parasitemia <- parasiteProportion_parasitemia[-2,]

cor(as.numeric(parasiteProportion_parasitemia[,4]), as.numeric(as.character(parasiteProportion_parasitemia[,2])))
# 0.51

para.lm <- lm(Parasite_percent ~ as.numeric(as.character(percentage.parasitemia)), data=parasiteProportion_parasitemia)

removeOutliers <- parasiteProportion_parasitemia[-c(11,21,42),] # as seen from plots
para.lm <- lm(Parasite_percent ~ as.numeric(as.character(percentage.parasitemia)), data=removeOutliers)

cor(as.numeric(removeOutliers[,4]), as.numeric(as.character(removeOutliers[,2])))
# 0.52 