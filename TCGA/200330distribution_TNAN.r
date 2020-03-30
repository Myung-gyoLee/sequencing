
library(dplyr)
setwd("H:TCGA/00ARCHs4_h5")
load(file="log2_quantile_normalization.RData")


metaT=read.csv("../01AllCancerTissue/200214TCGAmeta.csv")
metaG=read.csv("../01AllCancerTissue/200214GTExmeta.csv")

head(metaT)
head(metaG)


## TCGA metafile 
rm(metaT1)
metaT1<-data.frame(matrix(nrow=nrow(metaT), ncol=4))
colnames(metaT1)=c("Sampleid","Tissue","Sampletype","Condition")

metaT1["Sampleid"] <- metaT$Barcode
metaT1["Tissue"] <- metaT$Cancertype
metaT1["Sampletype"] <- metaT$Sampletype

#table(metaT1$Sampletype)
metaT1["Condition"]="Tumor"
metaT1[grep("Normal",metaT1$Sampletype),]["Condition"]="adj.Normal"
#table(metaT1$Condition)

head(metaT1)


## GTEx metafile 
metaG1<-data.frame(matrix(nrow=nrow(metaG), ncol=4))
colnames(metaG1)=c("Sampleid","Tissue","Sampletype","Condition")

metaG1["Sampleid"] <- metaG$Sampleid
metaG1["Tissue"] <- metaG$Tissue
metaG1["Sampletype"] <- metaG$Experiment

metaG1["Condition"]="Normal"

head(metaG1)

## rbind metafile of tcga gtex
metaall <- rbind(metaT1,metaG1) 
#write.csv(metaall, "../01AllCancerTissue/200330metaTCGA_GTEx.csv")


## expression file transpose and remove duplication
Texp_mod <-Texpression[,!duplicated(colnames(Texpression))]


TGexpression <- cbind(Texp_mod, Gexpression)
rTGexpression <- t(TGexpression)
# dim(rTGexpression) [1] 20852 25150


metaall <- metaall[!duplicated(metaall$Sampleid,metaall$Tissue),]
rownames(metaall) <- metaall$Sampleid
metaall <- metaall[,-1]
metaall %>% head


TCGA_GTEx <- cbind(metaall, rTGexpression, copy = FALSE)

save(list= c("metaall","rTGexpression","TGexpression","TCGA_GTEx","Texp_mod"), file="log2_metaall_T_Gexpression.RData")
#rm(list = c("metaG", "metaG1", "metaT", "metaT1"))

saveRDS(TCGA_GTEx, file = "200330TCGA_GTEx.RDS")
write.csv(TCGA_GTEx, file = "200330TCGA_GTEx.csv")

rlist = ls()
rlist = rlist[rlist != "TCGA_GTEx"]
rlist
#rm(list = c("Texpression", "Gexpression","TGexpression","Texp_mod"))
