

# 200213_extract_metrix_meta.r

setwd("H:TCGA")
load("matrixARCHs4.Rdata")


## read data tables

Texpression=read.csv("Bar200213TCGAAllcancer.csv")
metaT=read.csv("200214TCGAmeta.csv")
Gexpression=read.csv("200213GTEXall.csv")
metaG=read.csv("200214GTExmeta.csv")

#------------------------------------------------------------------#
## TCGA-BRCA Primary Tumor and Metastatic Expression Table
#colnames(metaT)
Tbreast=metaT[grep("Breast", metaT$Tissue),]
#table(Tbreast$Sampletype)

# normal sample remove
Tbreast_tumor=Tbreast[-grep("Normal", Tbreast$Sampletype),]
#nrow(Tbreast[grep("Normal", Tbreast$Sampletype),]) # nrow=112
#colnames(Tbreast_tumor)
#head(Tbreast_tumor)
#table(Tbreast_tumor$Sampletype)

# barcode of Primary tumor and metastatic-> indexing Expression file
Tbrca_tumor_bar=Tbreast_tumor$Sampleid
#ncol(Texpression[,Tbrca_tumor_bar]) # 1134
TBRCA_ex=Texpression[,Tbrca_tumor_bar]
#head(rownames(TBRCA_ex))

#rm(Tbrca_tumor_bar)

#write.csv(TBRCA_ex,"200214BRCA_Tumor_Metastatic_1134.csv")
#write.csv(Tbreast_tumor, "200214metaBRCA_Tumor.csv")

## TCGA-BRCA Adjacent normal Expression Table
Tbreast_normal=Tbreast[grep("Solid", Tbreast$Sampletype),]
normalbar=Tbreast_normal$Barcode
TBRCA_normalEx=Texpression[,normalbar]
#ncol(TBRCA_normalEx)
#nrow(normalbar)
#table(Tbreast_normal$Sampletype)

#write.csv(TBRCA_normalEx, "200214BRCA_Normal_112.csv")
#write.csv(Tbreast_normal,"200214metaBRCA_Normal.csv")

#------------------------------------------------------------------#
## GTEx Breast Tissue Expression Table
#table(metaG$Tissue)
Gbreast=metaG[grep("Breast", metaG$Tissue),]
#nrow(metaG[grep("Breast", metaG$Tissue),]) = 218
Gbreast_id=Gbreast$Sampleid

GBreast_ex=Gexpression[,Gbreast_id]
#head(rownames(GBreast_ex))
#length(colnames(GBreast_ex))

#write.csv(GBreast_ex, "200214GTEx_Breast_218.csv")
#write.csv(Gbreast, "200214metaBreast_GTEx.csv")

list.files()

#------------------------------------------------------------------#
library(dplyr)

## join 1) TCGA tumor vs TCGA normal
Tbreast # meta file
#colnames(Tbreast)
Tbar=Tbreast$Barcode
BRCA_TN=Texpression[,Tbar] # expression file

#dir.create("./02Breast/021merge")

#write.csv(metaTtum_nor, "./02Breast/021merge/metaBRCA_Tumor_normal.csv")
#write.csv(BRCA_TN, "./02Breast/021merge/ExBRCA_Tumor_normal.csv")

# meta file
metaTtum_nor=metaall[metaall$Project == "TCGA" & metaall$Tissue == "Breast Invasive Carcinoma",]


#nrow(metaTtum_nor)
#nrow(metaall)


## join 2) TCGA tumor vs GTEx normal

colnames(TBRCA_ex)
Ttumor_Gnormal=merge(TBRCA_ex,GBreast_ex, by = 'row.names',all=TRUE)
#nrow(Ttumor_Gnormal)
#nrow(TBRCA_ex)
#tail(colnames(Ttumor_Gnormal))
#write.csv(Ttumor_Gnormal, "./02Breast/021merge/ExBRCA_Tumor_GTExnormal.csv")

# meta file

metaTtum_Gnor=rbind(metaTtum_nor[metaTtum_nor$Condition == "Tumor",],metaall[metaall$Project == "GTEx" & metaall$Tissue == "Breast",])

#write.csv(metaTtum_Gnor, "./02Breast/021merge/metaBRCATum_GTExnormal.csv")

#nrow(metaall[metaall$Project == "GTEx" & metaall$Tissue == "Breast",])
#str(metaTtum_Gnor)
#table(metaTtum_Gnor$Condition)
#write.csv(meta1,file="200214TCGAmeta.csv")

## join 3) TCGA tumor vs GTEx normal+TCGA matched normal
#ls()
TCGA_Gtex=merge(BRCA_TN,GBreast_ex, by = 'row.names', all=TRUE)
metaTtum_TGnor=rbind(metaall,metaall[metaall$Project == "GTEx" & metaall$Tissue == "Breast",])
#write.csv(TCGA_Gtex, "./02Breast/021merge/ExBRCA_GTExnorm_matchnorm.csv")
#write.csv(metaTtum_TGnor, "./02Breast/021merge/metaBRCATum_TGnor.csv")

#save(Texpression,Gexpression,metaT, metaG, file = "matrixARCHs4.Rdata")
save(list=ls(), file='matrixARCHs4.Rdata')




