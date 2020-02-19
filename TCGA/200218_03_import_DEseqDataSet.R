setwd("H:TCGA")
load("matrixARCHs4.Rdata")

#------------------------------------------------------------------#
## 200214_02_de_input_matrix.r ##

## join 1) TCGA tumor vs TCGA normal
#write.csv(BRCA_TN, "./02Breast/021merge/ExBRCA_Tumor_normal.csv")
#write.csv(metaTtum_nor, "./02Breast/021merge/metaBRCA_Tumor_normal.csv")

## join 2) TCGA tumor vs GTEx normal
#write.csv(Ttumor_Gnormal, "./02Breast/021merge/ExBRCA_Tumor_GTExnormal.csv")
#write.csv(metaTtum_Gnor, "./02Breast/021merge/metaBRCATum_GTExnormal.csv")

## join 3) TCGA tumor vs GTEx normal+TCGA matched normal
#write.csv(TCGA_Gtex, "./02Breast/021merge/ExBRCA_GTExnorm_matchnorm.csv")
#write.csv(metaTtum_TGnor, "./02Breast/021merge/metaBRCATum_TGnor.csv")

#------------------------------------------------------------------#
## join 1) TCGA tumor vs TCGA normal
library(DESeq2)

# reread Expression Data and Meta data
BRCA_TN1=read.csv("./02Breast/021merge/ExBRCA_Tumor_normal.csv", header = T, row.names = 1)
metaTtum_nor1=read.csv("./02Breast/021merge/metaBRCA_Tumor_normal.csv", header = T, row.names = 1)

dds = DESeqDataSetFromMatrix(countData = BRCA_TN1,
                             colData = metaTtum_nor1[,c(2,3,5)],
                             design = ~ Condition)

dds
#View(metaTtum_nor1)
#rownames(metaTtum_nor) = seq(1:length(rownames(metaTtum_nor)))
#rownames(BRCA_TN) = seq(1:length(rownames(BRCA_TN)))

#------------------------------------------------------------------#
save(list=ls(), file='./02Breast/021merge/ddsBreast.Rdata')
