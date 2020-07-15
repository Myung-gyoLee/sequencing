setwd("H:/TCGA/TCGA_marker")

#### 1. make each cancer, lineage dataframe per tissue type
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
pancreas_lineage = c("AKR7A3", "ALG5", "AMY1A", "AMY1B", "AMY1C", "AMY2A", "AMY2B", "ANXA4", "AQP12A", "AQP12B", "ARPIN", "BACE1", "BBIP1", "BNIP3", "CASP9", "CATSPERB", "CEL", "CELA2A", "CELA2B", "CELA3A", "CELA3B", "CELP", "CFTR", "CHRNB3", "CLPS", "CLPSL1", "CLPSL2", "COMTD1", "CPA1", "CPA2", "CPB1", "CRH", "CTRB1", "CTRB2", "CTRC", "CTRL", "CUZD1", "DAP", "EIF4EBP1", "EPB41L4B", "ERO1B", "FAM159B", "FKBP11","FUT1", "G6PC2", "GAT ", "GCG", "GLP1R", "GNMT", "GP2", "GPR119", "GRB10", "IAPP", "IFRD1", "IL22RA1", "INS", "INS-IGF2", "KCNK16", "KIRREL2", "KLK1", "KSR1", "LFNG", "LGR4", "LINC00339", "LINC01251","LINC01625", "LMF2","LOC101927701", "LOC105370616", "LOC105374344", "LOC154449", "LOC158434", "METTL21B", "MIR208B", "MKNK1", "MNX1", "MNX1-AS1", "NAA16", "NEURL3", "NEUROD1","NOMO1", "NOMO2", "NOMO3", "NR5A2", "OR10G4", "OR10G9", "OR4D5", "OR6M1", "OR6T1", "OR8D4","DCD4", "PDIA2", "PDX1", "PGGHG", "PLA2G1B", "PLEKHH3", "PNLIP", "PNLIPRP1", "PNLIPRP2","PPP2R2D", "PPY", "PRDX4", "PRSS1", "PRSS2", "PRSS", "PSEN2", "PSMD6", "RAB3D", "RBPJL", "REG1A", "REG1B", "REG3G", "RNF214", "RRBP1", "SEL1L", "SERP1", "SLC25A45", "SLC38A5", "SLC39A5","SLC43A1", "SLC4A4", "SPINK1", "SRPRA", "SRPRB", "SYBU", "TBC1D30", "TC2N", "TECPR1", "TM7SF3", "TMED11P", "TMEM97", "TPRN", "TPST2", "TRPV6", "XPNPEP1", "ZNF710")

pancreas_lineage = read.table("clipboard")
# AKR7A3       ALG5         AMY1A        AMY1B        AMY1C        AMY2A       
# AMY2B        ANXA4        AQP12A       AQP12B       ARPIN        BACE1       
# BBIP1        BNIP3        CASP9        CATSPERB     CEL          CELA2A      
# CELA2B       CELA3A       CELA3B       CELP         CFTR         CHRNB3      
# CLPS         CLPSL1       CLPSL2       COMTD1       CPA1         CPA2        
# CPB1         CRH          CTRB1        CTRB2        CTRC         CTRL        
# CUZD1        DAP          EIF4EBP1     EPB41L4B     ERO1B        FAM159B     
# FKBP11       FUT1         G6PC2        GAT          GCG          GLP1R       
# GNMT         GP2          GPR119       GRB10        IAPP         IFRD1       
# IL22RA1      INS          INS-IGF2     KCNK16       KIRREL2      KLK1        
# KSR1         LFNG         LGR4         LINC00339    LINC01251    LINC01625   
# LMF2         LOC101927701 LOC105370616 LOC105374344 LOC154449    LOC158434   
# METTL21B     MIR208B      MKNK1        MNX1         MNX1-AS1     NAA16       
# NEURL3       NEUROD1      NOMO1        NOMO2        NOMO3        NR5A2       
# OR10G4       OR10G9       OR4D5        OR6M1        OR6T1        OR8D4       
# DCD4         PDIA2        PDX1         PGGHG        PLA2G1B      PLEKHH3     
# PNLIP        PNLIPRP1     PNLIPRP2     PPP2R2D      PPY          PRDX4       
# PRSS1        PRSS2        PRSS         PSEN2        PSMD6        RAB3D       
# RBPJL        REG1A        REG1B        REG3G        RNF214       RRBP1       
# SEL1L        SERP1        SLC25A45     SLC38A5      SLC39A5      SLC43A1     
# SLC4A4       SPINK1       SRPRA        SRPRB        SYBU         TBC1D30     
# TC2N         TECPR1       TM7SF3       TMED11P      TMEM97       TPRN        
# TPST2        TRPV6        XPNPEP1      ZNF710      
#pancreas_lin <- pancreas_lineage$V1
#str(pancreas_lin)
pancreas_cancer = c("ADAM28","ADRA2A","ANXA10","BRS3","CASR","CELA2A","CELA2B","CELA3A","CELA3B","CLPS","CTRB1","CTRB2","CTRC","FAM159B","FOXL1","G6PC2","GCG","GPR119","GPR142","HTRA3","IAPP","INS","KCNK16","LINC00443","LOC284865","PNLIP","PNLIPRP1","PPY","PRSS1","PRSS2","RBPJL","SLC30A8","SPINK1","SST","SYCN","TCN1","TFF2","UCN3")

lung_lineage = c("A2M", "ABCA3", "ADAMTS8", "AGER", "CCDC102B", "CHIAP2", "CYP2B7P", "EMP2", "EPAS1", "FCN3", "GGTLC1", "KCNA10", "LAMP3", "LINC01108", "LOC101927123", "LOC105373580", "MIR4635", "MS4A15", "MYO16-AS1", "NAPSA", "NCKAP5", "RGCC", "RHOBTB2", "ROS1", "SCGB1A1", "SCGB3A2", "SFTA1P", "SFTPA1", "SFTPD", "SLC6A4", "TMEM100", "VIPR1")

lung_cancer = c("BPIFA1","LOC105373373","LOC642131","MIR650","NAPSA","NDNF","ROS1","SCGB1A1","SCGB3A1","SFTPA1","SFTPA2","SFTPB","SFTPD","TBX4","TFAP2D","TREM1")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
columns(org.Hs.eg.db)
cols <- c("GENENAME","ALIAS")

#dt_pancreas_lineage

dt_pancreas_lineage = select(org.Hs.eg.db, keys=as.character(pancreas_lineage), columns=cols, keytype="SYMBOL")

colnames(dt_pancreas_lineage) = c("Gene_symbol","Gene_name")
colnames(dt_pancreas_lineage) = c("Gene_symbol","Gene_name", "GTEx_Lineage_marker")
dt_pancreas_lineage["Lineage_marker"] = "Pancreas"
head(dt_pancreas_lineage)

#dt_pancreas_cancer
dt_pancreas_cancer = select(org.Hs.eg.db, keys=pancreas_cancer, columns=cols, keytype="SYMBOL")

colnames(dt_pancreas_cancer) = c("Gene_symbol","Gene_name")
dt_pancreas_cancer["TCGA_Cancer_marker"] = "Pancreas"
head(dt_pancreas_cancer)

#dt_lung_lineage
dt_lung_lineage = select(org.Hs.eg.db, keys=lung_lineage, columns=cols, keytype="SYMBOL")

colnames(dt_lung_lineage) = c("Gene_symbol","Gene_name")

dt_lung_lineage["GTEx_Lineage_marker"] = "Lung"
head(dt_lung_lineage)

#dt_lung_cancer
dt_lung_cancer = select(org.Hs.eg.db, keys=lung_cancer, columns=cols, keytype="SYMBOL")

colnames(dt_lung_cancer) = c("Gene_symbol","Gene_name")
dt_lung_cancer["TCGA_Cancer_marker"] = "Lung"
head(dt_lung_cancer)

#### 2. write xlsx per tissue type ####

library(xlsx)
### Pancreas
# Write the first data set in a new workbook
write.xlsx(dt_pancreas_cancer, file = "Pancreas_TCGA_GTEx_Cancer_Lineage_markers.xlsx", sheetName="pancreas_cancer", append= FALSE, row.names = FALSE)

# Add a second data set in a new worksheet
write.xlsx(dt_pancreas_lineage, file = "Pancreas_TCGA_GTEx_Cancer_Lineage_markers.xlsx", sheetName = "pancreas_lineage", append = TRUE, row.names = FALSE)




### Lung
# Write the first data set in a new workbook
write.xlsx(dt_lung_cancer, file = "Lung_TCGA_GTEx_Cancer_Lineage_markers.xlsx", sheetName="lung_cancer", append=FALSE, row.names = FALSE)

# Add a second data set in a new worksheet
write.xlsx(dt_lung_lineage, file = "Lung_TCGA_GTEx_Cancer_Lineage_markers.xlsx", sheetName = "lung_lineage", append = TRUE, row.names = FALSE)



#### 3. outerjoin Cancer & lineage to Dataframe ####
library(dplyr)

### Pancreas
out_pancreas = full_join(dt_pancreas_cancer, dt_pancreas_lineage, by ="Gene_symbol")

head(out_pancreas)
colnames(out_pancreas)
# > colnames(out_pancreas)
# [1] "Gene_symbol"         "Gene_name.x"         "TCGA_Cancer_marker" 
# [4] "Gene_name.y"         "GTEx_Lineage_marker"

pancreas_merge = out_pancreas[,c("Gene_symbol", "TCGA_Cancer_marker", "GTEx_Lineage_marker")]
write.csv(pancreas_merge, "pancreas_join_cancer_lineage.csv", row.names = FALSE)


### Lung

out_lung = full_join(dt_lung_cancer, dt_lung_lineage, by ="Gene_symbol")
lung_merge = out_lung[,c("Gene_symbol", "TCGA_Cancer_marker", "GTEx_Lineage_marker")]

write.csv(lung_merge, "lung_join_cancer_lineage.csv", row.names = FALSE)


library(org.Hs.eg.db)
## dplyr error : #opancreas_anno = select(org.Hs.eg.db, keys=as.character(out_pancreas$Gene_symbol), columns=cols, keytype="SYMBOL")

#### 4. left_join 3. dataframe to input genelist ####
# pancreas foldchange >= 1.5
library(readxl)
filen = "./data/SEV003_B1Tissue_Conserved_Markers.xlsx"
clusterlist <- excel_sheets(filen)

d1=NULL
for(i in clusterlist){
  #cluster_num = i
  #sheetname = sprintf("Cluster_%s",cluster_num)
  sheetname = i
  d2 <- read_excel(filen, sheet = sheetname)
  d2 = subset(d2, blood1.avgFC >= 1.5 & tissue.avgFC >= 1.5)
  d1 <- rbind(d1,d2)
}

d1$gene

require(openxlsx)
list_of_datasets <- lapply(X = clusterlist, FUN = function(x){
  x <- read_excel(filen, sheet = x)
  x <- subset(x, blood1.avgFC >= 1.5 & tissue.avgFC >= 1.5)
  x <- x[,c(1,2,3,4,5)]
  x <- left_join(x, pancreas_merge, c("gene" = "Gene_symbol"))
})

names(list_of_datasets) = clusterlist
list_of_datasets

intersect(list_of_datasets$Cluster_5$gene,
          pancreas_merge$Gene_symbol)

intersect(d1$gene,
          pancreas_merge$Gene_symbol)

# list_of_datasets <- list("Name of DataSheet1" = dataframe1, "Name of Datasheet2" = dataframe2)
write.xlsx(list_of_datasets, file = "./output/SEV003_B1Tissue_Conserved_Markers_lineage_cancer_annotation.xlsx")

save.image("./intersect_conservedToPublicMarker.RData")
# left join
# > colnames(pancreas_merge)
# [1] "Gene_symbol"         "TCGA_Cancer_marker"  "GTEx_Lineage_marker"
load("H:/TCGA/TCGA_marker/intersect_conservedToPublicMarker.RData")


### Lung
#filen = "./data/SMC027_B1Tissue_Conserved_Markers.xlsx"
filen = "./data/SMC025_B1Tissue_Conserved_Markers.xlsx" #EPAS1
#filen = "./data/SMC024_B1Tissue_Conserved_Markers_Add.xlsx"
filen
clusterlist <- excel_sheets(filen)

d1=NULL
for(i in clusterlist){
  #cluster_num = i
  #sheetname = sprintf("Cluster_%s",cluster_num)
  sheetname = i
  d2 <- read_excel(filen, sheet = sheetname)
  d2["cluster"] = sheetname
  #d2["patient"] = sheetname
  d2 = subset(d2, blood1.avgFC >= 2 & tissue.avgFC >= 2)
  d1 <- rbind(d1,d2)
}

intersect(d1$gene,
          lung_merge$Gene_symbol)
head(d1$gene)

rm(list_of_lung)
list_of_lung <- lapply(X = clusterlist, FUN = function(x){
  x <- read_excel(filen, sheet = x)
  x <- subset(x, blood1.avgFC >= 2 & tissue.avgFC >= 2)
  x <- x[,c(1,2,3,4,5)]
  x <- left_join(x, lung_merge, c("gene" = "Gene_symbol"))
})

names(list_of_lung) = clusterlist
list_of_lung

write.xlsx(list_of_lung, file = "./output/SMC025_B1Tissue_Conserved_Markers_lineage_cancer_annotation.xlsx")

#### Add Gene alias ####
lung_merge = read.csv("lung_join_cancer_lineage.csv")
pancreas_merge = read.csv("pancreas_join_cancer_lineage.csv")

head(lung_merge)
head(pancreas_merge)

#lung_lineage, lung_cancer
anno_lung_lineage = select(org.Hs.eg.db, keys=lung_lineage, columns=cols, keytype="SYMBOL")
anno_lung_cancer = select(org.Hs.eg.db, keys=lung_cancer, columns=cols, keytype="SYMBOL")

lung_union = union(lung_lineage, lung_cancer)
anno_lung_union = select(org.Hs.eg.db, keys=lung_union, columns=cols, keytype="SYMBOL")


aggr_lung_union = aggregate(anno_lung_union$ALIAS, list(anno_lung_union$SYMBOL), FUN = paste, collapse =",")
head(aggr_lung_union)
tail(aggr_lung_union)

colnames(aggr_lung_union) = c("Gene", "Gene_alias")
library(dplyr)
lung_join_alias = left_join(lung_merge, aggr_lung_union, c("Gene_symbol" = "Gene"))

#pancreas_lineage,pancreas_cancer
anno_panc_lineage = select(org.Hs.eg.db, keys=pancreas_lineage, columns=cols, keytype="SYMBOL")
anno_panc_cancer = select(org.Hs.eg.db, keys=pancreas_cancer, columns=cols, keytype="SYMBOL")

pancreas_union = union(pancreas_lineage, pancreas_cancer)
anno_panc_union = select(org.Hs.eg.db, keys=pancreas_union, columns=cols, keytype="SYMBOL")
aggr_panc_union = aggregate(anno_panc_union$ALIAS, list(anno_panc_union$SYMBOL), FUN = paste, collapse=",")

head(aggr_panc_union)
colnames(aggr_panc_union) = c("Gene", "Gene_alias")
panc_join_alias = left_join(pancreas_merge, aggr_panc_union, c("Gene_symbol" = "Gene"))
head(panc_join_alias)

write.csv(lung_join_alias, "lung_join_cancer_lineage_alias.csv", row.names = FALSE)
write.csv(panc_join_alias, "pancreas_join_cancer_lineage_alias.csv", row.names = FALSE)
