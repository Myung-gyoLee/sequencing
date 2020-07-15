setwd("H:/TCGA/TCGA_marker")
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
columns(org.Hs.eg.db)

# read analysis result file

library(readxl)
filen = "./data/SMC024-025-027_SEV003_Marker_Annotation.xlsx"
clusterlist <- excel_sheets(filen)

clusterlist

sheetname = clusterlist[1]
d2 <- read_excel(filen, sheet = sheetname)
head(d2)
# Cluster Markers Gene  Gene.Name GOTERM_BP_DIRECT GOTERM_CC_DIRECT GOTERM_MF_DIRECT
# <dbl> <chr>   <chr> <chr>     <chr>            <chr>            <chr>           
#   1       0 Conser~ PPBP  pro-plat~ GO:0002523~leuk~ GO:0005576~extr~ GO:0005355~gluc~
#   2       0 Conser~ NRGN  neurogra~ GO:0007165~sign~ GO:0005634~nucl~ GO:0005516~calm~

d2$Gene



require(openxlsx)
list_of_lung <- lapply(X = clusterlist[1:3], FUN = function(x){
  x <- read_excel(filen, sheet = x)
  #x <- subset(x, blood1.avgFC >= 1.5 & tissue.avgFC >= 1.5)
  x <- x[,c(1,2,3,4)]
  x <- left_join(x, lung_merge, c("Gene" = "Gene_symbol"))
})

list_of_pancreas <- lapply(X = clusterlist[4], FUN = function(x){
  x <- read_excel(filen, sheet = x)
  #x <- subset(x, blood1.avgFC >= 1.5 & tissue.avgFC >= 1.5)
  x <- x[,c(1,2,3,4)]
  x <- left_join(x, pancreas_merge, c("Gene" = "Gene_symbol"))
})

list_of_lung
names(list_of_lung) = clusterlist[1:3]
names(list_of_pancreas) = clusterlist[4]

write.xlsx(list_of_lung, file = "./output/SMC024-025-027_B1Tissue_Marker_Annotation_add_lineage_cancer_annotation.xlsx")

write.xlsx(list_of_pancreas, file = "./output/SEV003_B1Tissue_Marker_Annotation_add_lineage_cancer_annotation.xlsx")

#
pancreas_union = union(pancreas_lineage,pancreas_cancer)
anno_panc_union = select(org.Hs.eg.db, keys=pancreas_union, columns=cols, keytype="SYMBOL")

View(anno_panc_union)

lung_union = union(lung_lineage, lung_cancer)
anno_lung_union = select(org.Hs.eg.db, keys=lung_union, columns=cols, keytype="SYMBOL")

View(anno_lung_union)
## confirm overlap

library(readxl)
filen = "./data/SMC024-025-027_SEV003_Marker_Annotation.xlsx"
clusterlist <- excel_sheets(filen)


d1=NULL
for(i in clusterlist[1:3]){
  #cluster_num = i
  #sheetname = sprintf("Cluster_%s",cluster_num)
  sheetname = i
  d2 <- read_excel(filen, sheet = sheetname)
  #d2 = subset(d2, blood1.avgFC >= 1.5 & tissue.avgFC >= 1.5)
  d1 <- rbind(d1,d2)
}


d1$Gene


inter_lung = intersect(d1$Gene,
          anno_lung_union$ALIAS)

intersect(inter_lung, lung_union)
setdiff(inter_lung, lung_union)


# pancreas

d1=NULL
for(i in clusterlist[4]){
  #cluster_num = i
  #sheetname = sprintf("Cluster_%s",cluster_num)
  sheetname = i
  d2 <- read_excel(filen, sheet = sheetname)
  #d2 = subset(d2, blood1.avgFC >= 1.5 & tissue.avgFC >= 1.5)
  d1 <- rbind(d1,d2)
}


d1$Gene


inter_pancreas = intersect(d1$Gene,
                       anno_panc_union$ALIAS)

intersect(inter_pancreas, pancreas_union)
setdiff(inter_pancreas, pancreas_union)
