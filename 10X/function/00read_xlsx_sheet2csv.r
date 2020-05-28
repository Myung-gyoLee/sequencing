library(readxl)
sheetname = sprintf("Cluster_%s",cluster_num)
d1 <- read_excel("SMC024_B1Tissue_Conserved_Markers.xlsx", sheet = sheetname)
#d <- read.csv("smc025b1t.csv")
#rm(d)
d1["cluster"] = cluster_num
head(d1)

d = subset(d1, blood1.avgFC >= 2 & tissue.avgFC >= 2)
d$gene
#cat(colnames(d))
#colnames(d)=c("gene","gene_name","blood1.avgFC","tissue.avgFC","avgVolume")
