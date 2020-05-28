###################################################
gene = c("JUN", "PLCG2", "ITGA2")
gene = c("FOS", "NCF2", "FCGR2A")
gene = c("ID1", "BMPR2", "ID3")
gene = c("RAMP3", "RAMP2", "CALCRL")
gene = c("ITGA5", "VIM", "BMPR2", "ZEB1")
gene = c("TNFSF10", "MX1", "STAT1")
gene = c("CXCL1", "JUN", "CXCL2")
gene = c("JUN", "PLCG2", "ITGA2")
#TNFSF10, MX1, STAT1
#hsa05132:Salmonella infection   CXCL1/JUN/CXCL2
#hsa05200:Pathways in cancer     JUN/PLCG2/ITGA2
###################################################
#library(org.Hs.eg.db)
#library(pathview)
gene.entrez <- bitr(gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
## Remove any Entrez duplicates
gene.entrez <- gene.entrez[which(duplicated(gene.entrez$ENTREZID) == F), ]
res_entrez <- gene.entrez$ENTREZID
res_entrez <- sort(res_entrez, decreasing = TRUE)
#res_entrez
cat(res_entrez)


egeneL = gene.entrez$ENTREZID
egeneL[1:length(egeneL)] = 1
#str(egeneL)
names(egeneL) = gene.entrez$ENTREZID


egeneL


pathview(gene.data  = gene.entrez$ENTREZID,
         pathway.id = "hsa05206",
         species    = "hsa",
         limit      = list(gene = 1, cpd = 1))

h = read.csv("cancerpathway.csv")
########################

library(org.Hs.eg.db)
library(DOSE)
library(pathview)
library(clusterProfiler)
library(AnnotationHub)
library(ensembldb)
library(tidyverse)


library(AnnotationHub)
hub <- AnnotationHub()
query(hub, "Homosapiens")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
### Gene mapping ###
gsymbol_input <- function(d,cluster_num, dfvalue){
  #df <- d[d$cluster == cluster_num,] # 
  df <- d
  cat("print head of data frame")
  cat("--------------------------------------------")
  head(df)
  table(df$cluster)
  geneList <- df[,dfvalue]
  cat("print all gene list")
  cat("--------------------------------------------")
  head(geneList)
  names(geneList) <- as.character(df[,1])
  geneList <- sort(geneList, decreasing = TRUE)
  return(geneList)
}
### heatplot 값 선택 
# df <- d[d$cluster == cluster_num,]
# geneList <- df[,5]
# colnames(df)
# geneList

### Entrez ID mapping ###
entrez_input <- function(gene){
  library(org.Hs.eg.db)
  gene.entrez <- bitr(gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
  ## Remove any Entrez duplicates
  gene.entrez <- gene.entrez[which(duplicated(gene.entrez$ENTREZID) == F), ]
  res_entrez <- gene.entrez$ENTREZID
  res_entrez <- sort(res_entrez, decreasing = TRUE)
  return(res_entrez)
}
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

setwd("H:SingleCell/function_smc025")

d <- read.csv("smc025b1t.csv")


#setwd("H:SingleCell/function_smc024")

#d <- read.csv("smc024c11c19.csv")
#d <- read.csv("smc024_conserved.csv", row.names = 1)
#rm(d)
head(d)
colnames(d)[1] = "gene"
cat(colnames(d))
# colnames(d)=c("gene","gene_name","blood1.avgFC","tissue.avgFC","avgVolume","blood1.data_p_val","blood1.data_avg_logFC","blood1.data_pct.1","blood1.data_pct.2",
#               "blood1.data_p_val_adj","tissue.data_p_val","tissue.data_avg_logFC","tissue.data_pct.1","tissue.data_pct.2","tissue.data_p_val_adj","max_pval","minimump_p_val","cluster")

### Get geneList << input (d, cluster_num) ###
# SMC024 2, 5, 7, 11, 14, 15, 17, 19, 22
# SMC025 7 12 14 18
unique(d$cluster)
cluster_num = 7
# SMC024 -cluster 14 SMC025-cluster7 :  no result
# cut off column !!!
dfvalue = 5
geneList <- gsymbol_input(d,cluster_num, dfvalue)

#gene <- names(geneList)[abs(geneList) >= 0.58]
gene <- names(geneList)

cat("print cut-off gene list\n")
cat("--------------------------------------------\n")
gene

### Entrez ID mapping ###
res_entrez <- entrez_input(gene)



for (i in res_entrez){
  cat(i, "red, black","\n"  )
}

# for (i in res_entrez){
#   cat(i,"\n"  )
# }

