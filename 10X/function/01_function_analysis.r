#https://yulab-smu.github.io/clusterProfiler-book/
#https://guangchuangyu.github.io/2016/01/go-analysis-using-clusterprofiler/
#Title: statistical analysis and visualization of functional profiles for genes and gene clusters
#Version: 3.14.3

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
  df <- d[d$cluster == cluster_num,] # df <- d
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
  bkgd.genes <- df[,1]
  cat("print length of background genes")
  cat("--------------------------------------------")
  length(bkgd.genes)
  bkgd.genes.entrez <- bitr(bkgd.genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
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


setwd("H:SingleCell/function_smc024")

#d <- read.csv("smc024c11c19.csv")
d <- read.csv("smc024_conserved.csv", row.names = 1)
#rm(d)
head(d)
colnames(d)[1] = "gene"
cat(colnames(d))
# colnames(d)=c("gene","gene_name","blood1.avgFC","tissue.avgFC","avgVolume","blood1.data_p_val","blood1.data_avg_logFC","blood1.data_pct.1","blood1.data_pct.2",
#               "blood1.data_p_val_adj","tissue.data_p_val","tissue.data_avg_logFC","tissue.data_pct.1","tissue.data_pct.2","tissue.data_p_val_adj","max_pval","minimump_p_val","cluster")

### Get geneList << input (d, cluster_num) ###
# 2, 5, 7, 14, 15, 17, 22
cluster_num = 22
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
if("cluster" %in% getwd()){
  nfol = sprintf("../cluster%s",cluster_num)
}else{nfol = sprintf("./cluster%s",cluster_num)}


dir.create(nfol)
setwd(nfol)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### Gene annotations
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
# gene<-rownames(results)
# gene.df <- bitr(gene, fromType = "SYMBOL",
#                 toType = c("ENSEMBL", "ENTREZID"),
#                 OrgDb = org.Hs.eg.db)

#results_anno <- geneAnnotations(input=results, keys=row.names(results), column=c("ENTREZID", "ENSEMBL"), keytype="SYMBOL", organism = "human")

ego <- enrichGO(gene         = gene,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'SYMBOL',
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)
head(summary(ego))


egobp <- enrichGO(gene         = gene,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = 'SYMBOL',
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)


egomf <- enrichGO(gene         = gene,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = 'SYMBOL',
                  ont           = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)

### MSigDb analysis
#BiocManager::install("msigdbr")
library(msigdbr)
msigdbr_show_species()

m_df <- msigdbr(species = "Homo sapiens")
head(m_df, 2) %>% as.data.frame

# C2

m_t2g <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, gene_symbol)
head(m_t2g)

em2 <- enricher(gene, TERM2GENE = m_t2g)
head(em2)

gs2 <- GSEA(geneList, TERM2GENE = m_t2g)
head(gs2)

# C5

m_t5g <- msigdbr(species = "Homo sapiens", category = "C5") %>% 
  dplyr::select(gs_name, gene_symbol)
head(m_t5g)

em5 <- enricher(gene, TERM2GENE = m_t5g)
head(em5)

gs5 <- GSEA(geneList, TERM2GENE = m_t5g)
head(gs5)



# C6: Oncogenic
#colnames(msigdbr(species = "Homo sapiens", category = "C6"))
m_t6g <- msigdbr(species = "Homo sapiens", category = "C6") %>% 
  dplyr::select(gs_name, gene_symbol)
head(m_t6g)


em6 <- enricher(gene, TERM2GENE = m_t6g)
head(em6)

gs6 <- GSEA(geneList, TERM2GENE = m_t6g)
head(em2)

#write.csv(em1@result,file = "enrichrresult.csv")
#write.csv(em2@result,file = "GSEA_oncogenic_result.csv")


# C7 immune

m_t7g <- msigdbr(species = "Homo sapiens", category = "C7") %>% 
  dplyr::select(gs_name, gene_symbol)
head(m_t7g)


em7 <- enricher(gene, TERM2GENE = m_t7g)
head(em7)

gs7 <- GSEA(geneList, TERM2GENE = m_t7g)
head(gs7)



## DO analysis
# https://yulab-smu.github.io/clusterProfiler-book/chapter4.html#enrichdo-function
# https://www.rdocumentation.org/packages/DOSE/versions/2.10.6/topics/enrichDO

### make Input 
# https://bioconductor.org/packages/release/data/annotation/manuals/org.Hs.eg.db/man/org.Hs.eg.db.pdf

head(res_entrez)
nrow(res_entrez)

#colnames(Id_Entrez)
#[1] "Gene"     "EntrezID"

### enrichDO
library(DOSE)

edo <- enrichDO(gene = res_entrez,
                ont = "DO",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                minGSSize     = 5,
                maxGSSize     = 500,
                qvalueCutoff  = 0.05,
                readable      = FALSE)


### enrichNCG : Network of Cancer Gene (NCG)(A. et al. 2016)
ncg <- enrichNCG(res_entrez)
head(ncg)
#write.csv(ncg@result,file = "NCG_result.csv")



### enrich DGN ; DisGeNET(Janet et al. 2015) gene-disease associations
dgn <- enrichDGN(res_entrez)
head(dgn)
#write.csv(dgn@result,file = "dgnresult.csv")


# unique(dgn$Description)
# Desc_Pav = data.frame(dgn$Description, dgn$p.adjust)
# 
# 
# #### sory by adjust P-value
# Desc_Pav = Desc_Pav[order(Desc_Pav$dgn.p.adjust),]
# head(Desc_Pav,10)

# #Order = grep("Breast",Desc_Pav$dgn.Description)
# #str(Order)
# #Disease = Desc_Pav$dgn.Description[Order]
# #P_value = Desc_Pav$dgn.p.adjust[Order]
# 
# #dgn_Breast = data.frame(Order, Disease)
# #dgn_Breast["P_adj"]=P_value
# 
# #head(dgn_Breast)



### kegg
# https://yulab-smu.github.io/clusterProfiler-book/chapter5.html


gkk <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(gkk)

ekk <- enrichKEGG(gene     = res_entrez,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)



ekk@result$Description
ekk@result$ID %>% sort
head(ekk)

filekg=sprintf("keggresult%s.csv",Sys.Date())
write.csv(ekk@result,file = filekg)

ekkR <- setReadable(ekk, 'org.Hs.eg.db', 'ENTREZID')
#ekkR
filekg=sprintf("keggresult%s.csv",Sys.Date())
write.csv(ekkR@result,file = filekg)

ekklist=ekk@result$ID[abs(ekk@result$p.adjust) < 0.05]
keggabs = sprintf("%s  %s  %s",ekklist, ekk@result$Description[abs(ekk@result$p.adjust) < 0.05], ekkR@result$geneID[abs(ekkR@result$p.adjust) < 0.05])
keggabs
fileakg=sprintf("kegg_adjpvaue%s.csv",Sys.Date())
write.csv(keggabs,file = fileakg)


ekk@result$p.adjust

egeneL = geneList
#str(egeneL)
egeneL[1:length(egeneL)] = 1
#egeneL

library("pathview")
for(i in ekklist){
  pathview(gene.data  = res_entrez,
           pathway.id = i,
           species    = "hsa",
           limit      = list(gene=max(abs(egeneL)), cpd=1))
}
# hsa04141 <- pathview(gene.data  = res_entrez,
#                       pathway.id = "hsa04141",
#                       species    = "hsa",
#                       limit      = list(gene=max(abs(geneList)), cpd=1))

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# save(list=ls(), file='SMC009.Rdata')

### Gene Set Enrichment Analysis
# ego -> ego Cellular Component
head(ego)
# Biological Process
head(egobp)
# Molecular Function
head(egomf)
head(em2)
head(gs2)
head(em5)
head(gs5)
head(em6)
head(gs6)
head(em7)
head(gs7)
head(edo)
head(ncg)
head(dgn)
head(gkk)
head(ekk)


## Cell Marker
#BiocManager::install("vroom")
library(dplyr)
library(vroom)
cell_markers <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt') %>%
  tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>% 
  dplyr::select(cellMarker, geneID) %>%
  dplyr::mutate(geneID = strsplit(geneID, ', '))
cell_markers

y <- enricher(res_entrez, TERM2GENE=cell_markers, minGSSize=1)
y <- setReadable(y, 'org.Hs.eg.db', 'ENTREZID')
DT::datatable(as.data.frame(y))


#### validation ####
# cst1 = "SPP1/APOE/CD74/SRGN/TYROBP/ALOX5AP/APOC1/FCER1G/C1QA/CTSB/FTL/CTSL/CCL3/FOS/C1QB/LAPTM5/CD14/AIF1/B2M/FCGR3A/CD68/GLUL/C1QC/CYBA/CTSS/HLA-B/FABP5/NPC2/CYBB/RGS1/GPR183/GRN/MRC1/CTSD/RNASET2/SGK1/FCGR2A/PFN1/ITGB2/LST1/CCL4L2/CCL4/LCP1/HLA-E/C5AR1/ASAH1/HLA-C/DUSP1/PTPRC/IGSF6/MNDA/NR4A2/JUNB/OLR1/CD83/BRI3/PLXDC2/PYCARD/TREM2/DAB2/RGS2/LIMS1/PSAP/PLEK/BTG2/SAT1/CD53/ARHGDIB/CD86/CTSC/IFNGR1/GNAI2/SERPINA1/TLN1/CSF1R"
# cstsp1 = strsplit(cst1, split = "/")
# intersect(cstsp1, gene)
# length(intersect(cstsp1, gene))
# str(gene)
# str(cstsp1)
# cstsp1=cstsp1[[1]]
