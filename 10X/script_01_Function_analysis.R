library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(dplyr)
library(vroom)
##----------------------------------------------------------------------------------##
## Preprocessing - read data
##----------------------------------------------------------------------------------##
setwd("H:TCGA")
setwd("./SingleCell/")
rm(d)
#d <- read.table(file = "SMC009.txt", sep = "\t")
d <- read.csv("SMC009.csv")
head(d)

### cl4
cl4 <- d[d$cluster == 4,]
head(cl4)
table(cl4$cluster) # n = 246

### cl5
cl5 <- d[d$cluster == 5,]
head(cl5)
table(cl5$cluster) # n = 541

### cl8
cl8 <- d[d$cluster == 8,]
head(cl8)
table(cl8$cluster) # n = 127


### cl14
cl14 <- d[d$cluster == 14,]
head(cl14)
table(cl14$cluster) # n = 255

##----------------------------------------------------------------------------------##
## read data
##----------------------------------------------------------------------------------##

### Set Working Directory ###

setwd("../")
getwd()
#dir.create("cluster4")
#dir.create("cluster5")
#dir.create("cluster8")
#dir.create("cluster14")
setwd("../")
#++++++++++++++++++++++++++#
setwd("./Funcion2/cluster4")
#++++++++++++++++++++++++++#
getwd()

rm(df)
rm(geneList)
rm(gene)

#### Set cluster dataframe ####
#### cl4 ; cl5 ; cl8 ; cl14 
#++++++++++++++++++++++++++#
df <- cl4

#++++++++++++++++++++++++++#
geneList <- df[,3]
head(geneList)

names(geneList) <- as.character(df[,1])

geneList <- sort(geneList, decreasing = TRUE)
#gene <-names(geneList)

#### choose Foldchange cut-off ####
gene <- names(geneList)[abs(geneList) >= 0.58]
gene

bkgd.genes <- df[,1]
length(bkgd.genes)

bkgd.genes.entrez <- bitr(bkgd.genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##
##up.genes <- cldf[cldf$avg_logFC > 1 & cldf$p_val_adj < 0.05, 1] 
#up.genes <- cldf[cldf$avg_logFC > 0.58 & cldf$p_val_adj < 0.05, 1] 
#length(up.genes)
#
##dn.genes <- cldf[cldf$avg_logFC < -1 & cldf$p_val_adj < 0.05, 1]
#dn.genes <- cldf[cldf$avg_logFC < -0.58 & cldf$p_val_adj < 0.05, 1]
#length(dn.genes)
#
#bkgd.genes <- cldf[,1]
#length(bkgd.genes)
#
#up.genes.entrez <- clusterProfiler::bitr(up.genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
#dn.genes.entrez <- bitr(dn.genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
#bkgd.genes.entrez <- bitr(bkgd.genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
#
#head(up.genes.entrez)
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##

### Entrez ID mapping 
library(org.Hs.eg.db)
gene.entrez <- bitr(gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)


## Remove any Entrez duplicates
gene.entrez <- gene.entrez[which(duplicated(gene.entrez$ENTREZID) == F), ]
res_entrez <- gene.entrez$ENTREZID
res_entrez <- sort(res_entrez, decreasing = TRUE)


##----------------------------------------------------------------------------------##
## Analysis Start
##----------------------------------------------------------------------------------##
getwd()
##----------------------------------------------------------------------------------##
### GO Analysis 
#### ego3, GeneSymbol Input ; Benjamini-Hotchberg; pvalue cut-off 0.01, q-value cut-off 0.05
dir.create("./enrichGO")
setwd("./enrichGO/")

ego3 <- enrichGO(gene         = names(geneList),
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'SYMBOL',
                 #universe = as.character(bkgd.genes),
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05)
head(ego3)
#write.csv(ego3@result,file = "enrichGO.csv")

filegp=sprintf("goplot%s.png",  Sys.Date())
png(filename = filegp, height=550, width=550, bg="white")
goplot(ego3)
dev.off()

gsego <-gseGO(geneList     = res_entrez,
              OrgDb        = org.Hs.eg.db,
              ont          = "CC",
              nPerm        = 1000,
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)

#write.csv(gsego@result,file = "gseGO.csv")

##----------------------------------------------------------------------------------##
### msigDB Analysis
#### C6 oncogenic :emC6, gsC6   C7 immune : emC7, gsC7
library(msigdbr)
dir.create("../MSigDb")
setwd("../MSigDb/")
msigdbr_show_species()
m_df <- msigdbr(species = "Homo sapiens")
head(m_df, 2) %>% as.data.frame


#### C2: curated gene sets

m_t2g <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, gene_symbol)
head(m_t2g)

emC2 <- enricher(gene, 
                 #universe = as.character(bkgd.genes), 
                 TERM2GENE = m_t2g, )
head(emC2)

gsC2 <- GSEA(geneList, TERM2GENE = m_t2g)
head(gsC2)

#write.csv(emC5@result,file = "enrich_msigGO_result.csv")
#write.csv(gsC5@result,file = "GSEA_msigGO_result.csv")




#### C5: GO gene sets

m_t5g <- msigdbr(species = "Homo sapiens", category = "C5") %>% 
  dplyr::select(gs_name, gene_symbol)
head(m_t5g)

emC5 <- enricher(gene, 
#                 universe = as.character(bkgd.genes), 
                 TERM2GENE = m_t5g, )
head(emC5)

gsC5 <- GSEA(geneList, TERM2GENE = m_t5g)
head(gsC5)

#write.csv(emC5@result,file = "enrich_msigGO_result.csv")
#write.csv(gsC5@result,file = "GSEA_msigGO_result.csv")


#### C6: Oncogenic
#colnames(msigdbr(species = "Homo sapiens", category = "C6"))
m_t6g <- msigdbr(species = "Homo sapiens", category = "C6") %>% 
  dplyr::select(gs_name, gene_symbol)
head(m_t6g)

emC6 <- enricher(gene, 
                 #universe = as.character(bkgd.genes), 
                 TERM2GENE = m_t6g)
head(emC6)

gsC6 <- GSEA(geneList, TERM2GENE = m_t6g)
head(gsC6)

#write.csv(emC6@result,file = "enrich_msigoncogenic_result.csv")
#write.csv(gsC6@result,file = "GSEA_msigoncogenic_result.csv")


#### C7: immune

m_t7g <- msigdbr(species = "Homo sapiens", category = "C7") %>% 
  dplyr::select(gs_name, gene_symbol)
head(m_t7g)

emC7 <- enricher(gene, TERM2GENE = m_t7g, universe = as.character(bkgd.genes))
head(emC7)

gsC7 <- GSEA(geneList, TERM2GENE = m_t7g)
head(gsC7)

#write.csv(emC7@result,file = "enrich_msigimmune_result.csv")
#write.csv(gsC7@result,file = "GSEA_msigimmune_result.csv")

##----------------------------------------------------------------------------------##
### WikiPathways Analysis
#### 

library(magrittr)
library(clusterProfiler)

dir.create("./wiki")
setwd("./wiki/")

wp.hs.gmt <- rWikiPathways::downloadPathwayArchive(organism="Homo sapiens", format = "gmt")

wp2gene <- clusterProfiler::read.gmt(wp.hs.gmt)
wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid,gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid,name) #TERM2NAME
wpid2gene
wpid2name


ewp.up <- clusterProfiler::enricher(
  res_entrez,
  #universe = bkgd.genes.entrez[[2]],
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.1, #p.adjust cutoff; relaxed for demo purposes
  TERM2GENE = wpid2gene,
  TERM2NAME = wpid2name)

head(ewp.up)

ewp.up <- DOSE::setReadable(ewp.up, org.Hs.eg.db, keytype = "ENTREZID")
head(ewp.up)

# ewp.dn <- enricher(
#   dn.genes.entrez[[2]],
#   #universe = bkgd.genes[[2]],  #hint: comment out to get any results for demo
#   pAdjustMethod = "fdr",
#   pvalueCutoff = 0.1, #p.adjust cutoff; relaxed for demo purposes
#   TERM2GENE = wpid2gene,
#   TERM2NAME = wpid2name)

# ewp.dn <- setReadable(ewp.dn, org.Hs.eg.db, keytype = "ENTREZID")
# head(ewp.dn)
# dotplot(ewp.dn, showCategory = 20)
# 
# 
# ### GSEA WikiPathways Analysis
# lung.expr$fcsign <- sign(lung.expr$log2FC)
# lung.expr$logfdr <- -log10(lung.expr$P.Value)
# lung.expr$sig <- lung.expr$logfdr/lung.expr$fcsign
# sig.lung.expr.entrez<-merge(lung.expr, bkgd.genes.entrez, by.x = "GeneID", by.y = "ENSEMBL")
# gsea.sig.lung.expr <- sig.lung.expr.entrez[,8]
# names(gsea.sig.lung.expr) <- as.character(sig.lung.expr.entrez[,9])
# gsea.sig.lung.expr <- sort(gsea.sig.lung.expr,decreasing = TRUE)
# 
# gwp.sig.lung.expr <- clusterProfiler::GSEA(
#   gsea.sig.lung.expr,
#   pAdjustMethod = "fdr",
#   pvalueCutoff = 0.05, #p.adjust cutoff
#   TERM2GENE = wpid2gene,
#   TERM2NAME = wpid2name)
# 
# gwp.sig.lung.expr.df = data.frame(ID=gwp.sig.lung.expr$ID,
#                                   Description=gwp.sig.lung.expr$Description,
#                                   enrichmentScore=gwp.sig.lung.expr$enrichmentScore,
#                                   NES=gwp.sig.lung.expr$NES,
#                                   pvalue=gwp.sig.lung.expr$pvalue,
#                                   p.adjust=gwp.sig.lung.expr$p.adjust,
#                                   rank=gwp.sig.lung.expr$rank,
#                                   leading_edge=gwp.sig.lung.expr$leading_edge
# )
# gwp.sig.lung.expr.df[which(gwp.sig.lung.expr.df$NES > 1),] #pathways enriched for upregulated lung cancer genes
# gwp.sig.lung.expr.df[which(gwp.sig.lung.expr.df$NES < -1),] #pathways enriched for downregulated lung cancer genes


##----------------------------------------------------------------------------------##
### enrichKEGG
kk1 <- enrichKEGG(gene     = res_entrez,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)
head(kk1)

kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)


library("pathview")
hsa04612 <- pathview(gene.data  = kk1@result$geneID,
                     pathway.id = "hsa04926",
                     species    = "hsa",
                     limit      = list(gene=max(abs(geneList)), cpd=1))

#hsa05323 https://www.genome.jp/kegg-bin/show_pathway?hsa05323 안될경우 직접 xml 다운로드 
##----------------------------------------------------------------------------------##
# ### MeSH Enrichment Analysis
# library(meshes)
# 
# emesh <- enrichMeSH(gene, MeSHDb = "MeSH.Hsa.eg.db", database='gene2pubmed', category = 'C')
# head(emesh)
# 
# y <- gseMeSH(geneList, MeSHDb = "MeSH.Hsa.eg.db", database = 'gene2pubmed', category = "G")
##----------------------------------------------------------------------------------##
### enrichDO Analysis
#### edo1


library(DOSE)

edo1 <- enrichDO(gene = res_entrez,
              ont = "DO",
              #universe = bkgd.genes.entrez[[2]],
              pvalueCutoff = 0.05,
              pAdjustMethod = "BH",
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.05,
              readable      = FALSE)

head(edo1)

dir.create("../enrichDO")
setwd("../enrichDO/")
getwd()

#write.csv(edo1@result,file = "enrichDOresult.csv")

##----------------------------------------------------------------------------------##
### ncg Analysis
#### ncg

#ncg <- enrichNCG(res_entrez, universe = bkgd.genes.entrez[[2]])
ncg <- enrichNCG(res_entrez)
head(ncg)
dir.create("../NCG")
setwd("../NCG/")
#write.csv(ncg@result,file = "NCG_result.csv")

##----------------------------------------------------------------------------------##
### dgn Analysis
#### dgn

#dgn <- enrichDGN(res_entrez, universe = bkgd.genes.entrez[[2]])
dgn <- enrichDGN(res_entrez)
head(dgn)
dir.create("../enrichDGN")
setwd("../enrichDGN/")


##----------------------------------------------------------------------------------##
### gseDGN Analysis

gsdgn <- gseDGN(res_entrez,
              nPerm         = 100,
              minGSSize     = 120,
              pvalueCutoff  = 0.2,
              pAdjustMethod = "BH",
              verbose       = FALSE)
#gsdgn <- setReadable(gsdgn, 'org.Hs.eg.db')
head(gsdgn, 3)




##----------------------------------------------------------------------------------##
### Cell Marker
setwd("../")
getwd()
library(vroom)
cell_markers <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt') %>%
  tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>% 
  dplyr::select(cellMarker, geneID) %>%
  dplyr::mutate(geneID = strsplit(geneID, ', '))

cell_markers

y <- enricher(res_entrez, 
              TERM2GENE=cell_markers, 
#              universe = bkgd.genes.entrez[[2]],
              minGSSize=1)

y <- setReadable(y, 'org.Hs.eg.db', 'ENTREZID')
DT::datatable(as.data.frame(y))



