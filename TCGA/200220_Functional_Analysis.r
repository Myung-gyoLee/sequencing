#https://yulab-smu.github.io/clusterProfiler-book/
#https://guangchuangyu.github.io/2016/01/go-analysis-using-clusterprofiler/


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

summary(res)
mcols(res)$description
# order results table by the smallest adjusted p value:
res <- res[order(res$padj),]

results = as.data.frame(dplyr::mutate(as.data.frame(res), sig=ifelse(res$padj<0.05, "FDR<0.05", "Not Sig")), row.names=rownames(res))
head(results)

DEgenes_DESeq <- results[which(abs(results$log2FoldChange) > log2(1.5) & results$padj < 0.05),]

### Volcano Plot
p = ggplot2::ggplot(results, ggplot2::aes(log2FoldChange, -log10(pvalue))) +
  ggplot2::geom_point(ggplot2::aes(col = sig)) +
  ggplot2::scale_color_manual(values = c("red", "black")) +
  ggplot2::ggtitle("Volcano Plot of DESeq2 analysis")

p + ggrepel::geom_text_repel(data=results[1:10, ], ggplot2::aes(label=rownames(results[1:10, ])))

# If there aren't too many DE genes:
#p + geom_text_repel(data = dplyr::filter(results, padj<0.05), aes(label = rownames(results[1:10, ])))

### MA-plot

DESeq2::plotMA(res, main="MA Plot", ylim=c(-2,2))

### Gene annotations
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
gene<-rownames(results)
gene.df <- bitr(gene, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)

#results_anno <- geneAnnotations(input=results, keys=row.names(results), column=c("ENTREZID", "ENSEMBL"), keytype="SYMBOL", organism = "human")

ego3 <- enrichGO(gene         = gene,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
head(summary(ego3))

### Visualization functions 

# Landscaped PDF size to 8 x 14
plotdir = "./02Breast/023report/2.GTExNormal/Function/"
fileeg=sprintf("%s/Dotplot%s.png", plotdir,Sys.Date())
png(filename = fileeg, height=800, width=500, bg="white")

dotplot(ego3, showCategory=30)
dev.off()

#categorySize="pvalue", showCategory = 5, foldChange=OE_foldchanges,vertex.label.font=6
# PDF size to 24 x 32

fileeg=sprintf("%s/cnetplot%s.png", plotdir,Sys.Date())
png(filename = fileeg, height=800, width=500, bg="white")


cnetplot(ego3, showCategory = 3)
dev.off()



## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
# PDF size to 24 x 32
emapplot(ego3, showCategory = 50)

plotGOgraph(ego3)

### Gene Set Enrichment Analysis
gsecc <- gseGO(geneList=geneList, ont="CC", OrgDb=org.Hs.eg.db, verbose=F)
head(summary(gsecc))
gseaplot(gsecc, geneSetID="GO:0031012")

## Cell Marker
#BiocManager::install("vroom")
library(vroom)
cell_markers <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt') %>%
  tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>% 
  dplyr::select(cellMarker, geneID) %>%
  dplyr::mutate(geneID = strsplit(geneID, ', '))
cell_markers

y <- enricher(gene, TERM2GENE=cell_markers, minGSSize=1)
DT::datatable(as.data.frame(y))

### MSigDb analysis
#BiocManager::install("msigdbr")
library(msigdbr)
msigdbr_show_species()

m_df <- msigdbr(species = "Homo sapiens")
head(m_df, 2) %>% as.data.frame

# C6: Oncogenic, C7 immune
#colnames(msigdbr(species = "Homo sapiens", category = "C6"))
m_t2g <- msigdbr(species = "Homo sapiens", category = "C7") %>% 
  dplyr::select(gs_name, gene_symbol)
head(m_t2g)


em1 <- enricher(gene, TERM2GENE = m_t2g)
head(em1)

head(results)
rm(geneList)
geneList <- results[,2]
names(geneList) <- as.character(rownames(results))
geneList <- sort(geneList, decreasing = TRUE)

em2 <- GSEA(geneList, TERM2GENE = m_t2g)
head(em2)

## DO analysis
# https://yulab-smu.github.io/clusterProfiler-book/chapter4.html#enrichdo-function
# https://www.rdocumentation.org/packages/DOSE/versions/2.10.6/topics/enrichDO

### make Input 
# https://bioconductor.org/packages/release/data/annotation/manuals/org.Hs.eg.db/man/org.Hs.eg.db.pdf
library(org.Hs.eg.db)

xx <- as.list(org.Hs.egALIAS2EG)
head(xx)

I_Entrez = xx[gene]

rm(Id_Entrez)
Id_Entrez = data.frame(matrix(nrow=length(gene)), ncol=NA)
colnames(Id_Entrez)=c("Gene","EntrezID")
Id_Entrez["Gene"] = as.character(gene)

head(Id_Entrez)


i=1
while(i < nrow(Id_Entrez)){
  if(is.null(I_Entrez[[Id_Entrez["Gene"][i,]]][1]) == FALSE){
    print(i)
    Id_Entrez["EntrezID"][i,] = I_Entrez[[Id_Entrez["Gene"][i,]]][1]
    i = i+1
  }else {
    print(i)
    Id_Entrez["EntrezID"][i,] <- NA
    i = i+1
    }
}

head(Id_Entrez)
tail(Id_Entrez)
sum(is.na(Id_Entrez))

Id_Entrez = na.omit(Id_Entrez)
sum(is.na(Id_Entrez))


## Remove any Entrez duplicates
res_entrez <- Id_Entrez[which(duplicated(Id_Entrez$EntrezID) == F), ]

head(res_entrez)
nrow(res_entrez)
nrow(Id_Entrez)
#colnames(Id_Entrez)
#[1] "Gene"     "EntrezID"

### enrichDO
library(DOSE)

x <- enrichDO(gene = Id_Entrez$EntrezID,
              ont = "DO",
              pvalueCutoff = 0.05,
              pAdjustMethod = "BH",
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.05,
              readable      = FALSE)

head(x)

### enrichNCG : Network of Cancer Gene (NCG)(A. et al. 2016)
ncg <- enrichNCG(Id_Entrez$EntrezID)
head(ncg)

### enrich DGN ; DisGeNET(Janet et al. 2015) gene-disease associations
dgn <- enrichDGN(Id_Entrez$EntrezID)
head(dgn)
unique(dgn$Description)

Desc_Pav = data.frame(dgn$Description, dgn$p.adjust)

#### sory by adjust P-value
Desc_Pav = Desc_Pav[order(Desc_Pav$dgn.p.adjust),]
head(Desc_Pav,10)

Order = grep("Breast",Desc_Pav$dgn.Description)
str(Order)
Disease = Desc_Pav$dgn.Description[Order]
P_value = Desc_Pav$dgn.p.adjust[Order]

dgn_Breast = data.frame(Order, Disease)
dgn_Breast["P_adj"]=P_value

head(dgn_Breast)


### gseGO
# https://yulab-smu.github.io/clusterProfiler-book/chapter5.html
head(results)
rm(geneList)
geneList <- results[,2]
names(geneList) <- as.character(rownames(results))
geneList <- geneList[res_entrez$Gene]
names(geneList) <- as.character(res_entrez$EntrezID)
geneList <- sort(geneList, decreasing = TRUE)


ego4 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "CC",
              nPerm        = 1000,
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)




geneList

kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)

filekg=sprintf("./02Breast/023report/2.GTExNormal/Function/keggresult%s.csv",Sys.Date())
write.csv(kk2@result,file = filekg)

library("pathview")
hsa04110 <- pathview(gene.data  = geneList,
                     pathway.id = "hsa05034",
                     species    = "hsa",
                     limit      = list(gene=max(abs(geneList)), cpd=1))


