library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
setwd("H:TCGA/SingleCell/")
list.files()


library(dplyr)
library(vroom)
cell_markers <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt') %>%
  tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>% 
  dplyr::select(cellMarker, geneID) %>%
  dplyr::mutate(geneID = strsplit(geneID, ', '))
cell_markers



## read data
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

### Analysis Start

#clt_name = "cluster5"
#dir.create("./Function/%s",clt_name)
setwd("H:/TCGA/SingleCell/Function/cluster8")
getwd()

rm(df)
rm(geneList)
rm(gene)
df <-cl8
geneList <- df[,3]
head(geneList)

names(geneList) <- as.character(df[,1])

geneList <- sort(geneList, decreasing = TRUE)
gene <-names(geneList)
gene <- names(geneList)[abs(geneList) >= 0.58]

gene


## 
rm(ego3)
ego3 <- enrichGO(gene         = names(geneList),
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)


#rm(plotdir)
#plotdir = "./Function/cluster4"
fileeg1=sprintf("enrichGo_Dotplot%s.png",Sys.Date())
png(filename = fileeg1, height=800, width=1200, bg="white")

dotplot(ego3, showCategory=30)
dev.off()

pdf("ego_dotplot.pdf", paper = "a4r")
dotplot(ego3, showCategory=30)
dev.off()



fileegh1=sprintf("enrichGO_heatplot%s.png",Sys.Date())
png(filename = fileegh1, height=400, width=800, bg="white")

heatplot(ego3, showCategory=30, foldChange = geneList)
dev.off()




#categorySize="pvalue", showCategory = 5, foldChange=OE_foldchanges,vertex.label.font=6
# PDF size to 24 x 32

fileeg2=sprintf("cnetplot%s.png", Sys.Date())
png(filename = fileeg2, height=500, width=500, bg="white")


cnetplot(ego3, showCategory = 5)
dev.off()



## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
# PDF size to 24 x 32
fileeg3=sprintf("emapplot%s.png", Sys.Date())
png(filename = fileeg3, height=800, width=600, bg="white")
emapplot(ego3, showCategory = 50)
dev.off()

pdf("emapplot.pdf", paper = "a4")
emapplot(ego3, showCategory = 50)
dev.off()

pdf("plotGograph1.pdf", paper = "a4")

plotGOgraph(ego3)
dev.off()




### MSigDb analysis
#BiocManager::install("msigdbr")
library(msigdbr)
msigdbr_show_species()

m_df <- msigdbr(species = "Homo sapiens")
head(m_df, 2) %>% as.data.frame

# C6: Oncogenic
#colnames(msigdbr(species = "Homo sapiens", category = "C6"))
m_t2g <- msigdbr(species = "Homo sapiens", category = "C6") %>% 
  dplyr::select(gs_name, gene_symbol)
head(m_t2g)


em1 <- enricher(gene, TERM2GENE = m_t2g)
head(em1)


em2 <- GSEA(geneList, TERM2GENE = m_t2g)
head(em2)

#write.csv(em1@result,file = "enrichrresult.csv")
write.csv(em2@result,file = "GSEA_oncogenic_result.csv")


# C7 immune
rm(m_t2g)
m_t2g <- msigdbr(species = "Homo sapiens", category = "C7") %>% 
  dplyr::select(gs_name, gene_symbol)
head(m_t2g)


em1 <- enricher(gene, TERM2GENE = m_t2g)
head(em1)


em2 <- GSEA(geneList, TERM2GENE = m_t2g)
head(em2)

#write.csv(em1@result,file = "enrichrresult.csv")
write.csv(em2@result,file = "GSEA_immune_result.csv")

## DO analysis
# https://yulab-smu.github.io/clusterProfiler-book/chapter4.html#enrichdo-function
# https://www.rdocumentation.org/packages/DOSE/versions/2.10.6/topics/enrichDO

### make Input 
# https://bioconductor.org/packages/release/data/annotation/manuals/org.Hs.eg.db/man/org.Hs.eg.db.pdf
library(org.Hs.eg.db)

xx <- as.list(org.Hs.egALIAS2EG)
head(xx)

rm(I_Entrez)
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
#Id_Entrez <- sort(Id_entrez, decreasing = TRUE)

head(res_entrez)
nrow(res_entrez)
nrow(Id_Entrez)
#colnames(Id_Entrez)
#[1] "Gene"     "EntrezID"

### enrichDO
library(DOSE)
rm(x)
x <- enrichDO(gene = Id_Entrez$EntrezID,
              ont = "DO",
              pvalueCutoff = 0.05,
              pAdjustMethod = "BH",
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.05,
              readable      = FALSE)
#x <- x[order(-x$Count),]
#str(x)
head(x)
fileDO=sprintf("DO_barplot%s.png", Sys.Date())
png(filename = fileDO, height=600, width=800, bg="white")
barplot(x, showCategory = 20)
dev.off()

### DO heatplot
rm(edox)
edox <- setReadable(x, 'org.Hs.eg.db', 'ENTREZID')
#edox <- edox[order(-edox$Count),]
head(edox)

fileheat=sprintf("DO_heatplot%s.png", Sys.Date())
png(filename = fileheat, height=350, width=600, bg="white")
heatplot(edox, foldChange=geneList)
dev.off()

### heat-bar DO plot
filehc=sprintf("DO_heat_bar%s.png", Sys.Date())
png(filename = filehc, height=600, width=1400, bg="white")
p1 <- barplot(x, showCategory = 20)
p2 <- heatplot(edox, showCategory = 20, foldChange = geneList)
cowplot::plot_grid(p1, p2, nrow = 1, labels=LETTERS[1:2], rel_widths = c(1, 2))
dev.off()



write.csv(x@result,file = "DOresult.csv")
### enrichNCG : Network of Cancer Gene (NCG)(A. et al. 2016)
ncg <- enrichNCG(Id_Entrez$EntrezID)
head(ncg)
write.csv(ncg@result,file = "NCG_result.csv")

filencg=sprintf("ncg_barplot%s.png", Sys.Date())
png(filename = filencg, height=350, width=600, bg="white")
barplot(ncg, foldChange=geneList)
dev.off()


### enrich DGN ; DisGeNET(Janet et al. 2015) gene-disease associations
dgn <- enrichDGN(Id_Entrez$EntrezID)
head(dgn)
write.csv(dgn@result,file = "dgnresult.csv")

### DGN barplot
library(enrichplot)
filedgn=sprintf("dgn_barplot%s.png", Sys.Date())
png(filename = filedgn, height=600, width=800, bg="white")
barplot(dgn, showCategory = 20)
dev.off()



unique(dgn$Description)

Desc_Pav = data.frame(dgn$Description, dgn$p.adjust)


#### sory by adjust P-value
Desc_Pav = Desc_Pav[order(Desc_Pav$dgn.p.adjust),]
head(Desc_Pav,10)

#Order = grep("Breast",Desc_Pav$dgn.Description)
#str(Order)
#Disease = Desc_Pav$dgn.Description[Order]
#P_value = Desc_Pav$dgn.p.adjust[Order]

#dgn_Breast = data.frame(Order, Disease)
#dgn_Breast["P_adj"]=P_value

#head(dgn_Breast)

### dgn heatplot
## convert gene ID to Symbol
rm(edox)
edox <- setReadable(dgn, 'org.Hs.eg.db', 'ENTREZID')

fileheat=sprintf("dgn_heatplot%s.png", Sys.Date())
png(filename = fileheat, height=900, width=800, bg="white")
heatplot(edox, foldChange = geneList)
dev.off()

### heat-bar DGN plot
filehcg=sprintf("DGN_heat_bar%s.png", Sys.Date())
png(filename = filehcg, height=400, width=1500, bg="white")
p1 <- barplot(dgn, showCategory = 20)
p2 <- heatplot(edox, showCategory = 20, foldChange = geneList)
cowplot::plot_grid(p1, p2, nrow = 1, labels=LETTERS[1:2], rel_widths = c(1, 3))
dev.off()




### gseGO
# https://yulab-smu.github.io/clusterProfiler-book/chapter5.html
head(results)
rm(geneList)

geneList <- geneList[res_entrez$Gene]
names(geneList) <- as.character(res_entrez$EntrezID)
geneList <- sort(geneList, decreasing = TRUE)

geneList

kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)

filekg=sprintf("keggresult%s.csv",Sys.Date())
write.csv(kk2@result,file = filekg)

library("pathview")
hsa04110 <- pathview(gene.data  = geneList,
                     pathway.id = "hsa05034",
                     species    = "hsa",
                     limit      = list(gene=max(abs(geneList)), cpd=1))

save(list=ls(), file='SMC009.Rdata')

### Gene Set Enrichment Analysis

gsecc <- gseGO(geneList=geneList, ont="CC", OrgDb=org.Hs.eg.db, verbose=F)
head(summary(gsecc))
gseaplot(gsecc, geneSetID="GO:0031012")

## Cell Marker
#BiocManager::install("vroom")

head(res_entrez)
y <- enricher(res_entrez$EntrezID, TERM2GENE=cell_markers, minGSSize=1)
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
