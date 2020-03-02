##----------------------------------------------------------------------------------##
## Plot Start
##----------------------------------------------------------------------------------##

#### set analysis type ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### GO Analysis
#### enrich GO : ego 
#### gsea GO : gsego

### msigDB Analysis 
#### C5 GO: emC5, gsC5 | C6 oncogenic :emC6, gsC6 | C7 immune : emC7, gsC7

### enrichDO Analysis #### edo1
### enrichGO Analysis #### ego3
### ncg Analysis #### ncg
### dgn Analysis #### dgn

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
getwd()
#setwd("../MSigDb/")
analysis_type = "enrichGO"
analysis1 <- ego3




### emapplot
fileeg3=sprintf("%s_emapplot%s.png", analysis_type, Sys.Date())
#png(filename = fileeg3, height=900, width=600, bg="white")
png(filename = fileeg3, height=500, width=800, bg="white")
emapplot(analysis1, showCategory = 30)
dev.off()

### plotGOgraph
fpdf = sprintf("%s_plotGograph%s.pdf", analysis_type, Sys.Date())
pdf(fpdf, paper = "a4")
plotGOgraph(analysis1, ontology())
dev.off()


### heatplot
rm(edox)
edox <- analysis1

## if EntrezID ### ego doesn't need
edox <- setReadable(analysis1, 'org.Hs.eg.db', 'ENTREZID')
head(edox)

filers=sprintf("%s_result%s.csv", analysis_type, Sys.Date())
write.csv(edox@result,file = filers)

#pheight <- length(analysis1@result$ID)*100

filehp=sprintf("%s_heatplot%s.png", analysis_type, Sys.Date())
png(filename = filehp, height=300, width=950, bg="white")
heatplot(edox, foldChange=geneList)
dev.off()


### dotplot
fileeg1=sprintf("%s_Dotplot%s.png", analysis_type, Sys.Date())
png(filename = fileeg1, height=250, width=650, bg="white")

dotplot(analysis1, showCategory=30, decreasing=FALSE)
dev.off()


### cnetplot
fileeg2=sprintf("%s_cnetplot%s.png", analysis_type, Sys.Date())
png(filename = fileeg2, height=500, width=400, bg="white")
cnetplot(edox, showCategory = 5)
dev.off()




