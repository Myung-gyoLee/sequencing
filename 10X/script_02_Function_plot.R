##----------------------------------------------------------------------------------##
## Plot Start
##----------------------------------------------------------------------------------##

#### set analysis type ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### GO Analysis
#### enrichGO : ego3 
#### gseaGO : gsego

### msigDB Analysis 
#### C2 curated: emC2, gsC2 |  C5 GO: emC5, gsC5 | 
#### C6 oncogenic :emC6, gsC6 | C7 immune : emC7, gsC7 |

### Wiki Pathways #### ewp.up
### kegg #### kk1

### enrichDO Analysis #### edo1
### ncg Analysis #### ncg
### dgn Analysis #### dgn

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

getwd()
dir.create('../../kegg')
setwd('../../kegg')

dir.create('../kegg')
setwd('../kegg')

rm(analysis1)
rm(analysis_type)
rm(edox)

analysis_type = "Wiki"
analysis1 <- ewp.up



### emapplot
fileeg3=sprintf("%s_emapplot%s.png", analysis_type, Sys.Date())
#png(filename = fileeg3, height=900, width=900, bg="white")
#png(filename = fileeg3, height=400, width=550, bg="white")
png(filename = fileeg3, height=500, width=350, bg="white")
emapplot(analysis1, showCategory = 30)
dev.off()

### plotGOgraph
# fpdf = sprintf("%s_plotGograph%s.pdf", analysis_type, Sys.Date())
# pdf(fpdf, paper = "a4")
# plotGOgraph(analysis1, ontology('CC'))
# dev.off()

head(analysis1)
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
png(filename = filehp, height=200, width=550, bg="white")
heatplot(edox, foldChange=geneList)
dev.off()
# 400

### dotplot
fileeg1=sprintf("%s_Dotplot%s.png", analysis_type, Sys.Date())
png(filename = fileeg1, height=300, width=550, bg="white")
dotplot(analysis1, showCategory=30, decreasing=FALSE)
dev.off()

# 400, 550

### cnetplot
fileeg2=sprintf("%s_cnetplot%s.png", analysis_type, Sys.Date())
png(filename = fileeg2, height=400, width=400, bg="white")
#png(filename = fileeg2, height=400, width=550, bg="white")
#png(filename = fileeg2, height=600, width=600, bg="white")
cnetplot(edox, showCategory = 5)
dev.off()




