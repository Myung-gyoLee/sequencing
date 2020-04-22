analysis_type = "gs7"
ego3 <- gs7
if (analysis_type %in% c("edo", "ncg", "dgn", "gkk", "ekk")){
  edox <- setReadable(ego3, 'org.Hs.eg.db', 'ENTREZID') 
}else {
  edox <- ego3
}

# ego
# em2
# gs2
# em5
# gs5
# em6
# gs6
# em7
# gs7
# edo
# ncg
# dgn
# gkk
# ekk

fileeg=sprintf("%s_Dotplot%s.png", analysis_type, Sys.Date())
png(filename = fileeg, height=800, width=1200, bg="white")
dotplot(ego3, showCategory=30)
dev.off()

fileegh1=sprintf("%s_heatplot%s.png",analysis_type,Sys.Date())
png(filename = fileegh1, height=400, width=1400, bg="white")
heatplot(edox, showCategory=30, foldChange = geneList)
dev.off()


#categorySize="pvalue", showCategory = 5, foldChange=OE_foldchanges,vertex.label.font=6
# PDF size to 24 x 32

fileeg2=sprintf("%s_cnetplot%s.png", analysis_type, Sys.Date())
png(filename = fileeg2, height=500, width=500, bg="white")
cnetplot(edox, showCategory = 5)
dev.off()


## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
# PDF size to 24 x 32
fileeg3=sprintf("%s_emapplot%s.png", analysis_type, Sys.Date())
png(filename = fileeg3, height=800, width=800, bg="white")
emapplot(ego3, showCategory = 50)
dev.off()

