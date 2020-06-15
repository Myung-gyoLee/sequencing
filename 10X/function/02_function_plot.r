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

#axis.text.x = element_text(size = 15 ... : x axis label size
#axis.text.y = element_text(size = 15 ... : y axis label size

##================= egobp =====================##
analysis_type1 = "egobp"
ego3 <- egobp
if (analysis_type1 %in% c("edo", "ncg", "dgn", "gkk", "ekk")){
  edox <- setReadable(ego3, 'org.Hs.eg.db', 'ENTREZID') 
}else {
  edox <- ego3
}

# ego
# egobp
# egomf
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

file1eg=sprintf("%s_Dotplot%s.png", analysis_type1, Sys.Date())
png(filename = file1eg, height=500, width=600, bg="white")
dotplot(ego3, showCategory=30) +theme(panel.grid.major = element_blank(),
                                      axis.text.y = element_text(size = 14, angle = 0))
dev.off()

file1egh1=sprintf("%s_heatplot%s.png",analysis_type1,Sys.Date())
png(filename = file1egh1, height=450, width=1200, bg="white")
heatplot(edox, showCategory=30, foldChange = geneList)  +
theme(panel.grid.major = element_blank(),
      axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1), 
      #axis.text.x = element_text(size = 15, angle = 60, hjust = 1), 
      axis.text.y = element_text(size = 14, angle = 0))
#heatplot.enrichResult(edox, showCategory=30, foldChange = geneList)
dev.off()


#categorySize="pvalue", showCategory = 5, foldChange=OE_foldchanges,vertex.label.font=6
# PDF size to 24 x 32

file1eg2=sprintf("%s_cnetplot%s.png", analysis_type1, Sys.Date())
png(filename = file1eg2, height=300, width=500, bg="white")
cnetplot(edox, showCategory = 5)
dev.off()


## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
# PDF size to 24 x 32
file1eg3=sprintf("%s_emapplot%s.png", analysis_type1, Sys.Date())
png(filename = file1eg3, height=400, width=400, bg="white")
emapplot(ego3, showCategory = 20)
dev.off()

##================= ekk =====================##

analysis_type2 = "ekk"
ego1 <- ekk
if (analysis_type2 %in% c("edo", "ncg", "dgn", "gkk", "ekk")){
  edox1 <- setReadable(ego1, 'org.Hs.eg.db', 'ENTREZID') 
}else {
  edox1 <- ego1
}

# ego
# egobp
# egomf
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

file2eg=sprintf("%s_Dotplot%s.png", analysis_type2, Sys.Date())
png(filename = file2eg, height=300, width=500, bg="white")
dotplot(ego1, showCategory=30) +theme(panel.grid.major = element_blank(),
                                      axis.text.y = element_text(size = 14, angle = 0))
dev.off()

file2egh1=sprintf("%s_heatplot%s.png",analysis_type2,Sys.Date())
png(filename = file2egh1, height=250, width=750, bg="white")
heatplot(edox1, showCategory=30, foldChange = geneList) +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1), 
        #axis.text.x = element_text(size = 15, angle = 60, hjust = 1), 
        axis.text.y = element_text(size = 14, angle = 0))
dev.off()


#categorySize="pvalue", showCategory = 5, foldChange=OE_foldchanges,vertex.label.font=6
# PDF size to 24 x 32

file2eg2=sprintf("%s_cnetplot%s.png", analysis_type2, Sys.Date())
png(filename = file2eg2, height=300, width=500, bg="white")
cnetplot(edox1, showCategory = 5)
dev.off()


## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
# PDF size to 24 x 32
file2eg3=sprintf("%s_emapplot%s.png", analysis_type2, Sys.Date())
png(filename = file2eg3, height=400, width=400, bg="white")
emapplot(ego1, showCategory = 20)
dev.off()


