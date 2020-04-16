
load(file = "200410dim12_b1b2t.RData")



### top10 heatmap

top10 <- smc024.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
filet10h=sprintf("12top10_heatmap_%s_%s.png", project_name, Sys.Date())
png(filename = filet10h, width = 1500, height = 1000)
DoHeatmap(smc024.merged, features = top10$gene) + NoLegend()
dev.off()



## top 10 heatmap
filet10h=sprintf("12top10_heatmap_%s_%s.png", project_name, Sys.Date())
png(filename = filet10h, width = 1500, height = 960)
ph1<-DoHeatmap(smc024.merged, features = top10$gene) + NoLegend()
dev.off()


## cluster dimplot
diN = 11
fileum=sprintf("09Dim_umap_%s_dim%s_%s.png", project_name, diN,Sys.Date())
png(filename = fileum, width = 1160, height = 400)
p1 <- DimPlot(smc024.merged, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(smc024.merged, reduction = "umap", label = TRUE)
p1+p2
dev.off()

fileums=sprintf("09Dim_split_umap_%s_dim%s_%s.png", project_name, diN,Sys.Date())
png(filename = fileums, width = 1160, height = 400)
DimPlot(smc024.merged, reduction = "umap", split.by = "orig.ident")
dev.off()



### nFeature RNA, nCount RNA
filef1=sprintf("11_1feat_nFRNA_%s_%s.png", project_name, Sys.Date())
png(filename = filef1, width = 560, height = 500)
FeaturePlot(smc024.merged, features=c("nFeature_RNA"))
dev.off()

filef2=sprintf("11_2feat_nCRNA_%s_%s.png", project_name, Sys.Date())
png(filename = filef2, width = 560, height = 500)
FeaturePlot(smc024.merged, features=c("nCount_RNA"))
dev.off()


filef1s=sprintf("11_1feat_split_nFRNA_%s_%s.png", project_name, Sys.Date())
png(filename = filef1s, width = 1660, height = 600)
FeaturePlot(smc024.merged, features=c("nFeature_RNA"), split.by = "orig.ident")
dev.off()


filef2s=sprintf("11_2feat_split_nCRNA_%s_%s.png", project_name, Sys.Date())
png(filename = filef2s, width = 1660, height = 600)
FeaturePlot(smc024.merged, features=c("nCount_RNA"), split.by = "orig.ident")
dev.off()

filef1=sprintf("11_1vln_nFRNA_%s_%s.png", project_name, Sys.Date())
png(filename = filef1, width = 660, height = 400)
VlnPlot(smc024.merged, features=c("nFeature_RNA"))
dev.off()

filef2=sprintf("11_2vln_nCRNA_%s_%s.png", project_name, Sys.Date())
png(filename = filef2, width = 660, height = 400)
VlnPlot(smc024.merged, features=c("nCount_RNA"))
dev.off()

## marker plotting
VlnPlot(smc024.merged, features=c("CEACAM5", "CEACAM6"))
VlnPlot(smc024.merged, features = c("CEACAM5", "CEACAM6"), slot = "counts", log = TRUE)
FeaturePlot(smc024.merged, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))


#FeaturePlot(smc024.merged, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))

fileseu=sprintf("13pbmc_feat_%s_%s.png", project_name, Sys.Date())
png(filename = fileseu, width = 1000, height = 600)
FeaturePlot(smc024.merged, features = c("MS4A1","GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))
dev.off()

## marker plotting
#"VIMENTIN"%in%rownames(smc024.merged[["RNA"]]) : FALSE

#CTC : EPCAM, CD44, KRT8, KRT19, Vimentin, MET, KRAS, NRAS
filectc1=sprintf("13ctc_vln_%s_%s.png", project_name, Sys.Date())
png(filename = filectc1, width = 1000, height = 600)
VlnPlot(smc024.merged, features=c("EPCAM","CD44","KRT8","KRT19","VIM","MET","KRAS","NRAS"), split.by = "orig.ident", pt.size = 0)
dev.off()


#VlnPlot(smc024.merged, features = c("EPCAM","CD44","KRT8","KRT19","VIM","MET","KRAS","NRAS"), slot = "counts", log = TRUE)
filectc2=sprintf("13ctc_feat_%s_%s.png", project_name, Sys.Date())
png(filename = filectc2, width = 800, height = 400)
FeaturePlot(smc024.merged, features = c("EPCAM","CD44","KRT8","KRT19"), split.by = "orig.ident", max.cutoff = 3)
dev.off()

filectc3=sprintf("13ctc_feat2_%s_%s.png", project_name, Sys.Date())
png(filename = filectc3, width = 800, height = 400)
FeaturePlot(smc024.merged, features = c("VIM","MET","KRAS","NRAS"), split.by = "orig.ident", max.cutoff = 3)
dev.off()


DoHeatmap(smc024.merged, features = c("EPCAM","CD44","KRT8","KRT19","VIM","MET","KRAS","NRAS"), slot="data")
#DoHeatmap(smc024.merged, features = c("EPCAM","CD44","KRT8","KRT19","VIM","MET","KRAS","NRAS"))


# CTC-SCLC neuroendocrine markers ## UCHL1, NCAM1, SYP, CHGA
# EMT signature ## ASCL1, NEUROD1, POU2F3
#"UCHL1", "NCAM1", "SYP", "CHGA", "ASCL1", "NEUROD1", "POU2F3"

filectc5=sprintf("13ctc5_vln_%s_%s.png", project_name, Sys.Date())
png(filename = filectc5, width = 1000, height = 600)
VlnPlot(smc024.merged, features=c("UCHL1", "NCAM1", "SYP", "CHGA", "ASCL1", "NEUROD1", "POU2F3"), split.by = "orig.ident", pt.size = 0)
dev.off()

#VlnPlot(smc024.merged, features = c("EPCAM","CD44","KRT8","KRT19","VIM","MET","KRAS","NRAS"), slot = "counts", log = TRUE)
filectc6=sprintf("13ctc6_feat_%s_%s.png", project_name, Sys.Date())
png(filename = filectc6, width = 800, height = 400)
FeaturePlot(smc024.merged, features = c("UCHL1", "NCAM1","SYP", "CHGA"), split.by = "orig.ident", max.cutoff = 3)
dev.off()


filectc7=sprintf("13ctc_feat7_%s_%s.png", project_name, Sys.Date())
png(filename = filectc7, width = 800, height = 400)
FeaturePlot(smc024.merged, features = c("ASCL1", "NEUROD1", "POU2F3","POU2F3"), split.by = "orig.ident", max.cutoff = 3)
dev.off()


# filectc8=sprintf("13ctc_feat8_%s_%s.png", project_name, Sys.Date())
# png(filename = filectc8, width = 600, height = 600)
# FeaturePlot(smc024.merged, features = c(), split.by = "orig.ident", max.cutoff = 3)
# dev.off()


#NSCLC : CEA, CYFRA 21-1, EGFR, HER2, BRAF, p53, VEGF, PI3k, mTOR, RAS, MEK, c-KIT
# VlnPlot(smc024.merged, features=c("CEACAM5", "CEACAM6", "KRT19", "EGFR", "ERBB2", "BRAF", "TP53", "VEGFA", "PIK3CA", "MTOR", "KRAS", "NRAS", "HRAS", "MAP2K1", "KIT"), y.max=4, same.y.lims=TRUE)
# FeaturePlot(smc024.merged, features = c("CEACAM5", "CEACAM6", "KRT19", "EGFR", "ERBB2", "BRAF", "TP53", "VEGFA", "PIK3CA", "MTOR", "KRAS", "NRAS", "HRAS", "MAP2K1", "KIT"))
# DoHeatmap(smc024.merged, features = c("CEACAM5", "CEACAM6", "KRT19", "EGFR", "ERBB2", "BRAF", "TP53", "VEGFA", "PIK3CA", "MTOR", "KRAS", "NRAS", "HRAS", "MAP2K1", "KIT"), slot="data")

filecan1=sprintf("13cancer_vln_%s_%s.png", project_name, Sys.Date())
png(filename = filecan1, width = 1100, height = 600)
VlnPlot(smc024.merged, features=c("ERBB2", "BRAF", "TP53", "VEGFA", "PIK3CA", "MTOR", "KRAS", "NRAS", "HRAS", "MAP2K1", "KIT"), split.by = "orig.ident", pt.size = 0)
dev.off()

filecan2=sprintf("13cancer_feat1_%s_%s.png", project_name, Sys.Date())
png(filename = filecan2, width = 1000, height = 600)
FeaturePlot(smc024.merged, features = c("ERBB2", "BRAF", "TP53"), split.by = "orig.ident", max.cutoff = 3)
dev.off()


filecan3=sprintf("13cancer_feat2_%s_%s.png", project_name, Sys.Date())
png(filename = filecan3, width = 1000, height = 600)
FeaturePlot(smc024.merged, features = c("VEGFA", "PIK3CA", "MTOR"), split.by = "orig.ident", max.cutoff = 3)
dev.off()


filecan4=sprintf("13cancer_feat3_%s_%s.png", project_name, Sys.Date())
png(filename = filecan4, width = 1000, height = 600)
FeaturePlot(smc024.merged, features = c("KRAS", "NRAS", "HRAS"), split.by = "orig.ident", max.cutoff = 3)
dev.off()


filecan5=sprintf("13cancer_feat4_%s_%s.png", project_name, Sys.Date())
png(filename = filecan5, width = 1000, height = 600)
FeaturePlot(smc024.merged, features = c("MAP2K1", "KIT", "ERBB2"), split.by = "orig.ident", max.cutoff = 3)
dev.off()


#T cell : CD3D, CD3E, CD3G, IGFLR1, RGS1, TCF7, NKIRAS2, IL7R, CD8A, KLRG1, GZMB, GZMH, PRF1
filetc1=sprintf("13tcell_vln_%s_%s.png", project_name, Sys.Date())
png(filename = filetc1, width = 1000, height = 600)
VlnPlot(smc024.merged, features=c("CD3D", "CD3E", "CD3G", "IGFLR1", "RGS1", "TCF7", "NKIRAS2", "IL7R", "CD8A", "KLRG1", "GZMB", "GZMH", "PRF1"), split.by = "orig.ident", pt.size = 0)
dev.off()

filetc2=sprintf("13tcell_feat_%s_%s.png", project_name, Sys.Date())
png(filename = filetc2, width = 1000, height = 600)
FeaturePlot(smc024.merged, features = c("CD3D", "CD3E", "CD3G"), split.by = "orig.ident", max.cutoff = 3)
dev.off()

filetc3=sprintf("13tcell_feat1_%s_%s.png", project_name, Sys.Date())
png(filename = filetc3, width = 1000, height = 600)
FeaturePlot(smc024.merged, features = c("IGFLR1", "RGS1", "TCF7"), split.by = "orig.ident", max.cutoff = 3)
dev.off()

filetc4=sprintf("13tcell_feat2_%s_%s.png", project_name, Sys.Date())
png(filename = filetc4, width = 1000, height = 600)
FeaturePlot(smc024.merged, features = c("NKIRAS2", "IL7R", "CD8A"), split.by = "orig.ident", max.cutoff = 3)
dev.off()

filetc5=sprintf("13tcell_feat3_%s_%s.png", project_name, Sys.Date())
png(filename = filetc5, width = 1000, height = 600)
FeaturePlot(smc024.merged, features = c("KLRG1", "GZMB", "GZMH"), split.by = "orig.ident", max.cutoff = 3)
dev.off()

filetc6=sprintf("13tcell_feat4_%s_%s.png", project_name, Sys.Date())
png(filename = filetc6, width = 1000, height = 600)
FeaturePlot(smc024.merged, features = c("GZMB", "GZMH", "PRF1"), split.by = "orig.ident", max.cutoff = 3)
dev.off()


#B cell : CD79A, CD63, IGHA1, IGHA2, IGHD, VPREB3
filebc1=sprintf("13bcell_vln_%s_%s.png", project_name, Sys.Date())
png(filename = filebc1, width = 1000, height = 600)
VlnPlot(smc024.merged, features=c("CD79A", "CD63", "IGHA1", "IGHA2", "IGHD", "VPREB3"), split.by = "orig.ident", pt.size = 0)
dev.off()

filebc2=sprintf("13bcell_feat1_%s_%s.png", project_name, Sys.Date())
png(filename = filebc2, width = 1000, height = 600)
FeaturePlot(smc024.merged, features = c("CD79A", "CD63", "VPREB3"), split.by = "orig.ident", max.cutoff = 3)
dev.off()

filebc3=sprintf("13bcell_feat2_%s_%s.png", project_name, Sys.Date())
png(filename = filebc3, width = 1000, height = 600)
FeaturePlot(smc024.merged, features = c("IGHA1", "IGHA2", "IGHD"), split.by = "orig.ident", max.cutoff = 3)
dev.off()

# CTC
filectca1=sprintf("13ctca1_vln_%s_%s.png", project_name, Sys.Date())
png(filename = filectca1, width = 800, height = 400)
VlnPlot(smc024.merged, features=c("AXL", "ERBB2", "ERBB3", "MET", "CSF1R", "KRT19"), split.by = "orig.ident", pt.size = 0)
dev.off()

filectca2=sprintf("13ctca2_feat1_%s_%s.png", project_name, Sys.Date())
png(filename = filectca2, width = 800, height = 400)
FeaturePlot(smc024.merged, features = c("AXL", "ERBB2", "ERBB3", "MET"), split.by = "orig.ident", max.cutoff = 3)
dev.off()

filectca3=sprintf("13ctca3_feat2_%s_%s.png", project_name, Sys.Date())
png(filename = filectca3, width = 400, height = 400)
FeaturePlot(smc024.merged, features = c("CSF1R", "KRT19"), split.by = "orig.ident", max.cutoff = 3)
dev.off()




