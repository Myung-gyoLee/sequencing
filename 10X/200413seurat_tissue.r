#!/usr/bin/env Rscript

library(Seurat)
library(dplyr)
library(ggplot2)
library(sctransform)
library(patchwork)

setwd('H:SingleCell/smc024')

setwd('./analysis_tissue')

# Load the dataset

SMC_024_blood1_data <- Read10X(data.dir = "smc024_blood1/filtered_feature_bc_matrix")

SMC_024_blood2_data <- Read10X(data.dir = "smc024_blood2/filtered_feature_bc_matrix")

SMC_024_tissue_data <- Read10X(data.dir = "smc024_tissue/filtered_feature_bc_matrix")


# SMC_024_blood1_data <- Read10X(data.dir = "/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/10X/HN00124804/HN00124804_result_10X/SMC_024_blood1/filtered_feature_bc_matrix")
# 
# SMC_024_blood2_data <- Read10X(data.dir = "/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/10X/HN00124804/HN00124804_result_10X/SMC_024_blood2/filtered_feature_bc_matrix")
# 
# SMC_024_tissue_data <- Read10X(data.dir = "/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/10X/HN00124804/HN00124804_result_10X/SMC_024_tissue/filtered_feature_bc_matrix")

#dir.create("./analysis_tissue")
setwd("./analysis_tissue")



# create seurat object

SMC_024_blood1 <- CreateSeuratObject(counts = SMC_024_blood1_data, project = "SMC_024_blood1")
SMC_024_blood2 <- CreateSeuratObject(counts = SMC_024_blood2_data, project = "SMC_024_blood2")
SMC_024_tissue <- CreateSeuratObject(counts = SMC_024_tissue_data, project = "SMC_024_tissue")


# https://satijalab.org/seurat/v3.1/merge_vignette.html
smc024.merged <- merge(SMC_024_blood1, y = c(SMC_024_blood2, SMC_024_tissue), add.cell.ids = c("SMC_024_blood1", "SMC_024_blood2", "SMC_024_tissue") , project = "SMC024")

saveRDS(smc024.merged, file = "200403merged.RDS")

smc024.merged[[]]$orig.ident %>% unique

# ==================================================================================================== #

# store mitochondrial percentage in object meta data
smc024.merged <- PercentageFeatureSet(smc024.merged, pattern = "^MT-", col.name = "percent.mt")
# Visualize QC metrics as a violin plot
project_name="SMC024_tissue"


filevp=sprintf("01QCvlnplot_%s_%s.png", project_name, Sys.Date())
png(filename = filevp, width = 660, height = 500)
VlnPlot(smc024.merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
#head(tenx009@meta.data, 5)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
filefc=sprintf("02feat_scatt_%s_%s.png", project_name, Sys.Date())
png(filename = filefc, width = 1020, height = 330)
plot1 <- FeatureScatter(smc024.merged, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(smc024.merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()


## filtering
smc024.merged <- subset(smc024.merged, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 30)

# # Normalizing the data
# smc024.merged <- NormalizeData(smc024.merged, normalization.method = "LogNormalize", scale.factor = 10000)
# 
# smc024.merged <- FindVariableFeatures(smc024.merged, selection.method = "vst", nfeatures = 2000)
# 
# top10 <- head(VariableFeatures(smc024.merged), 10)
# 
# filevfv=sprintf("02volcano_%s_%s.png", project_name, Sys.Date())
# png(filename = filevfv, width = 1020, height = 330)
# plot1 <- VariableFeaturePlot(smc024.merged)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# plot1+plot2
# dev.off()
# 
# all.genes <- rownames(smc024.merged)
# smc024.merged <- ScaleData(smc024.merged, features = all.genes)
# 
# smc024.merged <- RunPCA(smc024.merged, features = VariableFeatures(object = smc024.merged))
# 
# print(smc024.merged[["pca"]], dims = 1:5, nfeatures = 5)
# 

# run sctransform
smc024.merged <- SCTransform(smc024.merged, vars.to.regress = "percent.mt", conserve.memory = TRUE, verbose = FALSE)
save(smc024.merged,file='smc024.merged_sct.Rdata')

#tenx009 <- SCTransform(tenx009, vars.to.regress = "percent.mt", conserve.memory = TRUE, verbose = FALSE, return.only.var.genes = FALSE)
#pbmc[["SCT"]]@scale.data contains the residuals (normalized values), and is used directly as input to PCA. Please note that this matrix is non-sparse, and can therefore take up a lot of memory if stored for all genes. To save memory, we store these values only for variable genes, by setting the return.only.var.genes = TRUE by default in the SCTransform function call.
#The ?€˜corrected?€? UMI counts are stored in pbmc[["SCT"]]@counts. We store log-normalized versions of these corrected counts in pbmc[["SCT"]]@data, which are very helpful for visualization.
#You can use the corrected log-normalized counts for differential expression and integration. However, in principle, it would be most optimal to perform these calculations directly on the residuals (stored in the scale.data slot) themselves. This is not currently supported in Seurat v3, but will be soon.


smc024.merged <- RunPCA(smc024.merged, verbose = FALSE)

#load(file='smc024.merged_pca_umap.Rdata')
#table(smc024.merged, smc024.merged@meta.data$orig.ident)
filevd=sprintf("03vizdim_%s_%s.png", project_name, Sys.Date())
png(filename = filevd, width = 670, height = 600)
VizDimLoadings(smc024.merged, dims = 1:2, reduction = "pca")
dev.off()

filedp=sprintf("04dimpc_%s_%s.png", project_name, Sys.Date())
png(filename = filedp, width = 670, height = 600)
DimPlot(smc024.merged, reduction = "pca")
dev.off()

filedh1=sprintf("05dimheat_%s_%s.png", project_name, Sys.Date())
png(filename = filedh1, width = 670, height = 600)
DimHeatmap(smc024.merged, dims = 1, cells = 500, balanced = TRUE)
dev.off()

filedh2=sprintf("06dimheat15_%s_%s.png", project_name, Sys.Date())
png(filename = filedh2, width = 960, height = 1000)
DimHeatmap(smc024.merged, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()




smc024.merged <- JackStraw(smc024.merged, num.replicate = 100)
smc024.merged <- ScoreJackStraw(smc024.merged, dims = 1:20) #1:30<-error



filejs=sprintf("07JackStraw_d15_%s_%s.png", project_name, Sys.Date())
png(filename = filejs, width = 660, height = 600)
JackStrawPlot(smc024.merged, dims = 1:15)
dev.off()

fileel=sprintf("08elbow_%s_%s.png", project_name, Sys.Date())
png(filename = fileel, width = 660, height = 600)
ElbowPlot(smc024.merged)
dev.off()




smc024.merged <- FindNeighbors(object = smc024.merged, dims = 1:12, verbose = FALSE)
smc024.merged <- FindClusters(object = smc024.merged, resolution = 0.5, verbose = FALSE)

smc024.merged <- RunUMAP(smc024.merged, dims = 1:12, verbose = FALSE)

save(smc024.merged,file='smc024.merged_pca_umap.Rdata')
#load(file = 'smc024.merged_pca_umap.Rdata' )

# find all markers

smc024.markers <- FindAllMarkers(smc024.merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- smc024.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) 



save(smc024.markers, file='smc024.merged_umap20_markers.Rdata')

smc024.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) %>% select(gene) 

### top10 heatmap

top10 <- smc024.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
filet10h=sprintf("12top10_heatmap_%s_%s.png", project_name, Sys.Date())
png(filename = filet10h, width = 1500, height = 1000)
DoHeatmap(smc024.merged, features = top10$gene) + NoLegend()
dev.off()



write.csv(smc024.markers,file="DEG_SCTsampleType_mt30_rna7000_rmMT.csv")
cellpercluster <- smc024.merged@meta.data$seurat_clusters %>% table
write.csv(cellpercluster, file = "smc024_cell_per_cluster.csv")

smc024_ident_cluster <- smc024.merged@meta.data %>%  select(orig.ident, seurat_clusters) %>% group_by(orig.ident, seurat_clusters) %>% tally()
write.csv(smc024_ident_cluster, file = "smc024_ident_cluster.csv")


filet10h=sprintf("12top10_heatmap_%s_%s.png", project_name, Sys.Date())
png(filename = filet10h, width = 1500, height = 960)
ph1<-DoHeatmap(smc024.merged, features = top10$gene) + NoLegend()
dev.off()



#smc024.merged[[]]
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

save(smc024.merged, file='smc024.merged_umap20.Rdata')

#setwd("./dim11")

# # https://github.com/satijalab/seurat/issues/497
# # https://rdrr.io/cran/Seurat/man/FindConservedMarkers.html




filef1=sprintf("11_1feat_nFRNA_%s_%s.png", project_name, Sys.Date())
png(filename = filef1, width = 660, height = 600)
FeaturePlot(smc024.merged, features=c("nFeature_RNA"))
dev.off()

filef1s=sprintf("11_1feat_split_nFRNA_%s_%s.png", project_name, Sys.Date())
png(filename = filef1s, width = 1660, height = 600)
FeaturePlot(smc024.merged, features=c("nFeature_RNA"), split.by = "orig.ident")
dev.off()

filef2=sprintf("11_2feat_nCRNA_%s_%s.png", project_name, Sys.Date())
png(filename = filef2, width = 660, height = 600)
FeaturePlot(smc024.merged, features=c("nCount_RNA"))
dev.off()

filef2s=sprintf("11_2feat_split_nCRNA_%s_%s.png", project_name, Sys.Date())
png(filename = filef2s, width = 1660, height = 600)
FeaturePlot(smc024.merged, features=c("nCount_RNA"), split.by = "orig.ident")
dev.off()

filef1=sprintf("11_1vln_nFRNA_%s_%s.png", project_name, Sys.Date())
png(filename = filef1, width = 660, height = 600)
VlnPlot(smc024.merged, features=c("nFeature_RNA"))
dev.off()

filef2=sprintf("11_2vln_nCRNA_%s_%s.png", project_name, Sys.Date())
png(filename = filef2, width = 660, height = 600)
VlnPlot(smc024.merged, features=c("nCount_RNA"))
dev.off()



DimPlot(tenx009, group.by = c("tech", "celltype"), combine = FALSE)

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
png(filename = filectc2, width = 600, height = 600)
FeaturePlot(smc024.merged, features = c("EPCAM","CD44","KRT8"), split.by = "orig.ident", max.cutoff = 3)
dev.off()

filectc3=sprintf("13ctc_feat2_%s_%s.png", project_name, Sys.Date())
png(filename = filectc3, width = 600, height = 600)
FeaturePlot(smc024.merged, features = c("KRT19","VIM","MET"), split.by = "orig.ident", max.cutoff = 3)
dev.off()

filectc4=sprintf("13ctc_feat3_%s_%s.png", project_name, Sys.Date())
png(filename = filectc4, width = 600, height = 600)
FeaturePlot(smc024.merged, features = c("KRAS","NRAS","VIM"), split.by = "orig.ident", max.cutoff = 3)
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
png(filename = filectc6, width = 600, height = 600)
FeaturePlot(smc024.merged, features = c("UCHL1", "NCAM1"), split.by = "orig.ident", max.cutoff = 3)
dev.off()


filectc7=sprintf("13ctc_feat7_%s_%s.png", project_name, Sys.Date())
png(filename = filectc7, width = 600, height = 600)
FeaturePlot(smc024.merged, features = c("SYP", "CHGA"), split.by = "orig.ident", max.cutoff = 3)
dev.off()


filectc8=sprintf("13ctc_feat8_%s_%s.png", project_name, Sys.Date())
png(filename = filectc8, width = 600, height = 600)
FeaturePlot(smc024.merged, features = c("ASCL1", "NEUROD1", "POU2F3"), split.by = "orig.ident", max.cutoff = 3)
dev.off()



DoHeatmap(smc024.merged, features = c("EPCAM","CD44","KRT8","KRT19","VIM","MET","KRAS","NRAS"), slot="data")
#DoHeatmap(smc024.merged, features = c("EPCAM","CD44","KRT8","KRT19","VIM","MET","KRAS","NRAS"))



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



#DoHeatmap(smc024.merged, features = c("ERBB2", "BRAF", "TP53", "VEGFA", "PIK3CA", "MTOR", "KRAS", "NRAS", "HRAS", "MAP2K1", "KIT"), slot="data")


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


#DoHeatmap(smc024.merged, features = c("CD3D", "CD3E", "CD3G", "IGFLR1", "RGS1", "TCF7", "NKIRAS2", "IL7R", "CD8A", "KLRG1", "GZMB", "GZMH", "PRF1"), slot="data")

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

# DoHeatmap(smc024.merged, features = c("CD79A", "CD63", "IGHA1", "IGHA2", "IGHD", "VPREB3"), slot="data")

# CTC
filectca1=sprintf("13ctca1_vln_%s_%s.png", project_name, Sys.Date())
png(filename = filectca1, width = 800, height = 400)
VlnPlot(smc024.merged, features=c("AXL", "ERBB2", "ERBB3", "MET", "CSF1R", "KRT19"), split.by = "orig.ident", pt.size = 0)
dev.off()

filectca2=sprintf("13ctca2_feat1_%s_%s.png", project_name, Sys.Date())
png(filename = filectca2, width = 800, height = 600)
FeaturePlot(smc024.merged, features = c("AXL", "ERBB2", "ERBB3"), split.by = "orig.ident", max.cutoff = 3)
dev.off()

filectca3=sprintf("13ctca3_feat2_%s_%s.png", project_name, Sys.Date())
png(filename = filectca3, width = 800, height = 600)
FeaturePlot(smc024.merged, features = c("MET","CSF1R", "KRT19"), split.by = "orig.ident", max.cutoff = 3)
dev.off()



save.image(file = "200410dim12_b1b2t.RData")
load(file = "200410dim12_b1b2t.RData")
