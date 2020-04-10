#!/usr/bin/env Rscript
library(dplyr)
library(Seurat)
library(ggplot2)
library(sctransform)

setwd('/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/10X/HN00124804/SMC024/analysis')

#setwd("seurat_10X_SMC024/")

# Load the dataset

SMC_024_blood1_data <- Read10X(data.dir = "/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/10X/HN00124804/SMC024/run_count_10x_024_blood1/outs/filtered_feature_bc_matrix")

SMC_024_blood2_data <- Read10X(data.dir = "/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/10X/HN00124804/SMC024/run_count_10x_024_blood2/outs/filtered_feature_bc_matrix")

SMC_024_tissue_data <- Read10X(data.dir = "/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/10X/HN00124804/SMC024/run_count_10x_024_tissue/outs/filtered_feature_bc_matrix")


# SMC_024_blood1_data <- Read10X(data.dir = "/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/10X/HN00124804/HN00124804_result_10X/SMC_024_blood1/filtered_feature_bc_matrix")
# 
# SMC_024_blood2_data <- Read10X(data.dir = "/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/10X/HN00124804/HN00124804_result_10X/SMC_024_blood2/filtered_feature_bc_matrix")
# 
# SMC_024_tissue_data <- Read10X(data.dir = "/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/10X/HN00124804/HN00124804_result_10X/SMC_024_tissue/filtered_feature_bc_matrix")


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
project_name="SMC024"


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

smc024.merged <- subset(smc024.merged, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 30)


# run sctransform
smc024.merged <- SCTransform(smc024.merged, vars.to.regress = "percent.mt", conserve.memory = TRUE, verbose = FALSE)
save(smc024.merged,file='smc024.merged_sct.Rdata')

#tenx009 <- SCTransform(tenx009, vars.to.regress = "percent.mt", conserve.memory = TRUE, verbose = FALSE, return.only.var.genes = FALSE)
#pbmc[["SCT"]]@scale.data contains the residuals (normalized values), and is used directly as input to PCA. Please note that this matrix is non-sparse, and can therefore take up a lot of memory if stored for all genes. To save memory, we store these values only for variable genes, by setting the return.only.var.genes = TRUE by default in the SCTransform function call.
#The ‘corrected’ UMI counts are stored in pbmc[["SCT"]]@counts. We store log-normalized versions of these corrected counts in pbmc[["SCT"]]@data, which are very helpful for visualization.
#You can use the corrected log-normalized counts for differential expression and integration. However, in principle, it would be most optimal to perform these calculations directly on the residuals (stored in the scale.data slot) themselves. This is not currently supported in Seurat v3, but will be soon.




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

#====>
smc024.merged <- JackStraw(smc024.merged, num.replicate = 100)
smc024.merged <- ScoreJackStraw(smc024.merged, dims = 1:15) #1:30<-error

filejs=sprintf("07JackStraw_d15_%s_%s.png", project_name, Sys.Date())
png(filename = filejs, width = 660, height = 600)
JackStrawPlot(smc024.merged, dims = 1:15)
dev.off()

fileel=sprintf("08elbow_%s_%s.png", project_name, Sys.Date())
png(filename = fileel, width = 660, height = 600)
ElbowPlot(smc024.merged)
dev.off()


smc024.merged <- RunPCA(smc024.merged, verbose = FALSE)
smc024.merged <- RunUMAP(smc024.merged, dims = 1:20, verbose = FALSE)
smc024.merged <- FindNeighbors(object = smc024.merged, dims = 1:20, verbose = FALSE)
smc024.merged <- FindClusters(object = smc024.merged, verbose = FALSE)


save(smc024.merged,file='smc024.merged_pca_umap.Rdata')
load(file = 'smc024.merged_pca_umap.Rdata' )


#tenx009 <- FindNeighbors(tenx009, dims = 1:20, verbose = FALSE)
#tenx009 <- FindClusters(tenx009, verbose = FALSE)
#tenx009 <- RunUMAP(tenx009, dims = 1:20, verbose = FALSE)
#tenx009 <- FindNeighbors(tenx009, dims = 1:14, verbose = FALSE)
#tenx009 <- FindClusters(tenx009, verbose = FALSE)
#tenx009 <- RunUMAP(tenx009, dims = 1:14, verbose = FALSE)
#smc024.merged <- FindNeighbors(smc024.merged, dims = 1:10, verbose = FALSE)
#smc024.merged <- FindClusters(smc024.merged, verbose = FALSE)
#smc024.merged <- RunUMAP(smc024.merged, dims = 1:10, verbose = FALSE)
#tenx009 <- FindNeighbors(tenx009, dims = 1:5, verbose = FALSE)
#tenx009 <- FindClusters(tenx009, verbose = FALSE)
#tenx009 <- RunUMAP(tenx009, dims = 1:5, verbose = FALSE)

DimPlot(smc024.merged, label = TRUE) + NoLegend()
DimPlot(smc024.merged, reduction = "umap")

#smc024.merged[[]]
fileum=sprintf("09Dim_umap_%s_%s.png", project_name, Sys.Date())
png(filename = fileum, width = 1160, height = 400)
p1 <- DimPlot(smc024.merged, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(smc024.merged, reduction = "umap", label = TRUE)
p1+p2
dev.off()

save(smc024.merged, file='smc024.merged_umap20.Rdata')



# # https://github.com/satijalab/seurat/issues/497
# # https://rdrr.io/cran/Seurat/man/FindConservedMarkers.html
# smc024.conmarkers <- FindConservedMarkers(smc024.merged,  ident.1 = 0, ident.2 = 1, grouping.var = "orig.ident")
# smc024.conmarkers %>% group_by(orig.ident) %>% top_n(n = 2, wt = avg_logFC)
# 
# write.table(smc024.conmarkers,file="DEG_SCTsampleType_mt30_rna7000_rmMT.tsv",row.names = F,sep="\t")

# ===================>>>>>

#for (i in 0:20){}
filev1=sprintf("10_1violin_%s_%s.png", project_name, Sys.Date())
png(filename = filev1, width = 1060, height = 1000)
VlnPlot(smc024.merged, features = smc024.markers.top10$gene[1:9])
dev.off()


filev2=sprintf("10_2violin_%s_%s.png", project_name, Sys.Date())
png(filename = filev2, width = 1060, height = 1000)
VlnPlot(smc024.merged, features = smc024.markers.top10$gene[10:19])
dev.off()

filev3=sprintf("10_3violin_%s_%s.png", project_name, Sys.Date())
png(filename = filev3, width = 1060, height = 1000)
VlnPlot(smc024.merged, features = smc024.markers.top10$gene[20:29])
dev.off()

VlnPlot(smc024.merged, features=c("CEACAM5", "CEACAM6"))
VlnPlot(smc024.merged, features = c("CEACAM5", "CEACAM6"), slot = "counts", log = TRUE)
FeaturePlot(smc024.merged, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))

filef1=sprintf("11_1feat_nFRNA_%s_%s.png", project_name, Sys.Date())
png(filename = filef1, width = 660, height = 600)
FeaturePlot(smc024.merged, features=c("nFeature_RNA"))
dev.off()

filef2=sprintf("11_2feat_nCRNA_%s_%s.png", project_name, Sys.Date())
png(filename = filef2, width = 660, height = 600)
FeaturePlot(smc024.merged, features=c("nCount_RNA"))
dev.off()


top10 <- smc024.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
filet10h=sprintf("12top10_heatmap_%s_%s.png", project_name, Sys.Date())
png(filename = filet10h, width = 960, height = 1000)
DoHeatmap(smc024.merged) + NoLegend()
dev.off()

DimPlot(tenx009, group.by = c("tech", "celltype"), combine = FALSE)

## marker plotting
#"VIMENTIN"%in%rownames(tenx009[["RNA"]]) : FALSE

#CTC : EPCAM, CD44, KRT8, KRT19, Vimentin, MET, KRAS, NRAS
VlnPlot(tenx009, features=c("EPCAM","CD44","KRT8","KRT19","VIM","MET","KRAS","NRAS"), y.max=4, same.y.lims=TRUE)
VlnPlot(tenx009, features=c("CD44","VIM"), y.max=4, same.y.lims=TRUE)

#VlnPlot(tenx009, features = c("EPCAM","CD44","KRT8","KRT19","VIM","MET","KRAS","NRAS"), slot = "counts", log = TRUE)
FeaturePlot(tenx009, features = c("EPCAM","CD44","KRT8","KRT19","VIM","MET","KRAS","NRAS"))
DoHeatmap(tenx009, features = c("EPCAM","CD44","KRT8","KRT19","VIM","MET","KRAS","NRAS"), slot="data")
#DoHeatmap(tenx009, features = c("EPCAM","CD44","KRT8","KRT19","VIM","MET","KRAS","NRAS"))

#NSCLC : CEA, CYFRA 21-1, EGFR, HER2, BRAF, p53, VEGF, PI3k, mTOR, RAS, MEK, c-KIT
VlnPlot(tenx009, features=c("CEACAM5", "CEACAM6", "KRT19", "EGFR", "ERBB2", "BRAF", "TP53", "VEGFA", "PIK3CA", "MTOR", "KRAS", "NRAS", "HRAS", "MAP2K1", "KIT"), y.max=4, same.y.lims=TRUE)
FeaturePlot(tenx009, features = c("CEACAM5", "CEACAM6", "KRT19", "EGFR", "ERBB2", "BRAF", "TP53", "VEGFA", "PIK3CA", "MTOR", "KRAS", "NRAS", "HRAS", "MAP2K1", "KIT"))
DoHeatmap(tenx009, features = c("CEACAM5", "CEACAM6", "KRT19", "EGFR", "ERBB2", "BRAF", "TP53", "VEGFA", "PIK3CA", "MTOR", "KRAS", "NRAS", "HRAS", "MAP2K1", "KIT"), slot="data")

#T cell : CD3D, CD3E, CD3G, IGFLR1, RGS1, TCF7, NKIRAS2, IL7R, CD8A, KLRG1, GZMB, GZMH, PRF1
VlnPlot(tenx009, features=c("CD3D", "CD3E", "CD3G", "IGFLR1", "RGS1", "TCF7", "NKIRAS2", "IL7R", "CD8A", "KLRG1", "GZMB", "GZMH", "PRF1"), y.max=4, same.y.lims=TRUE)
FeaturePlot(tenx009, features = c("CD3D", "CD3E", "CD3G", "IGFLR1", "RGS1", "TCF7", "NKIRAS2", "IL7R", "CD8A", "KLRG1", "GZMB", "GZMH", "PRF1"))
DoHeatmap(tenx009, features = c("CD3D", "CD3E", "CD3G", "IGFLR1", "RGS1", "TCF7", "NKIRAS2", "IL7R", "CD8A", "KLRG1", "GZMB", "GZMH", "PRF1"), slot="data")

#B cell : CD79A, CD63, IGHA1, IGHA2, IGHD, VPREB3
VlnPlot(tenx009, features=c("CD79A", "CD63", "IGHA1", "IGHA2", "IGHD", "VPREB3"), y.max=5, same.y.lims=TRUE)
FeaturePlot(tenx009, features = c("CD79A", "CD63", "IGHA1", "IGHA2", "IGHD", "VPREB3"))
DoHeatmap(tenx009, features = c("CD79A", "CD63", "IGHA1", "IGHA2", "IGHD", "VPREB3"), slot="data")


######## code reference # https://github.com/satijalab/seurat/issues/497 ###########

# # Create a simulated grouping variable
# object@meta.data$groups <- sample(
#   x = c("g1", "g2", "g3"),
#   size = length(x = object@cell.names),
#   replace = TRUE
# )
# 
# # Create a simulated identity variable
# object@meta.data$ident <- sample(
#   x = c(1, 2),
#   size = length(x = object@cell.names),
#   replace = TRUE
# )
# 
# # Set cell identities to simulated identities
# object <- SetAllIdent(object = object, id = "ident")
# 
# # Find markers that are conserved between the groups for cell identity 1 compared to cell identity 2
# conserved.markers <- FindConservedMarkers(object = object, ident.1 = 1, ident.2 = 2, grouping.var = "groups")
# 
# # Select 'avg_logFC' columns indices
# avg_logFC_columns <- grep(pattern = "avg_logFC", x = colnames(x = conserved.markers))
# 
# # Compute mean 'avg_logFC'
# conserved.markers$mean_avg_logFC <- rowMeans(x = conserved.markers[avg_logFC_columns])
# 
# ## For upregulated genes
# # Order 'conserved.markers' by 'mean_avg_logFC'
# conserved.markers <- conserved.markers[order(conserved.markers$mean_avg_logFC, decreasing = TRUE), ]
# 
# # Get top 10 genes
# rownames(x = conserved.markers[1:10,])
# 
# ## For downregulated genes
# # Order 'conserved.markers' by 'mean_avg_logFC'
# conserved.markers <- conserved.markers[order(conserved.markers$mean_avg_logFC, decreasing = FALSE), ]
# 
# # Get top 10 genes
# rownames(x = conserved.markers[1:10,])
# 
#######################################################################################
