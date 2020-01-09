'''
conda activate scrna
'''

## load library
library(Seurat)
library(dplyr)
library(ggplot2)
library(sctransform)

## set working directory
setwd("/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/10X/HN00119076_10X/HN00119076_10X_RawData_Outs/10X_009/H72NHCCX2/seurat_10X_SMC009")


## Read 10x data

cellranger.data <- Read10X(data.dir = "/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/10X/HN00119076_10X/HN00119076_10X_RawData_Outs/10X_009/H72NHCCX2/run_count_10X_009/outs/filtered_feature_bc_matrix")


## create seurat object
SMC009 <- CreateSeuratObject(counts = cellranger.data, project = "SMC009", min.cells = 3, min.features = 200)

SMC009 <- PercentageFeatureSet(SMC009, pattern="^MT-", col.name = "percent.mt")

## Visualize QC metrics as a violin plot
png(filename="MT_violin_SMC009.png",width = 800, height=600)

VlnPlot(SMC009, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))

dev.off()

## Visualize QC metrics as a feature scatter plor

plot1 <- FeatureScatter(SMC009, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(SMC009, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

png(filename="QC_scatter_SMC009.png",width = 800, height=600)
CombinePlots(plots=list(plot1, plot2))
dev.off()
