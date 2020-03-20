```bash
conda activate scrna
```
## Cheatsheet Seurat
```r
pbmc.counts <- Read10X(data.dir = "~/Downloads/pbmc3k/filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.counts)
pbmc <- NormalizeData(object = pbmc)
pbmc <- FindVariableFeatures(object = pbmc)
pbmc <- ScaleData(object = pbmc)
pbmc <- RunPCA(object = pbmc)
pbmc <- FindNeighbors(object = pbmc)
pbmc <- FindClusters(object = pbmc)
pbmc <- RunTSNE(object = pbmc)
DimPlot(object = pbmc, reduction = "tsne")
```

## load library
```r
library(Seurat)
library(dplyr)
library(ggplot2)
library(sctransform)
```
## set working directory
```r
setwd("/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/10X/HN00119076_10X/HN00119076_10X_RawData_Outs/10X_009/H72NHCCX2/seurat_10X_SMC009")
```

## Read 10x data
```r
cellranger.data <- Read10X(data.dir = "/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/10X/HN00119076_10X/HN00119076_10X_RawData_Outs/10X_009/H72NHCCX2/run_count_10X_009/outs/filtered_feature_bc_matrix")
```

## create seurat object
```r
SMC009 <- CreateSeuratObject(counts = cellranger.data, project = "SMC009", min.cells = 3, min.features = 200)

SMC009 <- PercentageFeatureSet(SMC009, pattern="^MT-", col.name = "percent.mt")
```

## Visualize QC metrics as a violin plot
```r
png(filename="MT_violin_SMC009.png",width = 800, height=600)

VlnPlot(SMC009, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))

dev.off()
```

## Visualize QC metrics as a feature scatter plor
```r
plot1 <- FeatureScatter(SMC009, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(SMC009, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

png(filename="QC_scatter_SMC009.png",width = 800, height=600)
CombinePlots(plots=list(plot1, plot2))
dev.off()
```
### subset
```r
#SMC009 <- subset(SMC009, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt <5)
SMC009 <- subset(SMC009, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt <30)
```
### Normalizing the data
```r
SMC009 <- NormalizeData(SMC009, normalization.method = "LogNormalize", scale.factor = 10000)
```

## run sctransform
```r
SMC009 <- SCTransform(SMC009, vars.to.regress = "percent.mt", conserve.memory = TRUE, verbose = FALSE)
```
### save data
```r
save(SMC009, file = "CMC009_01_10X_Seurat_SCT.Rdata")
```

## run PCA
```r
SMC009 <- RunPCA(SMC009, verbose = FALSE)
```

## save data
```r
save(SMC009, file = "CMC009_02_10X_Seurat_PCA.Rdata")
```
