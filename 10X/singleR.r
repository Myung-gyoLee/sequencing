# https://static-content.springer.com/esm/art%3A10.1038%2Fs41590-018-0276-y/MediaObjects/41590_2018_276_MOESM1_ESM.pdf

## Obtaining data

################################# Example #######################################
# input : counts file
singler = CreateSinglerSeuratObject(counts.file, annot, project.name,
                                    min.genes = 500, technology, species = "Human" , citation,
                                    normalize.gene.length = F, min.cells = 2, npca = 10,
                                    regress.out = "nUMI", reduce.seurat.object = T)
### Combine multiple 10x data 
dirs = dir(paste0(path,"/10X"), full.names = T)
tenx = Combine.Multiple.10X.Datasets(dirs, random.sample=1000, min.genes=200)

#CreateSinglerSeuratObject() : species = "Human" or "Mouse"
singler = CreateSinglerSeuratObject(tenx$sc.data, tenx$orig.ident, 'SMC',
                                    variable.genes('de'), regress.out='nUMI',
                                    technology='10X', species = 'Human',
                                    citation='Zheng et al.', reduce.file.size = F,
                                    normalize.gene.length = F)
save(singler, file='10x (Zheng) - 10000cells.Rdata')

#################################################################################
library(Seurat)
#devtools::install_github('dviraran/SingleR')
library(SingleR)
library(SingleCellExperiment)
# set working directory

setwd("/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/10X/HN00119076_10X/HN00119076_10X_RawData_Outs/10X_009/H72NHCCX2/seurat_10X_SMC009")

# Read 10x data
cellranger.data <- Read10X(data.dir = "/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/10X/HN00119076_10X/HN00119076_10X_RawData_Outs/10X_009/H72NHCCX2/run_count_10X_009/outs/filtered_feature_bc_matrix")

# create seurat object
SMC009 <- CreateSeuratObject(counts = cellranger.data, project = "SMC009", min.cells = 3, min.features = 200)


singler = CreateSinglerSeuratObject(cellranger.data, annot = NULL, 'SMC009',
                                    min.genes = 200, technology, species = "Human" , citation,
                                    normalize.gene.length = F, min.cells = 3, npca = 10,
                                    regress.out = "nUMI", reduce.seurat.object = T, numCores = 2, 
                                    fine.tune = T, do.signatures = T)

save(singler, file='10x_SMC009_singler.RData')
################################# running Example #######################################
# [1] "Create SingleR object..."
# [1] "Dimensions of counts data: 33538x20461"
# [1] "Annotating data with HPCA..."
# [1] "Variable genes method: de"
# [1] "Number of DE genes:4337"
# [1] "Number of cells: 20461"
# [1] "Fine-tuning round on top cell types (using 2 CPU cores):"
################################# Example #######################################
## SingleR analysis
setwd("/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/10X/HN00119076_10X/HN00119076_10X_RawData_Outs/10X_009/H72NHCCX2/seurat_10X_SMC009")
load('10x_SMC009_singler.RData')

out = SingleR.PlotTsne(singler$singler[[1]]$SingleR.single,
                       singler$meta.data$xy, do.label = F,
                       do.letters = T, labels = singler$meta.data$orig.ident,
                       dot.size = 1.3, alpha=0.5, label.size = 6)
out$p



