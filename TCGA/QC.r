update.packages(ask = FALSE)
library(DEGreport)
library(DESeq2)
setwd("/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/TCGA")
dir.create("00script")

list.dirs()
data(humanGender)

idx <- c(1:10, 75:85)
dds <- DESeqDataSetFromMatrix(assays(humanGender)[[1]][1:1000, idx],
                              colData(humanGender)[idx,], design=~group)

colData(humanGender)
dds <- DESeq(dds)
res <- results(dds)


counts <- counts(dds, normalized = TRUE)
design <- as.data.frame(colData(dds))

degObj(counts, design, "degObj.rda")
library(shiny)
shiny::runGitHub("lpantano/shiny", subdir="expression")


ma = assay(rlog(dds))[row.names(res)[1:100],]
res <- degPatterns(ma, design, time = "group")

filter_count <- degFilter(counts(dds),
                          design, "group",
                          min=1, minreads = 50)
cat("gene in final count matrix", nrow(filter_count))

BiocManager::install('pheatmap')
library(ComplexHeatmap)
th <- HeatmapAnnotation(df = colData(dds),
                        col = degColors(colData(dds), TRUE))
Heatmap(log2(counts(dds) + 0.5)[1:10,],
        top_annotation = th)

degColors(colData(dds), TRUE)

colData(dds)

library(pheatmap)
pheatmap(log2(counts(dds) + 0.5)[1:10,], 
         annotation_col = as.data.frame(colData(dds))[,1:4],
         annotation_colors = degColors(colData(dds)[1:4],
                                       con_values = c("white","red")))
                                       