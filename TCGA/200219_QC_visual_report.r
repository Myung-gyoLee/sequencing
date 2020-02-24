#http://lpantano.github.io/DEGreport/articles/DEGreport.html
#http://lpantano.github.io/DEGreport/reference/degCovariates.html

library(DEGreport)
library(DESeq2)
library(BiocParallel)

setwd("H:TCGA")
list.dirs()
#dir.create('./02Breast/023report')


dds = DESeqDataSetFromMatrix(countData = BRCA_TN1,
                             colData = metaTtum_nor1[,c(2,3,5)],
                             design = ~ Condition)
load("./02Breast/021merge/Breast2dds_GTExnormal.Rdata")
register(SnowParam(4))
dds = DESeq(dds)
res = results(dds)
counts = counts(dds, normalized = TRUE)
design = as.data.frame(colData(dds))




## Size factor QC
degCheckFactors(counts[, 1:6])
degQC(counts, design[["Condition"]], pvalue = res[["pvalue"]])

resCov = degCovariates(log2(counts(dds)+0.5),
                       colData(dds), legacy = TRUE)
cor = degCorCov(colData(dds))
names(cor)

resCov$scatterPlot[[1]]

##QC report

createReport(colData(dds)[["Condition"]], counts(dds, normalized = TRUE),
             row.names(res)[1:20], res[["pvalue"]], path = './02Breast/023report')

resultsNames(dds)
#[1] "Intercept"                 "Condition_Tumor_vs_Normal"

degs <- degComps(dds, combs = "Condition",
                 contrast = list("Condition_Tumor_vs_Normal",
                                 c("Condition", "Normal", "Tumor")))


names(degs)
# [1] "Condition_Tumor_vs_Normal" "Condition_Normal_vs_Tumor"

deg(degs[[1]])

deg(degs[[1]], "raw", "tibble")

Brca_sig1=significants(degs[[1]], fc = 0.58, fdr = 0.05)

significants(degs, fc = 0, fdr = 0.05)
significants(degs, fc = 0, fdr = 0.05, full = TRUE)

### Since log2FoldChange are shrunken, the method for DEGSet class now can plot these changes as follow
degMA(degs[[1]], diff = 2, limit = 3)

### plot the original MA plot
degMA(degs[[1]], diff = 2, limit = 3, raw = TRUE)

###  the correlation between the original log2FoldChange and the new ones
degMA(degs[[1]], limit = 3, correlation = TRUE)


## Volcano plots
res[["id"]] <- row.names(res)
show <- as.data.frame(res[1:10, c("log2FoldChange", "padj", "id")])
degVolcano(res[,c("log2FoldChange", "padj")], plot_text = show)
degVolcano(degs[[1]], plot_text = show)


## Gene plots
### Plot top genes coloring by group. Very useful for experiments with nested groups. xs can be time or WT/KO, and group can be treated/untreated. Another classification can be added, like batch that will plot points with different shapes.
degPlot(dds = dds, res = res, n = 6, xs = "Condition")

### Another option for plotting genes in a wide format
degPlotWide(dds, rownames(dds)[1:5], group="Condition")
#rownames(dds)



## Full report
#https://rdrr.io/bioc/DEGreport/src/R/clustering.R

rld = vst(dds)

### Most significants, FDR< 0.05  and log2FC >  0.1 :  20278
dir.create('./02Breast/023report/2.GTExNormal/DEGreport')
resreport <- degResults(dds = dds, rlogMat = assay(rld),  name = "test", org = NULL,
                        FDR = 0.05, do_go = FALSE, FC = 0.1, group = "Condition", 
                        xs = "Condition",
                        path_results = NULL)

## Detect patterns of expression
ma = assay(vst(dds))[row.names(res)[1:100],]
res1 <- degPatterns(ma, design, time = "Condition")
head(res1[["normalized"]])

filept=sprintf("%s/patternDE_%s.png", plotdir,Sys.Date())
png(filename = filept, height=400, width=400, bg="white")
ggplot(res1[["normalized"]],
       aes(Condition, value, color = colored, fill = colored)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(dodge.width = 0.9)) +
  # change the method to make it smoother
  geom_smooth(aes(group=colored), method = "lm")

dev.off()

## Filter genes by group
filter_count <- degFilter(counts(dds),
                          design, "Condition",
                          min=1, minreads = 50)
cat("gene in final count matrix", nrow(filter_count))

## Generate colors for metadata variables
#BiocManager::install('ComplexHeatmap')


library(ComplexHeatmap)
library(circlize)

col_fun = colorRamp2(c(-1,13,25), c("blue","white","orange"))
col_fun(seq(-3,3))

df1 = as.data.frame(log2(counts(dds) + 0.5))
df_top20 = df1[Brca_sig1[1:20],]


Heatmap(df_top20, col = col_fun)


### heatmap with annotation(Condition, sizeFactor, replaceable)

th <- HeatmapAnnotation(df = colData(dds)[3:5],
                        col = degColors(colData(dds)[3:5], TRUE ))

fileht=sprintf("%s/Heatmap_top20_%s.png", plotdir,Sys.Date())
png(filename = fileht, height=1000, width=1000, bg="white")
Heatmap(df_top20, top_annotation = th)
dev.off()

### top 25
df_top25 = df1[Brca_sig1[1:25],]

topnostart=1
topnostop=topnostart+40-1
df_top40 = df1[Brca_sig1[topnostart:topnostop],]

### heatmap with annotation(Condition, sizeFactor, replaceable)

th <- HeatmapAnnotation(df = colData(dds)[3:5],
                        col = degColors(colData(dds)[3:5], TRUE ))

fileht=sprintf("%s/Heatmap%s_%s_%s.png", plotdir,topnostart, topnostop,Sys.Date())
png(filename = fileht, height=1000, width=1800, bg="white")
Heatmap(df_top40, top_annotation = th, name = "Matrix")
dev.off()


save(list=ls(), file='./02Breast/021merge/desBreast1.Rdata')
