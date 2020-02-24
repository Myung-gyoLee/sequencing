
setwd("H:TCGA")
library(DESeq2)


load('./02Breast/021merge/Breast2dds_GTExnormal.Rdata')

View(counts(dds))
#------------------------------------------------------------------#
## 02. Count Normalization with DESeq median of rations methods

# https://github.com/hbctraining/DGE_workshop/blob/master/lessons/02_DGE_count_normalization.md
# median of ratios method of normalization, DESeq2 has a single estimateSizeFactors() function


#dds <- estimateSizeFactors(dds)

sizeFactors(dds)

normalized_counts <- counts(dds, normalized=TRUE)

### create plotdir
dir.create("./02Breast/023report/2.GTExNormal")
plotdir = "./02Breast/023report/2.GTExNormal"

#dir.create("./02Breast/022normalized_counts")
#write.table(normalized_counts, file="./02Breast/022normalized_counts/200218normalized_counts.txt", sep="\t", quote=F, col.names=NA)
#------------------------------------------------------------------#
## 03. Quality Control
# https://github.com/hbctraining/DGE_workshop/blob/master/lessons/03_DGE_QC_analysis.md
# Transform normalized counts using the rlog transformation
# To improve the distances/clustering for the PCA and heirarchical clustering visualization methods, we need to moderate the variance across the mean by applying the rlog transformation to the normalized counts.

### Transform counts for data visualization
#rld <- rlog(dds, blind=TRUE) # take a long time

#rld = rlogTransformation(dds)
#rld = vst(dds)


# Input is a matrix of log transformed values
#rld <- rlog(dds, blind=T)
rld <- vst(dds)
### Extract the rlog matrix from the object
rld_mat <- assay(rld)    ## assay() is function from the "SummarizedExperiment" package that was loaded when you loaded DESeq2
pca <- prcomp(t(rld_mat))

### Plot PCA 
file1=sprintf("%s/plotPCA_BRCATN_%s.png", plotdir,Sys.Date())
png(filename = file1, height=400, width=400, bg="white")
plotPCA(rld, intgroup="Condition")
dev.off()


# Create data frame with metadata and PC3 and PC4 values for input to ggplot
df <- cbind(meta, pca$x)
ggplot(df) + geom_point(aes(x=PC3, y=PC4, color = Condition))


## Hierarchical Clustering


### Compute pairwise correlation values
rld_cor <- cor(rld_mat)    ## cor() is a base R function

View(rld_cor)   ## check the output of cor(), make note of the rownames and colnames


### Plot heatmap
library(pheatmap)
library(RColorBrewer)
pheatmap(rld_cor)

heat.colors <- brewer.pal(6, "Blues")
pheatmap(rld_cor, color = heat.colors, border_color=NA, fontsize = 10, 
         fontsize_row = 10, height=20)

#------------------------------------------------------------------#
# 04. DESeq2 analysis

#https://github.com/hbctraining/DGE_workshop/blob/master/lessons/04_DGE_DESeq2_analysis.md
## Step 1: Estimate size factors
## Check the size factors
sizeFactors(dds)

## Total number of raw counts per sample
colSums(counts(dds))

## Total number of normalized counts per sample
colSums(counts(dds, normalized=T))

## Step 2: Estimate gene-wise dispersion
#Dispersion is a measure of spread or variability in the data. Variance, standard deviation, IQR, among other measures, can all be used to measure dispersion.

## Plot dispersion estimates
dds <- estimateDispersions(dds)
file2=sprintf("%s/plotDisp_BRCATN_%s.png", plotdir, Sys.Date())
png(filename = file2, height=400, width=400, bg="white")
plotDispEsts(dds)
dev.off()

#------------------------------------------------------------------#
## 05. Differential Expression

#https://github.com/hbctraining/DGE_workshop/blob/master/lessons/05_DGE_DESeq2_analysis2.md
## Shrunken log2 foldchanges (LFC)




## Define contrasts, extract results table, and shrink the log2 fold changes
library(BiocParallel)
register(SnowParam(4))

dds = DESeq(dds)
# dds = DESeq(dds)

contrast_oe <- c("Condition", "Tumor", "Normal")

res_tableOE_unshrunken <- results(dds, contrast=contrast_oe, alpha = 0.05)

res_tableOE <- lfcShrink(dds, contrast=contrast_oe, res=res_tableOE_unshrunken)


# MA plot
# unshrunken
file3=sprintf("%s/plotMA_unshr_BRCATN_%s.png", plotdir, Sys.Date())
png(filename = file3, height=400, width=400, bg="white")
plotMA(res_tableOE_unshrunken, ylim=c(-2,2))
dev.off()

# shrunken
file4=sprintf("%s/plotMA_Shr_BRCATN_%s.png", plotdir, Sys.Date())
png(filename = file4, height=400, width=400, bg="white")
plotMA(res_tableOE, ylim=c(-2,2))
dev.off()

# Summarizing results
summary(res_tableOE)
file5=sprintf("%s/Summary_DE_%s.txt", plotdir, Sys.Date())
capture.output(summary(res_tableOE), 
               file = file5,
               append = TRUE)


#------------------------------------------------------------------#
## Extracting significant differentially expressed genes

### Set thresholds
# The lfc.cutoff is set to 0.58; remember that we are working with log2 fold changes so this translates to an actual fold change of 1.5 which is pretty reasonable.
library(dplyr)
library(tidyverse)
mcols(res_tableOE)

res_tableOE %>% data.frame() %>% View()


padj.cutoff <- 0.05
lfc.cutoff <- 0.58


res_tableOE_tb <- res_tableOE %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

sigOE <- res_tableOE_tb %>%
  subset(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

sigOE
#results(dds, contrast = contrast_oe, alpha = 0.05, lfcThreshold = 0.58)
fileOE=sprintf("%s/SigOE%s.csv", plotdir, Sys.Date())
write.csv(sigOE, file = fileOE)

#------------------------------------------------------------------#
## 06. Visualizing Results
# https://github.com/hbctraining/DGE_workshop/blob/master/lessons/06_DGE_visualizing_results.md
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(DEGreport)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)


######################### have to rewrite ###########################
DEGreport::degPlot(dds = dds, res = res_tableOE_tb, n = 20, xs = "type", group = "condition") # dds object is output from DESeq2

DEGreport::degVolcano(
  data.frame(res[,c("log2FoldChange","padj")]), # table - 2 columns
  plot_text = data.frame(res[1:10,c("log2FoldChange","padj","id")])) # table to add names

# Available in the newer version for R 3.4
DEGreport::degPlotWide(dds = dds, genes = row.names(res)[1:5], group = "condition")

# Create tibbles including row names
B1_meta <- metaTtum_nor1[,c(2,3,5)] %>% 
  rownames_to_column(var="samplename") %>% 
  as_tibble()

normalized_counts <- normalized_counts %>% 
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

# Using ggplot2 to plot multiple genes (e.g. top 20)
## Order results by padj values
top20_sigOE_genes <- res_tableOE_tb %>% 
  arrange(padj) %>% 	#Arrange rows by padj values
  pull(gene) %>% 		#Extract character vector of ordered genes
  head(n=20) 		#Extract the first 20 genes

## normalized counts for top 20 significant genes
top20_sigOE_norm <- normalized_counts %>%
  subset(gene %in% top20_sigOE_genes)

# Gathering the columns to have normalized counts to a single column
gathered_top20_sigOE <- top20_sigOE_norm %>%
  gather(colnames(top20_sigOE_norm)[2:9], key = "Condition", value = "normalized_counts")

## check the column header in the "gathered" data frame
View(gathered_top20_sigOE)

gathered_top20_sigOE <- inner_join(B1_meta, gathered_top20_sigOE)

## plot using ggplot2
ggplot(gathered_top20_sigOE) +
  geom_point(aes(x = gene, y = normalized_counts, color = Condition)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Top 20 Significant DE Genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))

#------------------------------------------------------------------#
#library(reshape)
library(reshape)
norm_counts = counts(dds, normalized = TRUE) # Extract the normalized counts
pseudoCount = log2(norm_counts + 1) # convert to log-scale for visualization
df = melt(pseudoCount) # transpose the matrix
