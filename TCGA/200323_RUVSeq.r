# https://bioconductor.org/packages/release/bioc/vignettes/RUVSeq/inst/doc/RUVSeq.pdf

library(RUVSeq)
setwd("H:TCGA")
load("matrixARCHs4.Rdata")


## read data tables

Texpression=read.csv("Bar200213TCGAAllcancer.csv")
metaT=read.csv("200214TCGAmeta.csv")
Gexpression=read.csv("200213GTEXall.csv")
metaG=read.csv("200214GTExmeta.csv")

### Breast data

#Ttumor_Gnormal=merge(TBRCA_ex,GBreast_ex, by = 'row.names',all=TRUE)


breastTG1 <- Ttumor_Gnormal[, !duplicated(colnames(Ttumor_Gnormal))]

### set row names
row.names(breastTG1) <- breastTG1[,1]
breastTG1 <- breastTG1[,-1]
filter <- apply(breastTG1, 1, function(x) length(x[x>5])>=2)
filtered <- breastTG1[filter,]
breastTG_sample <- colnames(filtered)


metaTtum_Gnor_mod <- unique(metaTtum_Gnor)
row.names(metaTtum_Gnor_mod) <- metaTtum_Gnor_mod$Sampleid

#metaTtum_Gnor_mod[breastTG_sample,]$Condition

phenobreast <- data.frame(condition=metaTtum_Gnor_mod[breastTG_sample,]$Condition, row.names = breastTG_sample)
#phenobreast <- data.frame(phenobreast[-1,])


set <- newSeqExpressionSet(as.matrix(filtered), phenoData = phenobreast )

set

library(RColorBrewer)
colors <- brewer.pal(3, "Set2")
x <- metaTtum_Gnor_mod[breastTG_sample,]$Condition
plotRLE(set, outline = FALSE, ylim=c(-4,4), col=colors[x])
plotPCA(set, col=colors[x])

### upper-quartile normalization from EDASeq
set <- betweenLaneNormalization(set, which = "upper")
plotRLE(set, outline=FALSE, ylim=c(-4,4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)


### CIDX genes as gtex gene
filter_gtex <- apply(GBreast_ex, 1, function(x) length(x[x>5])>=2)
filtered_gtex <- GBreast_ex[filter_gtex,]
breastGgenes <- rownames(filtered_gtex)

set1 <- RUVg(set, breastGgenes, k=1)

save(set, set1, file = "./03RUVSeq/200324breastset.RData")

pData(set1)
plotPCA(set1, col=colors[x], cex=1.2)


### Differential expression analysis with edgeR
# we are ready to look for differentially expressed genes, using the negative binomial 
# GLM approach implemented in edgeR

design <- model.matrix(~condition + W_1, data = pData(set1))
y <- DGEList(counts = counts(set1), group = x)
y <- calcNormFactors(y, method = "upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)

toptags_lrt <- topTags(lrt)
toptags_lrt

write.table(toptags_lrt, file = "./03RUVSeq/200324BRCAedgeR_toptags_lrt.txt", sep = "\t")
### Differential expression analysis with DESeq2
#### Wald test
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts(set1),
                              colData = pData(set1),
                              design = ~ W_1 + condition)

dds <- DESeq(dds)
res <- results(dds)
res
write.table(as.data.frame(res), file = "./03RUVSeq/200324BRCA_DESeq_Wald.txt", sep = "\t")

#### likelihood ratio test
ddslrt <- DESeq(dds, test = "LRT", reduced=as.formula("~ W_1"))
reslrt <- results(ddslrt)
reslrt

write.table(as.data.frame(reslrt), file = "./03RUVSeq/200324BRCA_DESeq_lrt.txt", sep = "\t")


save(list=c("set","set1","y","design","dds","res","ddslrt","reslrt"), file="./03RUVSeq/200324BRCA.RData")

# https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/

# extract log2 fold change
original_gene_list <- reslrt$log2FoldChange

# name the vector
names(original_gene_list) <- rownames(reslrt)

# omit any NA values
gene_list <- na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# Set the desired organism here
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
gse <- gseGO(geneList=gene_list,
             ont = "ALL",
             keyType = "SYMBOL",
             nPerm = 10000,
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = organism
             )

ego <- enrichGO(gene = names(gene_list),
             ont = "ALL",
             keyType = "SYMBOL",
             pvalueCutoff = 0.05,
             OrgDb = organism
)

dotplot(ego, showCategory=30)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# try tcga whole data

### We filter out non-expressed genes, by requiring more than 5 reads in at least two samples for each gene

filter_TCGA <- apply(Texpression, 1, function(x) length(x[x>5])>=2)
filtered_tcga <- Texpression[filter_TCGA,]
TCGAgenes <- rownames(filtered_tcga)

filter_gtex <- apply(Gexpression, 1, function(x) length(x[x>5])>=2)
filtered_gtex <- Gexpression[filter_gtex,]
gtexgenes <- rownames(filtered_gtex)

library(dplyr)
tcga_gtex <- full_join(as.data.frame(Texpression), as.data.frame(Gexpression), by = NULL, copy = FALSE)
tcga_gtex <- full_join(as.data.frame(filtered_tcga), as.data.frame(filtered_gtex))
length(colnames(filtered_tcga))
length(unique(colnames(filtered_tcga)))

tcga_gtex <-cbind(as.data.frame(filtered_tcga), as.data.frame(filtered_gtex))
saveRDS(tcga_gtex, "cbind_tcga_gtex.rds")

#tcga_gtex <- readRDS("cbind_tcga_gtex.rds")

set <- newSeqExpressionSet(as.matrix(unique(tcga_gtex)), phenoData = data.frame(unique(metaall)$Condition,row.names = unique(metaall)$Sampleid))
### file size error occur !!!

