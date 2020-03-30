#https://rpubs.com/tiagochst/TCGAworkshop
library(TCGAbiolinks)
library(MultiAssayExperiment)
library(maftools)
library(dplyr)
library(ComplexHeatmap)

### Clinical data
# Access indexed clinical data
clinical <- GDCquery_clinic("TCGA-BRCA")
head(clinical)

# Same case as figure above
clinical %<%
  dplyr::filter(submitter_id == "TCGA-AA-3562") %>%
  t %>%
  as.data.frame

query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Clinical",
                  data.type = "Clinical Supplement",
                  data.format = "BCR Biotab")

GDCdownload(query)
clinical.BCRtab.all <- GDCprepare(query)

names(clinical.BCRtab.all)

clinical.BCRtab.all$clinical_drug_brca %>%
  head %>%
  as.data.frame

### RNA-Seq data
query.exp.hg38 <- GDCquery(project = "TCGA-BRCA",
                           data.category = "Transcriptome Profiling",
                           data.type = "Gene Expression Quantification",
                           workflow.type = "HTSeq - Counts")

GDCdownload(query.exp.hg38)
raw.counts <- GDCprepare(query = query.exp.hg38, summarizedExperiment = FALSE)

head(raw.counts)

query.exp.hg38 <- GDCquery(project = "TCGA-BRCA",
                           data.category = "Transcriptome Profiling",
                           data.type = "Gene Expression Quantification",
                           workflow.type = "HTSeq - FPKM-UQ")

GDCdownload(query.exp.hg38)
fpkm.uq.counts <- GDCprepare(query = query.exp.hg38, summarizedExperiment = FALSE)


### Mutation

#GDCquery_Maf download the data from GDC
maf <- GDCquery_Maf("BRCA", pipelines = "muse")

maf %>% head %>% as.data.frame

# using maftools for data summary
maftools.input <- maf %>% read.maf

# Check summary
png("plotmaf.png", width = 1000, height = 1000)
plotmafSummary(maf = maftools.input,
               rmOutlier = TRUE,
               addStat = 'median',
               dashboard = TRUE)
dev.off()

png("oncoplot.png", width = 500, height = 500)
oncoplot(maf = maftools.input,
         top = 10,
         removeNonMutated = TRUE)
dev.off()


# classifies Single Nucleotide Variants into Transtions and Transversions
titv = titv(maf = maftools.input,
            plot = FALSE,
            useSyn = TRUE)
png("TiTv.png", width = 600, height = 600)
plotTiTv(res = titv)
dev.off()

getSampleSummary(maftools.input)


# Copy number alteration
rm(query)
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Copy Number Variation",
                  data.type = "Gene Level Copy Number Scores",
                  access = "open")
GDCdownload(query)
scores <- GDCprepare(query)
scores[1:5,1:5]



# Removes metadata from the first 3 columns
scores.matrix <- scores %>% 
  dplyr::select(-c(1:3)) %>% 
  as.matrix

rownames(scores.matrix) <- paste0(scores$'Gene Symbol',"_",scores$Cytoband)

# gain in more than 100 samples
gain.more.than.hundred.samples <- which(rowSums(scores.matrix == 1) >100)

# loss in more than 100 samples
loss.more.than.hunded.samples <- which(rowSums(scores.matrix == -1) > 100)

lines.selected <- c(gain.more.than.hundred.samples, loss.more.than.hunded.samples)

png("heatmap.png", width = 1000, height = 1000)
Heatmap(scores.matrix[lines.selected,],
        show_column_names = FALSE,
        show_row_names = TRUE,
        row_names_gp = gpar(fontsize = 8),
        col = circlize::colorRamp2(c(-1,0,1), colors = c("red","white","blue")))
dev.off()


### DNA methylation data

query_met.hg38 <- GDCquery(project = "TCGA-BRCA", 
                           data.category = "DNA Methylation",
                           platform = "Illumina Human Methylation 27")

GDCdownload(query_met.hg38)
data.hg38 <- GDCprepare(query_met.hg38, summarizedExperiment = TRUE)
data.hg38

data.hg38 %>% rowRanges %>% as.data.frame

data.hg38 %>% colData %>% as.data.frame

data.hg38 %>% assay %>% head %>% as.data.frame


data.hg38 %>% 
  assay %>% 
  rowVars %>% 
  order(decreasing = TRUE) %>% 
  head(10) -> idx

head(idx)

pal_methylation <- colorRampPalette(c("#000436",
                                      "#021EA9",
                                      "#1632FB",
                                      "#6E34FC",
                                      "#C732D5",
                                      "#FD619D",
                                      "#FF9965",
                                      "#FFD32B",
                                      "#FFFC5A"))(100)

png("heatmap1.png", width = 1000, height = 1000)
Heatmap(assay(data.hg38)[idx,],
        show_column_names = TRUE,
        show_row_names = TRUE,
        name = "Methylation Beta-value",
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8),
        col = circlize::colorRamp2(seq(0, 1, by = 1/99), pal_methylation))
dev.off()


save(list=ls(), file='TCGA_BRCA_multi.Rdata')

library("ELMER")
