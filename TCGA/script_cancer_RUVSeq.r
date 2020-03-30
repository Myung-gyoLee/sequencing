# https://bioconductor.org/packages/release/bioc/vignettes/RUVSeq/inst/doc/RUVSeq.pdf


setwd("H:TCGA/03RUVSeq/")
#setwd("../")

#save(list=c("Texpression","Gexpression","metaall","metaT","metaG"), file="matrixARCHs4.Rdata")
load("matrixARCHs4.Rdata")

organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

#------------------------------------------------------------------#
## TCGA-BRCA Primary Tumor and Metastatic Expression Table
#  , 
tissue="Kidney"
#setwd("../03RUVSeq/")
tissue_dir=sprintf("./%s",tissue)

dir.create(tissue_dir)
setwd(tissue_dir)

#------------------------------------------------------------------#
#Rprof("TCGA_de.out")
TCGA_de(Texpression, Gexpression, metaall, metaT, metaG, tissue)

#------------------------------------------------------------------#
TCGA_de <- function(Texpression, Gexpression, metaall, metaT, metaG, tissue){
  Tcancer=metaT[grep(tissue, metaT$Cancertype),]
  
  # normal sample remove
  Ttumor=Tcancer[-grep("Normal", Tcancer$Sampletype),]
  # barcode of Primary tumor and metastatic-> indexing Expression file
  Ttumor_bar=Ttumor$Barcode
  tcga_ex=Texpression[,Ttumor_bar]
  print("#------------------------------------------------------------------#")
  ## GTEx Breast Tissue Expression Table
  Gtissue=metaG[grep(tissue, metaG$Tissue),]
  Gtissue_id=Gtissue$Sampleid
  gtex_ex=Gexpression[,Gtissue_id]
  
  ## join 2) TCGA tumor vs GTEx normal
  Ttumor_Gnormal=merge(tcga_ex,gtex_ex, by = 'row.names',all=TRUE)
  sample_TGnor=rbind(Ttumor_bar,Gtissue_id)
  metaT_Gnor=rbind()
  metaT_Gnor=rbind(metaall[metaall$Condition == "Tumor" & metaall$Tissue %in% unique(Ttumor$Cancertype),]
                   ,metaall[metaall$Project == "GTEx" & metaall$Tissue == tissue,])
  
  print("## join 2) TCGA tumor vs GTEx normal done")
  print("#------------------------------------------------------------------#")
  
  ### Each cancer type data
  # Ttumor_Gnormal
  # metaTtum_Gnor
  
  library(RUVSeq)
  
  breastTG1 <- Ttumor_Gnormal[, !duplicated(colnames(Ttumor_Gnormal))]
  
  ### set row names
  row.names(breastTG1) <- breastTG1[,1]
  breastTG1 <- breastTG1[,-1]
  filter <- apply(breastTG1, 1, function(x) length(x[x>5])>=2)
  filtered <- breastTG1[filter,]
  breastTG_sample <- colnames(filtered)
  
  
  metaT_Gnor <- unique(metaT_Gnor)
  row.names(metaT_Gnor) <- metaT_Gnor$Sampleid
  
  phenobreast <- data.frame(condition=metaT_Gnor[breastTG_sample,]$Condition, row.names = breastTG_sample)
  
  print("#------------------------------------------------------------------#")
  print("set newSeqExpressionSet")
  
  set <- newSeqExpressionSet(as.matrix(filtered), phenoData = phenobreast )
  
  library(RColorBrewer)
  colors <- brewer.pal(3, "Set2")
  x <- metaT_Gnor[breastTG_sample,]$Condition
  
  print("upper-quartile normalization")
  ### upper-quartile normalization from EDASeq
  set <- betweenLaneNormalization(set, which = "upper")
  
  
  ### CIDX genes as gtex gene
  filter_gtex <- apply(gtex_ex, 1, function(x) length(x[x>5])>=2)
  filtered_gtex <- gtex_ex[filter_gtex,]
  breastGgenes <- rownames(filtered_gtex)
  
  print("analysis RUVg")
  set1 <- RUVg(set, breastGgenes, k=1)
  
  pData(set1)
  plotPCA(set1, col=colors[x], cex=0.4)
  
  
  ### Differential expression analysis with edgeR
  # we are ready to look for differentially expressed genes, using the negative binomial 
  # GLM approach implemented in edgeR
  print("#------------------------------------------------------------------#")
  print("edgeR start")
  design <- model.matrix(~condition + W_1, data = pData(set1))
  y <- DGEList(counts = counts(set1), group = x)
  y <- calcNormFactors(y, method = "upperquartile")
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit, coef=2)
  
  toptags_lrt <- topTags(lrt)
  
  edgerf = sprintf("./%s%s_edgeR.txt", Sys.Date(),tissue)
  write.table(toptags_lrt, file = edgerf, sep = "\t")
  ### Differential expression analysis with DESeq2
  #### Wald test
  print("#------------------------------------------------------------------#")
  print("DESeq2 start")
  library(DESeq2)
  dds <- DESeqDataSetFromMatrix(countData = counts(set1),
                                colData = pData(set1),
                                design = ~ W_1 + condition)
  dds <- DESeq(dds)
  res <- results(dds)
  
  waldf = sprintf("./%s%s_DESeq_Wald.txt", Sys.Date(),tissue)
  write.table(as.data.frame(res), file = waldf, sep = "\t")
  
  #### likelihood ratio test
  print("DESeq2 likelihood ratio test start")
  ddslrt <- DESeq(dds, test = "LRT", reduced=as.formula("~ W_1"))
  reslrt <- results(ddslrt)
  
  lrtf = sprintf("./%s%s_DESeq_lrt.txt", Sys.Date(),tissue)
  write.table(as.data.frame(reslrt), file = lrtf, sep = "\t")
  print("DESeq2 done")
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
  BiocManager::install(organism, character.only = TRUE, update = TRUE, ask = FALSE)
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
  print("gse done")
  ego <- enrichGO(gene = names(gene_list),
                  ont = "ALL",
                  keyType = "SYMBOL",
                  pvalueCutoff = 0.05,
                  OrgDb = organism
  )
  print("ego done")
  file_plot=sprintf("%sdotplot_%s.png", tissue, Sys.Date())
  png(filename = file_plot, height=500, width=600, bg="white")
  p <- dotplot(ego, showCategory=30)
  print(p)
  dev.off()
  
  rdf = sprintf("./%s%s.RData", Sys.Date(),tissue)
  save(list=c("set","set1","y","design","dds","res","ddslrt","reslrt", "gse", "ego"), file=rdf)
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#Rprof(NULL)
#summaryRprof("TCGA_de.out")$by.self

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
