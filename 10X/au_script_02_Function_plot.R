



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# > head(d)
# gene p_val avg_logFC   avg_FC pct.1 pct.2 p_val_adj cluster
# 1     LYZ     0  2.706244 6.526205 0.947 0.093         0       5
# 2    SPP1     0  2.652908 6.289336 0.702 0.252         0       5
# 3 HLA-DRA     0  2.355595 5.118052 0.960 0.311         0       5
# 4    APOE     0  2.290552 4.892431 0.693 0.187         0       5
# 5    CD74     0  2.240689 4.726226 0.960 0.269         0       5
# 6    SRGN     0  2.105911 4.304696 0.981 0.122         0       5
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++# 


#### set analysis type ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# Enter analysis_type "lower-case"
# "go","c2","c5","c6","c7","kegg","wiki","do","ncg","dgn","cellmarker"


### msigDB Analysis 
#### C2 curated |  C5 GO | C6 oncogenic | C7 immune |


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(dplyr)
library(vroom)
library(magrittr)


### Set Working Directory ###
setwd("H:TCGA")
setwd("./SingleCell/")

d <- read.csv("SMC009.csv")
head(d)

### Set Output Directory ###
dir.create("./test")
setwd("./test")

### Get geneList << input (d, cluster_num) ###
cluster_num = 4
geneList <- gsymbol_input(d,cluster_num)

#### choose Foldchange cut-off << input (geneList) ####
gene <- names(geneList)[abs(geneList) >= 0.58]
cat("print cut-off gene list\n")
cat("--------------------------------------------\n")
gene

### Entrez ID mapping ###
res_entrez <- entrez_input(gene)

analysis_type = "Wiki"

res_enrich <- anno_cluster_enrich(analysis_type, geneList, res_entrez)

res_gsea <- anno_cluster_gsea(analysis_type, geneList, res_entrez)


### create plot << input(analysis1 << res_enrich, res_gsea | plot_name << "total","emap","heat","dot","cnet")
analysis1 <- res_enrich
plot_name = "emap"


create_plot(analysis_type, analysis1, plot_name, 500, 500)
  




##----------------------------------------------------------------------------------##
## function
##----------------------------------------------------------------------------------##
## if cluster condition 필요 없으면 df <- d 
gsymbol_input <- function(d,cluster_num){
  df <- d[d$cluster == cluster_num,] # df <- d
  cat("print head of data frame")
  cat("--------------------------------------------")
  head(df)
  table(df$cluster)
  geneList <- df[,3]
  cat("print all gene list")
  cat("--------------------------------------------")
  head(geneList)
  names(geneList) <- as.character(df[,1])
  geneList <- sort(geneList, decreasing = TRUE)
  bkgd.genes <- df[,1]
  cat("print length of background genes")
  cat("--------------------------------------------")
  length(bkgd.genes)
  bkgd.genes.entrez <- bitr(bkgd.genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
  return(geneList)
}


### Entrez ID mapping ###
entrez_input <- function(gene){
  library(org.Hs.eg.db)
  gene.entrez <- bitr(gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
  ## Remove any Entrez duplicates
  gene.entrez <- gene.entrez[which(duplicated(gene.entrez$ENTREZID) == F), ]
  res_entrez <- gene.entrez$ENTREZID
  res_entrez <- sort(res_entrez, decreasing = TRUE)
  return(res_entrez)
}


# 반드시 analysis_type 지정 
anno_cluster_enrich <- function(analysis_type, geneList, res_entrez){
  gene <- names(geneList)[abs(geneList) >= 0.58]
  cat("print cut-off gene list")
  cat("--------------------------------------------")
  gene
  if(analysis_type == "go"){
    res_enrich <- enrichGO(gene         = names(geneList),
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'SYMBOL',
                     #universe = as.character(bkgd.genes),
                     ont           = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05)
    filegp=sprintf("goplot%s.png",  Sys.Date())
    png(filename = filegp, height=550, width=550, bg="white")
    goplot(res_enrich)
    dev.off()
  }else if(analysis_type == "c6"){
    library(msigdbr)
    msigdbr_show_species()
    m_df <- msigdbr(species = "Homo sapiens")
    head(m_df, 2) %>% as.data.frame
    m_t2g <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
      dplyr::select(gs_name, gene_symbol)
    res_enrich <- enricher(gene, 
                           #universe = as.character(bkgd.genes), 
                           TERM2GENE = m_t2g, )
  }else if(analysis_type == "c5"){
    library(msigdbr)
    msigdbr_show_species()
    m_df <- msigdbr(species = "Homo sapiens")
    head(m_df, 2) %>% as.data.frame
    m_t5g <- msigdbr(species = "Homo sapiens", category = "C5") %>% 
      dplyr::select(gs_name, gene_symbol)
    res_enrich <- enricher(gene, 
    #                 universe = as.character(bkgd.genes), 
                     TERM2GENE = m_t5g, )
  }else if(analysis_type == "c6"){
    library(msigdbr)
    msigdbr_show_species()
    m_df <- msigdbr(species = "Homo sapiens")
    head(m_df, 2) %>% as.data.frame
    m_t6g <- msigdbr(species = "Homo sapiens", category = "C6") %>% 
      dplyr::select(gs_name, gene_symbol)
    res_enrich <- enricher(gene, 
                           #universe = as.character(bkgd.genes), 
                           TERM2GENE = m_t6g)
  }else if(analysis_type == "c7"){
    m_t7g <- msigdbr(species = "Homo sapiens", category = "C7") %>% 
      dplyr::select(gs_name, gene_symbol)
    res_enrich <- enricher(gene, TERM2GENE = m_t7g, universe = as.character(bkgd.genes))
  }else if(analysis_type == "wiki"){
    library(magrittr)
    library(clusterProfiler)
    wp.hs.gmt <- rWikiPathways::downloadPathwayArchive(organism="Homo sapiens", format = "gmt")
    wp2gene <- clusterProfiler::read.gmt(wp.hs.gmt)
    wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
    wpid2gene <- wp2gene %>% dplyr::select(wpid,gene) #TERM2GENE
    wpid2name <- wp2gene %>% dplyr::select(wpid,name) #TERM2NAME
    wpid2gene
    wpid2name
    res_enrich <- clusterProfiler::enricher(
      res_entrez,
      #universe = bkgd.genes.entrez[[2]],
      pAdjustMethod = "fdr",
      pvalueCutoff = 0.1, #p.adjust cutoff; relaxed for demo purposes
      TERM2GENE = wpid2gene,
      TERM2NAME = wpid2name)
    #res_enrich <- DOSE::setReadable(res_enrich, org.Hs.eg.db, keytype = "ENTREZID")
  }else if(analysis_type == "kegg"){
    res_enrich <- enrichKEGG(gene     = res_entrez,
                      organism     = 'hsa',
                      pvalueCutoff = 0.05)
  }else if(analysis_type == "do"){
    library(DOSE)
    res_enrich <- enrichDO(gene = res_entrez,
                     ont = "DO",
                     #universe = bkgd.genes.entrez[[2]],
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "BH",
                     minGSSize     = 5,
                     maxGSSize     = 500,
                     qvalueCutoff  = 0.05,
                     readable      = FALSE)
  }else if(analysis_type == "ncg"){
    #ncg <- enrichNCG(res_entrez, universe = bkgd.genes.entrez[[2]])
    res_enrich <- enrichNCG(res_entrez)
    head(res_enrich)
  }else if(analysis_type == "dgn"){
    res_enrich <- enrichDGN(res_entrez)
  }else if(analysis_type == "cellmarker"){
    library(vroom)
    cell_markers <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt') %>%
      tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>% 
      dplyr::select(cellMarker, geneID) %>%
      dplyr::mutate(geneID = strsplit(geneID, ', '))
    cell_markers
    y <- enricher(res_entrez, 
                  TERM2GENE=cell_markers, 
                  #              universe = bkgd.genes.entrez[[2]],
                  minGSSize=1)
    y <- setReadable(y, 'org.Hs.eg.db', 'ENTREZID')
    DT::datatable(as.data.frame(y))
  }else{
    cat("wrong input of analysis_type!!\n")
  }
  cat("print enrichment %s\n",analysis_type)
  cat("--------------------------------------------\n")
  head(res_enrich)
  return(res_enrich)
}



anno_cluster_gsea <- function(analysis_type, geneList, res_entrez){
  gene <- names(geneList)[abs(geneList) >= 0.58]
  cat("print cut-off gene list")
  cat("--------------------------------------------")
  gene
  if(analysis_type == "go"){
    res_gsea <-gseGO(geneList     = res_entrez,
                     OrgDb        = org.Hs.eg.db,
                     ont          = "CC",
                     nPerm        = 1000,
                     minGSSize    = 100,
                     maxGSSize    = 500,
                     pvalueCutoff = 0.05,
                     verbose      = FALSE)
    filegp=sprintf("goplot%s.png",  Sys.Date())
    png(filename = filegp, height=550, width=550, bg="white")
    goplot(res_gsea)
    dev.off()
  }else if(analysis_type == "c6"){
    library(msigdbr)
    msigdbr_show_species()
    m_df <- msigdbr(species = "Homo sapiens")
    head(m_df, 2) %>% as.data.frame
    m_t2g <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
      dplyr::select(gs_name, gene_symbol)
    res_gsea <- GSEA(geneList, TERM2GENE = m_t2g)
  }else if(analysis_type == "c5"){
    library(msigdbr)
    msigdbr_show_species()
    m_df <- msigdbr(species = "Homo sapiens")
    head(m_df, 2) %>% as.data.frame
    m_t5g <- msigdbr(species = "Homo sapiens", category = "C5") %>% 
      dplyr::select(gs_name, gene_symbol)
    res_gsea <- GSEA(geneList, TERM2GENE = m_t5g)
  }else if(analysis_type == "c6"){
    library(msigdbr)
    msigdbr_show_species()
    m_df <- msigdbr(species = "Homo sapiens")
    head(m_df, 2) %>% as.data.frame
    m_t6g <- msigdbr(species = "Homo sapiens", category = "C6") %>% 
      dplyr::select(gs_name, gene_symbol)
    res_gsea <- GSEA(geneList, TERM2GENE = m_t6g)
  }else if(analysis_type == "c7"){
    m_t7g <- msigdbr(species = "Homo sapiens", category = "C7") %>% 
      dplyr::select(gs_name, gene_symbol)
    res_gsea <- GSEA(geneList, TERM2GENE = m_t7g)
  }else if(analysis_type == "kegg"){
    res_gsea <- gseKEGG(geneList     = geneList,
                        organism     = 'hsa',
                        nPerm        = 1000,
                        minGSSize    = 120,
                        pvalueCutoff = 0.05,
                        verbose      = FALSE)
  }else if(analysis_type == "dgn"){
    res_gsea <- gseDGN(res_entrez,
                       nPerm         = 100,
                       minGSSize     = 120,
                       pvalueCutoff  = 0.2,
                       pAdjustMethod = "BH",
                       verbose       = FALSE)
    #gsdgn <- setReadable(gsdgn, 'org.Hs.eg.db')
  }else{
    cat("wrong input of analysis_type!!\n")
  }
  cat("print GSEA %s",analysis_type)
  cat("--------------------------------------------\n")
  head(res_gsea)
  return(res_gsea)
}



create_plot <- function(analysis_type, analysis1, plot_name, width, height){
  file_plot=sprintf("%s_%s%s.png", plot_name, analysis_type, Sys.Date())
  width <- as.integer(width)
  height <- as.integer(height)
  if(plot_name == "total"){
    file_plot=sprintf("%s_%semap%s.png", plot_name, analysis_type, Sys.Date())
    png(filename = file_plot, height=height, width=width, bg="white")
    pe <- emapplot(analysis1, showCategory = 30)
    print(pe)
    dev.off()
    edox <- analysis1
    ## if EntrezID ### ego doesn't need
    edox <- setReadable(analysis1, 'org.Hs.eg.db', 'ENTREZID')
    head(edox)
    filers=sprintf("%s_result%s.csv", analysis_type, Sys.Date())
    write.csv(edox@result,file = filers)
    file_plot=sprintf("%s_%sheat%s.png", plot_name, analysis_type, Sys.Date())
    png(filename = file_plot, height=height, width=width, bg="white")
    ph <- heatplot(edox, foldChange=geneList)
    print(ph)
    dev.off()
    ### dotplot
    file_plot=sprintf("%s_%sdot%s.png", plot_name, analysis_type, Sys.Date())
    png(filename = file_plot, height=height, width=width, bg="white")
    pd <- dotplot(analysis1, showCategory=30, decreasing=FALSE)
    print(pd)
    dev.off()
    ### cnetplot
    file_plot=sprintf("%s_%scnet%s.png", plot_name, analysis_type, Sys.Date())
    png(filename = file_plot, height=height, width=width, bg="white")
    pc <- cnetplot(edox, showCategory = 5)
    print(pc)
    dev.off()
  }else if(plot_name == "emap"){
    ### emapplot
    png(filename = file_plot, height=height, width=width, bg="white")
    pe <- emapplot(analysis1, showCategory = 30)
    print(pe)
    dev.off()
  }else if(plot_name == "heat"){
    edox <- analysis1
    ## if EntrezID ### ego doesn't need
    edox <- setReadable(analysis1, 'org.Hs.eg.db', 'ENTREZID')
    head(edox)
    filers=sprintf("%s_result%s.csv", analysis_type, Sys.Date())
    write.csv(edox@result,file = filers)
    png(filename = file_plot, height=height, width=width, bg="white")
    ph <- heatplot(edox, foldChange=geneList)
    print(ph)
    dev.off()
  }else if(plot_name == "dot"){
    ### dotplot
    png(filename = file_plot, height=height, width=width, bg="white")
    pd <- dotplot(analysis1, showCategory=30, decreasing=FALSE)
    print(pd)
    dev.off()
  }else if(plot_name == "cnet"){
    ### cnetplot
    png(filename = file_plot, height=height, width=width, bg="white")
    pc <- cnetplot(edox, showCategory = 5)
    print(pc)
    dev.off()
  }else{
    cat("Wrong plot!!")
  }
}




