# --------------------------------------------------------------------------- #
# save(list= c("metaall","rTGexpression","TGexpression","TCGA_GTEx","Texp_mod"), file="log2_metaall_T_Gexpression.RData")
# #rm(list = c("metaG", "metaG1", "metaT", "metaT1"))
# 
# saveRDS(TCGA_GTEx, file = "200330TCGA_GTEx.RDS")
# write.csv(TCGA_GTEx, file = "200330TCGA_GTEx.csv")
# --------------------------------------------------------------------------- #
library(dplyr)
setwd("H:TCGA/00ARCHs4_h5")

TCGA_GTEx <- readRDS(file = "200330TCGA_GTEx.RDS")
colnames(TCGA_GTEx)[2] = "Tissue.raw"

stattissue <- function(Tfind, Tout){
  TCGA_GTEx[grep(Tfind,TCGA_GTEx$Tissue.raw),] %>% nrow
  TCGA_GTEx[grep(Tfind,TCGA_GTEx$Tissue.raw),]$Tissue.raw %>% table
}

getTissue <- function(TCGA_GTEx, Tfind, Tout){
  TCGA_GTEx[grep(Tfind,TCGA_GTEx$Tissue.raw),]["Tissue"]=Tout
  return(TCGA_GTEx)
}

stattissue("Bladder","")
dis_df["Tissue"] = ""
dis_df <- getTissue(dis_df,"Testis", "Testis")

table(dis_df$Tissue)

dis_df %>% select(Sampleid,Tissue.raw,Tissue) %>% filter(Tissue == "") %>% select(Tissue.raw) %>% table

dis_df %>% select(Sampleid,Tissue.raw,Tissue) %>% filter(Tissue.raw == "Breast") %>% head
saveRDS(dis_df, file = "200331TCGA_GTEx_match_tissue.RDS")
write.csv(dis_df, file = "200331TCGA_GTEx_match_tissue.csv")
# --------------------------------------------------------------------------- #
setwd("H:TCGA/00ARCHs4_h5")
dis_df <- readRDS("200331TCGA_GTEx_match_tissue.RDS")


# select gene
ingene = "EGFR"
#colnames(TCGA_GTEx) %>% head
select_gene = dis_df[,c("Sampleid", "Tissue", "Sampletype", "Condition", ingene)]
select_gene %>% head

## plot 

library(ggplot2)
library(tidyverse)
p <- ggplot(dis_df, aes(x=Tissue, y= EGFR, fill = Condition)) +
              geom_boxplot()
p +  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# just change label order
p <- ggplot(dis_df, aes(x=Tissue, y= EGFR, fill = Condition)) +
  geom_boxplot()
p +  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_fill_discrete(name = "Sample type",limits = c("Tumor", "adj.Normal", "Normal"))

# facet
p <- ggplot(dis_df, aes(x=Condition, y= EGFR))
p + geom_boxplot(aes(fill = Condition)) +
  facet_wrap(~ Tissue) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# facet grid
p <- ggplot(dis_df, aes(x=Condition, y= EGFR))
p + geom_boxplot(aes(fill = Condition)) +
  facet_grid(~ Tissue) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

