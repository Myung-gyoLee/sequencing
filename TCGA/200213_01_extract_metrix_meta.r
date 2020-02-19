
library(rhdf5)
library(preprocessCore)


install.packages(c("DESeq2","tidyverse","RColorBrewer","pheatmap","ggrepel","cowplot","clusterProfiler","DEGreport","org.Hs.eg.db","DOSE","pathview","tximport","AnnotationDbi","EnsDb.Hsapiens.v86","AnnotationHub","ensembldb"))

BiocManager::install("preprocessCore")
BiocManager::install("rhdf5")
BiocManager::install(c("clusterProfiler","DEGreport","org.Hs.eg.db","DOSE","tximport","AnnotationDbi","EnsDb.Hsapiens.v86","AnnotationHub","ensembldb"))

setwd('H:TCGA')

## TCGA h5 file

TCGAh5<-H5Fopen("tcga_matrix.h5")
Tgene=h5read(TCGAh5, "meta/genes")
analytes.submitter_id=h5read(TCGAh5, "/meta/gdc_cases.samples.portions.analytes.submitter_id")
Tsample_id=h5read(TCGAh5,"/meta/sampleid")
Texpression = h5read(TCGAh5, "data/expression")

Tcancertype=h5read(TCGAh5,"/meta/cancertype")
Tsample_type=h5read(TCGAh5,"/meta/gdc_cases.samples.sample_type")      

site_locations = which(Tprimarysite %in% "Breast")                     
analytes.submitter_id[site_locations]

rownames(Texpression) = Tgene
#colnames(Texpression) = Tsample_id
colnames(Texpression) = analytes.submitter_id

#write.csv(Texpression,file="200213TCGAAllcancer.csv")
#write.csv(Texpression,file="Bar200213TCGAAllcancer.csv")

## make TCGA meta file "Sampleid","Barcode","Cancertype","Sampletype"
metaT<-data.frame(matrix(nrow=length(analytes.submitter_id), ncol=4))
colnames(metaT)=c("CaseID","Sampleid","Cancertype","Sampletype")
#rownames(metaT)

metaT["Sampleid"]=analytes.submitter_id
metaT["Cancertype"]=Tcancertype
metaT["Sampletype"]=Tsample_type
metaT["CaseID"]=gdc_cases.case_id
metaT$Sampletype=paste0("TCGA ", metaT$Sampletype)
metaT["Project"]="TCGA"


#head(metaT)
#table(metaT$Sampletype)
#write.csv(metaT,file="200214TCGAmeta.csv")
#head(metaT[,c("Project","Barcode","Cancertype", "Sampletype")])

#metaT1=metaT[,c("Project","Barcode","Cancertype", "Sampletype")]
#colnames(metaT1)[2]="Sampleid"    
#colnames(metaT1)


#table(metaT1[grep("Metastatic",metaT1$Sampletype),]["Sampletype"])


metaT1[grep("Normal",metaT1$Sampletype),]["Condition"]="Normal"
metaT1["Condition"]="Tumor"
#table(metaT1$Condition)
#head(metaT1)



## gtex h5 file

#h5ls("gtex_matrix.h5", all=TRUE)
gtexh5<-H5Fopen("gtex_matrix.h5")
Gsample_id=h5read(gtexh5, "meta/sampid")
#head(Gsample_id)
Gexpression=h5read(gtexh5, "data/expression")
Gtissue=h5read(gtexh5, "meta/tissue")
Ggene=h5read(gtexh5, "meta/genes")
#head(Ggene)
Gplatform = h5read(gtexh5, "meta/smnabtcht")
Gbatch=h5read(gtexh5,"meta/smnabtch")

rownames(Gexpression) = Ggene
colnames(Gexpression) = Gsample_id

head(Gexpression)
write.csv(Gexpression,file="200213GTEXall.csv")


## make gtex meta file "Sampleid","Tissue","Batch", "Experiment"
metaG<-data.frame(matrix(nrow=length(Gsample_id), ncol=4))
colnames(metaG)=c("Sampleid","Tissue","Batch", "Experiment")
#rownames(meta1)

metaG["Sampleid"]=Gsample_id
metaG["Tissue"]=Gtissue
metaG["Batch"]=Gbatch
metaG["Experiment"]=Gplatform
metaG["Sampletype"]="GTEx Normal"
metaG["Project"]="GTEx"
metaG["Condition"]="Normal"

metaG1=metaG[,c("Project","Sampleid","Tissue", "Sampletype","Condition")]
#metaG1["Condition"]="Normal"
#head(metaG)
#write.csv(metaG,file="200214GTExmeta.csv")

metaall=rbind(metaT1, metaG1)
tail(metaall)
nrow(metaall)
nrow(metaT1) + nrow(metaG1)


h5closeAll()
