library(dplyr)
setwd("/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/oncomine/rm_head_full")

sample10_raw = read.csv("sample_10_full_rhead.tsv",sep = "\t", header = F) 
sample16_raw = read.csv("sample_16_full_rhead.tsv" ,sep = "\t", header = F) 
sample15_raw = read.csv("sample_15_full_rhead.tsv",sep = "\t", header = F) 
sample4_raw = read.csv("sample_4_full_rhead.tsv",sep = "\t", header = F) 
sample3_raw = read.csv("sample_3_full_rhead.tsv",sep = "\t", header = F) 
sample9_raw = read.csv("sample_9_full_rhead.tsv",sep = "\t", header = F) 

#colnames(sample10_raw)=c("locus","type","ref","length","genotype","filter","pvalue","cnv_pvalue","coverage","allele_coverage","allele_ratio","%_frequency","maf","hrun","iscn","confidence","precision","Subtype","Call","no_call_reason","genetranscript","location","function","codon","exon","protein","coding","sift","polyphen","grantham","normalizedAlt","5000Exomes","NamedVariants","clinvar","dbsnp","dgv","drugbank","exac","go","pfam","phylop","Oncomine Variant Annotator v2.5","MyVariantDefaultDb_hg19")
row_cal=c("locus","type","ref","length","genotype","filter","pvalue","cnv_pvalue","coverage","allele_coverage","allele_ratio","%_frequency","maf","hrun","iscn","confidence","precision","Subtype","Call","no_call_reason","gene","transcript","location","function","codon","exon","protein","coding","sift","polyphen","grantham","normalizedAlt","5000Exomes","NamedVariants","clinvar","dbsnp","dgv","drugbank","exac","go","pfam","phylop","Oncomine Variant Annotator v2.5","MyVariantDefaultDb_hg19")

colnames(sample10_raw) = row_cal
colnames(sample16_raw) = row_cal
colnames(sample15_raw) = row_cal
colnames(sample4_raw) = row_cal
colnames(sample3_raw) = row_cal
colnames(sample9_raw) = row_cal


sample03 = sample3_raw 
sample04 = sample4_raw
sample09 = sample9_raw
sample10 = sample10_raw 
sample15 = sample15_raw
sample16 = sample16_raw 

sample10["sampleNum"] = "Sample10"
sample16["sampleNum"] = "Sample16"
sample15["sampleNum"] = "Sample15"
sample04["sampleNum"] = "Sample04"
sample03["sampleNum"] = "Sample03"
sample09["sampleNum"] = "Sample09"


sample10_raw["sampleNum"]  %>% table
sample16_raw["sampleNum"] %>% table
sample15_raw["sampleNum"] %>% table
sample4_raw["sampleNum"] %>% table
sample3_raw["sampleNum"] %>% table
sample9_raw["sampleNum"] %>% table

# Sample10 # 89 
# Sample16 # 85 
# Sample15 # 68 
# Sample04 # 199 
# Sample03 # 46 
# Sample09 # 72 

#dir.create("../S03addsample")
setwd("../S03addsample")

write.table(sample03, file = "sample03_mod1.tsv",sep ="\t" ,row.names = F)
write.table(sample04, file = "sample04_mod1.tsv",sep ="\t" ,row.names = F)
write.table(sample09, file = "sample09_mod1.tsv",sep ="\t" ,row.names = F)
write.table(sample10, file = "sample10_mod1.tsv",sep ="\t" ,row.names = F)
write.table(sample15, file = "sample15_mod1.tsv",sep ="\t" ,row.names = F)
write.table(sample16, file = "sample16_mod1.tsv",sep ="\t" ,row.names = F)



# # rbind 2 samples
# rbind_03_04 = rbind(sample03, sample04) 
# rbind_09_10 = rbind(sample09, sample10)
# rbind_15_16 = rbind(sample15, sample16)
# 
# rbind_03_04 %>% select("sampleNum") %>% table
# rbind_09_10 %>% select("sampleNum") %>% table
# rbind_15_16 %>% select("sampleNum") %>% table
# 
# # > rbind_03_04 %>% select("sampleNum") %>% table 
# # Sample03 Sample04 
# # 46      199 
# # > rbind_09_10 %>% select("sampleNum") %>% table
# # Sample09 Sample10 
# # 72       89 
# # > rbind_15_16 %>% select("sampleNum") %>% table
# # Sample15 Sample16 
# # 68       85 
# 
# # rbind all
# rbind_all = rbind(rbind_03_04,rbind_09_10)
# rbind_all = rbind(rbind_all, rbind_15_16)
# 
# rbind_all %>% select(maf) %>% table
# rbind_all %>% select(sampleNum) %>% table
# 
# #dir.create("../S04rbind")
# setwd("../S04rbind")
# write.table(rbind_all, file = "rbindall_mod1.tsv",sep ="\t" ,row.names = F)
# write.table(rbind_03_04, file = "rbind_03_04_mod1.tsv",sep ="\t" ,row.names = F)
# write.table(rbind_09_10, file = "rbind_09_10_mod1.tsv",sep ="\t" ,row.names = F)
# write.table(rbind_15_16, file = "rbind_15_16_mod1.tsv",sep ="\t" ,row.names = F)
# 
# #save.image("raw_rbind200422.Rdata")
# 
# # aggregate
# 
# # CNV_0304 = rbind_03_04 %>% select(c(type, sampleNum, locus, location, coding)) %>%  arrange(type) %>% 
# #   filter(type=="CNV") %>% 
# #   select(c(type, locus, sampleNum)) %>% arrange(locus) 
# 
# aggregate(CNV_0304$sampleNum,list(CNV_0304$type,CNV_0304$locus),FUN = paste)
# 
# # locus
# attach(rbind_15_16)
# 
# res_aggr_0304 = aggregate(sampleNum,list(type,gene,locus),FUN = paste) %>% arrange(Group.1)
# colnames(res_aggr_0304) = c("type","gene","locus", "sampleNum")
# #res_aggr_0304 = sapply(res_aggr_0304, FUN = paste)
# res_aggr_0304 
# 
# #setwd("./locus/")
# #capture.output(cat(res_aggr_0304), "./res_aggr_0910.txt", append = TRUE)
# detach(rbind_15_16)
# #write.csv(as.data.frame(res_aggr_0304), file = "aggr_03_04_locus.tsv",row.names = F)
# 
# # gene
# 
# colnames(rbind_03_04)[45] = "blank"
# sel_0304 = rbind_03_04 %>% select(type, gene, sampleNum) %>% unique.data.frame()
# attach(sel_0304)
# res_aggr_0304 = aggregate(sampleNum,list(type,gene),FUN = paste) %>% arrange(Group.1)
# colnames(res_aggr_0304) = c("type","gene", "sampleNum")
# #res_aggr_0304 = sapply(res_aggr_0304, FUN = paste)
# res_aggr_0304 
# 
# #setwd("./locus/")
# #capture.output(cat(res_aggr_0304), "./res_aggr_0910.txt", append = TRUE)
# detach(rbind_03_04)
# 
# rbind_03_04
# rbind_09_10 
# rbind_15_16 
# 
# setwd("../../")
# save.image("aggr200422.Rdata")

# set1 <- letters[1:5]
# set2 <- letters[4:8]
# set3 <- letters[5:9]
# 
# install.packages("eulerr")
# ## S3 method for class'euler'
# plot(x,fills = TRUE,edges = TRUE,legend = FALSE,labels = identical(legend, FALSE),quantities = FALSE,strips = NULL,main = NULL,n = 200L,adjust_labels = TRUE,...)
# ## S3 method for class'venn'
# plot(x,fills = TRUE,edges = TRUE,legend = FALSE,labels = identical(legend, FALSE),quantities = TRUE,strips = NULL,main = NULL,n = 200L,adjust_labels = TRUE,)



sample03
sample04
sample09
sample10
sample15
sample16

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#-------------------------------------------------------------------#
### gene ####
######## Change variable!!########
library(eulerr)
project_name = "gene03vs04"
# "CNV", "INDEL", "SNV"
typeN = "INDEL"
typeNo = sprintf("%s,NOCALL" , typeN)
sampleA = sample03
sampleB = sample04
#-------------------------------------------------------------------#
sampleA$type %>% table
sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$gene %>% length
sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$gene %>% unique %>% length
sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$gene %>% unique
v_sampleA <- sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$gene 


sampleB$type %>% table
sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$gene %>% length
sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$gene %>% unique %>% length
sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$gene %>% unique
v_sampleB <- sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$gene 

intersect(unique(v_sampleA),unique(v_sampleB)) %>% sort()
setdiff(unique(v_sampleA),unique(v_sampleB)) %>% sort()
setdiff(unique(v_sampleB), unique(v_sampleA)) %>% sort()

titleN =  sprintf("%s %s" , typeN ,  project_name)
filev1=sprintf("venn_%s_%s_%s.png", typeN, project_name, Sys.Date())
png(filename = filev1, width = 300, height = 350)
plot(euler(list(A = unique(v_sampleA), B = unique(v_sampleB))), 
     labels = identical(legend, FALSE), quantities = TRUE, 
     main =  titleN,
     fills = list(fill = c("Orange", "steelblue4"), alpha = 0.5))
dev.off()
euler(list(A = unique(v_sampleA), B = unique(v_sampleB)))

rm(sampleA)
rm(sampleB)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#-------------------------------------------------------------------#
### locus ####
######## Change variable!!########
library(eulerr)
project_name = "locus03vs04"
# "CNV", "INDEL", "SNV"
typeN = "CNV"
typeNo = sprintf("%s,NOCALL" , typeN)
sampleA = sample03
sampleB = sample04
#-------------------------------------------------------------------#
sampleA$type %>% table
sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$locus %>% length
sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$locus %>% unique %>% length
sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$locus %>% unique
v_sampleA <- sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$locus 


v_sampleA_gene_locus <- paste(sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$gene, "~", sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$locus)


v_sampleA_locus_length <- paste(sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$gene, "~", 
                                sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$locus, "~", 
                                sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$length)

### sampleB

sampleB$type %>% table
sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$locus %>% length
sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$locus %>% unique %>% length
sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$locus %>% unique
v_sampleB <- sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$locus


v_sampleB_gene_locus <- paste(sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$gene , "~", 
                              sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$locus)

v_sampleB_locus_length <-paste(sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$gene , "~", 
                               sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$locus , "~", 
                               sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$length)

# gene only
setdiff(unique(v_sampleA),unique(v_sampleB)) %>% sort()
intersect(unique(v_sampleA),unique(v_sampleB)) %>% sort()
setdiff(unique(v_sampleB), unique(v_sampleA)) %>% sort()

# gene_locus
setdiff(unique(v_sampleA_gene_locus),unique(v_sampleB_gene_locus)) %>% sort()
intersect(unique(v_sampleA_gene_locus),unique(v_sampleB_gene_locus)) %>% sort()
setdiff(unique(v_sampleB_gene_locus), unique(v_sampleA_gene_locus)) %>% sort()



# gene_locus_length
setdiff(unique(v_sampleA_locus_length),unique(v_sampleB_locus_length)) %>% sort()
intersect(unique(v_sampleA_locus_length),unique(v_sampleB_locus_length)) %>% sort()
setdiff(unique(v_sampleB_locus_length), unique(v_sampleA_locus_length)) %>% sort()


### write csv
calN = "gene_locus_length"

Aonly <- setdiff(unique(v_sampleA_locus_length),unique(v_sampleB_locus_length)) %>% sort()
Aonly = data.frame(Aonly)
Aonly["sample"] = "Aonly"
colnames(Aonly)[1] = calN


AinterB <- intersect(unique(v_sampleA_locus_length),unique(v_sampleB_locus_length)) %>% sort()
AinterB = data.frame(AinterB)
AinterB["sample"] = "AinterB"
colnames(AinterB)[1] = calN

Bonly <- setdiff(unique(v_sampleB_locus_length), unique(v_sampleA_locus_length)) %>% sort()
Bonly = data.frame(Bonly)
Bonly["sample"] = "Bonly"
colnames(Bonly)[1] = calN


venndf_gene_locus_length = rbind(Aonly, AinterB, Bonly)
filedf=sprintf("vennlist_%s_%s_%s_%s.csv", typeN, project_name, colN, Sys.Date())

# sd1 <- setdiff(unique(v_sampleB), unique(v_sampleA)) %>% sort()
# write.csv(as.data.frame(sd1), "sample4_locus_indel.csv")


### draw venndiagram
titleN =  sprintf("%s %s" , typeN ,  project_name)
filev1=sprintf("venn_%s_%s_%s.png", typeN, project_name, Sys.Date())
png(filename = filev1, width = 300, height = 350)
plot(euler(list(A = unique(v_sampleA), B = unique(v_sampleB))), 
     labels = identical(legend, FALSE), quantities = TRUE, 
     main =  titleN,
     fills = list(fill = c("Orange", "steelblue4"), alpha = 0.5))
dev.off()
euler(list(A = unique(v_sampleA), B = unique(v_sampleB)))
intersect(unique(v_sampleA),unique(v_sampleB)) %>% sort()

rm(sampleA)
rm(sampleB)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#


