library(dplyr)
setwd("/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/oncomine/rm_head_full")

sample10_raw = read.csv("sample_10_full_rhead.tsv",sep = "\t", header = F) 
sample16_raw = read.csv("sample_16_full_rhead.tsv" ,sep = "\t", header = F) 
sample15_raw = read.csv("sample_15_full_rhead.tsv",sep = "\t", header = F) 
sample4_raw = read.csv("sample_4_full_rhead.tsv",sep = "\t", header = F) 
sample3_raw = read.csv("sample_3_full_rhead.tsv",sep = "\t", header = F) 
sample9_raw = read.csv("sample_9_full_rhead.tsv",sep = "\t", header = F) 

#colnames(sample10_raw)=c("locus","type","ref","length","genotype","filter","pvalue","cnv_pvalue","coverage","allele_coverage","allele_ratio","%_frequency","maf","hrun","iscn","confidence","precision","Subtype","Call","no_call_reason","genetranscript","location","function","codon","exon","protein","coding","sift","polyphen","grantham","normalizedAlt","5000Exomes","NamedVariants","clinvar","dbsnp","dgv","drugbank","exac","go","pfam","phylop","Oncomine Variant Annotator v2.5","MyVariantDefaultDb_hg19")
row_cal=c("locus","type","ref","length","genotype","filter","pvalue","cnv_pvalue","coverage","allele_coverage","allele_ratio","%_frequency","maf","hrun","iscn","confidence","precision","Subtype","Call","no_call_reason","genetranscript","location","function","codon","exon","protein","coding","sift","polyphen","grantham","normalizedAlt","5000Exomes","NamedVariants","clinvar","dbsnp","dgv","drugbank","exac","go","pfam","phylop","Oncomine Variant Annotator v2.5","MyVariantDefaultDb_hg19")

colnames(sample10_raw) = row_cal
colnames(sample16_raw) = row_cal
colnames(sample15_raw) = row_cal
colnames(sample4_raw) = row_cal
colnames(sample3_raw) = row_cal
colnames(sample9_raw) = row_cal

sample10_raw["sampleNum"] = "Sample10"
sample16_raw["sampleNum"] = "Sample16"
sample15_raw["sampleNum"] = "Sample15"
sample4_raw["sampleNum"] = "Sample04"
sample3_raw["sampleNum"] = "Sample03"
sample9_raw["sampleNum"] = "Sample09"


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

dir.create("../S03addsample")
setwd("../S03addsample")

write.table(sample3_raw, file = "sample03_mod1.tsv",sep ="\t" ,row.names = F)
write.table(sample4_raw, file = "sample04_mod1.tsv",sep ="\t" ,row.names = F)
write.table(sample9_raw, file = "sample09_mod1.tsv",sep ="\t" ,row.names = F)
write.table(sample10_raw, file = "sample10_mod1.tsv",sep ="\t" ,row.names = F)
write.table(sample15_raw, file = "sample15_mod1.tsv",sep ="\t" ,row.names = F)
write.table(sample16_raw, file = "sample16_mod1.tsv",sep ="\t" ,row.names = F)

```bash
cd /media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/oncomine/S03addsample
cut -f45,1,2,5,6,20,21,24,28,43 *.tsv |sort |uniq -c


awk '/CNV/' sample03_mod1.tsv |cut -f45,1,2,5,6,20,21,24,28
awk '/INDEL/' sample03_mod1.tsv |cut -f45,1,2,5,6,20,21,24,28
awk '/SNV/' sample03_mod1.tsv |cut -f45,1,2,5,6,20,21,24,28

awk '/CNV/' sample*_mod1.tsv |cut -f45,1,2,5,6,20,21,24,28
awk '/INDEL/' sample*_mod1.tsv |cut -f45,1,2,5,6,20,21,24,28
awk '/SNV/' sample*_mod1.tsv |cut -f45,1,2,5,6,20,21,24,28


/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/oncomine/S07twogroup/1.group030915/

/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/oncomine/S07twogroup/2.group041016/



awk '/CNV/' /media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/oncomine/S07twogroup/1.group030915/sample*_mod1.tsv |cut -f45,1,2,5,6,20,21,24,28

awk '/INDEL/' /media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/oncomine/S07twogroup/1.group030915/
sample*_mod1.tsv |cut -f45,1,2,5,6,20,21,24,28

awk '/SNV/' /media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/oncomine/S07twogroup/1.group030915/
sample*_mod1.tsv |cut -f45,1,2,5,6,20,21,24,28



awk '/CNV/' /media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/oncomine/S07twogroup/2.group041016/
sample*_mod1.tsv |cut -f45,1,2,5,6,20,21,24,28

awk '/INDEL/' /media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/oncomine/S07twogroup/2.group041016/
sample*_mod1.tsv |cut -f45,1,2,5,6,20,21,24,28

awk '/SNV/' /media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/oncomine/S07twogroup/2.group041016/
sample*_mod1.tsv |cut -f45,1,2,5,6,20,21,24,28



awk '/CNV/' /media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/oncomine/S07twogroup/1.group030915/sample*_mod1.tsv |cut -f1,21 |sort|uniq -d 

awk '/INDEL/' /media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/oncomine/S07twogroup/1.group030915/sample*_mod1.tsv |cut -f1,21 |sort|uniq -d 

awk '/SNV/' /media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/oncomine/S07twogroup/1.group030915/sample*_mod1.tsv |cut -f1,21 |sort|uniq -d 


awk '/CNV/' /media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/oncomine/S07twogroup/2.group041016/sample*_mod1.tsv |cut -f1,6,21 |sort|uniq -d 

awk '/INDEL/' /media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/oncomine/S07twogroup/2.group041016/sample*_mod1.tsv |cut -f1,21 |sort|uniq -d 

awk '/SNV/' /media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/oncomine/S07twogroup/2.group041016/sample*_mod1.tsv |cut -f1,21 |sort|uniq -d 


```

