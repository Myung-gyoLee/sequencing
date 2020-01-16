### read BRCA Clinical file
setwd
## Stemness
```r
clin_BRCA = read.table("/media/cytogenbi2/6eaf3ba8-a866-4e8a-97ef-23c61f7da612/BreastCancer/data/etc/TCGAStemness/BRCAClinmut.tsv", header = TRUE, sep = "\t")
```

## TCGA GDC
```r

txtdir="/media/cytogenbi2/6eaf3ba8-a866-4e8a-97ef-23c61f7da612/BreastCancer/data/etc/GDC_Harmonized/clinical_data/txt/"

setwd(txtdir)

drug_brca=read.table("nationwidechildrens.org_clinical_drug_brca.txt", header = TRUE, sep = "\t")

follow_up_v1.5_brca=read.table("nationwidechildrens.org_clinical_follow_up_v1.5_brca.txt", header = TRUE, sep = "\t")

follow_up_v2.1_brca=read.table("nationwidechildrens.org_clinical_follow_up_v2.1_brca.txt", header = TRUE, sep = "\t")

follow_up_v4.0_brca=read.table("nationwidechildrens.org_clinical_follow_up_v4.0_brca.txt", header = TRUE, sep = "\t")

follow_up_v4.0_nte_brca=read.table("nationwidechildrens.org_clinical_follow_up_v4.0_nte_brca.txt", header = TRUE, sep = "\t")

nte_brca=read.table("nationwidechildrens.org_clinical_nte_brca.txt", header = TRUE, sep = "\t")

omf_v4.0_brca=read.table("nationwidechildrens.org_clinical_omf_v4.0_brca.txt", header = TRUE, sep = "\t")

patient_brca=read.table("nationwidechildrens.org_clinical_patient_brca.txt", header = TRUE, sep = "\t")

patient_brca=read.table("clinical_patient_brca.csv", header = TRUE, sep = "\t", fill = TRUE)


radiation_brca=read.table("nationwidechildrens.org_clinical_radiation_brca.txt", header = TRUE, sep = "\t")
```
[R with SQL](http://the-r.kr/2018/07/03/r-%EC%97%90%EC%84%9C-%EB%8D%B0%EC%9D%B4%ED%84%B0%EB%B2%A0%EC%9D%B4%EC%8A%A4-%EC%9E%91%EC%97%85%ED%95%98%EA%B8%B0/)
```
clinical_TCGA_BRCA<-src_sqlite("TCGA_BRCA", create = TRUE)
copy_to(clinical_TCGA_BRCA,drug_brca,temporary=FALSE)
copy_to(clinical_TCGA_BRCA,follow_up_v1.5_brca,temporary=FALSE)
copy_to(clinical_TCGA_BRCA,follow_up_v2.1_brca,temporary=FALSE)
copy_to(clinical_TCGA_BRCA,follow_up_v4.0_brca,temporary=FALSE)
copy_to(clinical_TCGA_BRCA,follow_up_v4.0_nte_brca,temporary=FALSE)
copy_to(clinical_TCGA_BRCA,nte_brca,temporary=FALSE)
copy_to(clinical_TCGA_BRCA,omf_v4.0_brca,temporary=FALSE)
copy_to(clinical_TCGA_BRCA,patient_brca,temporary=FALSE)
copy_to(clinical_TCGA_BRCA,radiation_brca,temporary=FALSE)
```

## connect
```
db1=src_sqlite("TCGA_BRCA", create =  FALSE)
```
## table list of database
```
src_tbls(db1)
```

```r
library(dplyr)
library(psych)

unique($

unique(follow_up_v1.5_brca$bcr_patient_uuid)
fu15_bar = unique(follow_up_v1.5_brca$bcr_patient_barcode)
fu15_fubar = unique(follow_up_v1.5_brca$bcr_followup_barcode)
unique(follow_up_v1.5_brca$bcr_followup_uuid)

setdiff(fu15_bar, fu15_fubar) # 115
"TCGA-BH-A0RX"
setdiff(fu15_fubar, fu15_bar) # 115
"TCGA-BH-A0RX-F4875"

```

```r
attach(follow_up_v2.1_brca)

fu21_bar = bcr_patient_barcode
fu21_fubar = bcr_followup_barcode

detach(follow_up_v2.1_brca)
```

```r
attach(follow_up_v4.0_nte_brca )

fu40nte_bar = bcr_patient_barcode
fu40nte_fubar = bcr_followup_barcode

detach(follow_up_v4.0_nte_brca )

head(fu40_bar)
```


```r
attach(follow_up_v4.0_brca)

fu40_bar = bcr_patient_barcode
fu40_fubar = bcr_followup_barcode

detach(follow_up_v4.0_brca)

head(fu40_bar)
```

```r
attach(nte_brca)

nte_bar = bcr_patient_barcode
nte_fubar = bcr_followup_barcode

detach(nte_brca)

head(nte_bar)
```


```r
attach(patient_brca)

pat_bar = bcr_patient_barcode
pat_fubar = bcr_followup_barcode

detach(patient_brca)

head(pat_bar)
```
#---------------------------
## get length of barcode id per table
```r
length(

length(fu15_bar)
length(fu21_bar)
length(fu40_bar)
length(fu40nte_bar)
length(nte_bar)
length(pat_bar)

```

## variable of TCGA clinical data

```
$fu15_bar
$fu21_bar
$fu40_bar
$fu40nte_bar
$nte_bar
$pat_bar


fu15_fubar
fu21_fubar
fu40_fubar
fu40nte_fubar
nte_fubar
pat_fubar

#---------------------------
## subtype information table
follow_up_v1.5_brca
follow_up_v2.1_brca     
follow_up_v4.0_brca
follow_up_v4.0_nte_brca
nte_brca
patient_brca
```

#---------------------------
## Total BCR information table
```
drug_brca               
follow_up_v1.5_brca
follow_up_v2.1_brca     
follow_up_v4.0_brca
follow_up_v4.0_nte_brca
nte_brca
omf_v4.0_brca           
patient_brca
radiation_brca
```
#---------------------------
```
setdiff(fu15_bar,fu21_bar)
setdiff(fu15_bar,fu40_bar)
setdiff(fu15_bar,fu40nte_bar)
setdiff(fu15_bar,nte_bar)
setdiff(fu15_bar,pat_bar)


setdiff(fu21_bar,fu15_bar)
setdiff(fu21_bar,fu40_bar)
setdiff(fu21_bar,fu40nte_bar)
setdiff(fu21_bar,nte_bar)
setdiff(fu21_bar,pat_bar)


setdiff(fu40_bar,fu15_bar)
setdiff(fu40_bar,fu21_bar)
setdiff(fu40_bar,fu40nte_bar)
setdiff(fu40_bar,nte_bar)
setdiff(fu40_bar,pat_bar)


setdiff(fu40nte_bar,fu15_bar)
setdiff(fu40nte_bar,fu21_bar)
setdiff(fu40nte_bar,fu40_bar)
setdiff(fu40nte_bar,nte_bar)
setdiff(fu40nte_bar,pat_bar)


setdiff(nte_bar,fu15_bar)
setdiff(nte_bar,fu21_bar)
setdiff(nte_bar,fu40_bar)
setdiff(nte_bar,fu40nte_bar)
setdiff(nte_bar,pat_bar)


setdiff(pat_bar,fu15_bar)
setdiff(pat_bar,fu21_bar)
setdiff(pat_bar,fu40_bar)
setdiff(pat_bar,fu40nte_bar)
setdiff(pat_bar,nte_bar)

setdiff(pat_bar,fu15_bar)
setdiff(pat_bar,fu21_bar)
setdiff(pat_bar,fu40_bar)
setdiff(pat_bar,fu40nte_bar)
setdiff(pat_bar,nte_bar)
```

## intersect patient_brca vs other
```r

intersect(pat_bar,fu15_bar)
intersect(pat_bar,fu21_bar)
intersect(pat_bar,fu40_bar)
intersect(pat_bar,fu40nte_bar)
intersect(pat_bar,nte_bar)
```
#---------------------------
## setdiff() result
#---------------------------------------------------------------------------
> setdiff(fu15_bar,pat_bar) # fu15_bar = follow_up_v1.5_brca
character(0)

> setdiff(fu21_bar,pat_bar) # fu21_bar = follow_up_v2.1_brca
 [1] "TCGA-C8-A12L" "TCGA-C8-A12M" "TCGA-C8-A12N" "TCGA-C8-A12O" "TCGA-C8-A12P"
 [6] "TCGA-C8-A12Q" "TCGA-C8-A12U" "TCGA-C8-A12V" "TCGA-C8-A12W" "TCGA-C8-A12X"
[11] "TCGA-C8-A131" "TCGA-C8-A1HK" "TCGA-C8-A1HL" "TCGA-C8-A1HM" "TCGA-C8-A1HN"
[16] "TCGA-C8-A26Z" "TCGA-C8-A273" "TCGA-C8-A274"

> setdiff(fu40_bar,pat_bar) # fu40_bar = follow_up_v4.0_brca
 [1] "TCGA-C8-A12Y" "TCGA-C8-A12Z" "TCGA-C8-A131" "TCGA-C8-A1HM" "TCGA-C8-A1HN"
 [6] "TCGA-C8-A1HO" "TCGA-C8-A26V" "TCGA-C8-A26W" "TCGA-C8-A26X" "TCGA-C8-A26Y"
[11] "TCGA-C8-A26Z" "TCGA-C8-A273" "TCGA-C8-A274"

> setdiff(fu40nte_bar,pat_bar) # fu40nte_bar = follow_up_v4.0_nte_brca
character(0)

> setdiff(nte_bar,pat_bar) # nte_bar = nte_brca # pat_bar = patient_brca
[1] "TCGA-C8-A12L" "TCGA-C8-A12N"


#---------------------------------------------------------------------------
> table(clin_BRCA$metastatic_breast_carcinoma_immunohistochemistry_pr_pos_cell_score)

                [Not Available]
                            694

> table(clin_BRCA$metastatic_breast_carcinoma_her2_erbb_method_calculation_method_text)

                [Not Available]
              1             693

#### Extract column with sapply nlevels > 1
```r
library(dplyr)
BRCAmutMod=clin_BRCA[,sapply(clin_BRCA, nlevels)>1]
```
##### colnames(clin_BRCA) #57 BRCAmutMod column colnames(BRCAmutMod)  #53
