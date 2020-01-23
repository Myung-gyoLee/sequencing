/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/clinical/02cancer33_txt

## loading data
```r
txtdir = "/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/clinical"
setwd(txtdir)

cl_files=list.files("02cancer33_txt", pattern="*.txt", full.names=TRUE)


```

## Make column to Dictionary & count column
```python
import os, glob, pandas as pd
txtdir = "/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/clinical/02cancer33_txt"
txtfiles=glob.glob("%s/*.txt"%txtdir)
txtfiles


heads={}
columnl=[]
for tfile in txtfiles:
    tf=open(tfile)
    rline=[tf.readline().strip().split("\t")]
    rline2=tf.readline().strip().split("\t")
    print(rline)
    heads[tfile.split("/")[-1].split("_")[-1].split(".")[0]]=rline
    for rr in rline2:
        columnl.append(rr)
    tf.close()

columnl
len(columnl)

# count column at table
import collections
counter=collections.Counter(columnl)

# write file
with open("/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/clinical/statcolumn_mod.txt","w") as f:
    print("\n".join(str(counter).split(",")), file = f)

f.close()

col_uniq=list(set(columnl))
len(col_uniq)


with open("/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/clinical/cart/clinical.tsv") as inf:
    cart_head=inf.readline().strip().split("\t")
    print(cart_head)

inf.close()


with open("/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/clinical/clinical.tsv") as inf1:
    clin_head=inf1.readline().strip().split("\t")
    print(clin_head)

inf1.close()



# intersection of col_uniq , cart_head

inter_cart=list(set(col_uniq).intersection(set(cart_head)))
inter_cart
len(inter_cart)


with open("/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/clinical/column_cart.txt","w") as f:
    print("\n".join(cart_head), file = f)

f.close()
```
>>> len(columnl) 2549
>>> len(col_uniq) 700
>>> inter_cart
['days_to_birth', 'ethnicity', 'masaoka_stage', 'race', 'laterality', 'non_nodal_tumor_deposits',
'perineural_invasion_present', 'igcccg_stage', 'circumferential_resection_margin', 'mitotic_count'
, 'gender', 'vital_status', 'supratentorial_localization', 'days_to_death']
>>> len(inter_cart)
14

## intersection of clin_head , cart_head
```python
inter_clin=list(set(clin_head).intersection(set(cart_head)))
inter_clin
len(inter_clin)


with open("/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/clinical/column_inter_brca.txt","w") as f:
    print("\n".join(inter_clin), file = f)

f.close()


with open("/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/clinical/column_clin_brca.txt","w") as f:
    print("\n".join(clin_head), file = f)

f.close()


```

>>> len(inter_clin) 147
>>> len(clin_head) 161
>>> len(cart_head) 161

## TCGA clinical data load
### > length(colnames(clind)) [1] 161
### > length(colnames(clind1)) [1] 55


```r
library(dplyr)

#TCGA clinical data
cart="/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/clinical/cart/clinical.tsv"
clind=read.table(cart, header =  TRUE, sep="\t")

length(colnames(clind))

clind1 = clind[,sapply(clind, nlevels) >1 ]

attach(clind1)
table(treatment_or_therapy)
table(treatment_type)
table(prior_malignancy)
table(prior_treatment)

detach(clind1)

# archs4 data
archf="/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/clinical/200117tcga_norm_sample_info.tsv"
arch_tcga=read.table(archf, header=TRUE, sep="\t")

strsp=arch_tcga$ID
#strsp=strsplit(as.str(arch_tcga$ID),'-')

substr(strsp, 1,12)

> intersect(substr(strsp, 1,12),submitter_id)


> length(arch_tcga$ID) # ARCHs4 metadata
[1] 10406

> length(intersect(substr(strsp, 1,12),submitter_id)) # intersection of ARCHs4 and TCGA clinical cart

[1] 10169

```



## intersect
```bash

cut -f2 /media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/clinical/clinical.tsv |sort|uniq>brca_gdc_bar.txt

cut -f2  /media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/clinical/02cancer33_txt/nationwidechildrens.org_clinical_patient_brca.txt|sort|uniq>brca_nat_bar.txt

#-------------------------------------------------------------------------------#

(base) cytogenbi2@cytogenbi2-B365M-DS3H:/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/clinical$ diff brca_nat_bar.txt brca_gdc_bar.txt
1d0
< CDE_ID:2003301
1099c1098
< bcr_patient_barcode
---
> submitter_id




```

## add demographic information - join 2 table
```r

col_clind=clind[,c("submitter_id","age_at_index","ethnicity","gender","race","vital_status","ajcc_clinical_m","ajcc_clinical_n","ajcc_clinical_stage","ajcc_clinical_t","ajcc_pathologic_m","ajcc_pathologic_n","ajcc_pathologic_stage","ajcc_pathologic_t","ajcc_staging_system_edition","tumor_stage")]


substr(strsp, 1,12)

arch_tcga1=arch_tcga
arch_tcga["submitter_id"]= substr(strsp, 1,12)

arch_join = left_join(arch_tcga, col_clind, by=c("submitter_id"="submitter_id"))

> grep("-02B", strsp, value = TRUE)
[1] "TCGA-14-1034-02B" "TCGA-FG-5965-02B" "TCGA-TQ-A7RK-02B" "TCGA-DU-6407-02B"
[5] "TCGA-DU-6404-02B" "TCGA-DD-AACA-02B"

grep("TCGA-14-1034", strsp, value = TRUE)
grep("TCGA-FG-5965", strsp, value = TRUE)
grep("TCGA-TQ-A7RK", strsp, value = TRUE)
grep("TCGA-DU-6407", strsp, value = TRUE)
grep("TCGA-DU-6404", strsp, value = TRUE)
grep("TCGA-DD-AACA", strsp, value = TRUE)


A02=grep("-02A", strsp, value = TRUE)

for (i in A02){
    ii=substr(strsp, 1,12)
    pline=grep(ii, strsp, value = TRUE)
    print(pline)
}

write.table(arch_join, "/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/clinical/200122_tcga_sample_info_demographic_add.tsv", row.names = FALSE, quote=FALSE, sep = "\t")
```


## Join brca subtype table
```r
patient_brca=read.table("/media/cytogenbi2/6eaf3ba8-a866-4e8a-97ef-23c61f7da612/BreastCancer/data/etc/GDC_Harmonized/clinical_data/txt/nationwidechildrens.org_clinical_patient_brca.txt", header = TRUE, sep = "\t", fill = TRUE)

sub_brca = patient_brca[,c("bcr_patient_barcode","er_status_by_ihc","nte_er_status","pr_status_by_ihc","nte_pr_status_by_ihc","her2_status_by_ihc","nte_her2_status")]

# remove row 1,2
sub_brca=sub_brca[-c(1,2),]
grep("C50.9", sub_brca$bcr_patient_barcode)
sub_brca=sub_brca[-c(686),]

grep("C50.9", sub_brca[-c(686),]$bcr_patient_barcode)


# intersect 2 different barcode id

#> length(intersect(arch_tcga$submitter_id, sub_brca$bcr_patient_barcode))
#[1] 1066
#> Breast Invasive Carcinoma 1134

nrow(arch_tcga[arch_tcga$Group=="Breast Invasive Carcinoma",])
#> nrow(arch_tcga[arch_tcga$Group=="Breast Invasive Carcinoma",])    [1] 1134

# extract brca
arch_brca = arch_tcga[arch_tcga$Group=="Breast Invasive Carcinoma",]


> length(unique(arch_brca$submitter_id))
[1] 1093
> length(unique(sub_brca$bcr_patient_barcode))
[1] 1072


intersect(arch_brca$submitter_id, sub_brca$bcr_patient_barcode)
setdiff(arch_brca$submitter_id, sub_brca$bcr_patient_barcode)
setdiff(sub_brca$bcr_patient_barcode,arch_brca$submitter_id)
# > length(intersect(arch_brca$submitter_id, sub_brca$bcr_patient_barcode))
# [1] 1066

> setdiff(arch_brca$submitter_id, sub_brca$bcr_patient_barcode)
 [1] "TCGA-C8-A12Y" "TCGA-C8-A12Z" "TCGA-C8-A12L" "TCGA-C8-A26Y" "TCGA-C8-A1HK"
 [6] "TCGA-C8-A273" "TCGA-C8-A12X" "TCGA-C8-A12W" "TCGA-C8-A26X" "TCGA-C8-A12T"
[11] "TCGA-C8-A12V" "TCGA-C8-A12O" "TCGA-C8-A1HN" "TCGA-C8-A26V" "TCGA-C8-A1HM"
[16] "TCGA-C8-A12U" "TCGA-C8-A12Q" "TCGA-C8-A1HL" "TCGA-C8-A12P" "TCGA-C8-A1HO"
[21] "TCGA-C8-A12N" "TCGA-C8-A131" "TCGA-C8-A275" "TCGA-C8-A274" "TCGA-C8-A26Z"
[26] "TCGA-C8-A12M" "TCGA-C8-A26W"

> setdiff(sub_brca$bcr_patient_barcode,arch_brca$submitter_id)
[1] "TCGA-A7-A0DC"     "TCGA-AC-A5EI"     "TCGA-AR-A0U1"     "[Not Applicable]"
[5]  "TCGA-C8-A9FZ"



length(unique(arch_brca$submitter_id))
length(unique(sub_brca$bcr_patient_barcode))


```

## deal with not in patient_brca
```r
txtdir="/media/cytogenbi2/6eaf3ba8-a866-4e8a-97ef-23c61f7da612/BreastCancer/data/etc/GDC_Harmonized/clinical_data/txt/"

setwd(txtdir)
follow_up_v2.1_brca=read.table("nationwidechildrens.org_clinical_follow_up_v2.1_brca.txt", header = TRUE, sep = "\t")

follow_up_v4.0_brca=read.table("nationwidechildrens.org_clinical_follow_up_v4.0_brca.txt", header = TRUE, sep = "\t")


attach(follow_up_v2.1_brca)

fu21_bar = bcr_patient_barcode
fu21_fubar = bcr_followup_barcode

detach(follow_up_v2.1_brca)

attach(follow_up_v4.0_brca)

fu40_bar = bcr_patient_barcode
fu40_fubar = bcr_followup_barcode

detach(follow_up_v4.0_brca)

setdiff(diff_fu21, diff_fu40)

diff_fu21=setdiff(fu21_bar,pat_bar) # fu21_bar = follow_up_v2.1_brca
diff_fu40=setdiff(fu40_bar,pat_bar) # fu40_bar = follow_up_v4.0_brca
diff_nte= setdiff(nte_bar,pat_bar)

union1=union(diff_fu21, diff_fu40)

diff_asub=setdiff(arch_brca$submitter_id, sub_brca$bcr_patient_barcode)
setdiff(diff_asub, union1) #[1] "TCGA-C8-A12T" "TCGA-C8-A275"

```



## join
```r
# join
sub_join = left_join(arch_tcga, sub_brca, by=c("submitter_id"="bcr_patient_barcode"))

# write
write.table(arch_join, "/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/clinical/200123_tcga_sample_info_subtype_add.tsv", row.names = FALSE, quote=FALSE, sep = "\t")

```


##
```r

```


##
```r

```


##
```r

```
