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

##
```python



```


##
```python



```

##
```python



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


##
```r

```
