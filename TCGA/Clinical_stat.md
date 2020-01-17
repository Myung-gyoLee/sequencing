# python  tracking TCGA patient
## read txt file of TCGA BCR and readline()
```python
txtdir="/media/cytogenbi2/6eaf3ba8-a866-4e8a-97ef-23c61f7da612/BreastCancer/data/etc/GDC_Harmonized/clinical_data/txt/"

drugf=open("%snationwidechildrens.org_clinical_drug_brca.txt"%(txtdir))
drug_brca=drugf.readlines()
drugf.close()

follow15f=open("%snationwidechildrens.org_clinical_follow_up_v1.5_brca.txt"%txtdir)
follow15=follow15f.readlines()
follow15f.close()

follow21f=open("%snationwidechildrens.org_clinical_follow_up_v2.1_brca.txt"%txtdir)
follow21=follow21f.readlines()
follow21f.close()

follow40f=open("%snationwidechildrens.org_clinical_follow_up_v4.0_brca.txt"%txtdir)
follow40=follow40f.readlines()
follow40f.close()


follow40ntef=open("%snationwidechildrens.org_clinical_follow_up_v4.0_nte_brca.txt"%txtdir)
follow40nte=follow40_ntef.readlines()
follow40ntef.close()


ntef=open("%snationwidechildrens.org_clinical_nte_brca.txt"%txtdir)
nte=ntef.readlines()
ntef.close()

patientf=open("%snationwidechildrens.org_clinical_patient_brca.txt"%txtdir)
patient=patientf.readlines()
patientf.close()

```

## from patient_brca get key of Dictionary
>>> patient[0].split("\t")[1]
'bcr_patient_barcode'

```python
# function addFollowup
def addFollowup(patientD,follow,table_name):
    ExPatient=[];barcodel=[]
    for fline in follow:
        barcode=fline.split("\t")[1]
        if "TCGA" not in barcode:
            print(barcode)
        else :
            barcodel.append(barcode)
    print(barcodel)
    for key, value in patientD.items():
        if "TCGA" not in key:
            continue
        elif key not in barcodel:
            patientD[key]+="\t"
        elif key in barcodel:
            patientD[key]+="\t%s"%(table_name)    
    print(ExPatient)
    return patientD


# execute
patientD={pline.split("\t")[1]:"patient_BRCA" for pline in patient}

patientD1=addFollowup(patientD,follow15,"follow_up_v1.5_brca")

patientD2=addFollowup(patientD1,follow21,"follow_up_v2.1_brca")
patientD2

patientD3=addFollowup(patientD2,follow40,"follow_up_v4.0_brca")
patientD3

patientD4=addFollowup(patientD3,follow40nte,"follow_up_v4.0_nte_brca")
patientD4

# write txt file

ouf=open("/media/cytogenbi2/6eaf3ba8-a866-4e8a-97ef-23c61f7da612/BreastCancer/data/etc/GDC_Harmonized/clinical_data/stat/200116track_TCGA_BRCA_clinical_form.txt","w")

for wline in patientD4.keys():
    outline="%s\t%s\n"%(wline, str(patientD4[wline]))
    print(outline)
    ouf.write(outline)

ouf.close()

with open('/media/cytogenbi2/6eaf3ba8-a866-4e8a-97ef-23c61f7da612/BreastCancer/data/etc/GDC_Harmonized/clinical_data/stat/200116track_TCGA_BRCA_clinical.txt','w') as f:
    print(patientD4, file = f)
'''
# code test
patientD1={}
ExPatient=[];barcodel=[]
for fline in follow15:
    barcode=fline.split("\t")[1]
    if "TCGA" not in barcode:
        print(barcode)
    else :
        barcodel.append(barcode)

print(barcodel)
for key, value in patientD.items():
    if "TCGA" not in key:
        continue
    elif key not in barcodel:
        patientD[key]+="\t"
    elif key in barcodel:
        patientD[key]+="\t%s"%("follow_up_v1.5_brca")


print(ExPatient)
patientD
'''
'''
def addFollowup(patientD,follow,table_name):
    ExPatient=[]
    for fline in follow:
        barcode=fline.split("\t")[1]
        if "TCGA" not in barcode:
            print(barcode)
        elif barcode in patientD.keys():
            patientD[barcode]+="\t%s"%(table_name)
        else :
            ExPatient.append(barcode)
    print(ExPatient)
    return patientD

'''

### read BRCA Clinical file
setwd
## Stemness
```r
clin_BRCA = read.table("/media/cytogenbi2/6eaf3ba8-a866-4e8a-97ef-23c61f7da612/BreastCancer/data/etc/TCGAStemness/BRCAClinmut.tsv", header = TRUE, sep = "\t")
```

## read TCGA GDC
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

#patient_brca=read.table("clinical_patient_brca.csv", header = TRUE, sep = "\t", fill = TRUE)


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
## Make dataframe
### barcode |bcr_patient_barcode,er_status_by_ihc,nte_er_status,pr_status_by_ihc,nte_pr_status_by_ihc,her2_status_by_ihc,nte_her2_status
``` r
bar_sub_patient_brca=patient_brca[,c("bcr_patient_barcode","er_status_by_ihc","nte_er_status","pr_status_by_ihc","nte_pr_status_by_ihc","her2_status_by_ihc","nte_her2_status")]

write.csv(bar_sub_patient_brca, "/media/cytogenbi2/6eaf3ba8-a866-4e8a-97ef-23c61f7da612/BreastCancer/data/etc/GDC_Harmonized/clinical_data/stat/bar_sub_patient_brca.csv")

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
