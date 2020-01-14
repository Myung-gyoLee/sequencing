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

radiation_brca=read.table("nationwidechildrens.org_clinical_radiation_brca.txt", header = TRUE, sep = "\t")
```

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
