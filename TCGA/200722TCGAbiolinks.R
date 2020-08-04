#read.csv("200128_tcga_sample_info_subtype_add.tsv")

# https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/clinical.html

library(TCGAbiolinks)
library(dplyr)

query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Clinical",
                  data.type = "Clinical Supplement",
                  data.format = "BCR Biotab")

GDCdownload(query)
clinical.BCRtab.all <- GDCprepare(query)

# All available tables
names(clinical.BCRtab.all)

# colnames from clinical_patient_brca
tibble::tibble(sort(colnames(clinical.BCRtab.all$clinical_patient_brca)))

colnames(clinical.BCRtab.all$clinical_patient_brca)


# ER status count
plyr::count(clinical.BCRtab.all$clinical_patient_brca$er_status_by_ihc)


# ER content 

status.cols <- grep("status_by_ihc", colnames(clinical.BCRtab.all$clinical_patient_brca))
clinical.BCRtab.all$clinical_patient_brca[,c(2,status.cols)] %>% 
  DT::datatable(options = list(scrollX = TRUE))

# barcode + _status_by_ihc
ih_status <- clinical.BCRtab.all$clinical_patient_brca[,c(2,status.cols)]
View(ih_status)

### write csv "clinical_patient_brca_BCRtab_200727.csv"
## write.csv(ih_status, "clinical_patient_brca_BCRtab_200727.csv")

er.cols <- grep("^er",colnames(clinical.BCRtab.all$clinical_patient_brca))
clinical.BCRtab.all$clinical_patient_brca[,c(2,er.cols)] %>% 
  DT::datatable(options = list(scrollX = TRUE))

pr.cols <- grep("^pr",colnames(clinical.BCRtab.all$clinical_patient_brca))
clinical.BCRtab.all$clinical_patient_brca[,c(2,pr.cols)] %>% 
  DT::datatable(options = list(scrollX = TRUE))

her2.cols <- grep("^her2",colnames(clinical.BCRtab.all$clinical_patient_brca))
clinical.BCRtab.all$clinical_patient_brca[,c(2,her2.cols)] %>% 
  DT::datatable(options = list(scrollX = TRUE))



clinical <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")

clinical %>%
  head %>% 
  DT::datatable(filter = 'top', 
                options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),  
                rownames = FALSE)



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
ih_status = read.csv("clinical_patient_brca_BCRtab_200727.csv")
Tbreast = read.csv("TCGA_breast_ARCHS4_meta_200727.csv")

#Tbreast = Tbreast[,c(-1)]

# # normal sample remove
# Tbreast_tumor=Tbreast[-grep("Normal", Tbreast$Sampletype),]

# > colnames(ih_status)
# [1] "bcr_patient_barcode"  "er_status_by_ihc"     "pr_status_by_ihc"     "her2_status_by_ihc"   "nte_pr_status_by_ihc"

# > colnames(Tbreast)
# [1] "CaseID"     "Barcode"    "Cancertype" "Sampletype"

# > head(Tbreast$Barcode)
# [1] TCGA-E2-A1IU-01A-11R TCGA-A1-A0SB-01A-11R TCGA-A2-A04W-01A-31R TCGA-AN-A0AM-01A-11R TCGA-LL-A440-01A-11R TCGA-A7-A26E-01B-06R
# 11192 Levels: TCGA-02-0047-01A-01R TCGA-02-0055-01A-01R TCGA-02-2483-01A-01R TCGA-02-2485-01A-01R TCGA-02-2486-01A-01R ... TCGA-ZX-AA5X-01A-11R

# > head(ih_status$bcr_patient_barcode)
# [1] "bcr_patient_barcode" "CDE_ID:2003301"      "TCGA-3C-AAAU"        "TCGA-3C-AALI"        "TCGA-3C-AALJ"        "TCGA-3C-AALK"       

split_clinBar = substr(Tbreast$Barcode, 1, 12)
Tbreast["Barcode_mod"] = substr(Tbreast$Barcode, 1, 12)

# inspect intersect
intersect(split_clinBar, ih_status$bcr_patient_barcode) %>% length


join_meta_clin = left_join(Tbreast, ih_status, c("Barcode_mod" = "bcr_patient_barcode"))

write.csv(join_meta_clin, "./output/TCGA_BRCA_ER_PR_HER2.csv",row.names = FALSE)

colnames(join_meta_clin)
join_meta_clin[(join_meta_clin$er_status_by_ihc == "Negative") & (join_meta_clin$pr_status_by_ihc == "Negative") & (join_meta_clin$her2_status_by_ihc == "Negative"),] %>% nrow
# 131

join_meta_clin["er_status_by_ihc" == "Negative" , ] %>% length
