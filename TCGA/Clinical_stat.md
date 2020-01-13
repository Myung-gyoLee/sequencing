### read BRCA Clinical file
setwd

clin_BRCA = read.table("/media/cytogenbi2/6eaf3ba8-a866-4e8a-97ef-23c61f7da612/BreastCancer/data/etc/TCGAStemness/BRCAClinmut.tsv", header = TRUE, sep = "\t")

> table(clin_BRCA$metastatic_breast_carcinoma_immunohistochemistry_pr_pos_cell_score)

                [Not Available]
                            694

> table(clin_BRCA$metastatic_breast_carcinoma_her2_erbb_method_calculation_method_text)

                [Not Available]
              1             693

#### Extract column with sapply nlevels > 1
library(dplyr)
BRCAmutMod=clin_BRCA[,sapply(clin_BRCA, nlevels)>1]
##### colnames(clin_BRCA) #57 BRCAmutMod column colnames(BRCAmutMod)  #53
