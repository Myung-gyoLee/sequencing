setwd("H:SingleCell/smc024/compare")
b1b2_default_raw = read.csv("b1b2_default.csv")
b1b2_option_raw = read.csv("b1b2_option.csv")
b1b2b3_option_raw = read.csv("b1b2b3_option.csv")

head(b1b2_default)
head(b1b2_option)
library(dplyr)

#b1b2_default[,c(1,4,8)] %>% head

b1b2_default = b1b2_default[,c(1,4,8)]
#head(b1b2_default)

b1b2_option = b1b2_option[,c(1,4,9)]

colnames(b1b2_option) = c("gene", "avg_logFC", "cluster" )

over_gene = intersect(b1b2_default$gene, b1b2_option$gene)

b1b2_def_inter = b1b2_default[b1b2_default$gene %in% over_gene,]
b1b2_opt_inter = b1b2_option[b1b2_option$gene %in% over_gene,]

write.csv(b1b2_def_inter,"b1b2_default.csv")
write.csv(b1b2_opt_inter, "b1b2_option_inter.csv")

join_inter = inner_join(b1b2_def_inter, b1b2_opt_inter, by = c("gene" = "gene"))

names(join_inter) = c("gene","avg_logFC.def","cluster.def","avg_logFC.opt","cluster.opt")
write.csv(join_inter, "join_b1b2_def_opt.csv")

## scientific reports - CTC https://pubmed.ncbi.nlm.nih.gov/30068984/ -> non-CTC markers
## table s3
hcc_ctc = read.csv("non_CTC.csv")
head(hcc_ctc)

join_inter %>%left_join(hcc_ctc, by = c("gene" = "Marker.gene"))
join_inter %>%inner_join(hcc_ctc, by = c("gene" = "Marker.gene"))

join_inter %>%inner_join(hcc_ctc, by = c("gene" = "Marker.gene")) %>% tail(48)

aggr_hcc = aggregate(hcc_ctc$Cluster, list(hcc_ctc$Marker.gene), FUN = paste, collapse =",")
colnames(aggr_hcc) = c("Marker.gene", "cluster_ident")
head(aggr_hcc)


###
colnames(b1b2_option_raw)[1]="gene"
colnames(b1b2b3_option_raw)[1]="gene"


option_join_hcc = b1b2_option_raw %>%left_join(aggr_hcc, by = c("gene" = "Marker.gene"))
option_join_hcc = b1b2b3_option_raw %>%left_join(aggr_hcc, by = c("gene" = "Marker.gene"))

write.csv(option_join_hcc, file = "optionb1b2b3_join_hcc.csv")
