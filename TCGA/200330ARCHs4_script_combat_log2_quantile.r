library(rhdf5)
library(preprocessCore)
setwd("H:TCGA/00ARCHs4_h5")
load(file="log2_quantile_normalization.RData")

destination_file = "tcga_matrix.h5"
extracted_expression_file = "tcga_expression_log2.tsv"

# normalize samples and correct for differences in gene count distribution
expression = h5read(destination_file, "data/expression")
expression = log2(expression+1)
expression = normalize.quantiles(expression)

genes = h5read(destination_file, "meta/genes")
samples = h5read(destination_file, "/meta/gdc_cases.samples.portions.analytes.submitter_id")
series = h5read(destination_file, "meta/gdc_center.code")


rownames(expression) = genes
colnames(expression) = samples 

library(sva)
# correct batch effects in gene expression
batchid = match(series, unique(series))
Texpression <- ComBat(dat = expression, batch = batchid, par.prior = TRUE, prior.plots = FALSE)
#Texpression <- correctedExpression
#rm(correctedExpression)

#colnames(correctedExpression) = samples


# gtex data
gexpression = h5read("gtex_matrix.h5", "data/expression")
gexpression = log2(gexpression+1)
gexpression = normalize.quantiles(gexpression)

ggenes = h5read("gtex_matrix.h5", "meta/genes")
gsamples = h5read("gtex_matrix.h5", "meta/sampid")
gseries = h5read("gtex_matrix.h5", "meta/smnabtch")

rownames(gexpression) = ggenes
colnames(gexpression) = gsamples

gbatchid = match(gseries, unique(gseries))
Gexpression <- ComBat(dat = gexpression, batch = gbatchid, par.prior = TRUE, prior.plots = FALSE)

rm(expression)
rm(gexpression)
rm(batchid)
rm()

rlist = ls()

rlist = rlist[rlist != "Texpression"]
rlist = rlist[rlist != "Gexpression"]
              
save(Gexpression,Texpression, file="log2_quantile_normalization.RData")

rm(list = rlist)
