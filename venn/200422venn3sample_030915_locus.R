#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
### CNV
#-------------------------------------------------------------------#
### locus ####
######## Change variable!!########
library(eulerr)
project_name = "locus04vs10vs16"
rm(sampleA)
rm(sampleB)
rm(sampleC)
sampleA = sample04
sampleB = sample10
sampleC = sample16


# "CNV", "INDEL", "SNV"
typeN = "CNV"
typeNo = sprintf("%s,NOCALL" , typeN)


#-------------------------------------------------------------------#
### sampleA
sampleA$type %>% table
sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$locus %>% length
sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$locus %>% unique %>% length
sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$locus %>% unique
v_sampleA <- sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$locus 


v_sampleA_gene_locus <- paste(sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$gene, "~", sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$locus)


v_sampleA_locus_length <- paste(sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$gene, "~", 
                                sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$locus, "~", 
                                sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$length)

### sampleB

sampleB$type %>% table
sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$locus %>% length
sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$locus %>% unique %>% length
sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$locus %>% unique
v_sampleB <- sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$locus


v_sampleB_gene_locus <- paste(sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$gene , "~", 
                              sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$locus)

v_sampleB_locus_length <-paste(sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$gene , "~", 
                               sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$locus , "~", 
                               sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$length)



### sampleC

sampleC$type %>% table
sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$locus %>% length
sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$locus %>% unique %>% length
sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$locus %>% unique
v_sampleC <- sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$locus


v_sampleC_gene_locus <- paste(sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$gene , "~", 
                              sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$locus)

v_sampleC_locus_length <-paste(sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$gene , "~", 
                               sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$locus , "~", 
                               sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$length)


### draw venndiagram
titleN =  sprintf("%s %s" , typeN ,  project_name)
filev1=sprintf("venn_%s_%s_%s.png", typeN, project_name, Sys.Date())
png(filename = filev1, width = 300, height = 350)
plot(euler(list(A = unique(v_sampleA_locus_length), B = unique(v_sampleB_locus_length), C = unique(v_sampleC_locus_length))), 
     labels = identical(legend, FALSE), quantities = TRUE, 
     main =  titleN,
     fills = list(fill = c("Pink", "steelblue4", "Green"), alpha = 0.5))
dev.off()
euler(list(A = unique(v_sampleA_locus_length), B = unique(v_sampleB_locus_length), C = unique(v_sampleB_locus_length)))

# rm(sampleA)
# rm(sampleB)
# rm(sampleC)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#


for(i in v_sampleA_locus_length){cat(i);cat(";", sep = "\n")}
for(i in v_sampleB_locus_length){cat(i);cat(";", sep = "\n")}
for(i in v_sampleC_locus_length){cat(i);cat(";", sep = "\n")}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
### INDEL
#-------------------------------------------------------------------#
### locus ####
######## Change variable!!########
library(eulerr)

# "CNV", "INDEL", "SNV"
typeN = "INDEL"
typeNo = sprintf("%s,NOCALL" , typeN)

#-------------------------------------------------------------------#
### sampleA
sampleA$type %>% table
sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$locus %>% length
sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$locus %>% unique %>% length
sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$locus %>% unique
v_sampleA <- sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$locus 


v_sampleA_gene_locus <- paste(sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$gene, "~", sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$locus)


v_sampleA_locus_length <- paste(sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$gene, "~", 
                                sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$locus, "~", 
                                sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$length)

### sampleB

sampleB$type %>% table
sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$locus %>% length
sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$locus %>% unique %>% length
sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$locus %>% unique
v_sampleB <- sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$locus


v_sampleB_gene_locus <- paste(sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$gene , "~", 
                              sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$locus)

v_sampleB_locus_length <-paste(sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$gene , "~", 
                               sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$locus , "~", 
                               sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$length)



### sampleC

sampleC$type %>% table
sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$locus %>% length
sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$locus %>% unique %>% length
sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$locus %>% unique
v_sampleC <- sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$locus


v_sampleC_gene_locus <- paste(sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$gene , "~", 
                              sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$locus)

v_sampleC_locus_length <-paste(sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$gene , "~", 
                               sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$locus , "~", 
                               sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$length)


### draw venndiagram
titleN =  sprintf("%s %s" , typeN ,  project_name)
filev1=sprintf("venn_%s_%s_%s.png", typeN, project_name, Sys.Date())
png(filename = filev1, width = 300, height = 350)
plot(euler(list(A = unique(v_sampleA_locus_length), B = unique(v_sampleB_locus_length), C = unique(v_sampleC_locus_length))), 
     labels = identical(legend, FALSE), quantities = TRUE, 
     main =  titleN,
     fills = list(fill = c("Pink", "steelblue4", "Green"), alpha = 0.5))
dev.off()
euler(list(A = unique(v_sampleA_locus_length), B = unique(v_sampleB_locus_length), C = unique(v_sampleB_locus_length)))

# rm(sampleA)
# rm(sampleB)
# rm(sampleC)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#


for(i in v_sampleA_locus_length){cat(i);cat(";", sep = "\n")}
for(i in v_sampleB_locus_length){cat(i);cat(";", sep = "\n")}
for(i in v_sampleC_locus_length){cat(i);cat(";", sep = "\n")}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
### SNV
#-------------------------------------------------------------------#
### locus ####
######## Change variable!!########
library(eulerr)

# "CNV", "INDEL", "SNV"
typeN = "SNV"
typeNo = sprintf("%s,NOCALL" , typeN)

#-------------------------------------------------------------------#
### sampleA
sampleA$type %>% table
sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$locus %>% length
sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$locus %>% unique %>% length
sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$locus %>% unique
v_sampleA <- sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$locus 


v_sampleA_gene_locus <- paste(sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$gene, "~", sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$locus)


v_sampleA_locus_length <- paste(sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$gene, "~", 
                                sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$locus, "~", 
                                sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$length)

### sampleB

sampleB$type %>% table
sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$locus %>% length
sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$locus %>% unique %>% length
sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$locus %>% unique
v_sampleB <- sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$locus


v_sampleB_gene_locus <- paste(sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$gene , "~", 
                              sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$locus)

v_sampleB_locus_length <-paste(sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$gene , "~", 
                               sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$locus , "~", 
                               sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$length)



### sampleC

sampleC$type %>% table
sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$locus %>% length
sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$locus %>% unique %>% length
sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$locus %>% unique
v_sampleC <- sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$locus


v_sampleC_gene_locus <- paste(sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$gene , "~", 
                              sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$locus)

v_sampleC_locus_length <-paste(sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$gene , "~", 
                               sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$locus , "~", 
                               sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$length)


### draw venndiagram
titleN =  sprintf("%s %s" , typeN ,  project_name)
filev1=sprintf("venn_%s_%s_%s.png", typeN, project_name, Sys.Date())
png(filename = filev1, width = 300, height = 350)
plot(euler(list(A = unique(v_sampleA_locus_length), B = unique(v_sampleB_locus_length), C = unique(v_sampleC_locus_length))), 
     labels = identical(legend, FALSE), quantities = TRUE, 
     main =  titleN,
     fills = list(fill = c("Pink", "steelblue4", "Green"), alpha = 0.5))
dev.off()
euler(list(A = unique(v_sampleA_locus_length), B = unique(v_sampleB_locus_length), C = unique(v_sampleB_locus_length)))

# rm(sampleA)
# rm(sampleB)
# rm(sampleC)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#


for(i in v_sampleA_locus_length){cat(i);cat(";", sep = "\n")}
for(i in v_sampleB_locus_length){cat(i);cat(";", sep = "\n")}
for(i in v_sampleC_locus_length){cat(i);cat(";", sep = "\n")}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#



# gene only
setdiff(unique(v_sampleA),unique(v_sampleB)) %>% sort()
intersect(unique(v_sampleA),unique(v_sampleB)) %>% sort()
setdiff(unique(v_sampleB), unique(v_sampleA)) %>% sort()

# gene_locus
setdiff(unique(v_sampleA_gene_locus),unique(v_sampleB_gene_locus)) %>% sort()
intersect(unique(v_sampleA_gene_locus),unique(v_sampleB_gene_locus)) %>% sort()
setdiff(unique(v_sampleB_gene_locus), unique(v_sampleA_gene_locus)) %>% sort()


# gene_locus_length
setdiff(unique(v_sampleA_locus_length),unique(v_sampleB_locus_length)) %>% sort()
intersect(unique(v_sampleA_locus_length),unique(v_sampleB_locus_length)) %>% sort()
setdiff(unique(v_sampleB_locus_length), unique(v_sampleA_locus_length)) %>% sort()


### write csv
calN = "gene_locus_length"

Aonly <- setdiff(unique(v_sampleA_locus_length),unique(v_sampleB_locus_length)) %>% sort()
Aonly = data.frame(Aonly)
Aonly["sample"] = "Aonly"
colnames(Aonly)[1] = calN


AinterB <- intersect(unique(v_sampleA_locus_length),unique(v_sampleB_locus_length)) %>% sort()
AinterB = data.frame(AinterB)
AinterB["sample"] = "AinterB"
colnames(AinterB)[1] = calN

Bonly <- setdiff(unique(v_sampleB_locus_length), unique(v_sampleA_locus_length)) %>% sort()
Bonly = data.frame(Bonly)
Bonly["sample"] = "Bonly"
colnames(Bonly)[1] = calN


venndf_gene_locus_length = rbind(Aonly, AinterB, Bonly)
filedf=sprintf("vennlist_%s_%s_%s_%s.csv", typeN, project_name, colN, Sys.Date())

# sd1 <- setdiff(unique(v_sampleB), unique(v_sampleA)) %>% sort()
# write.csv(as.data.frame(sd1), "sample4_locus_indel.csv")


