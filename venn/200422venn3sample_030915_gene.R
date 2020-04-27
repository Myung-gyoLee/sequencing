#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
### CNV
#-------------------------------------------------------------------#
### gene ####
######## Change variable!!########
library(eulerr)
project_name = "gene04vs10vs16"
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
sampleA$type %>% table
sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$gene %>% length
sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$gene %>% unique %>% length
sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$gene %>% unique
v_sampleA <- sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$gene 


sampleB$type %>% table
sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$gene %>% length
sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$gene %>% unique %>% length
sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$gene %>% unique
v_sampleB <- sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$gene 


sampleC$type %>% table
sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$gene %>% length
sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$gene %>% unique %>% length
sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$gene %>% unique
v_sampleC <- sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$gene 



titleN =  sprintf("%s %s" , typeN ,  project_name)
filev1=sprintf("venn_%s_%s_%s.png", typeN, project_name, Sys.Date())
png(filename = filev1, width = 300, height = 350)
plot(euler(list(A = unique(v_sampleA), B = unique(v_sampleB), C = unique(v_sampleC))), 
     labels = identical(legend, FALSE), quantities = TRUE, 
     main =  titleN,
     fills = list(fill = c("Pink", "steelblue4", "Green"), alpha = 0.5))
dev.off()
euler(list(A = unique(v_sampleA), B = unique(v_sampleB), C = unique(v_sampleB)))


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#


for(i in v_sampleA){cat(i, sep = "\n")}
for(i in v_sampleB){cat(i, sep = "\n")}
for(i in v_sampleC){cat(i, sep = "\n")}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
### INDEL
#-------------------------------------------------------------------#
### gene ####
######## Change variable!!########
library(eulerr)

# "CNV", "INDEL", "SNV"
typeN = "INDEL"
typeNo = sprintf("%s,NOCALL" , typeN)

#-------------------------------------------------------------------#
sampleA$type %>% table
sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$gene %>% length
sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$gene %>% unique %>% length
sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$gene %>% unique
v_sampleA <- sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$gene 


sampleB$type %>% table
sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$gene %>% length
sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$gene %>% unique %>% length
sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$gene %>% unique
v_sampleB <- sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$gene 


sampleC$type %>% table
sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$gene %>% length
sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$gene %>% unique %>% length
sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$gene %>% unique
v_sampleC <- sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$gene 

Reduce(intersect, list(a,b,c))
intersect(intersect(unique(v_sampleA),unique(v_sampleB)), unique(v_sampleC)) %>% sort()
setdiff(unique(v_sampleA),unique(v_sampleB)) %>% sort()
setdiff(unique(v_sampleB), unique(v_sampleA)) %>% sort()

titleN =  sprintf("%s %s" , typeN ,  project_name)
filev1=sprintf("venn_%s_%s_%s.png", typeN, project_name, Sys.Date())
png(filename = filev1, width = 300, height = 350)
plot(euler(list(A = unique(v_sampleA), B = unique(v_sampleB), C = unique(v_sampleC))), 
     labels = identical(legend, FALSE), quantities = TRUE, 
     main =  titleN,
     fills = list(fill = c("Pink", "steelblue4", "Green"), alpha = 0.5))
dev.off()
euler(list(A = unique(v_sampleA), B = unique(v_sampleB), C = unique(v_sampleB)))


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#


for(i in v_sampleA){cat(i, sep = "\n")}
for(i in v_sampleB){cat(i, sep = "\n")}
for(i in v_sampleC){cat(i, sep = "\n")}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
### SNV
#-------------------------------------------------------------------#
### gene ####
######## Change variable!!########
library(eulerr)

# "CNV", "INDEL", "SNV"
typeN = "SNV"
typeNo = sprintf("%s,NOCALL" , typeN)

#-------------------------------------------------------------------#
sampleA$type %>% table
sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$gene %>% length
sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$gene %>% unique %>% length
sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$gene %>% unique
v_sampleA <- sampleA[sampleA$type == typeN | sampleA$type == typeNo,]$gene 


sampleB$type %>% table
sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$gene %>% length
sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$gene %>% unique %>% length
sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$gene %>% unique
v_sampleB <- sampleB[sampleB$type == typeN | sampleB$type == typeNo,]$gene 


sampleC$type %>% table
sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$gene %>% length
sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$gene %>% unique %>% length
sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$gene %>% unique
v_sampleC <- sampleC[sampleC$type == typeN | sampleC$type == typeNo,]$gene 

Reduce(intersect, list(a,b,c))
intersect(intersect(unique(v_sampleA),unique(v_sampleB)), unique(v_sampleC)) %>% sort()
setdiff(unique(v_sampleA),unique(v_sampleB)) %>% sort()
setdiff(unique(v_sampleB), unique(v_sampleA)) %>% sort()

titleN =  sprintf("%s %s" , typeN ,  project_name)
filev1=sprintf("venn_%s_%s_%s.png", typeN, project_name, Sys.Date())
png(filename = filev1, width = 300, height = 350)
plot(euler(list(A = unique(v_sampleA), B = unique(v_sampleB), C = unique(v_sampleC))), 
     labels = identical(legend, FALSE), quantities = TRUE, 
     main =  titleN,
     fills = list(fill = c("Pink", "steelblue4", "Green"), alpha = 0.5))
dev.off()
euler(list(A = unique(v_sampleA), B = unique(v_sampleB), C = unique(v_sampleB)))

#rm(sampleA)
#rm(sampleB)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#


for(i in v_sampleA){cat(i, sep = "\n")}
for(i in v_sampleB){cat(i, sep = "\n")}
for(i in v_sampleC){cat(i, sep = "\n")}


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

A&B&C
A&B_A&B&C
B&C_A&B&C
C&A_A&B&C

A_A&B&C_





intersect(intersect(unique(v_sampleA),unique(v_sampleB)), unique(v_sampleC)) %>% sort()
setdiff(unique(v_sampleA),unique(v_sampleB)) %>% sort()
setdiff(unique(v_sampleB), unique(v_sampleA)) %>% sort()