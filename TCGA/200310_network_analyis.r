setwd("H:TCGA")


## join 1) TCGA tumor vs TCGA normal
library(DESeq2)

# reread Expression Data and Meta data
BRCA_TN1=read.csv("./02Breast/021merge/ExBRCA_Tumor_normal.csv", header = T, row.names = 1)
metaTtum_nor1=read.csv("./02Breast/021merge/metaBRCA_Tumor_normal.csv", header = T, row.names = 1)

dds = DESeqDataSetFromMatrix(countData = BRCA_TN1[rowSums(BRCA_TN1)>0,],
                             colData = metaTtum_nor1[,c(2,3,5)],
                             design = ~ Condition)

dds
#View(metaTtum_nor1)
#rownames(metaTtum_nor) = seq(1:length(rownames(metaTtum_nor)))
#rownames(BRCA_TN) = seq(1:length(rownames(BRCA_TN)))

#------------------------------------------------------------------#

results <- results(DESeq(dds))

## Filter the data in "results" with the fold change log2FC>1, adjusted p-values<0.01. Store the data in "filter" object

res <- na.exclude(as.data.frame(results))
filter <- res[(abs(res$log2FoldChange)>1 & res$padj<0.01),]

write.table(filter,"filteredbreast200310.txt", quote = F, sep = "\t", col.names = NA)

BiocManager::install("EBSeq")
library(EBSeq)
NormData <- GetNormalizedMat(BRCA_TN1, MedianNorm(BRCA_TN1))
NormData.log <- log2(NormData+1)
Norm.interest <- NormData.log[rownames(filter),]

BiocManager::install("psych")
library("psych")

x <- t(Norm.interest)
Norm.interest.corr <- corr.test(x=t(Norm.interest), y=NULL, method = "pearson", ci = FALSE)

#------------------------------------------------------------------#
#------------------------------------------------------------------#
## Setting R session
setwd("H:TCGA")
getwd()
dir.create("./network")
setwd("./network")
# Load WGCNA package
library(WGCNA)
options(stringsAsFactors = FALSE)

## Simulation of expression and trait data
### Building the module structure

# Here are input parameters of the simulation model
# number of samples 
no.obs=50
# specify the true measures of eigengene significance
# recall that ESturquoise=cor(y,MEturquoise)
ESturquoise=0; ESbrown= -.6;
ESgreen=.6;ESyellow=0
# Note that we dont specify the eigengene significance of the blue module
# since it is highly correlated with the turquoise module.
ESvector=c(ESturquoise,ESbrown, ESgreen, ESyellow)
# number of genes
nGenes1=3000
# proportion of genes in the turquoise, blue, brown, green, and yellow module
simulateProportions1=c(0.2,0.15,0.08,0.06,0.04)
# Note that the proportions dont add up to 1. The remaining genes will be colored grey,
# ie the grey genes are non-module genes.
# set the seed of the random number generator. As a homework exercise change this seed.
set.seed(1)
#Step 1: simulate a module eigengene network.
# Training Data Set I
MEgreen=rnorm(no.obs)
scaledy=MEgreen*ESgreen+sqrt(1-ESgreen^2)*rnorm(no.obs)
y=ifelse(scaledy>median(scaledy),2,1)
MEturquoise= ESturquoise*scaledy+sqrt(1-ESturquoise^2)*rnorm(no.obs)
# we simulate a strong dependence between MEblue and MEturquoise
MEblue=.6*MEturquoise+sqrt(1-.6^2)*rnorm(no.obs)
MEbrown = ESbrown*scaledy+sqrt(1-ESbrown^2)*rnorm(no.obs)
MEyellow=ESyellow*scaledy+sqrt(1-.6^2)*rnorm(no.obs)
ModuleEigenegeneNetwork1=data.frame(y,MEturquoise, MEblue, MEbrown, MEgreen, MEyellow)

# The variable ModuleEigengeneNetwork1 contains the seed eigengens and a simulated clinical trait y. The eigengene network can be simply inspected  by

signif(cor(ModuleEigenegeneNetwork1, use="p"),2)

dat1=simulateDatExpr5Modules(MEturquoise = ModuleEigenegeneNetwork1$MEturquoise,
                             MEblue = ModuleEigenegeneNetwork1$MEblue,
                             MEbrown = ModuleEigenegeneNetwork1$MEbrown,
                             MEyellow=ModuleEigenegeneNetwork1$MEyellow,
                             MEgreen = ModuleEigenegeneNetwork1$MEgreen,
                             nGenes=nGenes1,
                             simulateProportions = simulateProportions1
                             )
names(dat1)

# attach the data "into the main search path" so we can use the component names directly, without referring to dat1

datExpr = dat1$datExpr
truemodule = dat1$truemodule
datME = dat1$datME
attach(ModuleEigenegeneNetwork1)

table(truemodule)
dim(datExpr)
table(truemodule)

datExpr=data.frame(datExpr)
ArrayName=paste("Sample",1:dim(datExpr)[[1]], sep="")
GeneName=paste("Gene",1:dim(datExpr)[[2]], sep="")
dimnames(datExpr)[[1]]=ArrayName
dimnames(datExpr)[[2]]=GeneName
rm(dat1) 
collectGarbage()
save.image("Simulated-dataSimulation.RData")


## Loading of expression and trait data
datGeneSummary=read.csv("./SimulatedData/GeneSummaryTutorial.csv")
datTraits=read.csv("./SimulatedData/TraitsTutorial.csv")
datMicroarrays=read.csv("./SimulatedData/MicroarrayDataTutorial.csv")

# This vector contains the microarray sample names
ArrayName= names(data.frame(datMicroarrays[,-1]))
# This vector contains the gene names
GeneName= datMicroarrays$GeneName
# We transpose the data so that the rows correspond to samples and the columns correspond to genes
# Since the first column contains the gene names, we exclude it.
datExpr=data.frame(t(datMicroarrays[,-1]))
names(datExpr)=datMicroarrays[,1]
dimnames(datExpr)[[1]]=names(data.frame(datMicroarrays[,-1]))
# true module color:
truemodule=datGeneSummary$truemodule
rm(datMicroarrays)
collectGarbage()
# First, make sure that the array names in the file datTraits line up with those in the microarray data
table(dimnames(datExpr)[[1]]==datTraits$ArrayName)
# Next, define the microarray sample trait
y=datTraits$y


## Identification of outlying samples
meanExpressionByArray=apply(datExpr, 1, mean, na.rm=T)
NumberMissingByArray=apply(is.na(data.frame(datExpr)),1,sum)
# simple way to examine the mean expression per array is to use
sizeGrWindow(9, 5)
barplot(meanExpressionByArray, 
        xlab = "Sample", ylab = "Mean expression",
        main = "Mean expression across samples",
        names.arg = c(1:50), cex.names = 0.7)

NumberMissingByArray

# Keep only arrays containint less than 500 missing entries
KeepArray= NumberMissingByArray<500
table(KeepArray)
datExpr=datExpr[KeepArray,]
y=y[KeepArray]
ArrayName[KeepArray]

## Handling missing data and zero variance in probe profiles
NumberMissingByGene = apply(is.na(data.frame(datExpr)),2,sum)
# one could do a barplot(NumberMissingByGene), but the barplot is empty in this case.
# It may be better to look at the numbers of missing samples using the summary method:
summary(NumberMissingByGene)
# Calculate the variances of the probes and the number of present entries
variancedatExpr=as.vector(apply(as.matrix(datExpr),2,var,na.rm=T))
no.presentdatExpr=as.vector(apply(!is.na(as.matrix(datExpr)),2,sum))
# Another way of summarizing the number of present entries
table(no.presentdatExpr)
# Keep only genes whose varianve is non-zero and have at least 4 present entries
KeepGenes= variancedatExpr>0 & no.presentdatExpr>=4
table(KeepGenes)
datExpr=datExpr[,KeepGenes]
GeneName=GeneName[KeepGenes]

# Rudimentary detection of outlier samples
sizeGrWindow(9,5)
plotClusterTreeSamples(datExpr = datExpr, y=y)

GS1= as.numeric(cor(y,datExpr, use = "p"))
