#===================================================================================
# Kallist Ensg with version annotation
#===================================================================================
inputTPM="/dt2/04rsem/06ENS3803_Merge191120/191120rsemTPMENSG.csv"
def ENSG2GeneName(inputTPM):
	import os
	file1="/media/cytogenbi1/D2C67EE7C67ECAED/BI/02ref/ensembl38.97/GeneE38_97EXversion.txt"
	file2="/media/cytogenbi2/6eaf3ba8-a866-4e8a-97ef-23c61f7da612/00ref/modGTF38/GeneE38_97EXversion.txt"
	if os.path.isfile(file1):
		inf1=open(file1)
		glines=inf1.readlines()
	elif os.path.isfile(file2):
		inf1=open(file2)
		glines=inf1.readlines()
	else :
		print("There's no gtf file!!")
	inf1.close()
	geneD={}
	for line in glines:
		geneD[line.split("\t")[1]]=line.split("\t")[3].strip("\n")
	res1=open("%s"%(str(inputTPM)))
	resline=res1.readlines()
	res1.close()
	outputFile=inputTPM.split(".csv")[0]+"_geneAnno.csv"
	outputFile2=inputTPM.split(".csv")[0]+"_gene_ENSGAnno.csv"
	ouf=open("%s"%(str(outputFile)),"w")
	ouf2=open("%s"%(str(outputFile2)),"w")
	for rline in resline:
		rkey=rline.split(",")[0].replace('"','')
		geneN = geneD.get(rkey)
		print(geneN,"_",rkey)
		geneNa = str(geneN)
		ouf.write('"%s",%s'%(geneNa, ",".join(rline.split(",")[1:])))
		gene_Na = str(geneN)+"_"+rkey
		ouf2.write('"%s",%s'%(gene_Na, ",".join(rline.split(",")[1:])))
	ouf.close()
	ouf2.close()
#===================================================================================

#===================================================================================
# Kallisto Ensg with version annotation
#===================================================================================
inputTPM="/media/cytogenbi1/D2C67EE7C67ECAED/BI/07kallisto/results/191114CMC11_19kallisto_raw.csv"
def ENSG2GeneName(inputTPM):
	import os
	file1="/media/cytogenbi1/D2C67EE7C67ECAED/BI/02ref/ensembl38.97/GeneE38_97Wthversion.txt"
	file2="/media/cytogenbi2/6eaf3ba8-a866-4e8a-97ef-23c61f7da612/00ref/modGTF38/GeneE38_97Wthversion.txt"
	if os.path.isfile(file1):
		inf1=open(file1)
		glines=inf1.readlines()
	elif os.path.isfile(file2):
		inf1=open(file2)
		glines=inf1.readlines()
	else :
		print("There's no gtf file!!")
	inf1.close()
	geneD={}
	for line in glines:
		geneD[line.split("\t")[1]]=line.split("\t")[3].strip("\n")
	res1=open("%s"%(str(inputTPM)))
	resline=res1.readlines()
	res1.close()
	outputFile=inputTPM.split(".csv")[0]+"_geneAnno.csv"
	ouf=open("%s"%(str(outputFile)),"w")
	for rline in resline:
		rkey=rline.split(",")[0].replace('"','')
		geneN = geneD.get(rkey)
		print(geneN,"_",rkey)
		geneNa = str(geneN)+"_"+rkey
		ouf.write('"%s",%s'%(geneNa, ",".join(rline.split(",")[1:])))
	ouf.close()
#===================================================================================
