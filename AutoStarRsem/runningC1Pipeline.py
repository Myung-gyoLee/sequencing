python ~/00script/BI_left/python/AutoStarRsem/C1Pipeline.py 
PBMC0102_42ea;08PBMC0102_191213

trim_galore /media/cytogenbi2/6eaf3ba8-a866-4e8a-97ef-23c61f7da612/01raw/PBMC0102_42ea/PBMC_2-59_1.fastq.gz /media/cytogenbi2/6eaf3ba8-a866-4e8a-97ef-23c61f7da612/01raw/PBMC0102_42ea/PBMC_2-59_2.fastq.gz --paired --phred33  -o /media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/02trim/08PBMC0102_191213_191213

fastqc -o /media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/01QC/08PBMC0102_191213_191213/beforeQC /media/cytogenbi2/6eaf3ba8-a866-4e8a-97ef-23c61f7da612/01raw/PBMC0102_42ea/PBMC_2-59_1.fastq.gz /media/cytogenbi2/6eaf3ba8-a866-4e8a-97ef-23c61f7da612/01raw/PBMC0102_42ea/PBMC_2-59_2.fastq.gz


STAR --runThreadN 2 --genomeDir /media/cytogenbi2/6eaf3ba8-a866-4e8a-97ef-23c61f7da612/00ref/STAR38 --sjdbGTFfile /media/cytogenbi2/6eaf3ba8-a866-4e8a-97ef-23c61f7da612/00ref/ensembl/GRCh38/Homo_sapiens.GRCh38.97.gtf --sjdbOverhang 150 --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.02 --alignIntronMin 20 --outFileNamePrefix /media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/03star/08PBMC0102_191213_191213/PBMC_2-59 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --readFilesIn /media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/02trim/08PBMC0102_191213_191213/PBMC_2-59_1_val_1.fq.gz,/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/02trim/08PBMC0102_191213_191213/PBMC_2-59_2_val_2.fq.gz --readFilesCommand zcat --quantMode TranscriptomeSAM

#========================================================================================================
# conda activate scrna

makebeforefastQC(project_name,rawFa,trimFa)
makeafterfastQC(project_name,rawFa,trimFa)


import os, glob, time
start=time.strftime("%y%m%d")
project_name="07SMC008_%s"%(start)

project_name="07SMC009_191209"
maindir="/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/"
maindir="/media/cytogenbi1/D2C67EE7C67ECAED/"

rawFa=glob.glob("/media/cytogenbi2/6eaf3ba8-a866-4e8a-97ef-23c61f7da612/01raw/SMC009_30ea/*_1.fastq.gz")
trimFa=[]
for rline in rawFa:
	trimdir="%s02trim/%s/"%(maindir,project_name)
	trimfile=rline.split("/")[-1].replace("_1.fastq.gz","_1_val_1.fq.gz")
	tline="%s%s"(trimdir,trimfile)
	trimFa.append(tline)
trimFa=[a.split("/")[-1].split(outExtension)[0] for a in rawFa]
trimFa=glob.glob("%s02trim/%s/*_1_val_1.fq.gz"%(maindir,project_name))

############
/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/02trim/07SMC009_191209/SMC_009-1_1_val_1.fq.gz
SMC_009-1_1.fastq.gz

#trimmed

trim1=glob.glob("%s/03trim/01CMC1113/*_1_val_1.fq.gz"%maindir)
print(len(trim1))
trim2=glob.glob("%s/03trim/02pbmc/*_1_val_1.fq.gz"%maindir)
print(len(trim2))
trim3=glob.glob("%s/03trim/03H69plate/*_1_val_1.fq.gz"%maindir)
print(len(trim3))
trim4=glob.glob("%s/03trim/04CMC1819/*_1_val_1.fq.gz"%maindir)
print(len(trim4))
trim5=glob.glob("%s/03trim/03H69C1/*_1_val_1.fq.gz"%maindir)
print(len(trim5))

trimFa=trim1+trim2+trim3+trim4+trim5

maindir="/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/"
trimFa=glob.glob("%s02trim/07SMC008_191212/*_1_val_1.fq.gz"%maindir)
print("trimFa length = ", len(trimFa))
print("sum of trim1...trimN = ", len(trim1)+len(trim2)+len(trim3)+len(trim4)+len(trim5))

#================================================================================================================
## Merge QC results
#================================================================================================================
# for file in `ls *.zip`; do unzip "${file}" -d "${file:0:-4}"; done


qcf =glob.glob("/media/desktop-bi-16/D2C67EE7C67ECAED/BI/03qc/*/*/summary.txt")
qcline=[]
for files in qcf:
	qcff=open(files)
	qcline.append(qcff.readlines())

failD={}
qctotal =[]
for aa in qcline:
	for bb in aa:
		qctotal.append(bb)
			if "fail" in bb:
			
		break
ouf=open("/media/desktop-bi-16/D2C67EE7C67ECAED/BI/03qc/resQC190917.txt","w")
for aaa in qctotal:
	if "FAIL" in aaa:
		wrine = "%s:%s\n"%(aaa.split("\t")[2].strip("\n"),aaa.split("\t")[1])
		print aaa.split("\t")[2].strip("\n"),":",aaa.split("\t")[1]
		ouf.write(wrine)
ouf.close()
		failD[aaa.split("\t")[2]]=aaa.split("\t")[1]
	



