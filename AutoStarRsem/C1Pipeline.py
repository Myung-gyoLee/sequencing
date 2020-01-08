input1 = input('Input rawfastq directory(ex:SMC009_30ea);project >')
 
print(input1)    

#================================================================================================================
### Commands of Make New Directory ###
#================================================================================================================


def makeNewDirectory(project_name):
	import os, glob
	if os.path.isdir("/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/") == 1:
		maindir="/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/"
	elif os.path.isdir("/media/cytogenbi1/D2C67EE7C67ECAED/") == 1:
		maindir="/media/cytogenbi1/D2C67EE7C67ECAED/"
	else : 
		print("Error: Change maindir path!!!!")	
	newTrimDir="%s02trim/%s"%(maindir,project_name)
	newQCdir="%s01QC/%s"%(maindir,project_name)
	newbeforeQCdir="%s01QC/%s/beforeQC"%(maindir,project_name)
	newafterQCdir="%s01QC/%s/afterQC"%(maindir,project_name)
	newStarDir="%s03star/%s"%(maindir,project_name)
	newQualiDir="%s03qualimap/%s"%(maindir,project_name)
	newRsemDir="%s04rsem/%s"%(maindir,project_name)
	newLogDir="%s04rsem/log_%s"%(maindir,project_name)
	if os.path.isdir(newTrimDir) == 0:
		os.system("mkdir %s"%(newTrimDir))
	if os.path.isdir(newQCdir) == 0:
		os.system("mkdir %s"%(newQCdir))
	if os.path.isdir(newbeforeQCdir) == 0:
		os.system("mkdir %s"%(newbeforeQCdir))
	if os.path.isdir(newafterQCdir) == 0:
		os.system("mkdir %s"%(newafterQCdir))
	if os.path.isdir(newStarDir) == 0:
		os.system("mkdir %s"%(newStarDir))
	if os.path.isdir(newQualiDir) == 0:
		os.system("mkdir %s"%(newQualiDir))
	if os.path.isdir(newRsemDir) == 0:
		os.system("mkdir %s"%(newRsemDir))
	if os.path.isdir(newLogDir) == 0:
		os.system("mkdir %s"%(newLogDir))
	return print(">>> Make new Directories Done!!!\n")


#================================================================================================================
### Commands of Trim_galore ###
#================================================================================================================
def makeTrimming(project_name,rawFa):
	if os.path.isdir("/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/") == 1:
		maindir="/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/"
	elif os.path.isdir("/media/cytogenbi1/D2C67EE7C67ECAED/") == 1:
		maindir="/media/cytogenbi1/D2C67EE7C67ECAED/"
	else : 
		print("Error: Change maindir path!!!!")	
	newTrimDir="%s02trim/%s"%(maindir,project_name)
	oufname="%s00script/run/01runtrim%s.sh"%(maindir,project_name)
	ouf=open(oufname,"w")
	lineNum=0
	for infile in rawFa:
		lineNum+=1
		#if lineNum == 1: continue
		#print("#",lineNum)	
		if "fastq.gz" in infile:
			fastq2= infile.split("_1.fastq.gz")[0]+"_2.fastq.gz"
		else :
			fastq2= infile.split("_1.fastq")[0]+"_2.fastq"
		cm_trim = "trim_galore %s %s --paired --phred33  -o %s\n\n\n"%(infile,fastq2,newTrimDir)
		#print(cm_trim)
		ouf.write("#%s\n"%lineNum)
		ouf.write(cm_trim)
	ouf.close()
	os.system("chmod 755 %s"%(oufname))
	return print("trim_galore:\n%s\n\n\n"%(oufname))	

#================================================================================================================
### Commands of fastQC ###
#================================================================================================================

def makebeforefastQC(project_name,rawFa,trimFa):
	import os, glob
	if os.path.isdir("/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/") == 1:
		maindir="/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/"
	elif os.path.isdir("/media/cytogenbi1/D2C67EE7C67ECAED/") == 1:
		maindir="/media/cytogenbi1/D2C67EE7C67ECAED/"
	else : 
		print("Error: Change maindir path!!!!")
	newQCdir="%s01QC/%s"%(maindir,project_name)
	newbeforeQCdir="%s01QC/%s/beforeQC"%(maindir,project_name)
	newafterQCdir="%s01QC/%s/afterQC"%(maindir,project_name)
	oufname="%s00script/run/02runbeforeQC%s.sh"%(maindir,project_name)
	ouf=open(oufname,"w")		 
	lineNum=0
	for tfile in rawFa:
		lineNum+=1
		#if lineNum == 1: continue
		#print("#",lineNum)	
		if "fastq.gz" in tfile:
			fastq2= tfile.split("_1.fastq.gz")[0]+"_2.fastq.gz"
		else :
			fastq2= tfile.split("_1.fastq")[0]+"_2.fastq"
		cmtr="fastqc -o %s %s %s\n\n"%(newbeforeQCdir, tfile, fastq2)
		#print(cmtr)
		ouf.write(cmtr)
	ouf.close()
	os.system("chmod 755 %s"%(oufname))
	return print("beforeQCfile:\n%s\n\n\n"%(oufname))

def makeafterfastQC(project_name,rawFa,trimFa):
	import os, glob
	if os.path.isdir("/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/") == 1:
		maindir="/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/"
	elif os.path.isdir("/media/cytogenbi1/D2C67EE7C67ECAED/") == 1:
		maindir="/media/cytogenbi1/D2C67EE7C67ECAED/"
	else : 
		print("Error: Change maindir path!!!!")
	newQCdir="%s01QC/%s"%(maindir,project_name)
	newbeforeQCdir="%s01QC/%s/beforeQC"%(maindir,project_name)
	newafterQCdir="%s01QC/%s/afterQC"%(maindir,project_name)
	oufname="%s00script/run/02runafterQC%s.sh"%(maindir,project_name)
	ouf=open(oufname,"w")
	lineNum=0
	for tfile in trimFa:
		lineNum+=1
		#if lineNum == 1: continue
		#print("#",lineNum)	
		fastq2 = tfile.replace("_1_val_1.fq.gz","_2_val_2.fq.gz")
		cmtr="fastqc -o %s %s %s\n\n"%(newafterQCdir, tfile, fastq2)
		#print(cmtr)
		ouf.write(cmtr)
	ouf.close()
	os.system("chmod 755 %s"%(oufname))
	return print("afterQCfile:\n%s\n\n\n"%(oufname))



#================================================================================================================
### Commands of Mapping with Star Qualimap Rsem ###
#================================================================================================================


def STAR2RSEM(project_name,trimFa):
	import os, glob
	if os.path.isdir("/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/") == 1:
		maindir="/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/"
	elif os.path.isdir("/media/cytogenbi1/D2C67EE7C67ECAED/") == 1:
		maindir="/media/cytogenbi1/D2C67EE7C67ECAED/"
	else : 
		print("Error: Change maindir path!!!!")
	newTrimDir="%s02trim/%s"%(maindir,project_name)
	newStarDir="%s03star/%s"%(maindir,project_name)
	newQualiDir="%s03qualimap/%s"%(maindir,project_name)
	newRsemDir="%s04rsem/%s"%(maindir,project_name)
	genomeDir="/media/cytogenbi2/6eaf3ba8-a866-4e8a-97ef-23c61f7da612/00ref/STAR38" 
	sjdbGTFfile="/media/cytogenbi2/6eaf3ba8-a866-4e8a-97ef-23c61f7da612/00ref/ensembl/GRCh38/Homo_sapiens.GRCh38.97.gtf"
	rsemrefDir="/media/cytogenbi2/6eaf3ba8-a866-4e8a-97ef-23c61f7da612/00ref/STAR38/RSEM"
	oufname="%s00script/run/03run%sstar.sh"%(maindir,project_name)
	ouf=open(oufname,"w")
	fline=[]
	lineNum=0
	for infile in trimFa:
		lineNum+=1
		#if lineNum == 1: continue
		#print("#",lineNum)	
		trfq2 = infile.replace("_1_val_1.fq.gz","_2_val_2.fq.gz")
		sample = infile.split("/")[-1].split("_1_val_1.fq.gz")[0]
		star_out_file="%s/%s"%(newStarDir, sample)
		quali_out_dir="%s/%s"%(newQualiDir, sample)
		rsem_out_file="%s/%s"%(newRsemDir, sample)
		cm_trimmed_star="STAR --runThreadN 2 --genomeDir %s --sjdbGTFfile %s --sjdbOverhang 150 --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.02 --alignIntronMin 20 --outFileNamePrefix %s --alignIntronMax 1000000 --alignMatesGapMax 1000000 --readFilesIn %s,%s --readFilesCommand zcat --quantMode TranscriptomeSAM\n\n\n"%(str(genomeDir), str(sjdbGTFfile), star_out_file, infile, trfq2)
		cm_qualimap="qualimap rnaseq -outdir %s -a proportional -bam %sAligned.out.sam -gtf %s --java-mem-size=8G \n\n\n"%(quali_out_dir, star_out_file,sjdbGTFfile)
		cm_rsem="rsem-calculate-expression --bam --no-bam-output --single-cell-prior -p 1 %sAligned.toTranscriptome.out.bam %s %sscQuant> %sscrsem.log\n\n\n"%(star_out_file, str(rsemrefDir), rsem_out_file, rsem_out_file)
		#print(cm_trimmed_star)
		#print(cm_qualimap)
		#print(cm_rsem)
		if lineNum ==1 : fline.append(cm_trimmed_star)
		ouf.write(cm_trimmed_star)
		ouf.write("echo %s %s star alignment done!\n"%(project_name, sample))
		ouf.write(cm_qualimap)
		ouf.write("echo %s %s Qualimap done!\n"%(project_name, sample))
		ouf.write(cm_rsem)
		ouf.write("echo %s %s rsem done!\n"%(project_name, sample))
	ouf.close()
	os.system("chmod 755 %s"%(oufname))
	return print("STAR_Qualimap_RSEM:\n%s\n\n\n"%(oufname))


#================================================================================================================
### Find missing RSEM results
#================================================================================================================

def compRaw(project_name,rawFa):
	import glob
	rawExtension="_1.fastq.gz"
	outExtension="scQuant.genes.results"
	outPath="%s04rsem/%s/"%(maindir,project_name)
	outFiles=glob.glob("%s*%s"%(outPath,outExtension))
	outsample=[a.split("/")[-1].split(outExtension)[0] for a in outFiles]
	rawsample=[r.split("/")[-1].split(rawExtension)[0] for r in rawFa]
	outset=set(outsample)
	missing=[x for x in rawsample if x not in outset]
	misPath=[line for line in missing if line in outFiles]
	print(",".join(missing))
	print("\n".join(misPath))
	return "\n".join(misPath)

#================================================================================================================
### Commands
#================================================================================================================
import os, glob, time
project=input1.split(";")[1]
rawdir=input1.split(";")[0] # rawdir example="SMC009_30ea"
start=time.strftime("%y%m%d")
project_name="%s_%s"%(project,start)

maindir="/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/"

rawFa=glob.glob("/media/cytogenbi2/6eaf3ba8-a866-4e8a-97ef-23c61f7da612/01raw/%s/*_1.fastq.gz"%(rawdir))
trimFa=[]
for rline in rawFa:
	trimdir="%s02trim/%s/"%(maindir,project_name)
	trimfile=rline.split("/")[-1].replace("_1.fastq.gz","_1_val_1.fq.gz")
	tline="%s%s"%(trimdir,trimfile)
	trimFa.append(tline)

makeNewDirectory(project_name)
makeTrimming(project_name,rawFa)
makebeforefastQC(project_name,rawFa,trimFa)
makeafterfastQC(project_name,rawFa,trimFa)
STAR2RSEM(project_name,trimFa)
#compRaw(project_name,rawFa)


