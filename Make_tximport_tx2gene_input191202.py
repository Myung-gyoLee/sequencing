#PC1

#inf1=open("/media/cytogenbi1/D2C67EE7C67ECAED/BI/02ref/ensembl38.97/Homo_sapiens.GRCh38.97.gtf")
#================================================================#
#PC2
ensembl gtf  38.90 path
/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/00ref/ensembl/GRCh38.90/Homo_sapiens.GRCh38.90.gtf
/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/00ref/ensembl/GRCh38.90/tx2gene3890EXversion.txt
/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/00ref/ensembl/GRCh38.90/tx2gene3890Wthversion.txt
#================================================================#


gtff="/media/cytogenbi2/6eaf3ba8-a866-4e8a-97ef-23c61f7da612/BreastCancer/data/etc/Homo_sapiens.GRCh38.90.gtf"
inf1=open(gtff)
lines=inf1.readlines()

inf1.close()


llist=[]
for num in range(5,len(lines)):
	line=lines[num]
	if "transcript_id" in line:
		chrnum=line.split("\t")[0]
		info=line.split("\t")[8]
		gene_id=info.split('gene_id "')[1].split('";')[0]
		gene_version=info.split('gene_version "')[1].split('";')[0]
		transcript_id=info.split('transcript_id "')[1].split('";')[0]
		transcript_version=info.split('transcript_version "')[1].split('";')[0]
		lline1="%s\t%s.%s\t%s.%s\n"%(chrnum, gene_id, gene_version, transcript_id,transcript_version)
		llist.append(lline1)
	else:
		continue


ouff="%s/tx2geneWthversion.txt"%("/".join(gtff.split("/")[:-1]))
ouf=open(ouff,"w")
ouf.write("".join(llist))


llist2=[]
for num in range(5,len(lines)):
	line=lines[num]
	if "transcript_id" in line:
		chrnum=line.split("\t")[0]
		info=line.split("\t")[8]
		gene_id=info.split('gene_id "')[1].split('";')[0]
		gene_version=info.split('gene_version "')[1].split('";')[0]
		transcript_id=info.split('transcript_id "')[1].split('";')[0]
		transcript_version=info.split('transcript_version "')[1].split('";')[0]
		lline1="%s\t%s\t%s\n"%(chrnum, gene_id, transcript_id)
		llist2.append(lline1)
	else:
		print(line)

ouff2="%s/tx2geneEXversion.txt"%("/".join(gtff.split("/")[:-1]))
ouf2=open(ouff2,"w")
ouf2.write("".join(llist2))




#=======================  test  ======================= #
num=2877401
line=lines[num]
chrnum=line.split("\t")[0]
print(chrnum)
info=line.split("\t")[8]
gene_id=info.split('gene_id "')[1].split('";')[0]
print(gene_id)
gene_version=info.split("gene_version "")[1].split("";")[0]
print(gene_version)
transcript_id=info.split("transcript_id "")[1].split("";")[0]
print(transcript_id)
transcript_version=info.split("transcript_version "")[1].split("";")[0]
print(transcript_version)
lline1="%s\t%s.%s\t%s.%s\n"%(chrnum, gene_id, gene_version, transcript_id,transcript_version)

#======================================================= #

