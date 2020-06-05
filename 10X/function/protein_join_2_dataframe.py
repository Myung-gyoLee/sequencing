# read excel file per sheet
import numpy as np
import pandas as pd
import os
import time


#==========>
def extract_protein(projectN,clusterN,inputgene):
	import re
	import pandas as pd
	oufname = "/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/singlecell_protein/%s/%s.txt"%(projectN,clusterN)
	ouf = open(oufname,"w")
	col = proteinD["Gene"].split("\t")
	colist = ["Gene",col[6],col[12],col[15],col[20],col[23],col[28],col[31],"Cancer prognostic","Cancer prognostic favourable","Cancer prognostic unfavourable","\n"]
	ouf.write("\t".join(colist))
	print(colist)
	noresult=[]
	synonym=[]
	for cgene in inputgene:
		print(cgene)
		try:
			p = re.compile("[^0-9.:]")
			val_P=proteinD[cgene].split("\t")
			indices = ["%s_%s"%(i,s.split("(")[0].replace('"','').rstrip()) for i, s in enumerate(val_P[55:72]) if 'favourable' in s]
			tissuespec = "".join(p.findall(val_P[15])).replace("; ",",")
			cancerfpkm = "".join(p.findall(val_P[23])).replace("; ",",")
			bloodspec = "".join(p.findall(val_P[31])).replace("; ",",")
			bloodlineage = "".join(p.findall(val_P[35])).replace("; ",",")
			if len(indices) >= 1:
				prol=[]
				unprol=[]
				for a in indices:
					if "un" in a:
						unprol.append(a)
					else:
						prol.append(a)
				prolval=["%s(%s)"%(prognosticD[str(a.split("_")[0])],a.split("_")[1]) for a in prol]
				unprolval=["%s(%s)"%(prognosticD[str(a.split("_")[0])],a.split("_")[1]) for a in unprol]
				print('"%s"'%cgene,val_P[6],val_P[12],tissuespec,val_P[20],cancerfpkm,val_P[28],bloodspec,'"%s"'%(",".join(cancerl)))
				wlinel = ['"%s"'%cgene,val_P[6],val_P[12],tissuespec,val_P[20],cancerfpkm,val_P[28],bloodspec,'"prognostic"','"%s"'%(",".join(prolval)),'"%s"'%(",".join(unprolval)),"\n"]
				writeline = "\t".join(wlinel)
			else : 
				wlinel = ['"%s"'%cgene,val_P[6],val_P[12],tissuespec,val_P[20],cancerfpkm,val_P[28],bloodspec,"not prognostic","","","\n"]
				writeline = "\t".join(wlinel)
		except KeyError:
			if cgene in proteinD.values():
				print(cgene,"gene synonym")
				synonym.append(cgene)
				writeline = "%s\n"%cgene
			else:
				print("no results")
				noresult.append(cgene)
				writeline = "%s\n"%cgene
		ouf.write(writeline)
	ouf.close()

#sheetname = "s025con_Cluster_7"
def protein_annotation(sheetname,writer):
	import pandas as pd
	sheetlname = "%s_%s"%(projectN,sheetname)
	vars()[sheetlname] = smc025conf.parse(sheetname)
	anfilename = "/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/singlecell_protein/%s/%s.txt"%(projectN, sheetlname)
	resname = "%sres"%(sheetlname)
	vars()[resname] = pd.read_csv(anfilename, sep='\t')
	#print(vars()[resname].head())
	#print(vars()[resname].columns)
	newcol = [col.replace(" ","_") for col in vars()[resname].columns]
	vars()[resname].columns = newcol
	#print(vars()[resname].columns)
	#print(vars()[sheetlname].columns)
	joinname = "%sjoin"%(sheetlname)
	vars()[joinname] = pd.merge(vars()[sheetlname],vars()[resname],
	how = 'left',
	left_on = 'gene',
	right_on = 'Gene', 
	sort = False)
	#print(vars()[joinname].columns)
	print(vars()[joinname].head())
	vars()[joinname].to_excel(writer, sheet_name= sheetname)


# extract gene name
with open("/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/singlecell_protein/raw/proteinatlas.tsv") as f:
    proteinl = f.readlines()

proteinD = {}
for pline in proteinl:
    proteinD[pline.split("\t")[0]]="\t".join(pline.split("\t")[1:])


prognosticD = {"0":"Breast cancer","1":"Cervical cancer","2":"Colorectal cancer","3":"Endometrial cancer","4":"Glioma","5":"Head and neck cancer","6":"Liver cancer","7":"Lung cancer","8":"Melanoma","9":"Ovarian cancer","10":"Pancreatic cancer","11":"Prostate cancer","12":"Renal cancer","13":"Stomach cancer","14":"Testis cancer","15":"Thyroid cancer","16":"Urothelial cancer"}

# 1. read raw input per sheet

"""cytogenbi2@cytogenbi2-B365M-DS3H:/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/singlecell_protein$ ls
SMC024_B1Tissue_Conserved_Markers.xlsx
SMC024_B1Tissue_Differential_Expressed_Markers.xlsx
SMC025_B1Tissue_Conserved_Markers.xlsx
SMC025_B1Tissue_Differential_Expressed_Markers.xlsx

"""
infname = '/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/singlecell_protein/raw/SMC025_B1Tissue_Differential_Expressed_Markers.xlsx'
smc025conf = pd.ExcelFile(infname)
projectN = infname.split('/')[-1].split('.')[0]
directory = "/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/singlecell_protein/%s"%(projectN)

import os
try:
	if not os.path.exists(directory):
		os.makedirs(directory)
except OSError:
	print('Error: Creating directory. '+ directory)

s025con_sheetl = smc025conf.sheet_names
parsesheetl = []

for sheetname in s025con_sheetl:
	sheetlname = "%s_%s"%(projectN,sheetname)
	parsesheetl.append(sheetlname)
	vars()[sheetlname] = smc025conf.parse(sheetname)
	inputgene = vars()[sheetlname].gene
	print(inputgene.head())
	extract_protein(projectN, sheetlname,inputgene)


mergefile = '/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/singlecell_protein/merged/%s_protein_atlas.xlsx'%(projectN)
writer = pd.ExcelWriter(mergefile, engine='xlsxwriter')
for sheetname in s025con_sheetl:
	protein_annotation(sheetname, writer)

writer.save()

#smc025conf = pd.ExcelFile('/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/singlecell_protein/SMC025_B1Tissue_Conserved_xcell_200527.xlsx')


for ll in parsesheetl:
	print(ll)

'''
s025con_Cluster_7
s025con_Cluster_12
s025con_Cluster_14
s025con_Cluster_18

'''
# 2. read annotation result file 
''' # test code
n=0
for sheetname in s025con_sheetl:
	n+=1
	sheetlname = "%s_%s"%("s025con",sheetname)
	vars()[sheetlname] = smc025conf.parse(sheetname)
	anfilename = "/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/singlecell_protein/smc025con/%s.txt"%sheetlname
	resname = "%sres"%(sheetlname)
	vars()[resname] = pd.read_csv(anfilename, sep='\t')
	#print(vars()[resname].head())
	#print(vars()[resname].columns)
	newcol = [col.replace(" ","_") for col in vars()[resname].columns]
	vars()[resname].columns = newcol
	#print(vars()[resname].columns)
	#print(vars()[sheetlname].columns)
	joinname = "%sjoin"%(sheetlname)
	vars()[joinname] = pd.merge(vars()[sheetlname],vars()[resname],
	how = 'left',
	left_on = 'gene',
	right_on = 'Gene', 
	sort = False)
	#print(vars()[joinname].columns)
	print(vars()[joinname].head())
	if n > 1:break 
'''
