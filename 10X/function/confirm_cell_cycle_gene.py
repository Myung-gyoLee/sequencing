# https://github.com/theislab/single-cell-tutorial/blob/master/Macosko_cell_cycle_genes.txt
Macosko EZ, Basu A, Satija R, Nemesh J, Shekhar K, Goldman M, Tirosh I,
Bialas AR, Kamitaki N, Martersteck EM et al (2015) Highly parallel genome-wide expression profiling of individual cells using nanoliter droplets. Cell 161: 1202 – 1214

import pandas as pd

#read cell cycle gene list file
cellcyclefile = "/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/cellcycle/Macosko_cell_cycle_genes.txt"

cycle_geneL = pd.read_csv(cellcyclefile, sep='\t')

'''>>> cycle_geneL.head()
     IG1.S         S      G2.M        M    M.G1  Unnamed: 5
0      ACD     ABCC5      ANLN     AHI1   AGFG1         NaN
1    ACYP1    ABHD10     AP3D1  AKIRIN2  AGPAT3         NaN
2  ADAMTS1  ANKRD18A  ARHGAP19  ANKRD40  AKAP13         NaN
3  ANKRD10     ASF1B     ARL4A     ANLN    AMD1         NaN
4    APEX2     ATAD2     ARMC1   ANP32B  ANP32E         NaN
'''

#read input excel file
infname = "/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/singlecell_protein/20200605_SMC027_result/SMC024-025-027_Marker_Annotation.xlsx"

resExcel = pd.ExcelFile(infname)


sheetlist = []
for sheetname in resExcel.sheet_names:
	sheetlname = "%s"%(sheetname) 
	#sheetlname = "%s_%s"%(projectN,sheetname)
	vars()[sheetlname] = resExcel.parse(sheetname)
	sheetlist.append(vars()[sheetlname])
	print(sheetlname) # sheetlname : SMC-024
	print(vars()[sheetlname].columns) #Index(['Cluster', 'Markers', 'Gene',],dtype='object')
	#inputgene = vars()[sheetlname].gene
	inputgene = vars()[sheetlname].Gene
	print(inputgene.head())



      
# intersection
cycle_geneL.columns = ['IG1S', 'S', 'G2M', 'M', 'MG1', 'Unnamed: 5']


sheetlist[0].Gene[sheetlist[0].Gene in cycle_geneL.IG1S]

set(sheetlist[0].Gene)&set(cycle_geneL.IG1S)
set(sheetlist[0].Gene)&set(cycle_geneL.S)
set(sheetlist[0].Gene)&set(cycle_geneL.G2M)
set(sheetlist[0].Gene)&set(cycle_geneL.M)
set(sheetlist[0].Gene)&set(cycle_geneL.MG1)




set(sheetlist[1].Gene)&set(cycle_geneL.IG1S)
set(sheetlist[1].Gene)&set(cycle_geneL.S)
set(sheetlist[1].Gene)&set(cycle_geneL.G2M)
set(sheetlist[1].Gene)&set(cycle_geneL.M)
set(sheetlist[1].Gene)&set(cycle_geneL.MG1)


set(sheetlist[2].Gene)&set(cycle_geneL.IG1S)
set(sheetlist[2].Gene)&set(cycle_geneL.S)
set(sheetlist[2].Gene)&set(cycle_geneL.G2M)
set(sheetlist[2].Gene)&set(cycle_geneL.M)
set(sheetlist[2].Gene)&set(cycle_geneL.MG1)


