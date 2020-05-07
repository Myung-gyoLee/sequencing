#!/bin/bash
import os
import glob
import time
from datetime import date
from itertools import chain

#inf1 = list(filter(os.path.isdir, glob.glob("/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/10X/HN00126855/HN00126855_10X_RawData_Outs/*")))

inf1 = list(filter(os.path.isdir, glob.glob("/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/10X/HN00126855/HN00126855_10X_RawData_Outs/*/H*")))

### except specific id
inf2 = [line for line in inf1 if "SMC_025" not in line]

### extract running id
idlist = list(set("_".join(line.split("/")[-2].split("_")[:2]) for line in inf2))
#idlist = list(set("_".join(line.split("/")[-1].split("_")[:2]) for line in inf2))

ouf = open("/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/10X/00script/run10x%s_%s.sh"%(str("_".join(idlist)),str(date.today())),"w")

#fastqdir = list(chain(*[glob.glob("%s/*"%idline) for idline in inf2]))

for idline in inf2:
	cdline = "cd %s\n\n"%("/".join(idline.split("/")[:-1]))
	print(cdline)
	ouf.write(cdline)
	cellcountline = """cellranger count --id=run_10x_%s \\
--fastqs=%s \\
--sample %s \\
--transcriptome=/home/cytogenbi2/singlecell/refdata-cellranger-GRCh38-3.0.0 \\
--localcores=4 \\
--localmem=32 
#--force-cells=16000
#--expect-cells=35000
"""%(idline.split("/")[-2], idline, idline.split("/")[-2])
	print(cellcountline)
	ouf.write(cellcountline)
	sleepline = "\nsleep 5m\n\n\n"
	print(sleepline)
	ouf.write(sleepline)

ouf.write("shutdown -h 0")

ouf.close()
