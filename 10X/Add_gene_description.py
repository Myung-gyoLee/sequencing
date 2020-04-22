genef = open("H:/SingleCell/smc024/gene_description.txt")
geneline = genef.readlines()
genef.close()

gdict = {}
for gline in geneline:
    gdict[gline.split("\t")[0]] = gline.split("\t")[1].strip()

#with open("/media/cytogenbi2/8e7f6c8b-bc45-4c58-816f-a062fd95b91a/10X/HN00124804/SMC024/blood1_2/02DEG_SCTclustermt30_rna7000_rmMT.csv") as inf:
#    infline = inf.readlines()
    
with open("H:/SingleCell/smc024/analysis_tissue/deg_b1b2t.csv") as inf:
    infline = inf.readlines()
    
    
outline=[]
for iline in infline:
    iline1 = iline.split(',')[0].split('"')[1]
    if len(iline1)<1:
        continue
    print(iline1)
    try:
        print(gdict[iline1])
        oline='%s,"%s",%s'%(iline.split(',')[0], gdict[iline1], ','.join(iline.split(',')[1:]))
        outline.append(oline)
    except KeyError:
        print(iline1)
        oline='%s,"",%s'%(iline.split(',')[0], ','.join(iline.split(',')[1:]))
        outline.append(oline)
        
ouf = open("H:/SingleCell/smc024/analysis_tissue/deg_b1b2t_gene.csv","w")
ouf.write('"gene","gene_name","p_val","avg_logFC","pct.1","pct.2","p_val_adj","cluster","gene"\n')
ouf.write(''.join(outline))
ouf.close()