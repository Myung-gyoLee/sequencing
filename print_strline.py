import sys, os

with open(sys.argv[1]) as f:
	fline = f.readlines()

f.close()


llist=[]
llist.append(fline[0])

for line in fline:
	if sys.argv[2] in line:
		llist.append(line)

jlen=len(fline[0].split("\t"))
ldic={}


if len(llist) == 1:
	print("No result")
elif len(llist) == 2:
	linenum=0
	for i in range(0,len(llist)):
		linenum+=1
		for j in range(0,jlen):
			print(linenum,":",fline[i].split("\t")[j], "|", fline[i+1].split("\t")[j])
else : 
	for j in range(0,jlen):
		for i in [1,len(llist)-1, len(llist)]:
			if i == 1:
				ldic[fline[0].split("\t")[j].strip()]=fline[i].split("\t")[j].strip()
			else :
				try:
					ldic[llist[0].split("\t")[j].strip()]+="|%s"%llist[i].split("\t")[j].strip()
				except IndexError:
					print('::::','IndexError:line number = ',i,'::::')
					continue

print("Matched line number = %s"%(len(llist)))

linenum2=0
for key, value in ldic.items():
	linenum2+=1
	print(linenum2, ":", key,value)


