import sys, os

print("~/read_column.py [txtfile] [delimiter]")
with open(sys.argv[1]) as f:
        fline = f.readline().split(sys.argv[2])

linenum=0
for line in fline:
        linenum+=1
        if '"' in line:print(linenum, ":", line.replace('"',''))
        else:print(linenum, ":", line)
