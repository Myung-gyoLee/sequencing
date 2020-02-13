import sys, os

with open(sys.argv[1]) as f:
	fline = f.readline().split("\t")

linenum=0
for line in fline:
	linenum+=1
	print(linenum, ":", line)


