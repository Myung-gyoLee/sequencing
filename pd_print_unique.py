from pandas import DataFrame
import pandas as pd
import sys

df = pd.read_csv(sys.argv[1], sep = '\t')
collist = df.columns.tolist()

linenum = 0
for col in collist:
	linenum+=1
	print('::::',linenum,'-',col,'::::')
	print(df[col].unique())

