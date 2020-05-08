#import GZip
#feature_bc = h5read("H:\\SingleCell\\smc009\\filtered_feature_bc_matrix.h5")

using MatrixMarket
mat = MatrixMarket.mmread("H:\\SingleCell\\smc009\\unzip\\matrix.mtx")

add DataFrames
using DataFrames

features = readtable("H:\\SingleCell\\smc009\\unzip\\features.tsv", header = false)
feature_ids = features[1] # ENSG id
gene_names = features[2]
barcodes = readtable("H:\\SingleCell\\smc009\\unzip\\barcodes.tsv", header = false)
barcode = barcodes[1]


df = convert(DataFrame, mat)
names!(df, Symbol.(barcode))
###
using Clustering
R = kmeans(mat, 20; maxiter=200, display=:iter)

@assert nclusters(R) == 20

a = assignments(R)
c = counts(R)
M = R.centers

## add gene name 
df[:feature] = gene_names 

using CSV

CSV.write("H:\\SingleCell\\smc009\\unzip\\smc00910x_flt_barcode.csv", df)
