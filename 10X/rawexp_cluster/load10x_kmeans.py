#!bin/python3
import csv
import gzip
import os
import scipy.io

matrix_dir = "H:\\SingleCell\\smc009\\filtered_feature_bc_matrix\\"
mat = scipy.io.mmread(os.path.join(matrix_dir, "matrix.mtx.gz"))

features_path = os.path.join(matrix_dir, "features.tsv.gz")

feature_ids = [row[0] for row in csv.reader(gzip.open(features_path, mode = 'rt'), delimiter="\t")]
gene_names = [row[1] for row in csv.reader(gzip.open(features_path, mode = 'rt'), delimiter="\t")]
feature_types = [row[2] for row in csv.reader(gzip.open(features_path, mode = 'rt'), delimiter="\t")]
barcodes_path = os.path.join(matrix_dir, "barcodes.tsv.gz")
barcodes = [row[0] for row in csv.reader(gzip.open(barcodes_path), mode = 'rt', delimiter="\t")]


mat1 = pd.DataFrame.sparse.from_spmatrix(mat)

mat1.columns = barcodes
mat1.index = gene_names


mat1_transpose = mat1.transpose()
mat1_transpose.loc[["cluster_id"]]
mat1_t_del = mat1_transpose.drop('cluster_id', axis=0)

kmeans1 = KMeans(n_clusters = 30, random_state=0).fit(mat1_t_del)
mat1_t_del['cluster_id'] = kmeans1.labels_


mat1_t_del[mat1_t_del['cluster_id'] == 1].loc[:, mat1_t_del.columns != 'cluster_id'].sum(axis=1)

mat1_t_del[mat1_t_del['cluster_id'] == 1].loc[:, mat1_t_del.columns != 'cluster_id'].sum(axis=0)

cl1sum = mat1_t_del[mat1_t_del['cluster_id'] == 1].loc[:, mat1_t_del.columns != 'cluster_id'].sum(axis=0)

n=0
for clsnum in range(len(mat1_t_del['cluster_id'].unique())-1):
    n+=1
    clsname = "cl%ssum"%(n-1)
    print(clsname, n)
    vars()[clsname] = mat1_t_del[mat1_t_del['cluster_id'] == n].loc[:, mat1_t_del.columns != 'cluster_id'].sum(axis=0)
    
n=0
for clsnum in range(len(mat1_t_del['cluster_id'].unique())-1):
    n+=1
    clsname = "cl%ssum"%(n-1)
    print("%s.sort_values(ascending = False).head(50)"%clsname)


cl0sum.sort_values(ascending = False).head(50)
cl1sum.sort_values(ascending = False).head(50)
cl2sum.sort_values(ascending = False).head(50)
cl3sum.sort_values(ascending = False).head(50)
cl4sum.sort_values(ascending = False).head(50)
cl5sum.sort_values(ascending = False).head(50)
cl6sum.sort_values(ascending = False).head(50)
cl7sum.sort_values(ascending = False).head(50)
cl8sum.sort_values(ascending = False).head(50)
cl9sum.sort_values(ascending = False).head(50)
cl10sum.sort_values(ascending = False).head(50)
cl11sum.sort_values(ascending = False).head(50)
cl12sum.sort_values(ascending = False).head(50)
cl13sum.sort_values(ascending = False).head(50)
cl14sum.sort_values(ascending = False).head(50)
cl15sum.sort_values(ascending = False).head(50)
cl16sum.sort_values(ascending = False).head(50)
cl17sum.sort_values(ascending = False).head(50)
cl18sum.sort_values(ascending = False).head(50)
cl19sum.sort_values(ascending = False).head(50)
cl20sum.sort_values(ascending = False).head(50)
cl21sum.sort_values(ascending = False).head(50)
cl22sum.sort_values(ascending = False).head(50)
cl23sum.sort_values(ascending = False).head(50)
cl24sum.sort_values(ascending = False).head(50)
cl25sum.sort_values(ascending = False).head(50)
cl26sum.sort_values(ascending = False).head(50)
cl27sum.sort_values(ascending = False).head(50)
cl28sum.sort_values(ascending = False).head(50)


mat1.to_csv("H:\\SingleCell\\smc009\\smc009matrix_featurebarcode.csv")

"""
features_file = csv.reader(gzip.open(features_path, mode = 'rt'), delimiter="\t")

feature_ids = [row[0] for row in features_file]
gene_names = [row[1] for row in features_file]
feature_types = [row[2] for row in features_file]
barcodes_path = os.path.join(matrix_dir, "barcodes.tsv.gz")
barcodes = [row[0] for row in csv.reader(gzip.open(barcodes_path, mode = 'rt'), delimiter="\t")]

feature_ids = [row[0] for row in csv.reader(gzip.open(features_path), delimiter="\t")]
gene_names = [row[1] for row in csv.reader(gzip.open(features_path), delimiter="\t")]
feature_types = [row[2] for row in csv.reader(gzip.open(features_path), delimiter="\t")]
barcodes_path = os.path.join(matrix_dir, "barcodes.tsv.gz")
barcodes = [row[0] for row in csv.reader(gzip.open(barcodes_path), delimiter="\t")]
"""
