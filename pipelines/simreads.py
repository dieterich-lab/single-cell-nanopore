from scipy.io import mmread
import pandas as pd
import numpy as np
import gzip
import sys

mtx_file = gzip.open(sys.argv[3], 'r')
mat = (mmread(mtx_file))
ind = mat.nonzero()
row = ind[0]
col = ind[1]
v = len(ind[0])
n = int(sys.argv[4])
l = int(sys.argv[5])
s = '*'

name_file = gzip.open(sys.argv[1], 'rt')
cell_names = name_file.read()
if type(cell_names) is bytes: cell_names = cell_names.decode('utf-8')
cell_names = cell_names.splitlines()
cell_names = np.array(['CB:Z:'+name for name in cell_names])

feature_file = gzip.open(sys.argv[2], 'rt')
features = feature_file.read()
if type(features) is bytes: features = features.decode('utf-8')
features = features.splitlines()
features = np.array(['GN:Z:'+(feature.split('\t')[1] if (len(feature.split('\t')) > 1) else feature) for feature in features])

i = np.random.choice(np.arange(v), n, replace=False, p=mat.data/sum(mat.data))
c = np.array('ATGC', dtype='c')
d = np.random.rand(n,l)*4
umi = [a.tostring().decode("utf-8") for a in c[d.astype(int)]]
print('\t'.join(['@SQ','SN:15','LN:101991189']))
pd.DataFrame({'QNAME':list(range(n)), 'FLAG':[0]*n, 'RNAME':[15]*n, 'POS':[69453703]*n, 'MAPQ':[255]*n, 'CIGAR':[str(len(s))+'M']*n, 'RNEXT':['*']*n, 'PNEXT':[0]*n, 'TLEN':[0]*n, 'SEQ':[s]*n, 'QUAL':[s]*n, 'GN':[a for a in features[row[i]]],'CB':[a for a in cell_names[col[i]]],'UB':['UB:Z:'+(a) for a in umi]}).to_csv(sys.stdout, header=False, index=False, sep="\t")
