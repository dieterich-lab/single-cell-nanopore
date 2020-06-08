from scipy.io import mmread
import sys

mtx_file = gzip.open(sys.argv[3], 'r')
mat = (mmread(mtx_file))
ind = mat.nonzero()
row = ind[0]
col = ind[1]
v = len(ind[0])
n = int(sys.argv[4])
l = int(sys.argv[5])

name_file = gzip.open(sys.argv[1], 'rt')
cell_names = name_file.read().decode('utf-8').splitlines()
cell_names = np.array(['CB:Z:'+name for name in cell_names])

feature_file = gzip.open(sys.argv[2], 'rt')
features = feature_file.read().decode('utf-8').splitlines()
features = np.array(['GN:Z:'+(feature.split('\t')[1] if (len(feature.split('\t')) > 1) else feature) for feature in features])

i = np.random.choice(np.arange(v), n, replace=False, p=mat.data/sum(mat.data))
c = np.array('ATGC', dtype='c')
d = np.random.rand(n,l)*4
umi = [a.tostring().decode("utf-8") for a in c[d.astype(int)]]
pd.DataFrame({'GN':[a for a in features[row[i]]],'CB':['CB:Z:'+(a) for a in cell_names[col[i]]],'UB':['UB:Z:'+(a) for a in umi]}).to_csv(sys.argv[6], header=False, index=False, sep="\t")
