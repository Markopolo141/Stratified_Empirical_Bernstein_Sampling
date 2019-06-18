weights = [45, 41, 27, 26, 25, 21, 13, 13, 12, 12, 11, 11, 10, 10, 10]
weights = [w*1.0/50 for w in weights]
def v(p):
	return (sum([weights[pp] for pp in p])**2) % 1.0

import sys
sys.path.append('../')
from shapley_methods import *

file_root = "data_out_{}_{}.csv"
N=len(weights)
for m in [N*N*10,N*N*50,N*N*100,N*N*500,N*N*1000]:
	print "calculating Castro for {}".format(m)
	sample_error(file_root, castro, 40, N,m,v)
for m in [N*N*10,N*N*50,N*N*100,N*N*500,N*N*1000]:
	print "calculating Maleki for {}".format(m)
	sample_error(file_root, maleki, 40, N,m,v)
for m in [N*N*10,N*N*50,N*N*100,N*N*500,N*N*1000]:
	print "calculating Simple for {}".format(m)
	sample_error(file_root, simple, 40, N,m,v)
for m in [N*N*10,N*N*50,N*N*100,N*N*500,N*N*1000]:
	print "calculating ApproShapley for {}".format(m)
	sample_error(file_root, approshapley, 40, N,m,v)

