weights = [1,3,3,6,12,16,17,19,19,19,21,22,23,24,29]
sum_weights = sum(weights)
def v(p):
	if (sum([weights[pp] for pp in p]) > sum_weights/2):
		return 1.0
	return 0.0

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

