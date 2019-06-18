weights = [1,1,2,2,2,3,4,5,5,5,7,8,8,8,10]
def v(p):
	return 1.0*max([0]+[weights[pp] for pp in p])

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

