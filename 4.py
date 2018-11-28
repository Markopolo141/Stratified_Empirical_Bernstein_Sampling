from random import shuffle,random,randint,betavariate
import pdb
from copy import deepcopy as copy
from tqdm import tqdm

from bounds import *


import sys

d=1.0
m = 300
ps=int(sys.argv[1])


def bernoulli(p):
	if random()<p:
		return 1
	return 0

print ps
with open("data2_{}_{}.csv".format(m,ps),"w") as f:
	f.write("burgess_ideal,burgess,castro\n");
	vals = []
	N=2
	length = 1000
	vals.append([random() for i in range(length)])
	vals.append([0]*length)
	for i in range(ps):
		vals[1][randint(0,length-1)]=1
	for trial in tqdm(range(20000)):
		vals[0] = [random() for i in range(length)]
		vals[1] = [0]*length
		for i in range(ps):
			vals[1][randint(0,length-1)]=1
		
		collected_vals = sum(vals,[])
		mean = sum(collected_vals)*1.0/len(collected_vals)
		cvals = copy(vals)
		f.write("{},".format(abs(mean-burgess_ideal(cvals,m,d))))
		cvals = copy(vals)
		f.write("{},".format(abs(mean-burgess(cvals,m,d))))
		cvals = copy(vals)
		f.write("{},".format(abs(mean-super_castro(cvals,m))))
		f.write("\n")



