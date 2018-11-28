from random import shuffle,random,randint,betavariate
import pdb
from copy import deepcopy as copy
from tqdm import tqdm
import sys

from bounds import *


d=1.0
mperN=int(sys.argv[1])

print mperN


with open("data{}.csv".format(mperN),"w") as f:
	f.write("burgess_ideal,burgess,burgess_small,simple,simple_small,super_castro,super_castro_small\n");
	for trial in tqdm(range(5000)):
		vals = []
		N = randint(5,21)
		Ni = [randint(10,201) for i in range(N)]
		m = mperN * N
		while (sum(Ni)<m):
			N = randint(5,21)
			Ni = [randint(10,201) for i in range(N)]
			m = mperN * N
		for i in range(N):
			alpha = random()*4
			beta = random()*4
			vals.append([betavariate(alpha,beta) for ii in range(Ni[i])])
		collected_vals = sum(vals,[])
		mean = sum(collected_vals)*1.0/len(collected_vals)
		cvals = copy(vals)
		f.write("{},".format(abs(mean-burgess_ideal(cvals,m,d))))
		cvals = copy(vals)
		f.write("{},".format(abs(mean-burgess(cvals,m,d))))
		cvals = copy(vals)
		f.write("{},".format(abs(mean-burgess_small(cvals,m,d))))
		cvals = copy(vals)
		f.write("{},".format(abs(mean-simple(cvals,m))))
		cvals = copy(vals)
		f.write("{},".format(abs(mean-simple_small(cvals,m))))
		cvals = copy(vals)
		f.write("{},".format(abs(mean-super_castro(cvals,m))))
		cvals = copy(vals)
		f.write("{}\n".format(abs(mean-super_castro_small(cvals,m))))






