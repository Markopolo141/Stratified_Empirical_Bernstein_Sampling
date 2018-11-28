
import sys
import pdb
from statistics import mean

filename = sys.argv[1]
n = sys.argv[2]
v = sys.argv[3]

with open(filename,'r') as f:
	data = f.readlines()
data = [d.strip().split(",") for d in data]
data = [[float(dd) for dd in d if len(dd)>0] for d in data[1:]]
data = [[dd[i] for dd in data] for i in range(len(data[0]))]
data = [mean(dd) for dd in data]
print "({}, {:.6})".format(n,data[int(v)])
