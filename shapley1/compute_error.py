
from computed_shapley_value import *
#averages = [0.0666666667, 0.0666666667, 0.143571418, 0.1434954697, 0.143722395, 0.243242812, 0.354821483, 0.4802750499, 0.4789936291, 0.4796110297, 0.8797155552, 1.1296792262, 1.1296860895, 1.1337698345, 3.1274355565]

import sys
import pdb

filename = sys.argv[1]

with open(filename,'r') as f:
	data = f.readlines()
data = [d.strip().split(",") for d in data]
data = [[float(dd) for dd in d if len(dd)>0] for d in data]
data = [[abs(dd-averages[i]) for i,dd in enumerate(d)] for d in data]
data = sum(sum(data,[]))/(40*15)
print "{:.5}".format(10000*data)

with open("../Shapley1_Error_output.txt","a") as f:
	f.write(filename)
	f.write("\n{:.5}\n\n".format(10000*data))
