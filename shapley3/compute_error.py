
from computed_shapley_value import *
#averages = [0.6457383309, 0.5883115375, 0.387419688, 0.373121619, 0.3587359562, 0.3013723492, 0.1865498819, 0.186543054, 0.1721918502, 0.1722102474, 0.1578253315, 0.1578405863, 0.1435026585, 0.1434831816, 0.143495217]

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

with open("../Shapley3_Error_output.txt","a") as f:
	f.write(filename)
	f.write("\n{:.5}\n\n".format(10000*data))
