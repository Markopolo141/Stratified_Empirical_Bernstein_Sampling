
averages = [0.1033533858, 0.0961702681, 0.097771176, 0.0847447809, 0.0780227754, 0.0820221354, 0.05838364, 0.0594630579, 0.0535092258, 0.0547554271, 0.0462093208, 0.045666904, 0.02935292, 0.0288863418, 0.0287756776]

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
