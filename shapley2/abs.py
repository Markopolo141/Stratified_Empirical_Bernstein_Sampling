
averages = [0.0040299833, 0.0118521563, 0.0116650689, 0.0243446502, 0.053535758, 0.0674256858, 0.0717579693, 0.080559573, 0.0802304956, 0.0808407057, 0.0897705788, 0.0946289242, 0.0993384707, 0.1037406436, 0.1258651961]

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
