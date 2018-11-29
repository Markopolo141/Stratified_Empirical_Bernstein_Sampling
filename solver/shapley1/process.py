import sys
import pdb

filename = sys.argv[1]
with open(filename, "r") as f:
	data = f.readlines()
data = [d.strip().split(',') for d in data]
data = [[dd for dd in d if len(dd)>0] for d in data]
data = [[int(dd) if i==0 else float(dd) for i,dd in enumerate(d)] for d in data]
data_store = {}


for d in data:
	if d[0] not in data_store.keys():
		data_store[d[0]] = []
	data_store[d[0]].append(d[1:])

refined_data_store = []
for i in range(max(data_store.keys())+1):
	data = data_store[i]
	length = len(data[0])
	temp_data = [0.0 for i in range(length)]
	for d in data:
		for ii in range(length):
			temp_data[ii] += d[ii]
	for ii in range(length):
		temp_data[ii] /= len(data)
	refined_data_store.append(temp_data)

length = len(refined_data_store[0])
shapley_values = [0.0 for i in range(length)]
for i in range(len(refined_data_store)):
	for ii in range(length):
		shapley_values[ii] += refined_data_store[i][ii]

for ii in range(length):
	shapley_values[ii] /= len(refined_data_store)+0;

print shapley_values
print sum(shapley_values)

