#Experiment: Bernoulli and Uniform
#---------------------------------
#
# computes the average error for mean estimation in the 
# bernoulli and uniform stratified data case using different methods 
# with a sample budget between 1 and 20, outputs data to csv file


#do imports
from random import shuffle,random,randint,betavariate
from copy import deepcopy as copy
from methods import *
import sys
try:
	from tqdm import tqdm
	tqdm_enabled=True
except ImportError:
	tqdm_enabled=False


# parameters
d = 1.0 #data width is one in bernoulli and uniform
m = 300 #sample budget of 300
length = 1000 #strata have 1000 data points


print "Computing Bernoulli and Uniform Experiment"
print " for bernoulli successes 1 to 20 ---"
with open("Bernoulli_and_Uniform.csv","w") as f:
	f.write("successes,SEBM*,SEBM,Ney\n");
	for ps in range(1,21): # for bernoulli successes 1 to 20
		print ps

		# structures for strata population data points
		vals = [[],[]]
		# errors achived by the three methods
		error_values = [[],[],[]]

		iterator = range(20000)
		if tqdm_enabled:
			iterator = tqdm(iterator)
		for trial in iterator: #iterate a large number of times
			
			#setup all the data points
			vals[0] = [random() for i in range(length)]
			vals[1] = [0]*length
			for i in range(ps):
				vals[1][randint(0,length-1)]=1
			
			#calculate true population mean
			collected_vals = sum(vals,[])
			mean = sum(collected_vals)*1.0/len(collected_vals)
			
			#calculate error achieved using SEBM*
			cvals = copy(vals)
			error_values[0].append(abs(mean-burgess_ideal(cvals,m,d)))

			#calculate error achieved using SEBM
			cvals = copy(vals)
			error_values[1].append(abs(mean-burgess(cvals,m,d)))

			#calculate error achieved using Neyman sampling
			cvals = copy(vals)
			error_values[2].append(abs(mean-super_castro(cvals,m)))

		#average the errors achieved by each method and output
		error_values = [sum(errors)*1.0/len(errors) for errors in error_values]
		f.write("{},{},{},{}\n".format(ps,error_values[0],error_values[1],error_values[2]))

print "Finished Experiment"



