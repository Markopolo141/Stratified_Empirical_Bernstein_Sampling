#Experiment: Synthetic Data
#--------------------------
#
# computes the average error for mean estimation in the case
# of randomly generated strata numbers of beta distributed data
# for sample budgets between 10 and 150, outputs data to csv file


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

d=1.0 #data width is one for beta distribution

print "Computing Beta data Experiment"
print " for sample budget per strata of [10,50,100,150] ---"

#open csv file
with open("Beta_Synthetic.csv","w") as f:

	#for each sample budget
	for mperN in [10,50,100,150]:
		print mperN

		#data holding for the errors
		sampling_errors = [[],[],[],[],[],[],[]]
		
		iterator = range(5000)
		if tqdm_enabled:
			iterator = tqdm(iterator)
		for trial in iterator: #iterate a large number of times

			# generate strata dimensions
			N = randint(5,21)							#number of strata
			Ni = [randint(10,201) for i in range(N)]	#strata sizes
			m = mperN * N								#sample budget
			while (sum(Ni)<m):							#regenerate while sample budget is too large
				N = randint(5,21)
				Ni = [randint(10,201) for i in range(N)]
				m = mperN * N

			#generate population data values for the strata
			vals = []
			for i in range(N):
				alpha = random()*4
				beta = random()*4
				vals.append([betavariate(alpha,beta) for ii in range(Ni[i])])

			#calculate actual population mean
			collected_vals = sum(vals,[])
			mean = sum(collected_vals)*1.0/len(collected_vals)

			#calculate error in using SEBM* method
			cvals = copy(vals)
			sampling_errors[0].append(abs(mean-burgess_ideal(cvals,m,d)))

			#calculate error in using SEBM method
			cvals = copy(vals)
			sampling_errors[1].append(abs(mean-burgess(cvals,m,d)))

			#calculate error in using SEBM method with replacement
			cvals = copy(vals)
			sampling_errors[2].append(abs(mean-burgess_small(cvals,m,d)))

			#calculate error in using simple sampling method without replacement
			cvals = copy(vals)
			sampling_errors[3].append(abs(mean-simple(cvals,m)))

			#calculate error in using simple sampling method with replacement
			cvals = copy(vals)
			sampling_errors[4].append(abs(mean-simple_small(cvals,m)))

			#calculate error in using Neyman sampling method with replacement
			cvals = copy(vals)
			sampling_errors[5].append(abs(mean-super_castro(cvals,m)))

			#calculate error in using Neyman sampling method without replacement
			cvals = copy(vals)
			sampling_errors[6].append(abs(mean-super_castro_small(cvals,m)))

		#sorting the sampling errors for quartile reporting
		sampling_errors = [sorted(s) for s in sampling_errors]

		#write csv entry, using rough quartile reporting
		f.write("{} sample budget per strata\n".format(mperN))
		f.write(",SEBM*,SEBM,SEBM-W,simple,simple-w,Ney,Ney-W\n")
		for p in [0.09,0.25,0.5,0.75,0.91]:
			pp = [p]+[s[int(len(s)*p)] for s in sampling_errors]
			f.write("{} percentile,{},{},{},{},{},{},{}\n".format(*pp))

print "Finished Experiment"



