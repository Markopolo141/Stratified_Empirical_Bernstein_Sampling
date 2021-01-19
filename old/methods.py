#  Methods.py
#  ---------
#
# The Implementation of all the sampling methods used for the synthetic data experiments
# and some utility functions

from random import shuffle,choice
from math import floor,factorial,log,sqrt

'''given an array of numbers a, scale thoes numbers to closest integers such as to sum to integer v'''
def scale_array(a,v):
	if sum(a) != 0:
		vsuma = v*1.0/sum(a)
		a = [vsuma*aa for aa in a]
		overflow = [aa-floor(aa) for aa in a]
		a = [int(floor(aa)) for aa in a]
		while sum(a)<v:
			i = overflow.index(max(overflow))
			overflow[i]-=1
			a[i]+=1
		return a
	else:
		a = [int(v*1.0/len(a)) for aa in a]
		i=-1
		while sum(a)<=v:
			i += 1
			a[i] += 1
		a[i] -= 1
		return a


'''
SIMPLE sampling method - without replacement
given the stratified population data values, jumbles them together and takes the average of the first m
'''
def simple(vals,m):
	vals = sum(vals,[])
	shuffle(vals)
	return sum(vals[0:m])*1.0/m


'''
SIMPLE sampling method - with replacement
given the stratified population data values, jumbles them together and takes the average of m of them chosen randomly with replacement
'''
def simple_small(vals,m):
	vals = sum(vals,[])
	return sum([choice(vals) for i in range(m)])*1.0/m


'''
NEYMAN sampling method - with replacement
given the stratified data, computes the variance of the strata, and allocates samples per strata (above atleast 1 sample) according to the variance
computes the mean by sampling with replacement of thoes strata per thoes allocations
'''
def super_castro_small(vals,m):
	lengths = [len(v) for v in vals]
	sumN = sum([len(v) for v in vals])
	meanvals = [sum(vals[i])*1.0/len(vals[i]) for i in range(len(vals))]
	varvals = [sum([(vals[i][o]-meanvals[i])**2 for o in range(len(vals[i]))])*1.0/len(vals[i]) for i in range(len(vals))]
	samples = [1 for i in range(len(vals))]
	new_samples = scale_array([sqrt(varvals[i])*lengths[i] for i in range(len(vals))],m-sum(samples))
	for i in range(len(samples)):
		samples[i] += new_samples[i]
	return sum([sum([choice(vals[i]) for ii in range(ss)])*len(vals[i])*1.0/ss for i,ss in enumerate(samples)])/sumN


'''
NEYMAN sampling method - without replacement
given the stratified data, computes the variance of the strata, and tries to allocates the sample budget across strata (above atleast 1 sample each) according to the variance
if the any strata would be oversampled, it clips the samples of that strata to its maximum and reallocates remaining samplings among the remainder
computes the mean by sampling without replacement of thoes strata per thoes allocations
'''
def super_castro(vals,m):
	lengths = [len(v) for v in vals]
	sumN = sum([len(v) for v in vals])
	meanvals = [sum(vals[i])*1.0/len(vals[i]) for i in range(len(vals))]
	varvals = [sum([(vals[i][o]-meanvals[i])**2 for o in range(len(vals[i]))])*1.0/len(vals[i]) + 0.00000001 for i in range(len(vals))]
	samples = [1 for i in range(len(varvals))]
	allocated = sum(samples)

	while allocated < m:
		new_samples = scale_array([sqrt(varvals[i])*lengths[i] for i in range(len(vals))],m-allocated)
		for i in range(len(samples)):
			if samples[i]+new_samples[i] >= len(vals[i]):
				varvals[i]=0
				additional = len(vals[i])-samples[i]
			else:
				additional = new_samples[i]
			samples[i] += additional
			allocated += additional
	return sum([sum(vals[i][0:ss])*len(vals[i])*1.0/ss for i,ss in enumerate(samples)])/sumN


'''
Computes the SEBB for sampling without replacement:
N = integer, the number of strata
ni = array of integers, how much each strata is allready sampled
Ni = array of integers, the size of each strata
var = array of floats, the sample variance of each strata as sampled
d = float, the width of the population datas
r = float, the confidence level of the bound
'''
#OmegaBig = lambda n,N: sum([1.0/(k**2) for k in range(n,N)])
#PsiBig = lambda n,N: N*sum([1.0/(k**2*(k+1)) for k in range(n,N)])
OmegaBig = lambda n,N: (n+1)*(1-n*1.0/N)*1.0/(n**2)
PsiBig = lambda n,N: (N+1.0-n)/(n**2)
OmegaSmall = lambda n,N: 1.0/n
PsiSmall = lambda n,N: 1.0/n
def sebm_bound(N,ni,Ni,var,d,r):
	sumN = sum(Ni)
	onN = [0 for i in range(2)]
	max1 = [0 for i in range(2)]
	var1 = [0 for i in range(2)]
	d1 = [0 for i in range(2)]
	log6r = log(6/r)
	log3r = log(3/r)
	log2r = log(2/r)
	log6Nr = log(6*N/r)
	d2 = d*d;
	for i in range(N):
		tau = Ni[i]*1.0/sumN
		OB = OmegaBig(ni[i],Ni[i])
		OS = OmegaSmall(ni[i],Ni[i])
		PB = PsiBig(ni[i],Ni[i])
		PS = PsiSmall(ni[i],Ni[i])
		onN[0] += PB*min(OB,OS)*tau**2
		onN[1] += PS*min(OB,OS)*tau**2
		max1[0] = max(max1[0],PB*min(PB,PS)*tau**2)
		max1[1] = max(max1[1],PS*min(PB,PS)*tau**2)
		var1[0] += PB*((ni[i]-1)*var[i]*1.0/ni[i])*tau**2
		var1[1] += PS*((ni[i]-1)*var[i]*1.0/ni[i])*tau**2
		d1[0] += OB*tau**2
		d1[1] += OS*tau**2
	A = [0,0]
	A[0] = (d2*4.0/(17))*log6r*d1[0] + log6r*(sqrt(2*var1[0] + log6Nr*d2*onN[0] + log3r*d2*max1[0]) + sqrt(log3r*d2*max1[0]))**2
	A[1] = (d2*4.0/(17))*log6r*d1[1] + log6r*(sqrt(2*var1[1] + log6Nr*d2*onN[1] + log3r*d2*max1[1]) + sqrt(log3r*d2*max1[1]))**2
	return sqrt(min(A))

'''
Computes the SEBB for sampling with replacement:
N = integer, the number of strata
ni = array of integers, how much each strata is allready sampled
Ni = array of integers, the size of each strata
var = array of floats, the sample variance of each strata as sampled
d = float, the width of the population datas
r = float, the confidence level of the bound
'''
def sebm_bound_small(N,ni,Ni,var,d,r):
	sumN = sum(Ni)
	onN = 0
	max1 = 0
	var1 = 0
	d1 = 0
	log6r = log(6/r)
	log3r = log(3/r)
	log2r = log(2/r)
	logN = log(N)
	log2 = log(2)
	d2 = d*d;
	for i in range(N):
		tau = Ni[i]*1.0/sumN
		OS = OmegaSmall(ni[i],Ni[i])
		PS = PsiSmall(ni[i],Ni[i])
		onN += PS*OS*tau**2
		max1 = max(max1,PS*PS*tau**2)
		var1 += PS*((ni[i]-1)*var[i]*1.0/ni[i])*tau**2
		d1 += OS*tau**2
	A = [0]
	A[0] = (d2*4.0/(17))*log6r*d1 + log6r*(sqrt(2*var1 + (log6r+logN)*d2*onN/log2 + log3r*d2*max1) + sqrt(log3r*d2*max1))**2
	return sqrt(min(A))


'''
Computes SEBB* for sampling without replacement:
N = integer, the number of strata
ni = array of integers, how much each strata is allready sampled
Ni = array of integers, the size of each strata
var = array of floats, the sample variance of each strata as sampled
d = float, the width of the population datas
r = float, the confidence level of the bound
'''
def sebm_bound_ideal(N,Ni,ni,var,d):
	oversumN = 1.0/sum(Ni)
	v = 0
	d2 = d*d
	o17 = 1.0/17
	for i in range(N):
		a = (OmegaBig(ni[i],Ni[i])*d2*o17 + PsiBig(ni[i],Ni[i])*var[i]*0.5)*(Ni[i]*oversumN)**2
		b = (OmegaSmall(ni[i],Ni[i])*d2*o17 + PsiSmall(ni[i],Ni[i])*var[i]*0.5)*(Ni[i]*oversumN)**2
		#c = (OmegaBig(ni[i],Ni[i])*d2*1.0/2)*(Ni[i]*1.0/sumN)**2
		#d = (OmegaSmall(ni[i],Ni[i])*d2*1.0/2)*(Ni[i]*1.0/sumN)**2
		v += min(a,b)
		#v += min(a,b,c,d)
	return v


'''
Computes SEBB* for sampling with replacement:
N = integer, the number of strata
ni = array of integers, how much each strata is allready sampled
Ni = array of integers, the size of each strata
var = array of floats, the sample variance of each strata as sampled
d = float, the width of the population datas
r = float, the confidence level of the bound
'''
def sebm_bound_ideal_small(N,Ni,ni,var,d):
	oversumN = 1.0/sum(Ni)
	v = 0
	d2 = d*d
	o17 = 1.0/17
	for i in range(N):
		b = (OmegaSmall(ni[i],Ni[i])*d2*o17 + PsiSmall(ni[i],Ni[i])*var[i]*0.5)*(Ni[i]*oversumN)**2
		#d = (OmegaSmall(ni[i],Ni[i])*d2*1.0/2)*(Ni[i]*1.0/sumN)**2
		v += b
	return v


'''
The SEBM* method without replacement:
calculates the variance of the strata, and then iteratively allocates the budget to iteratively minimise SEBB*
'''
def sebm_ideal(vals,m,d):
	Ni = [len(v) for v in vals]
	N = len(Ni)
	ni = [1 for i in range(N)]
	mean = [sum(vals[i])*1.0/Ni[i] for i in range(len(vals))]
	var = [sum([(vals[i][o]-mean[i])**2 for o in range(Ni[i])])*1.0/Ni[i] for i in range(len(vals))]
	allocated = sum(ni)
	while allocated < m:
		min_i = -1
		min_advantage = float('inf')
		for i in range(N):
			if ni[i]<Ni[i]:
				ni[i] += 1
				adv = sebm_bound_ideal(N,Ni,ni,var,d)
				ni[i] -= 1
				if adv<min_advantage:
					min_advantage = adv
					min_i = i
		ni[min_i] += 1
		allocated += 1
	return sum([sum(vals[i][0:ss])*Ni[i]*1.0/ss for i,ss in enumerate(ni)])/sum(Ni)
		
'''
The SEBM* method with replacement:
calculates the variance of the strata, and then iteratively allocates the budget to iteratively minimise SEBB*
'''
def sebm_ideal_small(vals,m,d):
	Ni = [len(v) for v in vals]
	N = len(Ni)
	ni = [1 for i in range(N)]
	mean = [sum(vals[i])*1.0/Ni[i] for i in range(len(vals))]
	var = [sum([(vals[i][o]-mean[i])**2 for o in range(Ni[i])])*1.0/Ni[i] for i in range(len(vals))]
	allocated = sum(ni)
	while allocated < m:
		min_i = -1
		min_advantage = float('inf')
		for i in range(N):
			ni[i] += 1
			adv = sebm_bound_ideal_small(N,Ni,ni,var,d)
			ni[i] -= 1
			if adv<min_advantage:
				min_advantage = adv
				min_i = i
		ni[min_i] += 1
		allocated += 1
	return sum([sum([choice(vals[i]) for ii in range(ss)])*Ni[i]*1.0/ss for i,ss in enumerate(ni)])/sum(Ni)




'''
The SEBM method with replacement:
calculates the variance of the strata, and then iteratively allocates the budget to iteratively minimise SEBB
'''
def sebm(vals,m,d,r=0.5):
	Ni = [len(v) for v in vals]
	N = len(Ni)
	ni = [0 for i in range(N)]
	s = [0.0 for i in range(N)]
	s2 = [0.0 for i in range(N)]
	var = [0.0 for i in range(N)]
	samples = 0
	# seed with minimum initial two samples (if possible) for all strata
	for i in range(N): #strata i
		for p in range(2):
			if ni[i]<Ni[i] or Ni[i]==-1:
				v = vals[i].pop()
				s[i]+=v
				s2[i]+=v*v
				ni[i]+=1
				if ni[i]>1:
					var[i] = (s2[i] - s[i]**2*1.0/ni[i])/(ni[i]-1)
				samples+=1
	advantage = [0.0 for i in range(N)]
	while samples < m:
		#calculate the bound as it exists:
		bound = sebm_bound(N,ni,Ni,var,d,r)
		#calculate the advantages possible
		for i in range(N):
			if ni[i]<Ni[i] or Ni[i]==-1:
				ni[i]+=1
				advantage[i] = bound-sebm_bound(N,ni,Ni,var,d,r)
				ni[i]-=1
			else:
				advantage[i]=-float("inf")
		#detect the sample that maximises advantage
		maxi=0
		maxadvantage=-float("inf")
		for i in range(N):
			if advantage[i] > maxadvantage:
				maxi = i
				maxadvantage = advantage[i]
		#take the best sample
		v = vals[maxi].pop()
		s[maxi]+=v
		s2[maxi]+=v*v
		ni[maxi]+=1
		var[maxi] = (s2[maxi] - s[maxi]**2*1.0/ni[maxi])/(ni[maxi]-1)
		samples += 1
	#return the stratafied average from the samples
	ss = 0.0
	for i in range(N):
		ss += s[i]*Ni[i]/ni[i]
	ss /= sum(Ni)
	return ss


'''
The SEBM method without replacement:
calculates the variance of the strata, and then iteratively allocates the budget to iteratively minimise SEBB
'''
def sebm_small(vals,m,d,r=0.5):
	Ni = [len(v) for v in vals]
	N = len(Ni)
	ni = [0 for i in range(N)]
	s = [0.0 for i in range(N)]
	s2 = [0.0 for i in range(N)]
	var = [0.0 for i in range(N)]
	samples = 0
	# seed with minimum initial two samples (if possible) for all strata
	for i in range(N): #strata i
		for p in range(2):
			v = choice(vals[i])
			s[i]+=v
			s2[i]+=v*v
			ni[i]+=1
			if ni[i]>1:
				var[i] = (s2[i] - s[i]**2*1.0/ni[i])/(ni[i]-1)
			samples+=1
	advantage = [0.0 for i in range(N)]
	while samples < m:
		#calculate the bound as it exists:
		bound = sebm_bound_small(N,ni,Ni,var,d,r)
		#calculate the advantages possible
		for i in range(N):
			ni[i]+=1
			advantage[i] = bound-sebm_bound_small(N,ni,Ni,var,d,r)
			ni[i]-=1
		#detect the sample that maximises advantage
		maxi=0
		maxadvantage=-float("inf")
		for i in range(N):
			if advantage[i] > maxadvantage:
				maxi = i
				maxadvantage = advantage[i]
		#take the best sample
		v = choice(vals[maxi])
		s[maxi]+=v
		s2[maxi]+=v*v
		ni[maxi]+=1
		var[maxi] = (s2[maxi] - s[maxi]**2*1.0/ni[maxi])/(ni[maxi]-1)
		samples += 1
	#return the stratafied average from the samples
	ss = 0.0
	for i in range(N):
		ss += s[i]*Ni[i]/ni[i]
	ss /= sum(Ni)
	return ss



