from random import shuffle,random,randint
from math import floor,factorial,log,sqrt
from itertools import combinations
import pdb
from copy import deepcopy as copy

minus_floor = lambda x:x-floor(x)
int_floor = lambda x:int(floor(x))

def scale_array(a,v):
	# given an array of numbers, scale thoes numbers to closest integers such as to sum to integer v
	vsuma = v*1.0/sum(a)
	a = [vsuma*aa for aa in a]
	overflow = [aa-floor(aa) for aa in a]
	a = [int(floor(aa)) for aa in a]
	while sum(a)<v:
		i = overflow.index(max(overflow))
		overflow[i]-=1
		a[i]+=1
	return a

def scale_array2d(a,v):
	# given a wrapper2d array of numbers, scale thoes numbers to closest integers such as to sum to integer v
	vsuma = v*1.0/sum(a)
	a = copy(a)
	a*vsuma
	overflow = copy(a)
	overflow.map(minus_floor)
	a.map(int_floor)
	while sum(a)<v:
		i = overflow.random_index(max(overflow))
		overflow[i]-=1
		a[i]+=1
	return a

class wrapper2d(object):
	def __init__(self,o,d1,d2):
		self.o = o
		self.d1 = d1
		self.d2 = d2
	def __getitem__(self,i):
		ii = i/self.d1
		return self.o[ii][i-ii*self.d1]
	def __setitem__(self,i,v):
		ii = i/self.d1
		self.o[ii][i-ii*self.d1] = v
	def __mul__(self,m):
		for x in range(self.d1):
			for y in range(self.d2):
				self.o[x][y] *= m
	def __rsub__(self,v):
		for x in range(self.d1):
			for y in range(self.d2):
				self.o[x][y] -= v
	def __radd__(self,v):
		for x in range(self.d1):
			for y in range(self.d2):
				self.o[x][y] += v
	def map(self,operator):
		for x in range(self.d1):
			for y in range(self.d2):
				self.o[x][y] = operator(self.o[x][y])
	def index(self,v):
		i = 0
		for x in range(self.d1):
			for y in range(self.d2):
				if self.o[x][y] == v:
					return i
				i += 1
		raise ValueError("{} is not in list".format(v))
	def random_index(self,v):
		i = 0
		indices = []
		for x in range(self.d1):
			for y in range(self.d2):
				if self.o[x][y] == v:
					indices.append(i)
				i += 1
		if len(indices)==0:
			raise ValueError("{} is not in list".format(v))
		shuffle(indices)
		return indices[0]


weights = [1,1,2,2,2,3,4,5,5,5,7,8,8,8,10]
def v(p):
	return 1.0*max([0]+[weights[pp] for pp in p])


def gen1(l,i,N):
	r = list(range(N))
	r.remove(i)
	shuffle(r)
	return r[0:l]
def gen2(l,i,N):
	g = gen1(l,i,N)
	return v(g+[i])-v(g)


def gen3(l,i,N,gen_list):
	r = list(range(N))
	r.remove(i)
	shuffle(r)
	for c in combinations(r,l):
		cs = set(c)
		if cs not in gen_list:
			gen_list.append(cs)
			break
	else:
		raise Exception("cannot create novel coalition")
	return v(list(c)+[i])-v(c)
	


def castro(N,m,v):
	m1exp = int(floor(m*1.0/(2*N**2)))
	s = [[0.0 for i in range(N)] for i in range(N)]
	var = [[0.0 for i in range(N)] for i in range(N)]
	sumvar = 0
	for l in range(N):
		for i in range(N):
			for c in range(0,m1exp):
				x0 = gen2(l,i,N)
				s[l][i] += x0
				var[l][i] += x0**2
			ss = (var[l][i]-s[l][i]**2/m1exp)/(m1exp-1)
			var[l][i] = ss
			sumvar += ss
	m2 = m - N*N*m1exp
	for l in range(N):
		for i in range(N):
			if sumvar==0:
				var[l][i] = m1exp
			else:
				var[l][i] = m2*var[l][i]/sumvar - m1exp
				var[l][i] = max(0,var[l][i])
	var = scale_array2d(wrapper2d(var,N,N),m2).o
	for l in range(N):
		for i in range(N):
			for c in range(var[l][i]):
				x0 = gen2(l,i,N)
				s[l][i] += x0
			s[l][i] = s[l][i]*1.0/(var[l][i]+m1exp)

	'''for o in range(N):
		for i in range(N):
			print "{}\t".format(m1exp+var[i][o]),
		print ""'''

	ss = [0.0 for i in range(N)]
	for l in range(N):
		for i in range(N):
			ss[i] += s[l][i]
	for i in range(N):
		ss[i] /= N
	return ss

def maleki(N,m,v):
	mm = [[int(floor(m*pow(i+1,2.0/3)/(N*sum([pow(j+1,2.0/3) for j in range(N)])) )) for i in range(N)] for ii in range(N)]
	summ = sum([sum(aa) for aa in mm])
	i = 0
	while summ<m:
		ii = 0
		while summ<m and ii < N:
			mm[ii][i] += 1
			summ+=1
			ii+=1
		i = (i+1)%N
	s = [[0.0 for i in range(N)] for i in range(N)]
	for l in range(N):
		for i in range(N):
			for c in range(mm[l][i]):
				x0 = gen2(l,i,N)
				s[l][i] += x0
			s[l][i] = s[l][i]*1.0/mm[l][i]
	ss = [0.0 for i in range(N)]
	for l in range(N):
		for i in range(N):
			ss[i] += s[l][i]
	for i in range(N):
		ss[i] /= N
	return ss


def simple(N,m,v):
	mm = [[0 for i in range(N)] for ii in range(N)]
	ss = [[0 for i in range(N)] for ii in range(N)]
	for i in range(0,m-N,N):
		vector = list(range(N))
		shuffle(vector)
		for ii,vv in enumerate(vector): #size ii, player vv
			mm[ii][vv] += v(vector[0:ii+1])-v(vector[0:ii])
			ss[ii][vv] += 1
	vector = list(range(N))
	shuffle(vector)
	indices = list(range(N))
	for i in range(m%N):
		ind = indices.pop(int(random()*(len(indices))))
		vector_index = vector.index(ind)
		mm[vector_index][ind] += v(vector[0:vector_index+1])-v(vector[0:vector_index])
		ss[vector_index][ind] += 1
	for i in range(N):
		for ii in range(N):
			mm[i][ii] = mm[i][ii]/ss[i][ii] if ss[i][ii]>0 else 0
	mm = [sum([mm[o][i] for o in range(N)])/N for i in range(N)]
	return mm



def simple_simple(N,m,v):
	s = [[0.0 for i in range(N)] for i in range(N)]
	n = [[1 for i in range(N)] for i in range(N)]
	for l in range(N):
		for i in range(N):
			s[l][i] += gen2(l,i,N)
	for ii in range(m-N*N):
		i = randint(0,N-1)
		l = randint(0,N-1)
		s[l][i] += gen2(l,i,N)
		n[l][i] += 1

	'''for o in range(N):
		for i in range(N):
			print "{}\t".format(n[i][o]),
		print ""'''

	ss = [0.0 for i in range(N)]
	for l in range(N):
		for i in range(N):
			ss[i] += s[l][i]/n[l][i]
	for i in range(N):
		ss[i] /= N
	return ss




#OmegaBig = lambda n,N: sum([1.0/(k**2) for k in range(n,N)])
#PsiBig = lambda n,N: N*sum([1.0/(k**2*(k+1)) for k in range(n,N)])
OmegaBig = lambda n,N: (n+1)*(1-n*1.0/N)*1.0/(n**2)
PsiBig = lambda n,N: (N+1.0-n)/(n**2)
OmegaSmall = lambda n,N: 1.0/n
PsiSmall = lambda n,N: 1.0/n
def burgess_bound(N,ni,Ni,var,d,r):
	onN = [[0 for i in range(N)] for i in range(2)]
	max1 = [[0 for i in range(N)] for i in range(2)]
	var1 = [[0 for i in range(N)] for i in range(2)]
	d1 = [[0 for i in range(N)] for i in range(2)]
	log6r = log(6/r)
	log3r = log(3/r)
	log2r = log(2/r)
	logN = log(N)
	d2 = d*d;
	N2 = N*N;
	N4 = N2*N2;
	for i in range(N):
		for o in range(N):
			OB = OmegaBig(ni[i][o],Ni[o])
			OS = OmegaSmall(ni[i][o],Ni[o])
			PB = PsiBig(ni[i][o],Ni[o])
			PS = PsiSmall(ni[i][o],Ni[o])
			onN[0][i] += PB**2*min(OB,OS)/Ni[o]
			onN[1][i] += PS**2*min(OB,OS)/Ni[o]
			max1[0][i] = max(max1[0][i],PB*min(PB,PS))
			max1[1][i] = max(max1[1][i],PS*min(PB,PS))
			var1[0][i] += PB*(ni[i][o]-1)*var[i][o]/ni[i][o]
			var1[1][i] += PS*(ni[i][o]-1)*var[i][o]/ni[i][o]
			d1[0][i] += OB
			d1[1][i] += OS
	A = [[0 for i in range(N)],[0 for i in range(N)]]
	for i in range(N):
		A[0][i] = sqrt(min((d2*4.0/(17*N2))*log6r*d1[0][i] + 4*log6r*(sqrt((1.0/(2*N2))*var1[0][i] + (log6r+logN)*d2/(8*N4)*onN[0][i] + log3r*d2*max1[0][i]/(4*N2)) + sqrt(log3r*d2*max1[0][i]/(4*N2)))**2, log2r*d2*d1[0][i]/(2*N2)))
		A[1][i] = sqrt(min((d2*4.0/(17*N2))*log6r*d1[1][i] + 4*log6r*(sqrt((1.0/(2*N2))*var1[1][i] + (log6r+logN)*d2/(8*N4)*onN[1][i] + log3r*d2*max1[1][i]/(4*N2)) + sqrt(log3r*d2*max1[1][i]/(4*N2)))**2, log2r*d2*d1[1][i]/(2*N2)))
	return sum([min(A[0][i],A[1][i]) for i in range(N)])


def burgess(N,m,v,d,r=0.5):
	ni = [[0 for i in range(N)] for i in range(N)]
	Ni = [factorial(N-1)/(factorial(N-1-i)*factorial(i)) for i in range(N)]
	listsi = [[] for i in range(N)]
	s = [[0.0 for i in range(N)] for i in range(N)]
	s2 = [[0.0 for i in range(N)] for i in range(N)]
	var = [[0.0 for i in range(N)] for i in range(N)]
	samples = 0
	# seed with minimum initial two samples (if possible) for all strata
	for i in range(N): #player i
		for o in range(N): #coalition size o
			for p in range(2):
				if ni[i][o]<Ni[o] or Ni[o]==-1:
					v = gen3(o,i,N,listsi[i])
					s[i][o]+=v
					s2[i][o]+=v*v
					ni[i][o]+=1
					if ni[i][o]>1:
						var[i][o] = (s2[i][o] - s[i][o]**2*1.0/ni[i][o])/(ni[i][o]-1)
					samples+=1
	advantage = [[0.0 for i in range(N)] for i in range(N)]
	while samples < m:
		print "upto: {}%".format(samples*100.0/m)
		#calculate the bound as it exists:
		bound = burgess_bound(N,ni,Ni,var,d,r)
		#calculate the advantages possible
		for i in range(N):
			for o in range(N):
				if ni[i][o]<Ni[o] or Ni[o]==-1:
					ni[i][o]+=1
					advantage[i][o] = bound-burgess_bound(N,ni,Ni,var,d,r)
					ni[i][o]-=1
				else:
					advantage[i][o]=0
		#detect the sample that maximises advantage
		maxi=0
		maxo=0
		maxadvantage=0
		for i in range(N):
			for o in range(N):
				if advantage[i][o] > maxadvantage:
					maxi = i
					maxo = o
					maxadvantage = advantage[i][o]
		#take the best sample		
		v = gen3(maxo,maxi,N,listsi[maxi])
		s[maxi][maxo]+=v
		s2[maxi][maxo]+=v*v
		ni[maxi][maxo]+=1
		var[maxi][maxo] = (s2[maxi][maxo] - s[maxi][maxo]**2*1.0/ni[maxi][maxo])/(ni[maxi][maxo]-1)
		samples += 1
	#return the calculated shapley value from the samples
	ss = [0.0 for i in range(N)]
	for o in range(N):
		for i in range(N):
			ss[i] += s[i][o]/ni[i][o]
	for i in range(N):
		ss[i] /= N
	return ss,sqrt(bound*log(2.0/r)*0.5)



from tqdm import tqdm
def sample_average_error(f,n,N,m,v):
	data = []
	for i in tqdm(range(n)):
		data.append(f(N,m,v))
	data2 = []
	for i in range(len(data[0])):
		s =0
		for ii in range(len(data)):
			s += data[ii][i]
		data2.append(s*1.0/len(data))
	for i in range(len(data[0])):
		for ii in range(len(data)):
			data[ii][i] = abs(data[ii][i]-data2[i])
	data3 = []
	for i in range(len(data[0])):
		s =0
		for ii in range(len(data)):
			s += data[ii][i]
		data3.append(s*1.0/len(data))
	return 1000*sum(data3)*1.0/len(data3)


def sample_error(f,n,N,m,v):
	data = []
	for i in tqdm(range(n)):
		data.append(f(N,m,v))
	with open("data_out_{}_{}.csv".format(f.func_name,m),"w") as f:
		for d in data:
			f.write("{}\n".format(",".join([str(dd) for dd in d])))





N=len(weights)
for m in [N*N*10,N*N*50,N*N*100,N*N*500,N*N*1000]:
	print "calculating Castro for {}".format(m)
	print sample_error(castro, 40, N,m,v)
for m in [N*N*10,N*N*50,N*N*100,N*N*500,N*N*1000]:
	print "calculating Maleki for {}".format(m)
	print sample_error(maleki, 40, N,m,v)
for m in [N*N*10,N*N*50,N*N*100,N*N*500,N*N*1000]:
	print "calculating Simple for {}".format(m)
	print sample_error(simple, 40, N,m,v)
for m in [N*N*10,N*N*50,N*N*100,N*N*500,N*N*1000]:
	print "calculating SimpleSimple for {}".format(m)
	print sample_error(simple_simple, 40, N,m,v)

#print maleki(N,N*N*N,v)
#print castro(N,N*N*N,v)
#print burgess(N,N*N*4,v,1)
#print simple(N,N*N*N,v)

#m = N*N*10
#print sample_average_error(castro, 40, N,m,v)
#print sample_average_error(maleki, 40, N,m,v)
#print sample_average_error(simple, 40, N,m,v)


