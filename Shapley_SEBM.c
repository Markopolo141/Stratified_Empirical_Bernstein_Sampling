#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include <time.h>

#include "float.h"
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/ipc.h>
#include <sys/shm.h>

#include "utils.c"

//========== SIMULATION DATAS:
int m=10*N*N;
double r=0.5;
int max_listl = 500000;
int counter = 0;
int max_counter=5;//40;//5;

//========== IPC DATAS:
long cores;
int* workspace;
double* resultspace;
int child;
int parent;

//========== PARENT DATAS:
double* xs;
double* xsum;
double* x2sum;


// the combinations function.
unsigned long long nChoosek( unsigned n, unsigned k ) {
	if (k > n) return 0;
	if (k * 2 > n) k = n-k;
	if (k == 0) return 1;
	unsigned long long result = n;
	for( int i = 2; i <= k; ++i ) {
		result *= (n-i+1);
		result /= i;
	}
	return result;
}

// return 2^b
inline long shift(int b) {
	if (b==0)
		return (long)1;
	return ((long)2)<<(b-1);
}


// function returning a long composed of 'ticks' number of bits randomly flipped in the first 'length' of bits.
// since longs count as coalitions, this function generates a random coalition of size ticks
int* tick_data;
unsigned long long gen_comb(int ticks, int length) {
	for (int i=0; i < ticks; i++) {
		tick_data[i] = rand()%(length-i);
		for (int o=0; o < i; o++) {
			if (tick_data[i]>=tick_data[o]) {
				tick_data[i] += 1;
			}
		}
		for (int o=i-1; o>-1; o--) {
			if (tick_data[o+1] < tick_data[o]) {
				int t = tick_data[o+1];
				tick_data[o+1] = tick_data[o];
				tick_data[o] = t;
			} else {
				break;
			}
		}
	}
	unsigned long long ret = 0;
	for (int i=0; i < ticks; i++) {
		ret += shift(tick_data[i]);
	}
	return ret;
}




// given a combination, find the 'next' combination in iterated and wrapped sequence. done by bit flipping
// WARNING: technically destroyed purity of random sampling of combinations, only to be used sparingly when beginning to sample exhaustively without replacement
unsigned long long iterate_comb(unsigned long long comb, int length) {
	bool trigger = false;
	unsigned long long comb2 = comb;
	unsigned long long ze = 0;
	for (int i=0; i< length; i++) {
		if (!trigger) {
			if ((comb2 & 1)==0) {
				trigger = true;
			} else {
				ze += 1;
			}
		} else {
			if ((comb2 & 1)==1) {
				return comb ^ (((long)3)<<(i-1)) ^ ((((ze>>2)^comb)<<(sizeof(unsigned long long)*8-i))>>(sizeof(unsigned long long)*8-i));
			}
		}
		ze = ze << 1;
		comb2 = comb2 >> 1;
	}
	return ze>>1;
}



// tries to generate a unique coalition that hasnt been seen before
// after 10 tries it reverts to 'iterate_comb' for expediency
unsigned long long gen_new_comb(int ticks, int length, unsigned long long* combs, int* combs_length, int combs_max_length) {
	unsigned long long v;
	bool in_list= true;
	v = gen_comb(ticks, length);
	int tries=0;
	while (in_list) {
		in_list= false;
		for (int i = 0; i < *combs_length; i++) {
			if (combs[i]==v) {
				in_list= true;
				break;
			}
		}
		if (tries<10) {
			if (in_list) {
				v = gen_comb(ticks, length);
			}
			tries++;
		} else {
			if (in_list) {
				v = iterate_comb(v,length);
			}
		}
	}
	if (*combs_length < combs_max_length) {
		combs[*combs_length] = v;
		*combs_length += 1;
	}
	return v;
}



// samples a random coalition of size l without replacement and 
// calculates the marginal contribution of adding player i
double gen3(int l, int i, int N, unsigned long long* combs, int* combs_length, int combs_max_length) {
	unsigned long long vv = gen_new_comb(l, N-1, combs, combs_length, combs_max_length);
	unsigned long long v1 = vv%shift(i);
	v1 += 2*(vv-v1);
	return v(v1+shift(i))-v(v1);
}






// the SEBB (assuming sampling without replacement). the bound which the SEBM seeks to iteratively minimize.
// N is the number of strata
// ni is the amount that each strata is sampled
// Ni is the sizes of the strata
// var is the sample variances of the strata
// d is the population data width
// r is the confidence level
inline double OmegaBig(long n, long N) {
	return (n+1)*(1-n*1.0/N)*1.0/(n*n);
}
inline double PsiBig(long n, long N) {
	return (N+1.0-n)/(n*n);
}
inline double OmegaSmall(long n, long N) {
	return 1.0/n;
}
inline double PsiSmall(long n, long N) {
	return 1.0/n;
}
double* onN;
double* max1;
double* var1;
double* d1;
inline double burgess_bound(int N, long* ni, long* Ni, double* var, double d, double r) {
	for (int i=0; i<N*2; i++) {
		onN[i] = 0;
		max1[i] = 0;
		var1[i] = 0;
		d1[i] = 0;
	}
	double d2 = d*d;
	double tau;
	double log6N2ro2 = log(6*N*(N-2)/r);
	double log3r = log(3/r);
	double temp,PB,PS,OB,OS;
	for (int i=0; i<N; i++) {
		for (int o=0; o<N; o++) {
			tau = 1.0/(N*(N-2));
			PS = PsiSmall(ni[i*N+o],Ni[o]);
			OS = OmegaSmall(ni[i*N+o],Ni[o]);
			PB = PsiBig(ni[i*N+o],Ni[o]);
			OB = OmegaBig(ni[i*N+o],Ni[o]);
			onN[0*N+i] += PB*fmin(OB,OS)*tau;
			onN[1*N+i] += PS*fmin(OB,OS)*tau;
			temp = PB*fmin(PB,PS)*tau;
			if (temp > max1[0*N+i]) max1[0*N+i]=temp;
			temp = PS*fmin(PB,PS)*tau;
			if (temp > max1[1*N+i]) max1[1*N+i]=temp;
			var1[0*N+i] += PB*(ni[i*N+o]-1)*var[i*N+o]*tau/ni[i*N+o];
			var1[1*N+i] += PS*(ni[i*N+o]-1)*var[i*N+o]*tau/ni[i*N+o];
			d1[0*N+i] += OB*tau;
			d1[1*N+i] += OS*tau;
		}
	}
	double A = 0;
	for (int i=0; i<N; i++) {
		A += fmin(
			sqrt(
				(d2*4.0/(17))*d1[0*N+i] + 
				pow(
					sqrt(
						2*var1[0*N+i] + log6N2ro2*d2*onN[0*N+i] + log3r*d2*max1[0*N+i]
					) + 
					sqrt(log3r*d2*max1[0*N+i])
				,2)
			),
			sqrt(
				(d2*4.0/(17))*d1[1*N+i] + 
				pow(
					sqrt(
						2*var1[1*N+i] + log6N2ro2*d2*onN[1*N+i] + log3r*d2*max1[1*N+i]
					) + 
					sqrt(log3r*d2*max1[1*N+i])
				,2)
			)
		);
	}
	return A;
}

// the SEBM (assuming sampling without replacement)
// N is the number of strata
// m is the sample budget
// d is the data width
// r is the confidence level
double* burgess(int N, int m, double d, double r) {
	long* ni = (long*)calloc(sizeof(long),N*N);
	long* Ni = (long*)calloc(sizeof(long),N);
	for (int i=0; i<N; i++) {
		Ni[i] = nChoosek(N-1,i);
	}
	unsigned long long**listsi = (unsigned long long**)malloc(sizeof(unsigned long long*)*N*N);
	for (int i=0; i<N*N; i++)
		listsi[i] = (unsigned long long*)malloc(sizeof(unsigned long long)*max_listl);
	int* listsil = (int*)calloc(sizeof(int),N*N);
	double* s = (double*)calloc(sizeof(double),N*N);
	double* s2 = (double*)calloc(sizeof(double),N*N);
	double* var = (double*)calloc(sizeof(double),N*N);
	double* advantage = (double*)calloc(sizeof(double),N*N);
	double bound;

	int samples = 0;
	// seed with minimum initial two samples (if possible) for all strata
	for (int i=0; i<N; i++) { //player i
		for (int o=0; o<N; o++) { //coalition size o
			for (int p=0; p<2; p++) {
				if (ni[i*N+o]<Ni[o]) {
					double v = gen3(o,i,N,listsi[i*N+o],&(listsil[i*N+o]),max_listl);
					s[i*N+o]+=v;
					s2[i*N+o]+=v*v;
					ni[i*N+o]+=1;
					if (ni[i*N+o]>1) {
						var[i*N+o] = (s2[i*N+o] - s[i*N+o]*s[i*N+o]*1.0/ni[i*N+o])/(ni[i*N+o]-1);
					}
					samples+=1;
				}
			}
		}
	}
	while (samples < m) {
		//calculate the bound as it exists:
		bound = burgess_bound(N,ni,Ni,var,d,r);
		//calculate the advantages possible
		for (int i=0; i<N; i++) {
			for (int o=0; o<N; o++) {
				if (ni[i*N+o]<Ni[o]) {
					ni[i*N+o]+=1;
					advantage[i*N+o] = bound-burgess_bound(N,ni,Ni,var,d,r);
					ni[i*N+o]-=1;
					//advantage[i*N+o] = 1.0/ni[i*N+o];
				} else {
					advantage[i*N+o]=-INFINITY;
				}
			}
		}
		//detect the sample that maximises advantage
		int maxi=0;
		int maxo=0;
		double maxadvantage=-INFINITY;
		for (int i=0; i<N; i++) {
			for (int o=0; o<N; o++) {
				if (advantage[i*N+o] > maxadvantage) {
					maxi = i;
					maxo = o;
					maxadvantage = advantage[i*N+o];
				}
			}
		}
		if ((samples^((samples >>6)<<6))==0) printf("%i\t%i\t%i\n",samples,maxi,maxo);
		//take the best sample
		double v = gen3(maxo,maxi,N,listsi[maxi*N+maxo],&(listsil[maxi*N+maxo]),max_listl);
		s[maxi*N+maxo]+=v;
		s2[maxi*N+maxo]+=v*v;
		ni[maxi*N+maxo]+=1;
		var[maxi*N+maxo] = (s2[maxi*N+maxo] - s[maxi*N+maxo]*s[maxi*N+maxo]*1.0/ni[maxi*N+maxo])/(ni[maxi*N+maxo]-1);
		samples += 1;
	}
	//return the calculated shapley value from the samples
	double* ss = (double*)calloc(sizeof(double),N);
	for (int o=0; o<N; o++)
		for (int i=0; i<N; i++)
			ss[i] += s[i*N+o]/ni[i*N+o];
	for (int i=0; i<N; i++)
		ss[i] /= N;

	//free resources
	free(ni);
	free(Ni);
	for (int i=0; i<N*N; i++)
		free(listsi[i]);
	free(listsi);
	free(listsil);
	free(s);
	free(s2);
	free(var);
	free(advantage);
	
	return ss;
}



// output computed data to file
void file_output_data() {
	char* filename;
	filename = (char*)calloc(sizeof(char),150);
	sprintf(filename, "data_out_Burgess_%i.csv", m);
	FILE* f = fopen(filename,"w");
	free(filename);
	
	int n = max_counter*cores;
	for (int i=0; i<n; i++) {
		for (int o=0; o<N; o++) {
			fprintf(f,"%f, ",xs[i*N+o]);
		}
		fprintf(f,"\n");
	}
	fclose(f);
}


// output computed data to screen
void output_data() {
	int n = max_counter*cores;
	double v=0;
	double avg_abs_err=0;
	for (int o=0; o<N; o++) {
		for (int i=0; i<n; i++) {
			avg_abs_err += fabs(xs[i*N+o]-xsum[o]/n);
		}
	}
	printf("average result:\n");
	for (int o=0; o<N; o++) {
		printf("%f ",xsum[o]/n);
		v += (x2sum[o] - xsum[o]*xsum[o]*1.0/n)/(n-1);
	}
	printf("\n");
	printf("variance in results: %f\n",v);
	printf("avg err: %f\n", avg_abs_err/(N*n));
	printf("%i total samples\n", n);
}


// the main loop executed by each of the program multiprocesses
// the parent coordinates the children
void loop() {
	int* command = &workspace[0];
	int* state = &workspace[child+1];
	int all_done;
	while (*command != -1) {
		if ((*state == 0) && (*command == 1)) {
			*state = 1;
			double* v = burgess(N, m, d, r);
			for (int i=0; i<N; i++) {
				resultspace[child*N+i] = v[i];
			}
			free(v);
			*state = 2;
		}
		if (parent) {
			all_done = 1;
			for (int i = 0; i < cores; i++) {
				if (workspace[i+1]!=2)
					all_done=0;
			}
			if (all_done==1) {
				*command=0;
				for (int i=0; i<cores; i++) {
					for (int o=0; o<N; o++) {
						xs[(counter*cores+i)*N+o] = resultspace[i*N+o];
						xsum[o] += resultspace[i*N+o];
						x2sum[o] += resultspace[i*N+o]*resultspace[i*N+o];
					}
				}
				for (int i=1; i<cores+1; i++)
					workspace[i]=0;
				counter++;
				printf("%.2f%% complete\n",counter*100.0/max_counter);
				if (counter<max_counter) {
					*command=1;
				} else {
					*command=-1;
				}
			}
		}
	}
	if (parent) {
		file_output_data();
		while(wait(NULL)>0);
	}
}


// allocate needed memory
void memory_allocate(char* key) {
	tick_data = (int*)calloc(sizeof(int),N);
	onN = (double*)calloc(sizeof(double),N*2);
	max1 = (double*)calloc(sizeof(double),N*2);
	var1 = (double*)calloc(sizeof(double),N*2);
	d1 = (double*)calloc(sizeof(double),N*2);
	cores = sysconf(_SC_NPROCESSORS_ONLN);//1;
	workspace = (int*)salloc(key,'z',sizeof(int)*(cores+1),0666|IPC_CREAT);
	resultspace = (double*)salloc(key,'z',sizeof(double)*cores*N,0666|IPC_CREAT);
	workspace[0]=1;
}
// free all memory
void memory_free() {
	free(tick_data);
	free(onN);
	free(max1);
	free(var1);
	free(d1);
	sfree(workspace);
	sfree(resultspace);
}
// additional allocation of memory for parent
void parent_memory_allocate() {
	xs = (double*)calloc(sizeof(double),N*cores*max_counter);
	xsum = (double*)calloc(sizeof(double),N);
	x2sum = (double*)calloc(sizeof(double),N);
}
// additional frees of memory for parent
void parent_memory_free() {
	free(xs);
	free(xsum);
	free(x2sum);
}

// main function, forks as many children as there are processor cores and executes main loop
int main(int argc, char *argv[]) {
	memory_allocate(argv[0]);
	m = N*N*atoi(argv[1]);
	pid_t p = 999;
	child=0;
	while ((child < cores-1) && (p > 0)) {
		p = fork();
		if (p==-1) {
			perror("fork");
			exit(1);
		} else if (p>0) {
			child++;
		}
	}
	parent = (p==0) ? 0 : 1;
	clock_t t;
	if (parent) {
		t = clock();
		parent_memory_allocate();
	}
	srand (time(NULL)*(child+1));
	loop();
	memory_free();
	if (parent) {
		parent_memory_free();
		printf("Execution time: %fs\n",(1.0*clock()-t)/CLOCKS_PER_SEC);
	}
	return 0;
}


