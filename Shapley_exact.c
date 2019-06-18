#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include <time.h>

// return 2^b
inline long shift(int b) {
	if (b==0)
		return (long)1;
	return ((long)2)<<(b-1);
}


// calculate marginal contribution of player i, to coalition identified by bits in vv
double marginal(int i, long vv) {
	unsigned long long v1 = vv%shift(i);
	v1 += 2*(vv-v1);
	return v(v1+shift(i))-v(v1);
}

// generate and output all marginal contributions in the game to csv
int main(int argc, char *argv[]) {
	clock_t t;
	t = clock();
	FILE* f = fopen("marginals.csv","w");
	long len = shift(N-1);
	long lend = len/100;
	for (long i=0; i<len; i++) {
		if (i%lend==0)
			printf("percent complete:\t%f\n",i*100.0/len);
		int o=0;
		long ii=i;
		while (ii) {
			o += ii&1;
			ii >>= 1;
		}
		fprintf(f,"%i,",o);
		for (int oo = 0; oo < N; oo++) {
			fprintf(f,"%f,",marginal(oo,i));
		}
		fprintf(f,"\n");
	}
	fclose(f);
	printf("Execution time: %fs\n",(1.0*clock()-t)/CLOCKS_PER_SEC);
	return 0;
}


