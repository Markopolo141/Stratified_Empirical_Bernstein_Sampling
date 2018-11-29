#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include <time.h>


void printbin(long a) {
	while (a) {
		printf("%i",a&1);
		a >>= 1;
	}
	printf("\n");
}


inline long shift(int b) {
	if (b==0)
		return (long)1;
	return ((long)2)<<(b-1);
}

int N=15;
int weights[15] = {45, 41, 27, 26, 25, 21, 13, 13, 12, 12, 11, 11, 10, 10, 10};
double v(unsigned long long a) {
	double vv =0;
	for (int i=0; i < N; i++) {
		vv += weights[i]*(a&1)*1.0/50;
		a = a>>1;
	}
	vv = vv*vv;
	vv = vv - (int)vv;
	return vv;
}

double gen3(int i, long vv) {
	//printbin(vv);
	unsigned long long v1 = vv%shift(i);
	//printbin(v1);
	v1 += 2*(vv-v1);
	//printbin(v1);
	//printbin(v1+shift(i));
	return v(v1+shift(i))-v(v1);
}





void loop() {
	FILE* f = fopen("data_out.csv","w");
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
			fprintf(f,"%f,",gen3(oo,i));
		}
		fprintf(f,"\n");
	}
	fclose(f);
	return;
}


void memory_allocate(char* key) {}
void memory_free() {}

int main(int argc, char *argv[]) {
	memory_allocate(argv[0]);
	clock_t t;
	t = clock();
	srand (time(NULL));
	loop();
	//printf("%f\n", gen3(1,16384));
	memory_free();
	printf("Execution time: %fs\n",(1.0*clock()-t)/CLOCKS_PER_SEC);
	return 0;
}


