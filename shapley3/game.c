int N=15;
int weights[15] = {45, 41, 27, 26, 25, 21, 13, 13, 12, 12, 11, 11, 10, 10, 10};
double v(unsigned long long a) {
	double vv =0;
	for (int i=0; i < N; i++) {
		vv += weights[i]*(a&1)*1.0/100;
		a = a>>1;
	}
	return vv*vv/2.0;
}
double d=1.19025;

