int N=15;
int weights[15] = {1,1,2,2,2,3,4,5,5,5,7,8,8,8,10};
double v(unsigned long long a) {
	double vv =0;
	for (int i=0; i < N; i++) {
		int t = weights[i]*(a&1);
		if (t>vv)
			vv=t;
		a = a>>1;
	}
	return vv;
}
double d=10.0;

