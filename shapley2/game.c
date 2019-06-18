int N=15;
int weights[15] = {1,3,3,6,12,16,17,19,19,19,21,22,23,24,29};
int sum_weights = 234;
double v(unsigned long long a) {
	int vv =0;
	for (int i=0; i < N; i++) {
		vv += weights[i]*(a&1);
		a = a>>1;
	}
	if (vv>(sum_weights/2))
		return 1.0;
	return 0.0;
}
double d=1.0;

