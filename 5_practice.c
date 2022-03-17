#include <stdio.h>	// for printf function
#include <stdlib.h>	// for memory allocation
#include <time.h>	// for time calculation
#include <math.h>	// for sine and cosine functions
int main() {
	// Declare all the variables
	int k, n, N;
	double *x, *yr, *yi, a, b, an, bn, temp, P, c, s;
	time_t t;

	// Input the number N
	scanf("%d", &N);
	
	// Locate the memory for x, yr, yi;
	x = (double*)malloc(N*sizeof(double));
	yr = (double*)malloc(N*sizeof(double));
	yi = (double*)malloc(N*sizeof(double));

	// Initial setting for x, for example, x[k] = k
	for(k = 0; k < N; k++) x[k] = k;
	
	t = clock();	
	// yr[n]+i*yi[n] = sum(exp(-i 2 Pi k n / N)*x[k], k=0..N-1), n=0..N-1 
	a = cos(2*M_PI/N); b = -sin(2*M_PI/N);
	for(n = 0; n < N; n++){
		yr[n] = 0.0;
		yi[n] = 0.0;
		P = 2*M_PI/N;
		an = 1; bn = 0;
		c = 1; s = 0;
		for(k = 0; k < N; k++){
			yr[n] += c*x[k];
			yi[n] -= s*x[k];
			temp = c * an - s * bn;
			s = c * bn + s * an;
			c = temp;
//			double theta = (-2*M_PI*k*n)/N;
//			yr[n] = yr[n] + cos(theta)*x[k];
//			yi[n] = yi[n] + sin(theta)*x[k];
		}
		temp =  an * a - bn * b;
		bn = an * b + bn * a;
		an = temp; 
	}

	// output the results
	t = clock() - t;
	printf("%d ms for discrete Fourier Transform of %d elements\n", t, N);

	// free the memory located by x, yr, yi
	free(x);
	free(yr);
	free(yi);
	
	
	return 100;
}

