#include <stdio.h>	// for printf function
#include <stdlib.h>	// for memory allocation
#include <time.h>	// for time calculation
#include <math.h>	// for sine and cosine functions
int main() {
	// Declare all the variables
	int k, n, N;
	double P, theta, yf=0;
	double c, s, a, b, an, bn, temp;
	double *x, *yr, *yi;
	time_t t;

	// Input the number N
	printf("Pleas input N :\n");
	scanf("%d", &N);
	
	// Locate the memory for x, yr, yi;
	x = (double *) malloc(N*sizeof(double));
	yr = (double *) malloc(N*sizeof(double));
	yi = (double *) malloc(N*sizeof(double));
	// Initial setting for x, for example, x[k] = k
	for(n=0; n<N; n++){
		x[n] = n;
	}
	printf("x is :\n");
	for(n=0; n<N; n++){
		printf("%f ", x[n]);
	}
	printf("\n");
	
	t = clock();	
	// yr[n]+i*yi[n] = sum(exp(-i 2 Pi k n / N)*x[k], k=0..N-1), n=0..N-1 
	//                 (cos(2 pi k n /N) - i sin(2 pi k n /N))
	a = cos(2*M_PI/N); b = -sin(2*M_PI/N);
	for(n=0; n<N; n++) {
		yr[n] = 0.0;
		yi[n] = 0.0;
		an = 1; // cos(2*M_PI*n/N) 
		bn = 0; // sin(2*M_PI*n/N)
		c = 1; s = 0; // cos and sin
		for(k=0; k<N; ++k) {
			yr[n] += c * x[k];
			yi[n] -= s * x[k];
			temp = c * an - s * bn;
			s = c * an + s * an ;
			c = temp;
		}
		an = an * a - bn * b;
		bn = an * b + bn * a;
		an = temp;
	}





	// output the results
	t = clock() - t;
	for(n=0; n<N; n++) {
		printf("idx:%d, real:%f, image:%f \n",n, yr[n], yi[n]);
	}

	printf("%d ms for discrete Fourier Transform of %d elements\n", t, N);

	// free the memory located by x, yr, yi
	free(x);
	free(yr);
	free(yi);
	
	
	return 100;
}

