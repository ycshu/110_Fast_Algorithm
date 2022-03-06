#include <stdio.h>	// for printf function
#include <stdlib.h>	// for memory allocation
#include <time.h>	// for time calculation
#include <math.h>	// for sine and cosine functions
int main() {
	// Declare all the variables
	int k, n, N;
	double *x, *yr, *yi;
	time_t t;

	// Input the number N
	printf("Please input an integer : ");
	scanf("%d",&N);
	
	// Locate the memory for x, yr, yi;
	x = (double *) malloc(N*sizeof(double));
	yr = (double *) malloc(N*sizeof(double));
	yi = (double *) malloc(N*sizeof(double));
	
	// Initial setting for x, for example, x[k] = k
	for(k=0;k<N;k++){
		x[k] = k;
	}
	
	t = clock();	
	// yr[n]+i*yi[n] = sum(exp(-i 2 Pi k n / N)*x[k], k=0..N-1), n=0..N-1 
	for(n=0;n<N;n++){
		yr[n] = cos(0)*x[0];
		yi[n] = (-1)*sin(0)*x[0];
		for(k=1;k<N;k++){
			yr[n] = yr[n] + cos(2*M_PI*n*k/N)*x[k];
			yi[n] = yi[n] - sin(2*M_PI*n*k/N)*x[k];
		}
	}

	for(n=0;n<N;n++){
		printf("n = %d, yr[%d]+i*yi[%d] = %f + i * %f\n",n, n, n, yr[n], yi[n]);
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

