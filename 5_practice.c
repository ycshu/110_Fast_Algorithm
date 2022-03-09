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
	printf("Enter N = ");
	scanf("%d",&N);
	
	// Locate the memory for x, yr, yi;
	x = (double *) malloc(N*sizeof(double));
	yr = (double *) malloc(N*sizeof(double));
	yi = (double *) malloc(N*sizeof(double));
	
	// Initial setting for x, for example, x[k] = k
	for(int i = 0; i<N; i++){
		x[i] = i;
		yr[i] = 0.0;
		yi[i] = 0.0;
	}
	
	t = clock();	
	// yr[n]+i*yi[n] = sum(exp(-i 2 Pi k n / N)*x[k], k=0..N-1), n=0..N-1 
	//					cos(2 PI k n / N)		i sin(2 PI k n / N)
	double tmp = 2*M_PI/N;
	double theta;
	for (int i = 1; i<N; i++){
 		for(int j = 1; j<N; j++){
			theta = tmp*i*j;
			yr[i] += cos(theta)*x[i];
			yi[i] -= sin(theta)*x[i];
		}
	}

	// output the results
	t = clock() - t;
	printf("%d ms for discrete Fourier Transform of %d elements\n", t, N);

	// free the memory located by x, yr, yi
	free(yi);
	free(yr);
	free(x);

	return 0;
}

