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
	printf("Input the number N:\n");
	scanf("%d",&N);

	// Locate the memory for x, yr, yi;
	x = (double *)malloc(N * sizeof(double));
	yr = (double *)malloc(N * sizeof(double)); 
	yi = (double *)malloc(N * sizeof(double));

	// Initial setting for x, for example, x[k] = k
	for ( k = 0; k < N; k++)
	{
		x[k] = k;
		yr[k] = 0;
		yi[k] = 0;
	}
	
	
	
	
	t = clock();	
	// yr[n]+i*yi[n] = sum(exp(-i 2 Pi k n / N)*x[k], k=0..N-1), n=0..N-1 
	for ( n = 0; n < N; n++)
	{
		for ( k = 0; k < N; k++)
		{
			yr[n] += cos(2 * M_PI * k * n / N ) * x[k];
			yi[n] -= sin(2 * M_PI * k * n / N ) * x[k];
		}
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

