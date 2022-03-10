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
	scanf("%d", &N);
	
	// Locate the memory for x, yr, yi;

	// Initial setting for x, for example, x[k] = k
	
	
	
	
	t = clock();	
	// yr[n]+i*yi[n] = sum(exp(-i 2 Pi k n / N)*x[k], k=0..N-1), n=0..N-1 







	// output the results
	t = clock() - t;
	printf("%d ms for discrete Fourier Transform of %d elements\n", t, N);

	// free the memory located by x, yr, yi
	
	
	
	return 100;
}

