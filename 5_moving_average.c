#include <stdio.h>	// for printf function
#include <stdlib.h>	// for memory allocation
#include <time.h>	// for time calculation
#include <math.h>	// for sine and cosine functions
int main() {
	// Declare all the variables
	int k, n, M, N;
	double *x, *y;
	time_t t;

	// Input the number N
	
	
	// Locate the memory for x, y;

	
	// Initial setting for x[k] = sin(2*Pi*k/N), k = 0..N-1
	
	
	
	t = clock();	
	// y[n] = sum(x[k], k=n-M..n+M)/(2M+1), n=0..N-1 
	// zero padding for x[k] if k<0 or k>=N


	// Output the results
	t = clock() - t;
	printf("\n");

	// free the memory located by x, y
	
	
	
	return 100;
}
