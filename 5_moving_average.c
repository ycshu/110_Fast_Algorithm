#include <stdio.h>	// for printf function
#include <stdlib.h>	// for memory allocation
#include <time.h>	// for time calculation
#include <math.h>	// for sine and cosine functions
int main() {
	// Declare all the variables
	int k, n, M, N;
	double *x, *y;
	time_t t;

	// Input the number N, M
	printf("N M= ");	
	scanf("%d %d", &N, &M);

	// Locate the memory for x, y;
	x = (double *) malloc(N*sizeof(double));
	y = (double *) malloc(N*sizeof(double));
	// Initial setting for x[k] = sin(2*Pi*k/N), k = 0..N-1
	for(k=0;k<N;++k) {
		x[k] = sin(2*M_PI*k/N);
	}
	
	t = clock();	
	// y[n] = sum(x[k], k=n-M..n+M)/(2M+1), n=0..N-1 
	// zero padding for x[k] if k<0 or k>=N
	for(n=0;n<N;n++){
		y[n] = 0.0;
		for(k=n-M;k<=n+M;++k){
			if(k>=0 || k<N)
				y[n] = y[n] + x[k];
		}
		y[n] = y[n]/(2*M+1);
	}

	// Output the results
	t = clock() - t;
	printf("y[0]=%f, y[1]=%f\n",y[0],y[1]);

	// free the memory located by x, y
	
	
	
	return 100;
}
