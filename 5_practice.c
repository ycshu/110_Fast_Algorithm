#include <stdio.h>	// for printf function
#include <stdlib.h>	// for memory allocation
#include <time.h>	// for time calculation
#include <math.h>	// for sine and cosine functions
#include <string.h> // for memset


int main() {
	// Declare all the variables
	int N;
	double *x, *yr, *yi;
	clock_t t;

	// Input the number N
    printf("Please input N = ");
    scanf("%d", &N);    	
	
	// Locate the memory for x, yr, yi;
    x = (double*)malloc(sizeof(double)*N);
    yr = (double*)malloc(sizeof(double)*N);
    yi = (double*)malloc(sizeof(double)*N);

	// Initial setting for x, for example, x[k] = k
    for(int i = 0; i < N; ++i) x[i] = i;	
    
    memset(yr, 0, sizeof(double)*N);
    memset(yi, 0, sizeof(double)*N);
    
    double pi_2_n_div_N, theta;
	t = clock();	
	// yr[n]+i*yi[n] = sum(exp(-i 2 Pi k n / N)*x[k], k=0..N-1), n=0..N-1 
    for(int n = 0; n < N; ++n){
        pi_2_n_div_N = 2*M_PI*n / (double)N;
        for(int k = 0; k < N; ++k){
            theta = pi_2_n_div_N * k;
            yr[n] += cos(theta) * x[k];
            yi[n] -= sin(theta) * x[k];
        }
    }


	// output the results
	t = clock() - t;
	printf("%.0f ms for discrete Fourier Transform of %d elements\n", 1000*(double)t/CLOCKS_PER_SEC, N);

	// free the memory located by x, yr, yi
    free(x);
    free(yr);
    free(yi);
	
	
	return 100;
}

