#include <stdio.h>	// for printf function
#include <stdlib.h>	// for memory allocation
#include <time.h>	// for time calculation
#include <math.h>	// for sine and cosine functions

int main() {
	// Declare all the variables
	int k, n, N, i;
	double *x, *yr, *yi;
	double *yr1,*yi1;
	double an,bn,a0,b0,a,b;
	double temp,dr,di;
	time_t t,t1;
	
	// Input the number N
	printf("Input the number N :");
	scanf("%d",&N);

	// Locate the memory for x, yr, yi;
	x = (double *) malloc(N*sizeof(double));
	yr = (double *) malloc(N*sizeof(double));
	yi = (double *) malloc(N*sizeof(double));
	yr1 = (double *) malloc(N*sizeof(double));
	yi1 = (double *) malloc(N*sizeof(double));
	
	// Initial setting for x, for example, x[k] = k
	for(k=0;k<N;k++){
		x[k] = k;
	}
	
	t = clock();	
	// yr[n]+i*yi[n] = sum(exp(-i 2 Pi k n / N)*x[k], k=0..N-1), n=0..N-1 
	// Method 1 : ª½±µ­pºâ 
	for(n=0;n<N;n++){
		yr[n] = 0.0;
		yi[n] = 0.0;
		for(k=0;k<N;k++){
			yr[n] += cos(2.0*M_PI*k*n/N)*x[k];
			yi[n] -= sin(2.0*M_PI*k*n/N)*x[k];
		} 
		
	}
	// output the results
	t = clock() - t;
	
	// Method 2
	a = cos(2.0*M_PI*k/N);
	b = -sin(2.0*M_PI*k/N);
	an = 1.0;
	bn = 0.0;
	t1 = clock();
	for(n=0;n<N;n++){
		yr1[n] = 0.0;
		yi1[n] = 0.0;
		a0 = 1.0;
		b0 = 0.0;
		for(k=0;k<N;k++){
			yr1[n] += a0*x[k];
			yi1[n] -= b0*x[k];
			temp = an*a0-bn*b0;
			b0 = an*b0+bn*a0;
			a0 = temp;
		}
		temp = an*a - bn*b;
		bn = an*b + bn*a;  
		an = temp;
	}
	// output the results
	t1 = clock() - t1;
	printf("Method 1 : %e ms for discrete Fourier Transform of %d elements\n", t, N);
	printf("Method 2 : %e ms for discrete Fourier Transform of %d elements\n", t1, N);
	dr = 0.0;
	di = 0.0;
	for(i=0;i<N;i++){
		dr += pow(yr1[i]-yr[i],2);
		di += pow(yi1[i]-yi[i],2);
	}
	printf("the difference of real number for two methods: %e\n",sqrt(dr));
	printf("the difference of imaginary numberfor two methods: %e\n",sqrt(di));
	// free the memory located by x, yr, yi	
	free(yr);
	free(yi);
	free(yr1);
	free(yi1);
	free(x);
	
	return 100;
}

