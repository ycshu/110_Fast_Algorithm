#include <stdio.h>	// for printf function
#include <stdlib.h>	// for memory allocation
#include <time.h>	// for time calculation
#include <math.h>	// for sine and cosine functions
int main() {
	// Declare all the variables
	int k, m, n, N;
	double *x, *y, z, p;
	time_t t;

	// Input the number N
	printf("Input N: ");
	scanf("%d",&N);
	
	// Locate the memory for x and y;
	x = (double *) malloc(N*sizeof(double));
	y = (double *) malloc(N*sizeof(double));

	// Initial setting for x, for example, x[k] = 1.0*rand()/RAND_MAX
	srand( time(NULL) );
	for(k=0;k<N;++k){
		x[k] = y[k] = 1.0*rand()/RAND_MAX;
	}
	
	t = clock();	
	// sorting x (bubble sort) from small to large
	for(n=0;n<N;++n) {
		for(k=n+1;k<N;++k) {
			if (x[n] > x[k]) {
				z = x[n];
				x[n] = x[k];
				x[k] = z;
			}
		}
	}
	t = clock() - t;

	// print y, x, and time
	printf("Sorting %d elements: %f s\n", N, 1.0*t/CLOCKS_PER_SEC);
	if(N<10) {
		printf("y \t\t x\n");
		for(k=0;k<N;++k) {
			printf("%f\t%f\n",y[k],x[k]);
		}
	}

	// reset x
	for(k=0;k<N;++k){
		x[k] = y[k];
	}

	// quick sort
	t = clock();
	quick_sort(x, 0, N);
	t = clock() - t;	
	
	// free the memory located by x, y
	free(x);
	free(y);	
	
	return 0;
}

int quick_sort(double *x, int L, int R){
	//sort from L to R-1

	//set the break situation
	if (R<L) return 0; 

	//declare variables for pivot, tmp(z), index k n
	double p, z;
	int k, n;

	//set pivot and starting index, then start the process
	p = X[R-1]; n = L;
	for (k = L; k < R-1; k++){
		if(x[k] < p){
			z = x[n];
			x[n] = x[k];
			x[k] = z;
			n++;
		}
	}

	//swap the pivot to the correct location
	z = x[n];
	x[n] = x[R-1];
	x[R-1] = z;

	//do the quick sort for subsequence
	quick_sort(x,L,n);
	quick_sort(x,n+1,R);
	
	return 0;
}