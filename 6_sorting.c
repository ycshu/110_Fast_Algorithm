#include <stdio.h>	// for printf function
#include <stdlib.h>	// for memory allocation
#include <time.h>	// for time calculation
#include <math.h>	// for sine and cosine functions
int quick_sort(double *x, int L, int R);
double quick_sort_median(double *x1, int L, int R, int M1,int bl);
int main() {
	// Declare all the variables
	int k, m, n, N;
	double *x, *y,*x1, z, p , med, med1;
	time_t t;

	// Input the number N
	printf("Input N: ");
	scanf("%d",&N);
	
	// Locate the memory for x and y;
	x = (double *) malloc(N*sizeof(double));
	y = (double *) malloc(N*sizeof(double));
	x1 = (double *) malloc(N*sizeof(double));
	// Initial setting for x, for example, x[k] = 1.0*rand()/RAND_MAX
	srand( time(NULL) );
	for(k=0;k<N;k++){
		x[k] = y[k] = x1[k] = 1.0*rand()/RAND_MAX;
	}
	
	t = clock();	
	// sorting x;  
	for(n=0;n<N;n++) {
		for(k=n+1;k<N;k++) {
			if (x[n] > x[k]) {
				z = x[n];
				x[n] = x[k];
				x[k] = z;
			}
		}
	}
	t = clock() - t;
	// print y, x, and time
	printf("Sorting %d elements: %e s\n", N, 1.0*t/CLOCKS_PER_SEC);
	if(N<=10) {
		printf("y \t\t x\n");
		for(k=0;k<N;++k) {
			printf("%f\t%f\n",y[k],x[k]);
		}
	}

	for(k=0;k<N;++k){
		x[k] = y[k];
	}
	
	// quick sort x 
	t = clock();
	quick_sort(x,0,N);
	t = clock() - t;
	// x[0] is at the correct location!
	printf("Quick Sorting %d elements: %e s\n", N, 1.0*t/CLOCKS_PER_SEC);
	if(N<=10) {
		printf("y \t\t x\n");
		for(k=0;k<N;++k) {
			printf("%f\t%f\n",y[k],x[k]);
		}
	}
	
	t = clock();
	if( N % 2 != 0){
		med1 = quick_sort_median(x1,0,N,(N-1)/2,1);
	}
	else{
		med1 = quick_sort_median(x1,0,N,(N-1)/2,0);
	}
	t = clock() - t;
	if( N % 2 != 0){
		med = x[N/2];
	}
	else{
		med = (x[(N/2)-1]+x[N/2])/2.0;
	}
	printf("quick_sort_median spends : %e s\n", N, 1.0*t/CLOCKS_PER_SEC);
	printf("The error is %f \n",med-med1);

	// free the memory located by x, y
	free(x);
	free(y);	
	free(x1);
	return 100;
}

int quick_sort(double *x, int L, int R){
	int k, n;
	double z, p;
	if(L >= R){
		return 0;
	}
	p = x[R-1];
	n = L;
	
	for(k=L;k<R-1;k++) {
		// put all elements smaller than p from the beginning
		if (x[k] < p) {
			z = x[n];
			x[n] = x[k];
			x[k] = z;
			n++;
		}
	}
	z = x[n];
	x[n] = x[R-1];
	x[R-1] = z;

	quick_sort(x,L,n);
	quick_sort(x,n+1,R);
	return 0;
}

double quick_sort_median(double *x1, int L, int R, int M1,int bl){
	int k, n,m;
	double z, p;
	if(L >= R){
		return 0;
	}

	p = x1[R-1];
	n = L;
	for(k=L;k<R-1;k++) {
		// put all elements smaller than p from the beginning
		if (x1[k] < p) {
			z = x1[n];
			x1[n] = x1[k];
			x1[k] = z;
			n++;
		}
	}
	z = x1[n];
	x1[n] = x1[R-1];
	x1[R-1] = z;

	if(n > M1){
		quick_sort_median(x1,L,n,M1,bl);
	}
	else if(n < M1){
		quick_sort_median(x1,n+1,R,M1,bl);
	}
	else{
		if(bl){
			return x1[n];
		}
		else{
			quick_sort_median(x1,n+1,R,M1+1,bl);
			return (x1[n]+x1[n+1])/2.0;
		}
	}
}
