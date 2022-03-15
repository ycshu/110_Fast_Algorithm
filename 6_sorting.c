#include <stdio.h>	// for printf function
#include <stdlib.h>	// for memory allocation
#include <time.h>	// for time calculation
#include <math.h>	// for sine and cosine functions
int quick_sort(double *x, int L, int R);
double quick_sort_median(double *x, int L, int R, int M1,int bl);
double quick_sort_median1(double *x, int L, int R, int M1);
int main() {
	// Declare all the variables
	int k, m, n, N;
	double *x, *y, z, p , med, med1, med2,q1,q3;
	time_t t,t1,t2;

	// Input the number N
	printf("Input N: ");
	scanf("%d",&N);
	
	// Locate the memory for x and y;
	x = (double *) malloc(N*sizeof(double));
	y = (double *) malloc(N*sizeof(double));
	// Initial setting for x, for example, x[k] = 1.0*rand()/RAND_MAX

	srand( time(NULL) );
	for(k=0;k<N;k++){
		x[k] = y[k] = 1.0*rand()/RAND_MAX;
	}	
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
	for(k=0;k<N;++k){
		x[k] = y[k];
	}
	
	quick_sort(x,0,N);
	
	if(N<=20) {
		printf("y \t\t x\n");
		for(k=0;k<N;++k) {
			printf("%d %f\t%f\n",k,y[k],x[k]);
		}
	}
	
	if( N % 2 != 0){
		med = x[N/2];
	}
	else{
		med = (x[(N/2)-1]+x[N/2])/2.0;
	}

	for(k=0;k<N;++k){
		x[k] = y[k];
	}
	t = clock();
	if( N % 2 != 0){
		med1 = quick_sort_median(x,0,N,(N-1)/2,1);
	}
	else{
		med1 = quick_sort_median(x,0,N,(N-1)/2,0);
	}
	t = clock() - t;
	
	t2 = clock();
	if( (((N-1)/2)%2 != 0 && N%2 != 0) || (N % 2 ==0 && (N/2)%2 != 0)){
		q1 = quick_sort_median(x,0,N/2,((N/2)-1)/2,1); 
		q3 = quick_sort_median(x,N/2+1,N,3*N/4,1); 
	}
	else{
		q1 = quick_sort_median(x,0,N/2,((N/2)-1)/2,0);
		q3 = quick_sort_median(x,(N-1)/2+1,N,(3*N-1)/4,0);
	}
	t2 = t2-clock();
	
	for(k=0;k<N;++k){
		x[k] = y[k];
	}
	
	t1 = clock();
	if( N % 2 != 0){
	med2 = quick_sort_median1(x,0,N,(N-1)/2);
	  }
	else{
	  	med2 = quick_sort_median1(x,0,N,(N-1)/2);
	   	med2 = med2 + quick_sort_median1(x,0,N,N/2);
	   	med2 = med2/2.0;
	}
	t1 = t1-clock();

	printf("quick_sort_median spends : %e s\n", N, 1.0*t/CLOCKS_PER_SEC);
	printf("quick_sort_median1 spends : %e s\n", N, 1.0*t1/CLOCKS_PER_SEC);
	printf("The error with quick_sort_median is %f \n",med-med1);
	printf("The error with quick_sort_median1 is %f \n",med-med2);
	printf("q1 = %f q3 = %f time : %e \n", q1,q3,1.0*t2/CLOCKS_PER_SEC);
	// free the memory located by x, y
	free(x);
	free(y);	

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

double quick_sort_median(double *x, int L, int R, int M1,int bl){
	int k, n,m;
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

	if(n > M1){
		quick_sort_median(x,L,n,M1,bl);
	}
	else if(n < M1){
		quick_sort_median(x,n+1,R,M1,bl);
	}
	else{
		if(bl){
			return x[n];
		}
		else{
			quick_sort_median(x,n+1,R,M1+1,bl);
			return (x[n]+x[n+1])/2.0;
		}
	}
}

double quick_sort_median1(double *x, int L, int R, int M1){
	int k, n,m;
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

	if(n > M1){
		quick_sort_median1(x,L,n,M1);
	}
	else if(n < M1){
		quick_sort_median1(x,n+1,R,M1);
	}
	else{
		return x[n];
	}
	
}
