#include <stdio.h>	// for printf function
#include <stdlib.h>	// for memory allocation
#include <time.h>	// for time calculation
#include <math.h>	// for sine and cosine functions



void quick_sort(double*, int, int);
int Partition(double*, int, int);
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
	// sorting x;  
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
	
	printf("bubble sort : \n");
	for(k=0;k<N;++k) {
		printf("%f\n",x[k]);
	}
	

	for(k=0;k<N;++k){
		x[k] = y[k];
	}
	printf("------------------------\n");
	t = clock();
	quick_sort(x, 0, N - 1);
	printf("quick sort : \n");
	for(k = 0; k < N; k++){
		printf("%f\n",x[k]);
	}
	// put x[N-1] at the right location
//	p = x[N-1]; n = 0;
//	for(k=0;k<N-1;k++) {
//		// put all elements smaller than p from the beginning
//		if (x[k] < p) {
//			printf("swap %d <-> %d\n",k,n);
//			z = x[n];
//			x[n] = x[k];
//			x[k] = z;
//			n++;
//		}
//	}
//	z = x[n];
//	x[n] = x[N-1];
//	x[N-1] = z;
	
	// x[0] is at the correct location!
//	if(N<10) {
//		printf("y \t\t x\n");
//		for(k=0;k<N;++k) {
//			printf("%f\t%f\n",y[k],x[k]);
//		}
//	}
	// sorting


	t = clock() - t;	
	


	// free the memory located by x, y
	free(x);
	free(y);	
	
	return 100;
}

void quick_sort(double *x, int P, int r){
	double temp;
	if(P < r){
		int q = Partition(x, P, r);
		quick_sort(x, P, q - 1);
		quick_sort(x, q + 1, r);
	}
}

int Partition(double *x, int P, int r){
	double a = x[r];
	int i = P - 1, k;
	for(k = P; k < r; k++){
		if(x[k] <= a){
			i++;
			double temp;
			temp = x[i];
			x[i] = x[k];
			x[k] = temp;
		}
	}
	double temp;
	temp = x[i + 1];
	x[i + 1] = x[r];
	x[r] = temp;
	return (i+1);
}






















