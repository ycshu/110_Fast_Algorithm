#include <stdio.h>	// for printf function
#include <stdlib.h>	// for memory allocation
#include <time.h>	// for time calculation
#include <math.h>	// for sine and cosine functions



int quick_sort(double*, int, int);
void SWAP(double, double);

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
	if(N<10) {
		printf("y \t\t x\n");
		for(k=0;k<N;++k) {
			printf("%f\t%f\n",y[k],x[k]);
		}
	}

	for(k=0;k<N;++k){
		x[k] = y[k];
	}
	t = clock();
		
	// put x[N-1] at the right location
	p = x[N-1]; n = 0;
	for(k=0;k<N-1;k++) {
		// put all elements smaller than p from the beginning
		if (x[k] < p) {
			printf("swap %d <-> %d\n",k,n);
			z = x[n];
			x[n] = x[k];
			x[k] = z;
			n++;
		}
	}
	z = x[n];
	x[n] = x[N-1];
	x[N-1] = z;
	
	// x[0] is at the correct location!
	if(N<10) {
		printf("y \t\t x\n");
		for(k=0;k<N;++k) {
			printf("%f\t%f\n",y[k],x[k]);
		}
	}
	// sorting


	t = clock() - t;	
	


	// free the memory located by x, y
	free(x);
	free(y);	
	
	return 100;
}

int quick_sort(double *x, int L, int R){
	double temp;
	if(L < R){
		double pk = x[L];
		int i = L, j = R + 1;
		while(i < j){
			while(x[i] < pk) i++;
			while(x[j] > pk) j--;
			if(i < j){
				temp = x[i];
				x[i] = x[j];
				x[j] = temp;
			}
		}
		temp = x[L];
		x[L] = x[j];
		x[j] = temp; 
//		SWAP(x[L], x[j]);
		quick_sort(x, L, j - 1);
		quick_sort(x, j + 1, R);
	}
}

//void SWAP(double X, double Y){
//	X ^= Y;
//	Y ^= X;
//	X ^= Y;
//}





















