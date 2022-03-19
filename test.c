#include <stdio.h>	// for printf function
#include <stdlib.h>	// for memory allocation
#include <time.h>	// for time calculation
#include <math.h>	// for sine and cosine functions



void quick_sort(double*, int, int);
int Partition(double*, int, int);
double Select(double*, int, int, int);
double Find_Med(double*, int, int);
double QD(double*, int, int);
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
	printf("origin : \n");
	for(k = 0; k < N; k++){
		printf("%f\n", x[k]);
	}
	printf("----------------------\n");
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
	printf("Bubble sort %d elements: %f s\n", N, 1.0*t/CLOCKS_PER_SEC);
	
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
	// sorting
	t = clock() - t;
	
	
	printf("Quick sort %d elements: %f s\n", N, 1.0*t/CLOCKS_PER_SEC); 
	
	for(k=0;k<N;++k){
		x[k] = y[k];
	}
	
	t = clock();
	double Med = Find_Med(x, 0, N - 1);
	printf("-----------------------\n");
	printf("Median = %f\n", Med); 
	t = clock() - t;	
	printf("Find median in %d elements: %f s\n", N, 1.0*t/CLOCKS_PER_SEC);
	
	printf("-----------------------\n");
	t = clock();
	printf("Quartile Deviation = %f\n", QD(x, 0, N - 1));
	t = clock() - t;
	printf("Find Quartile Deviation in %d elements: %f s\n", N, 1.0*t/CLOCKS_PER_SEC);
		


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
	return (i + 1);
}

double Select(double *x, int P, int r, int i){
	if(P < r){
		int q = Partition(x, P, r);
		int k = q - P + 1;
		if(i == k) return x[q];
		else if(i < k) return Select(x, P, q - 1, i);
		else return Select(x, q + 1, r, i - k);
	}
	else{
		return x[r];
	}
}

double Find_Med(double *x, int P, int r){
	int number = r - P + 1;
	if(number % 2 == 0){
		double right = Select(x, P, r, number/2 + 1);
		double left = Select(x, P, r, number/2);
		return (right + left)/2;
	}
	else{
		return Select(x, P, r, number/2 + 1);
	}
}

double QD(double *x, int P, int r){		// ¥|¤À¦ì®t  Quartile Deviation, QD 
	int number = r - P + 1;
	if(number/2 == 0){
		double Q1 = Find_Med(x, P, number/2 - 1);
		printf("Q1 = %f\n", Q1);
		double Q3 = Find_Med(x, number/2, r);
		printf("Q3 = %f\n", Q3);
		return Q3 - Q1;	
	}
	else{
		double Q1 = Find_Med(x, P, number/2 - 1);
		printf("Q1 = %f\n", Q1);
		double Q3 = Find_Med(x, (number + 1)/2, r);
		printf("Q3 = %f\n", Q3);
		return Q3 - Q1;	
	}		
}
















