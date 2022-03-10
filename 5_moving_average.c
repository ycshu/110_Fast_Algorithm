#include <stdio.h>	// for printf function
#include <stdlib.h>	// for memory allocation
#include <time.h>	// for time calculation
#include <math.h>	// for sine and cosine functions
#include <string.h>

// #define TEST

#ifdef TEST
#include <assert.h>
#endif

int main(int argc, const char *argv[]) {
	// Declare all the variables
	int k, n, M, N;
	double *x, *y;
	time_t t;

	// Input the number N, M
    printf("N = ");
    scanf("%d", &N);	
    printf("M = ");
    scanf("%d", &M);
	
    double size = 2*M + 1;
	// Locate the memory for x, y;
    x = (double*)malloc(sizeof(double)*N);
    y = (double*)malloc(sizeof(double)*N);


#ifdef TEST
    double *y_test = (double*)malloc(sizeof(double)*N);
    memset(y, 0, sizeof(double)*N);
#endif
	
	// Initial setting for x[k] = sin(2*Pi*k/N), k = 0..N-1
    for(int i = 0; i < N; ++i){
        x[i] = sin(2*M_PI*i/(double)N);
    }
    memset(y, 0, sizeof(double)*N);

#ifdef TEST
    // use O(NM) moving average to perform test
    for(int n = 0; n < N; ++n){
        for(int k = n-M, upper_bound = n+M; k <= upper_bound; ++k){
            if(k >= 0 && k < N){
                y_test[n] += x[k];
            }
        }
        y_test[n] /= size;
    }
#endif
	
	t = clock();	
	// y[n] = sum(x[k], k=n-M..n+M)/(2M+1), n=0..N-1 
	// zero padding for x[k] if k<0 or k>=N
    //initialize y[0]
    for(int i = 0; i <= M; ++i) y[0] += x[i];
    

    // Use branchless method to reduce if-else
    for(int i = 1, lower = - M, upper = 1 + M; i < N; ++i, ++lower, ++upper){
        y[i] = (y[i - 1] - x[lower]) * (lower >= 0) + y[i-1]*(lower < 0);
        y[i] += x[i + M] * (upper < N);
    }
    
    for(int i = 0; i < N; ++i) y[i] /= size;


	// Output the results
	t = clock() - t;
    for(int i = 0; i < N; ++i){
#ifdef TEST
        assert(y[i] - y_test[i] < 0.001); // 0.001 for floating number bias
#endif
        printf("%.3f\n", y[i]);
    }
	printf("\n");

	// free the memory located by x, y
    free(x);
    free(y);

#ifdef TEST
    free(y_test);
#endif
	
	return 100;
}
