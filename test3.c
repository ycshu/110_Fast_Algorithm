#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <string.h>
#include  <assert.h>



int main(){
	int i, j, k, N, M;
	double **A, *x, *u, **U, *b, **F, *r, *y, *z;
	time_t t;
	
	N = 4;
	
	b = (double*) malloc(M*sizeof(double));
	F = (double**)malloc((N-1)*sizeof(double*));
	F[0] = b;
	for(i = 1; i < N-1; ++i){
		F[i] = F[i-1] + N - 1;
	}
	
	F[1][0] = 26.0;
	printf("%f\n", b[N-1]);
	
	return 0;
}
