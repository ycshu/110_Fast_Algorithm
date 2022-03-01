#include <stdio.h>
#include <stdlib.h>	// for memory allocation
#include <time.h>	// for random seed

int main() {
	int i, j, N=4000;
	//double A[100][100], x[100], b[100];		// only for small integer!
	double **A, *x, *b;
	time_t t1, t2;

	x = (double *) malloc(N*sizeof(double));
	b = (double *) malloc(N*sizeof(double));	
	A = (double **) malloc(N*sizeof(double*));
	for(i=0;i<N;++i) {
		A[i] = (double *) malloc(N*sizeof(double));
	}
	t1 = clock();
	srand( time(NULL) );
	for(i=0;i<N;i++) {
		x[i] = 1.0*rand()/RAND_MAX;
		for(j=0;j<N;++j) {
			A[i][j] = 1.0*rand()/RAND_MAX;
		}
	}
	t1 = clock()-t1;

	t2 = clock();
	for(i=0;i<N;i++) {
		b[i] = 0.0;
		for(j=0;j<N;++j) {
			b[i] += A[i][j]*x[j];
		}
	}
	t2 = clock()-t2;
	printf("Generate %dx%d random matrix: %d ms\nCompute b=Ax: %d ms\n",N,N,t1,t2);
	free(x);
	free(b);
	return 0;
}
