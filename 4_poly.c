#include <stdio.h>
#include <stdlib.h>	// for memory allocation
#include <time.h>	// for random seed
#include <math.h>	// for pow
int main() {
	int j, k, M = 10000000, N =20;
	double *c, x, p1, p2;
	time_t t1, t2;

	srand( time(NULL) );
	c = (double *) malloc(N*sizeof(double));
	x = 1.0*rand()/RAND_MAX;
	for(k=0;k<N;++k) {
		c[k] = 1.0*rand()/RAND_MAX;
	}
	t1 = clock();
	for(j=0;j<M;++j) {
		p1 = 0.0;
		for(k=0;k<N;++k) {
			p1 += c[k]*pow(x,k);
		}
	}
	t1 = clock()-t1;
	t2 = clock();
	for(j=0;j<M;++j) {
		p2 = c[N-1];	
		for(k=N-2;k>=0;--k) {
			p2 = p2*x+c[k];
		}
	}
	t2 = clock()-t2;
	printf("M1: %d ms, M2: %d ms\n", t1,t2);
	printf("Polynomial p(%f): M1: %f M2: %f,\ntheir difference is %e,\n", x, p1, p2, p1-p2);
	printf("where p(x)=\n");
	for(k=N-1;k>=0;k--) {
		printf("     %f * x^%d \n",c[k],k);
	}
	free(c);
	return 100;
}

