#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
time_t timer;
int main() 
{
	int j, k, p = 20, n, L, M, N = 1<<p; 
	complex *y, *x, t, w, u;
	
	printf("N = %d\n",N);
	x = (complex *) malloc(N*sizeof(complex));
	y = (complex *) malloc(N*sizeof(complex));
	
	for(k=0;k<N;++k){
		x[k] = k;
	}
	
	timer = clock();
	// bit-reverse swap
	M = 1 << (p-1);
	j = 0;
	for(k=1;k<N-1;++k) {
		L = M;
		while(j>=L) {
			j -= L;
			L /= 2;
		}
		j += L;
		if(j>k) {
			t = x[k];
			x[k] = x[j];
			x[j] = t;
		}
	}

	// Butterfly Radix
	j = 1; 
	while(j<N) {
		w = 1.0;
		u = cexp(-I*M_PI/j);
		for(n=0;n<j;n++) {
			for(k=0;k<N;k+=2*j) {
				// k+n+j --> multiply w
				x[k+n+j] *= w;
				// FFT2 for k+n and k+n+j		
				t = x[k+n];
				x[k+n  ] = x[k+n]+x[k+n+j];
				x[k+n+j] =      t-x[k+n+j];
			}
			w *= u;
		}
		j = j << 1;
	}
	timer = clock() - timer;
//	for(k=0;k<N;++k){
//		printf("%d:%f+%f i\n",k,creal(x[k]),cimag(x[k]));
//	}
	
	
	printf("cost time = %d ms\n", timer);
	return 0;
}

