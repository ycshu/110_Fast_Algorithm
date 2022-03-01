#include <stdio.h>
#include <time.h>
#include <math.h>

int main() {
	time_t t;
	unsigned int j, k, M=1<<28;
	double a = 0.5*M_PI-1e-9, r, b=M_PI, c, T;
	r = 1+1e-9;
	t = clock();
	for(k=0;k<M;k++) {
		a = sin(a) ;
	}
	t = clock()-t;
	T = 1.0*t/CLOCKS_PER_SEC;
	printf("a=%f, multiply: %d, %e %f\n",a,t,T/M, M/T);
	return 0;
}
