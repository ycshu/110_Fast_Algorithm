#include <time.h>
#include "fft.h"

int main(){
    int k, n, v, i;
	int L = 1<<20;
    complex *x, *y;

    x = (complex *) malloc( L * sizeof(complex) );
    y = (complex *) malloc( L * sizeof(complex) );

    for(k = 0;k<L;k++){
        x[k] = k;
        y[k] = 0.0;
    }
    
    time_t T = clock();
    printf("start FFT\n");
    FFT(y, x, L);
    printf("end FFT\n");
    T = clock() - T;
    
    printf("Use %d ms for doing FFT.\n",T);
    printf("Some results: \n");
    for(i = 0; i < 10; i++){
    	printf("y[%d] = %f + %f i  \n", i, creal(y[i]),cimag(y[i]));
	}
	
	free(y);
	free(x);

    return 0;
}