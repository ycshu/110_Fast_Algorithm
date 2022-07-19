#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "fft.h"

//declare the function FFT
int FFT(double *yr, double *yi, double *xr, double *xi, int N);

int main(){

    //declare variables
    //Suppose N is the total number, pqr are square of 235 respec.
    //Suppose 
    int k, n, N, p, q, r, ctrl;
    double *xr, *xi, *yr, *yi, t, a, b, an, bn, c, s;
    time_t T;

    xr = (double *) malloc(N*sizeof(double));
    xi = (double *) malloc(N*sizeof(double));
    yr = (double *) malloc(N*sizeof(double));
    yi = (double *) malloc(N*sizeof(double));

    printf("Memory Done! \n");

    for(k = 0;k<5;++k){
        printf("%d: %f \n",k,xr[k]);
    }

    T = clock();

    a = cos(2 * M_PI / N); b = sin(2 * M_PI / N);
    an = 1; bn = 0;
    for(n = 0; n < N; n++){
        yr[n] = 0.0;
        yi[n] = 0.0;
        c = 1; s = 0;
        for(k = 0; k<N; ++k){
            yr[n] += c*xr[k];
            yi[n] -= s*xr[k];
            t = c * an - bn * b;
            s = c * bn + s * an;
            c = t;
        }
        t = an * a - bn * b;
        bn = an * b + bn * a;
        an = t;
    }

    T = clock() - T;
    printf("%d ms for DFT of %d elements\n",T,N);

    free(xr);free(xi);free(yr);free(yi);

    return 0;
}

int FFT(double *yr, double *yi, double *xr, double *xi, int N){
    if(N == 2){
        yr[0] = xr[0]+xr[1];
        yi[0] = xi[0]+xi[1];
        yr[1] = xr[0]-xr[1];
        yr[1] = xi[0]-xi[1];
    } else {
        int k;
        double *yEr, *yEi, *yOr, *yOi;
        double *xEr, *xEi, *xOr, *xOi;
        double c, s;

        xEr = (double *) malloc((N/2)*sizeof(double));
        xEi = (double *) malloc((N/2)*sizeof(double));
        yEr = (double *) malloc((N/2)*sizeof(double));
        yEi = (double *) malloc((N/2)*sizeof(double));
        xOr = (double *) malloc((N/2)*sizeof(double));
        xOi = (double *) malloc((N/2)*sizeof(double));
        yEr = (double *) malloc((N/2)*sizeof(double));
        yOi = (double *) malloc((N/2)*sizeof(double));
        
        for(k = 0; k<N/2; ++k){
            xEr[k] = xr[2*k];
            xEi[k] = xi[2*k];
            xOr[k] = xr[2*k+1];
            xOi[k] = xi[2*k+1];

        }
        FFT(yEr, yEi, xEr, xEi, N/2);
        FFT(yOr, yOi, xOr, xOi, N/2);
        for(k = 0; k < N/2; ++k){
            c = cos(2*M_PI*k/N);
            s = -sin(2*M_PI*k/N);
            yr[k] = yEr[k] + (yOr[k]*c - yOi[k]*s);
            yi[k] = yEi[k] + (yOr[k]*s + yOi[k]*c);
            yr[N/2+k] = yEr[k] - (yOr[k]*c - yOi[k]*s);
            yi[N/2+k] = yEi[k] - (yOr[k]*s + yOi[k]*c);
        }
        
        free(xEr);free(xEi);
        free(xOr);free(xOi);
        free(yEr);free(yEi);
        free(yOr);free(yOi);
    }
    return 0;
}