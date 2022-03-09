//
//  DFT
//
//  Created by Liu Yuting on 2022/3/9.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int main(int argc, const char * argv[]) {
    // Declare variables.
    int N=7;
    int n, k;
    int index;
   
    double *x, *yr, *yi, *cos_array, *sin_array;
    x=(double *) malloc(N*sizeof(double));
    yr=(double *) malloc(N*sizeof(double));
    yi=(double *) malloc(N*sizeof(double));
    cos_array=(double *) malloc(N*sizeof(double));
    sin_array=(double *) malloc(N*sizeof(double));
    for(n=0; n<N; n++) {
        //Initialize x.
        x[n]=n;
    }
    time_t t;
   
    t=clock();
   
    cos_array[0]=1.0;
    sin_array[0]=0.0;
    cos_array[1]=cos(2*M_PI/N);
    sin_array[1]=sin(2*M_PI/N);
   
    for(n=2; n<N; n++) {
        //Compute table of cos and sin.
        cos_array[n]=cos_array[1]*cos_array[n-1]-sin_array[1]*sin_array[n-1];
        sin_array[n]=cos_array[1]*sin_array[n-1]+sin_array[1]*cos_array[n-1];
    }
   
    for(k=0; k<N; k++) {
        //Compute DFT.
        yr[k]=0.0;
        yi[k]=0.0;
        for(n=0; n<N; n++) {
            //Compute yr[k], yi[k].
            index=(n*k)%N;
            yr[k]+=x[n]*cos_array[index];
            yi[k]+=x[n]*sin_array[index];
        }
    }
   
    t=clock()-t;
   
    printf("x:");
    for(n=0; n<N; n++) {
        //Print values of x.
        printf(" %f", x[n]);
    }
    printf("\n");
   
    printf("yr:");
    for(k=0; k<N; k++) {
        //Print values of yr.
        printf(" %f", yr[k]);
    }
    printf("\n");
   
    printf("yi:");
    for(k=0; k<N; k++) {
        //Print values of yi.
        printf(" %f", yi[k]);
    }
    printf("\n");
   
    printf("Time: %f s\n", 1.0*t/CLOCKS_PER_SEC);
   
    free(yi);
    free(yr);
    free(x);

    return 0;
}
