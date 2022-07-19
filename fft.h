#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

int roundup(int a, int b){
    int c;
    if(a%b){
        a += (b-a%b);
    }
    c = a/b;
    return c;
}

complex FFT(complex *y, complex *x, int N){
    //declare some variables 
    int i, j, k, tmp;     //i,j,k for index using

    //prime factorization
    //delcare some variables
    int nextp, nextN, psize;
    int *prime, *ptime;     //record prime and square number
    int Ntmp = N, index = 0, t, ctrl, final=0;

    //memory allocation for prime and its square number
    psize = roundup(N,2);
    prime = (int *) malloc( psize * sizeof(int));
    ptime = (int *) malloc( psize * sizeof(int));

    for(i = 2; i < N; i++){
        t = 0;
        while(Ntmp%i == 0){
            Ntmp /= i;
            t++;                //sum up square number
            ctrl = 1;           //switch the position of the array
        }

        if(ctrl == 1){
            prime[index] = i;   //prime        
            ptime[index] = t;   //square number of the prime resp.
            index++;
        }
        ctrl = 0;
    }
    prime[index] = N;
    ptime[index] = 1;

    if((prime[0] == N && ptime[0] == 1)||(!index)){  //as N = p, do the final FFT
        nextp = N;
        final = 1;
        //printf("1. N is %d, nextp is %d \n",N,nextp);
    }else if(ptime[index-1] >= 1){
        nextp = prime[index-1];
        //printf("2. N is %d, nextp is %d \n",N,nextp);
    }else{
        nextp = prime[index-2];
    }

    //FFT
    complex wN = cexp(-2*M_PI*I/N);

    complex *w;
    w = (complex *) malloc( N * sizeof(complex) );
    for( j = 0; j < N; j++){
            w[j] = cpow(wN,j);
    }
    if(final){      //final subseq, p = N
        for(i = 0; i < N; i++){
            for(j = 0; j < N; j++){
                y[i] += x[j] * w[(i*j)%N];
            } 
        }
    }
    else{

        complex *xsub, *ysub;
        xsub = (complex *) malloc( N * sizeof(complex) );
        ysub = (complex *) malloc( N * sizeof(complex) );

        //memory allocation for subseq
        nextN = N/nextp;

        //assign x subseq and set y zero seq
        for(i = 0; i < nextp; i++){
            tmp = nextN*i;
            for(j = 0; j < nextN; j++){
                index = j + tmp;
                xsub[index] = x[nextp*j + i];
                ysub[index] = 0;
            }
        }
        tmp = 0;index = 0;

        //do the subseq FFT
        for(i = 0; i < nextp; i++){
            FFT( ysub+i*nextN, xsub+i*nextN, nextN);
        }

        for(i = 0; i < nextp; i++){           //in how many parts
            tmp = nextN*i;
            for( j = 0; j < nextN; j++){      //seq index
                index = tmp + j;
                for(k = 0; k < nextp; k++){   //subseq index
                    y[index] += ysub[nextN*k+j] * w[(k*index)%N];
                }
            }
        }
        index = 0;
        free(ysub);
        free(xsub);
    }
    free(w);
    free(ptime);
    free(prime);

    return 0;
}

int max_prime(int N){
    int i, p = 1, ctrl = 0, Ntmp = N;
    for(i = 2; i < roundup(N,2); i++){
        while( Ntmp%i == 0){
            Ntmp /= i;
            ctrl = 1;
            p = i;
        }
    }
    if(ctrl == 0){
        p = N;
    }
    return p;
}