#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "complex.h"

int roundup(int a, int b){
    int c;
    if(a%b){
        a += (b-a%b);
    }
    c = a/b;
    return c;
}

complex FFTP(complex *y, complex *x, int N, int p){
    
    //declare some variables 
    int i, j, k, Np = N/p;     //i,j,k for index using

    //prime factorization
    //delcare some variables
    int nextp, nextN, psize;
    int *prime, *ptime;     //record prime and square number
    int Ntmp = N, index = 0, t, ctrl, final;

    //memory allocation for prime and its square number
    psize = roundup(N,2);
    prime = (int *) malloc( psize * sizeof(int));
    ptime = (int *) malloc( psize * sizeof(int));
        
    for(i = 2; i < psize; i++){
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
    if(!index){
        nextp = N;
        final = 1;              //as N = p, do the final FFT
        prime[index] = N;
        ptime[index] = 1;
    }else if(ptime[index-1] > 1){
        nextp = prime[index-1];
    }else{
        nextp = prime[index-2];
    }

    //FFTP
    if(Np == 1 || final){
        p = N;
        complex *w;
        w = (complex *) malloc( p * sizeof(complex) );
        for(i = 0; i < p; i++){
            w[i] = wtrans(N,i);
        }
        for(i = 0; i < p; i++){
            y[i].Re = 0;
            y[i].Im = 0;
            for(j = 0; j < p; j++){
                y[i] = Cadd(y[i], Cmulti( x[j] ,w[(i*j)%N] ) );
            } 
        }
        free(w);
    }else{

        //memory allocation for subseq
        nextN = N / nextp;
        complex **xsub, **ysub;
        xsub = (complex **) malloc( nextp * sizeof(complex) );
        ysub = (complex **) malloc( nextp * sizeof(complex) );
        for(i = 0; i < nextp; i++){
            xsub[i] = (complex *) malloc( nextN * sizeof(complex) );
            ysub[i] = (complex *) malloc( nextN * sizeof(complex) );
        }

        //assign x subseq and y
        for(i = 0; i < nextN; i++){
            for(j = 0; j < nextp; j++){
                xsub[i][j] = x[nextp*j + i];
                y[i*nextN + j].Re = 0;
                y[i*nextN + j].Im = 0;
            }
        }

        //do the subseq FFT
        for(i = 0;i < nextp; i++){
            FFTP( *(ysub+i), *(xsub+i), nextN, nextp);
        }

        //combine the subseq 
        complex *w;
        w = (complex *) malloc( nextp * sizeof(complex));
        for( j = 0; j < nextp; j++){
            w[j] = wtrans(N,i);
        }
        for(i = 0; i < nextp; i++){
            for( j = 0; j < nextN; j++){
                for(k = 0; k < nextp; k++){
                    y[ nextp * i + j ] = Cadd( y[ nextp * i + j ] , \
					Cmulti(ysub[k][j], w[ ( k* (i*N/p+j) ) %N] ) );
                }
            }
        }
        //free every memory
        free(w);
        for(i = 0; i < nextp; i++){
            free(xsub[i]);
            free(ysub[i]);
        }
        free(xsub);free(ysub);
        
    }

    free(prime);free(ptime);

    //end
    complex c;
    c.Re = 0;
    c.Im = 0;
    return c;
}
